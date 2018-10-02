#!/bin/bash

###############################################################################
#                                   Input
###############################################################################
image=$1
ext=$2
imgType=$3

scriptpath=`dirname $0`
echo $scriptpath
###############################################################################
#                                   Parse argument
###############################################################################

export AFNI_NIFTI_TYPE_WARN=NO

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 0. Parse the arguments.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"

if [ "${1}" = "-h" ]; then
    printf "extract_vessels.sh filename ext type(TOF or SWI or OTHER) [-m --mask mask_filename] [-p --phase phase_filename] [-s --std denoised_STD_value] [-f --factor upsampleFactor] [-c --centerline]"
    exit
fi

if [ -z "$image" ]; then
    printf "Input image name is required"
    exit
fi

if [ -z "$ext" ]; then
    printf "Extension is required"
    exit
fi

if [ -z "$imgType" ]; then
    printf "Image type {TOF|SWI} is required"
    exit
else
    if [ "${imgType}" = "TOF" ]; then
        printf "Time-of-flight is selected. \n"
    elif [ "${imgType}" = "SWI" ]; then
        printf "SWI is selected. \n"
    elif [ "${imgType}" = "OTHER" ]; then
    printf "SWI is selected. \n"
    else
        printf "ImageType Error, must be {TOF|SWI}"
        exit
    fi
fi

############################### Optional #####################################

function checkArg {
    if [ -z "$1" ]; then
        printf "Must provide a value after key.\n"
        exit
    fi
    case $1 in
        -m|--mask)    
        printf "Value must not be a key.\n"    
        exit 
        ;;
        -p|--phase)
        printf "Value must not be a key.\n"    
        exit 
        ;;
        -f|--factor)
        printf "Value must not be a key.\n"    
        exit 
        ;;
        -c|--centerline)
        printf "Value must not be a key.\n"    
        exit 
        ;;        
        -d|--diameters)
        printf "Value must not be a key.\n"    
        exit 
        ;;
        *)
        ;;
    esac
}

if(( $(echo "$# > 3" | bc -l) )); then

    # Pass the 3 required arguments
    shift; shift; shift;

    while [[ $# -gt 0 ]]
    do
    key="$1"
    case $key in
        -m|--mask)        
        mask="$2"
        checkArg ${image}_mask
        shift # past argument
        ;;
        -p|--phase)
        phase="$2"
        checkArg ${phase}
        shift # past argument
        ;;
        -f|--factor)
        upsampleFactor="$2"
        checkArg ${upsampleFactor}
        shift # past argument
        ;;
        -c|--centerline)
        getCenterline=true
        ;;
        -d|--diameters)
        getDiameters=true
        ;;
        *)
        echo "unknow args!"
        exit
        ;;
    esac
    shift # past argument or value
    done
    

else
    printf "No optional args. \n"
fi

if [ -z "$upsampleFactor" ]; then
    printf "No upsample factor has been passed as argument.\nData won't be affect except a simple resample to isotropic for diameter extraction.\n"
else 
    if(( $(echo "${upsampleFactor} < 2" | bc -l) )); then
        printf "Error: Your upsample factor must be greater than 1.\n"
        exit
    fi
fi

echo INPUT FILE          = "${image}"
echo EXTENSION           = "${ext}"
echo IMAGE TYPE          = "${imgType}"
echo FILE MASK           = "${image}_mask"
echo FILE PHASE          = "${phase}"
echo STD                 = "${stdDenoised}"
echo UPSAMPLED FACTOR    = "${upsampleFactor}"
echo GENERATE DIAMETERS = "${getDiameters}"
echo GENERATE CENTERLINE = "${getCenterline}"

###############################################################################
#                               Execution
###############################################################################


printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 0. Autobox. \n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_autobox.${ext} ]; then
   if [ "${imgType}" = "OTHER" ]; then
      cp -rf ${image}.${ext} ${image}_autobox.${ext}
   else
        3dAutobox -overwrite -prefix ${image}_autobox.${ext} -npad 3 ${image}.${ext}
   fi

else
    printf "Autobox already exists for this subject.\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 1. SkullStripping. \n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_mask.${ext} ]; then
    if [ "${imgType}" = "TOF" ]; then
        #3dSkullStrip -overwrite -mask_vol -input ${image}_autobox.${ext} -prefix ${image}_mask.${ext} -blur_fwhm 2 -use_edge -avoid_vent -touchup
        #3dcalc -overwrite -a ${image}_autobox.${ext} -b ${image}_mask.${ext} -expr "step(b)*abs(a)" -prefix ${image}_ss.${ext} -datum short
        bet ${image}_autobox.${ext} ${image}_ss.${ext} -R -m -f 0.95
        mv ${image}_ss_mask.${ext} ${image}_mask.${ext}  
    elif [ "${imgType}" = "SWI" ]; then
        if [ ! -f ${image}_mask.${ext} ]; then
            if [ -f ${phase}.${ext} ] ; then
                #read -p "Skullstripping ${image}.${ext}"
                min=`3dBrickStat -slow -min ${phase}.${ext}  | awk '{print $1}'`
                max=`3dBrickStat -slow -max ${phase}.${ext} | awk '{print $1}'`
                3dresample -input ${phase}.${ext} -master ${image}_autobox.${ext} -prefix ${phase}_autobox.${ext}
                3dcalc -overwrite -a ${phase}_autobox.${ext} -expr "(a+${min})*1000" -prefix ${phase}_weird.${ext}
                bet ${image}_autobox.${ext} ${image}_ss.${ext} -A2 ${phase}_weird.${ext} -R -m -f 0.3  
                mv ${image}_ss_mask.${ext} ${image}_mask.${ext}    
            else
                bet ${image}_autobox.${ext} ${image}_ss.${ext} -R -m -f 0.65
                mv ${image}_ss_mask.${ext} ${image}_mask.${ext}  
            fi
        else
            printf "Mask already exists. Putting in same space as magnitude.\n"
            3dresample -master ${image}_autobox.${ext} -input ${image}_mask.${ext} -prefix ${image}_mask.${ext} -overwrite
        fi
    elif [ "${imgType}" = "OTHER" ]; then
        3dcalc -a ${image}_autobox.${ext} -expr "step(a)" -prefix ${image}_mask.nii.gz
        cp -rf ${image}_autobox.${ext} ${image}_ss.${ext}
    fi
else
    printf "Mask already exists for this subject.\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 2. Denoising with nl means and refit to space ORIG\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"

if [ ! -f ${image}_std${stdDenoised}_denoised.nii.gz ]; then
    ${scriptpath}/dipy_nlmeans.py -std 20 -mask ${image}_mask.${ext} -o ${image}_std${stdDenoised} ${image}_autobox.${ext}
else
    printf "Denoise file already exists for this subject.\n"
fi



printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 3. Vessel Extracting with VED.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"

if [ "${imgType}" = "SWI" ]; then
    if [ ! -f ${image}_Contrasted.nii.gz ]; then        
        ${scriptpath}/InvertContrast.py ${image}_std${stdDenoised}_denoised.nii.gz ${image}_mask.${ext} ${image}_Contrasted.nii.gz        
    else
        printf "Inverted contrast file already exists for this subject.\n"
    fi
else
    printf "No needs to invert the constrast for TOF.\n"
fi

# Checking for smallest pix dim.
smalldim=99999.0

xdim=`3dinfo ${image}_autobox.${ext} | grep "[R]" | grep "\-step\-" | awk '{print $9}'` ;  echo $xdim
ydim=`3dinfo ${image}_autobox.${ext} | grep "[A]" | grep "\-step\-" | awk '{print $9}'` ;  echo $ydim
zdim=`3dinfo ${image}_autobox.${ext} | grep "[I]" | grep "\-step\-" | awk '{print $9}'` ;  echo $zdim

if(( $(echo "${xdim} < $smalldim" | bc -l) )); then
    smalldim=$xdim
fi
 
if(( $(echo "${ydim} < $smalldim" | bc -l) )); then
    smalldim=$ydim
fi

if(( $(echo "${zdim} < $smalldim" | bc -l) )); then
    smalldim=$zdim
fi


printf "Smallest dim : %s \n" ${smalldim}

# Hack to change smallest diameter value (added to original)
# smalldim=0.599
a=0.5
b=`echo "scale=3; ${smalldim} * 1.2" | bc`
small_scale=$(( a > b ? a : b ))
large_scale=`echo "scale=3; ${small_scale} * 3" | bc`
large_scale_clarity=`echo "scale=3; ${small_scale} * 1.5" | bc`
 
if [ ! -f ${image}_Ved.${ext} ]; then
    if [ "${imgType}" = "TOF" ]; then
        3dresample -overwrite -dxyz ${smalldim} ${smalldim} ${smalldim} -rmode Cu -prefix ${image}_upsampled.${ext} -inset ${image}_std${stdDenoised}_denoised.nii.gz
        ${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${small_scale} -O -M ${large_scale} -t 0 -n 15 -s 2 -w 90 -I --out_folder "./${image}_iterations" 
        #${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${smalldim} -M 6 -t 18 -n 10 -s 5 -w 25 #--generate_scale -D 'scales'
    elif [ "${imgType}" = "SWI" ]; then
        3dresample -overwrite -dxyz ${smalldim} ${smalldim} ${smalldim} -rmode Cu -prefix ${image}_upsampled.${ext} -inset ${image}_Contrasted.nii.gz
        ${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${small_scale} -O -M ${large_scale} -t 0 -n 15 -s 2 -w 90 -I --out_folder "./${image}_iterations" 
        #${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${smalldim} -M 6 -t 18 -n 10 -s 5 -w 25
    elif [ "${imgType}" = "OTHER" ]; then
        3dresample -overwrite -dxyz ${smalldim} ${smalldim} ${smalldim} -rmode Cu -prefix ${image}_upsampled.${ext} -inset ${image}_std${stdDenoised}_denoised.nii.gz
        ${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${small_scale} -O -M ${large_scale_clarity} -t 0 -n 15 -s 2 -w 90 -I --out_folder "./${image}_iterations" 
    fi
    rm -rf ./${image}_iterations
else
    printf "Vessel enhancing diffusion file already exists for this subject.\n"
fi

# Hack that stops the process to look at the VED results from parameter change (if parameters were changed)
# read -p "look at your images now, or press any key to pursue script"

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 4. VED post-processing \n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_newVed_Thresh.${ext} ]; then        
    #post-processing of Ved 
    
    if [ ! -d ${image}_Ved_scales_nofilters ]; then
        mkdir ${image}_Ved_scales_nofilters 
    else 
        rm -rf ${image}_Ved_scales_nofilters/*
    fi
    mv -f Scale_NOW* ${image}_Ved_scales_nofilters/
    3dTcat -prefix ${image}_Ved_scales_nofilters/${image}_Ved_scales_nofilters.nii.gz  ${image}_Ved_scales_nofilters/Scale_NOW*
    #rm -rf Scale_NOW*
    mv -f ${image}_Ved_scales_nofilters/${image}_Ved_scales_nofilters.nii.gz ./
    3dTstat -max -prefix ${image}_Ved_nofilters_fromscales.nii.gz ${image}_Ved_scales_nofilters.nii.gz

    if [ ! -d ${image}_Ved_scales ]; then
        mkdir ${image}_Ved_scales 
    else 
        rm -rf ${image}_Ved_scales/*
    fi
    mv -f Scale_* ${image}_Ved_scales/
    3dTcat -prefix ${image}_Ved_scales/${image}_Ved_scales.nii.gz ${image}_Ved_scales/Scale*
    #rm -rf Scale_NOW*
    mv -f ${image}_Ved_scales/${image}_Ved_scales.nii.gz ./
    3dTstat -max -prefix ${image}_newVed.nii.gz ${image}_Ved_scales.nii.gz


    3dresample -overwrite -dxyz ${smalldim} ${smalldim} ${smalldim} -rmode Cu -prefix ${image}_mask_up.${ext} -inset ${image}_mask.${ext}
    #3dmask_tool -overwrite -input ${image}_mask_up.${ext} -prefix ${image}_mask_ero.${ext} -dilate_input -5
    
    ### Tricky thresh
    3dcalc -overwrite -a ${image}_newVed.${ext} -b ${image}_mask_up.${ext} -expr "step(b)*(a)" -prefix ${image}_newVed_masked.${ext}
    3dcalc -overwrite -a ${image}_newVed_masked.${ext} -expr "astep(a,200)*(1-astep(a, 1000))*a + astep(a, 1000)*1000" -prefix ${image}_newVed_stretched.${ext} -datum short
    3dcalc -overwrite -a ${image}_newVed_masked.${ext} -expr "astep(a, 50)" -prefix ${image}_newVed_50Thresh.${ext} -datum short
    3dcalc -overwrite -a ${image}_newVed_masked.${ext} -expr "astep(a, 100)" -prefix ${image}_newVed_100Thresh.${ext} -datum short
    3dcalc -overwrite -a ${image}_newVed_masked.${ext} -expr "astep(a, 200)" -prefix ${image}_newVed_200Thresh.${ext} -datum short
    3dcalc -overwrite -a ${image}_newVed_masked.${ext} -expr "astep(a, 300)" -prefix ${image}_newVed_300Thresh.${ext} -datum short
    3dcalc -overwrite -a ${image}_newVed_masked.${ext} -expr "astep(a, 400)" -prefix ${image}_newVed_400Thresh.${ext} -datum short
    3dcalc -overwrite -a ${image}_newVed_masked.${ext} -expr "astep(a, 500)" -prefix ${image}_newVed_500Thresh.${ext} -datum short
    3dcalc -overwrite -a ${image}_newVed_masked.${ext} -expr "astep(a, 600)" -prefix ${image}_newVed_600Thresh.${ext} -datum short
    3dcalc -overwrite -a ${image}_newVed_masked.${ext} -expr "astep(a, 700)" -prefix ${image}_newVed_700Thresh.${ext} -datum short
    3dcalc -overwrite -a ${image}_newVed_masked.${ext} -expr "astep(a, 800)" -prefix ${image}_newVed_800Thresh.${ext} -datum short
    3dcalc -overwrite -a ${image}_newVed_masked.${ext} -expr "astep(a, 100)" -prefix ${image}_newVed_Thresh.${ext} -datum short
    
    
else
    printf "Processed VED files already exists for this subject.\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 5. Diameters extraction.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_diameters.${ext} ]; then        
    echo "Diameter Extraction"
    halfsize=`echo "scale=3; ${smalldim} / 2.0" | bc`
    3dresample -overwrite -dxyz ${halfsize} ${halfsize} ${halfsize} -rmode Cu -prefix ${image}_Ved_Thresh_HALF.${ext} -inset ${image}_newVed_Thresh.${ext}
    ${scriptpath}/ExtractDiameter.py ${image}_Ved_Thresh_HALF.${ext} ${image}_diameters_HALF.${ext}
    3dresample -overwrite -master ${image}_newVed_Thresh.${ext} -rmode Cu -prefix ${image}_diameters.${ext} -inset ${image}_diameters_HALF.${ext}
    rm -rf *HALF*
else
    printf "Diameters file already exists for this subject.\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 6. Centerlines extraction.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_centerdia.${ext} ]; then        
    echo "Centerlines Extraction"
    ${scriptpath}/ExtractCenterline.py ${image}_newVed_Thresh.${ext} ${image}_skel.${ext} 
    echo "Diameters in centerlines only"
    3dcalc -a ${image}_skel.nii.gz -b ${image}_diameters.${ext}  -expr "step(a)*b" -prefix ${image}_centerdia.${ext} -datum float
else
    printf "Centerline file already exists for this subject.\n"
fi

printf "Pipeline process completed.\n"
unset AFNI_NIFTI_TYPE_WARN
