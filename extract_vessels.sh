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
    printf "extract_vessels.sh filename ext type(TOF or SWI) [-m --mask mask_filename] [-p --phase phase_filename] [-s --std denoised_STD_value] [-f --factor upsampleFactor] [-c --centerline]"
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
printf "Step 1. SkullStripping. \n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_autobox.${ext} ]; then
    3dAutobox -overwrite -prefix ${image}_autobox.${ext} -npad 3 ${image}.${ext}
    if [ "${imgType}" = "TOF" ]; then
        3dSkullStrip -overwrite -mask_vol -input ${image}_autobox.${ext} -prefix ${image}_mask.${ext} -blur_fwhm 2 -use_edge -avoid_vent -touchup
        3dcalc -overwrite -a ${image}_autobox.${ext} -b ${image}_mask.${ext} -expr "step(b)*abs(a)" -prefix ${image}_ss.${ext} -datum short
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
    fi
else
    printf "Autobox already exists for this subject.\n"
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

if [ ! -f ${image}_Ved.${ext} ]; then
    if [ "${imgType}" = "TOF" ]; then
        3dresample -overwrite -dxyz ${smalldim} ${smalldim} ${smalldim} -rmode Cu -prefix ${image}_upsampled.${ext} -inset ${image}_std${stdDenoised}_denoised.nii.gz
        ${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${smalldim} -M 6 -t 18 -n 10 -s 5 -w 25 --generate_scale -D 'scales'
    elif [ "${imgType}" = "SWI" ]; then
        3dresample -overwrite -dxyz ${smalldim} ${smalldim} ${smalldim} -rmode Cu -prefix ${image}_upsampled.${ext} -inset ${image}_Contrasted.nii.gz
        #3dLocalstat -nbhd 'RECT(1,1,1)' -stat max -prefix ${image}_mip.${ext} -overwrite ${image}_Contrasted.nii.gz 
        ${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${smalldim} -M 6 -t 18 -n 10 -s 5 -w 25
    fi
else
    printf "Vessel enhancing diffusion file already exists for this subject.\n"
fi

# Hack that stops the process to look at the VED results from parameter change (if parameters were changed)
# read -p "look at your images now, or press any key to pursue script"

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 4. VED post-processing \n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_Ved_Thresh.${ext} ]; then        
    #post-processing of Ved 
    dilatation_factor=`echo "scale=4; ${smalldim}*1.5" | bc `
    #3dmask_tool -input ${image}_Ved.${ext} -overwrite -prefix ${image}_ero.${ext} -dilate_input -1     
    3dresample -overwrite -dxyz ${smalldim} ${smalldim} ${smalldim} -rmode Cu -prefix ${image}_mask_up.${ext} -inset ${image}_mask.${ext}
    3dmask_tool -overwrite -input ${image}_mask_up.${ext} -prefix ${image}_mask_ero.${ext} -dilate_input -2
    3dcalc -overwrite -a ${image}_Ved.${ext} -b ${image}_mask_ero.${ext} -expr "step(b)*(1000000*a)" -prefix ${image}_Ved_boost.${ext}
    3dcalc -overwrite -a ${image}_Ved_boost.${ext} -expr "astep(a,50)*(1-astep(a, 1000))*a + astep(a, 1000)*1000" -prefix ${image}_Ved_stretched.${ext} -datum short
    fslmaths ${image}_Ved_stretched.${ext} -kernel sphere ${dilatation_factor} -dilF -eroF ${image}_Ved_filt.${ext}
    if [ "${imgType}" = "TOF" ]; then
        3dcalc -overwrite -a ${image}_Ved_filt.${ext} -expr "astep(a, 200)" -prefix ${image}_Ved_Thresh.${ext} -datum short
    elif [ "${imgType}" = "SWI" ]; then
        3dcalc -overwrite -a ${image}_Ved_filt.${ext} -expr "astep(a, 500)" -prefix ${image}_Ved_Thresh.${ext} -datum short
    fi
else
    printf "Processed VED files already exists for this subject.\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 5. Diameters extraction.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_diameters.${ext} ]; then        
    echo "Diameter Extraction"
    halfsize=`echo "scale=3;${smalldim} / 2.0" | bc`
    3dresample -overwite -dxyz ${halfsize} ${halfsize} ${halfsize} -rmode Cu -prefix ${image}_Ved_Thresh_HALF.${ext} -inset ${image}_Ved_Thresh.${ext}
    ${scriptpath}/ExtractDiameter.py ${image}_Ved_Thresh_HALF.${ext} ${image}_diameters_HALF.${ext}
    3dresample -overwrite -master ${image}_Ved_Thresh.${ext} -rmode Cu -prefix ${image}_diameters.${ext} -inset ${image}_diameters_HALF.${ext}
    rm -rf *HALF*
else
    printf "Diameters file already exists for this subject.\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 6. Centerlines extraction.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_centerdia.${ext} ]; then        
    echo "Centerlines Extraction"
    ${scriptpath}/ExtractCenterline.py ${image}_Ved_Thresh.${ext} ${image}_skel.${ext} 
    echo "Diameters in centerlines only"
    3dcalc -a ${image}_skel.nii.gz -b ${image}_diameters.${ext}  -expr "step(a)*b" -prefix ${image}_centerdia.${ext} -datum float
else
    printf "Centerline file already exists for this subject.\n"
fi

printf "Pipeline process completed.\n"
unset AFNI_NIFTI_TYPE_WARN
