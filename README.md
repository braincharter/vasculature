# Braincharter MRI Vasculature Extraction
Region- and voxel-based probabilistic mapping of whole-brain veins and arteries (densities and diameter) in the MNI space. The initial release of this work is based on the following paper:

Bernier, M., Cunnane, S. C., & Whittingstall, K. (2018). The morphology of the human cerebrovascular system. Human Brain Mapping. https://doi.org/10.1002/hbm.24337

*The script is being currently updated for more stable results across scanners and image types.*

## Template

The initial release of the arterial and vascular template (based on 41 healthy time-of-flight and susceptibility weighted imaging) are available to download in the '/Template_MNI_2018' folder. template are in their early stage and will be updated regularly (details to follow).

The folders contains, for SWI and ToF:
- mean_data_VED.nii.gz: The output of the extraction itself, unthresholded, averaged across all subjects
- mean_data_VED_thresh.nii.gz: The output of the extraction itself, thresholded, averaged across all subjects
- CAT_data_VED.nii.gz: The output of the extraction itself, unthresholded, where all subjects are concatenated in a 4D file
- CAT_data_VED_thresh.nii.gz: The output of the extraction itself, thresholded, where all subjects are concatenated in a 4D file

TODO: Add the computed MNI regional templates 

TODO: Add the script to create the template

## Installation

Developement has been first inspired from the VMTK toolbox, but as of today the script has and is currently beiing actively revamped. It is currently used on Linux and MACOS. The dependences are:

- ITK >= 4.9, < 5.0 (do not forget to set ITK_DIR)
- python 2.7 (not tested with 3.5 and higher) with numpy, dipy, scikit-image (non exhaustive list)
- ANTs (https://github.com/ANTsX/ANTs). ITK 5.0 will be installed during ANTs compilation; ignore this version
- cmake + cmake-gui
- MACOS: realpath

After cloning the repository, go to the cplusplus folder, make a build directory, cd in build and "ccmake ../". Generate and make, then copy and replace the itkVedMain in the main directory by the one generated in the 'build' folder.

Put the main directory in you're PATH, and you're good to go!

## Running the script

To call the process:
- Param1: prefix of the data (ex: ToF.nii.gz -> ToF)
- Param2: suffix of the data (ex: ToF.nii.gz -> .nii.gz)
- Param3: Type of the data (ex: TOF or SWI or OTHER) -> SWI means dark vessels, ToF means bright vessels, OTHER won't skullstrip and skip some process and assume bright blood.

ex: " extract_vessels.sh ToF 'nii.gz' TOF "

This script is highly configurable. For now, unfortunately, the parameters need to be changed manually in that script. The most important one are when "ComputeVED.py" is called; This script most important parameter is the 'scale' parameter (-M), which should be put at half the maximum expected vessel size. Further down is the threshold parameters (postprocessing) that can be modified. It was automatic before, but was not consistant across image types.

### Output
(TODO)
