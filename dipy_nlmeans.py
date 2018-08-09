#! /usr/bin/env python2.7
# Wrapper for the nlmeans python wrapper function

from __future__ import print_function

import nibabel as nib
import numpy as np
import os
import argparse

from dipy.denoise.non_local_means import non_local_means as nlmeans
from dipy.denoise.noise_estimate import estimate_sigma, piesno

DESCRIPTION = """
    Loads an image and launch dipy NLMEANS
    """


def buildArgsParser():

    p = argparse.ArgumentParser(description=DESCRIPTION)
    p.add_argument('input', action='store', metavar='DIR',
                   help='path of the image file to denoise')

    p.add_argument('-std', action='store', dest='std',
                   metavar='std', required=False, default=0, type=int,
                   help='standard deviation of the noise')

    p.add_argument('-mask', action='store', dest='mask', 
                   metavar='mask', required=False, default=None, type=str,
                   help='Mask to accelerate computation.')  

    p.add_argument('-o', action='store', dest='savename',
                   metavar='savename', required=False, default=None, type=str,
                   help='path and prefix for the saved denoised file. The name is always appended with _denoised.nii.gz')

    p.add_argument('-avg', action='store', dest='avg',
                   metavar='averaging', required=False, default=1, type=int,
                   help='Set to 0 to use Standard weighting. Set to 1 to use Rician weighting')



    return p


def main():
    parser = buildArgsParser()
    args = parser.parse_args()

    if args.savename is None:
        temp, ext = str.split(os.path.basename(args.input), '.', 1)
        filename = os.path.dirname(os.path.realpath(args.input)) + '/' + temp

    else:
        filename = args.savename
   

    # If the file already exists, check if the user want to overwrite it. If not, simply exit.
    filename += '_denoised.nii.gz'

    print ("Now denoising", os.path.realpath(args.input))
    
    vol = nib.load(args.input)
    affine = vol.get_affine()
    data = vol.get_data()

    if args.std is 0:
        print ("Compute sigma with piesno")
        sigma = piesno(data, N=1, return_mask=False)
        print ("sigma is: ", sigma)
        sigma = np.mean(sigma)
        print ("sigma is: ", sigma)
        
    else:
        print ("Sigma is: ", args.std) 
        sigma = args.std

    if args.mask is None:
        print ("No mask.. you should")
    else:
        print ("Mask is: ", args.mask)
        img_mask = nib.load(args.mask)
        mask_t1 = img_mask.get_data()

    img_denoised = nlmeans(data, sigma, mask=mask_t1)

    nib.save(nib.Nifti1Image(img_denoised, affine), filename)
    print ("Denoised file was saved as", filename)


if __name__ == "__main__":
    main()
