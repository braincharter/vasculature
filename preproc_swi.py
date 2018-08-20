#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Combine the phase image to the magnitude image for SWI 

Usage:
  preproc_swi.py <input_magnitude> <input_phase>
"""

import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)

import os
import argparse
import nibabel as nib
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.signal import convolve2d as conv2d
from scipy.signal import hann as hann
from numpy.fft import ifftn, fftn, fftshift, ifftshift

class SWIPhaseBoost:
    def __init__(self, 
                 input_mag_filename=None,
                 input_phase_filename=None,
                 from_cmd=False, ):

        self._input_mag = input_mag_filename
        self._input_phase = input_phase_filename

        if not from_cmd:
            self.valid_arg()

    def parse_arg(self, args):
        self._input_mag = args.input_mag
        self._input_phase = args.input_phase

    def valid_arg(self):

        if not self._input_mag:
            raise Exception('SWI magnitude input missing.')

        if not self._input_phase:
            raise Exception('SWI phase input missing.')

  
    #% Phase unwrap (deformation induced due to cyclic Fourier thingy) 
    def LaplacianUnwrap(self, Original_phase, VoxelSizes, padding_pixels):
        # Local Variables: AcqSpacing, padding_pixels, SliceSpacing, k, m, n, k2, Original_phase, size_data
        AcqSpacing = VoxelSizes[0]
        SliceSpacing = VoxelSizes[2]
        size_data_original = Original_phase.shape
        
        #% Zero-padding of the images
        Original_phase = np.lib.pad(Original_phase, ((padding_pixels,padding_pixels),(padding_pixels,padding_pixels),(padding_pixels,padding_pixels)), 'wrap')
        
        #%Mask=padarray(Mask,[padding_pixels/2 padding_pixels/2 padding_pixels/2]);
        size_data = Original_phase.shape
        
        #% Calculation of the k**2 matrix (we assume an isotropic inplane resolution)
        k2 = np.zeros(size_data)
        for k in np.arange(0., size_data[0]):
            for m in np.arange(0., size_data[1]):
                for n in np.arange(0., size_data[2]):
                    k2[int(k),int(m),int(n)]=(np.double(k)-np.double(np.floor(size_data[0]/2))-1.0)**2 + (np.double(m)-np.double(np.floor(size_data[1]/2))-1.0)**2 + ((np.double(n)-np.double(np.floor(size_data[2]/2))-1.0)*(SliceSpacing/AcqSpacing))**2 

        k2 = np.array(k2)

        #% Equation 13 from the paper from Li et al.
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3947438/
        Original_phase = fftn((np.cos(Original_phase)*ifftn((ifftshift(k2)*fftn(np.sin(Original_phase))))-np.sin(Original_phase)*ifftn((ifftshift(k2)*fftn(np.cos(Original_phase))))))/ifftshift(k2)
      
        #% To prevent errors arising from the k=0 point
        Original_phase[np.isnan(Original_phase)] = 0.
        Original_phase[np.isinf(Original_phase)] = 0.
       
        #% Back to the image domain
       
        #%Original_phase=ifftn(Original_phase)*Mask;
        Original_phase = ifftn(Original_phase)
       
        #% We remove the zero-padding
        Original_phase = Original_phase[[slice(padding_pixels, -padding_pixels) for _ in Original_phase.shape]]

        if Original_phase.shape != size_data_original:
            raise Exception('Padding is weird')


        print "Laplacian unwrap done"
        return np.array(np.real(Original_phase))
   
    #% Homodyne phase filtering (Haacke et al. 2005)
    def Homodyne2D(self, OriginalImage, FilterSize):
        ImageSize = OriginalImage.shape
        print "Image size: " + str(ImageSize)
        
        if FilterSize > ImageSize[0]:
            FilterSize = ImageSize[0]
        FilterSize_minus = int(FilterSize-1)
        Image_size_x = int(ImageSize[0])
        Image_size_y = int(ImageSize[1])
        Image_size_z = int(ImageSize[2])
        #% Create Hann filter
        hann_filt = np.ones((FilterSize_minus, FilterSize_minus, Image_size_z))
        #print "Filter size: " + str(hann_filt.shape)

        for axis, axis_size in enumerate(hann_filt.shape):
            filter_shape = [1, ] * hann_filt.ndim
            filter_shape[axis] = axis_size
            window = np.hanning(axis_size).reshape(filter_shape)
            hann_filt *= window
        
        mask = np.zeros((Image_size_x, Image_size_y, Image_size_z))
        
        index_x_minus = int(np.round(((Image_size_x-FilterSize)/2.+1.)))
        index_y_minus = int(np.round(((Image_size_y-FilterSize)/2.+1.)))
        index_x_plus = int(np.round(((Image_size_x+FilterSize)/2.)))
        index_y_plus = int(np.round(((Image_size_y+FilterSize)/2.)))
        
        mask[ index_x_minus:index_x_plus , index_y_minus:index_y_plus, : ] = hann_filt

        print "Filter size: " + str(hann_filt.shape)
        
        SmoothedImage = OriginalImage*mask

        #% Back to image space
        OriginalImage = ifftn(ifftshift(OriginalImage))
        SmoothedImage = ifftn(ifftshift(SmoothedImage))
        CorrectedImage = OriginalImage/SmoothedImage
        return CorrectedImage, OriginalImage, SmoothedImage
   
    #% PADRE phase/magnitude reconstruction
    def PADRE(self, MagnData,PhaseData,alpha,beta,sigma):  
        ImageSize = MagnData.shape
        Mask = np.ones((ImageSize))
        #condlist = [PhaseData>=sigma*np.pi, PhaseData<-sigma*np.pi, PhaseData>0]
        #choicelist = [np.exp(-alpha*(abs(PhaseData[PhaseData>=sigma*np.pi])-(sigma*np.pi))**beta), np.exp(-alpha*(abs(PhaseData[PhaseData<-sigma*np.pi])-(sigma*np.pi))**beta), 1]
        #Mask = np.select(condlist, choicelist)
        Mask[PhaseData>=sigma*np.pi] = np.exp(-alpha*(abs(PhaseData[PhaseData>=sigma*np.pi])-(sigma*np.pi))**beta)
        Mask[PhaseData<-sigma*np.pi] = np.exp(-alpha*(abs(PhaseData[PhaseData<-sigma*np.pi])-(sigma*np.pi))**beta)
        Mask[PhaseData>0] = 1
        Mask = Mask.real.astype(float)
        MagnData = MagnData.astype(np.float)
        return Mask*Mask*Mask*Mask*MagnData, Mask*Mask*MagnData, Mask

    def run(self):

        img_mag = nib.load(self._input_mag)
        img_phase = nib.load(self._input_phase)
        
        MagnData = np.abs(img_mag.get_data().astype(np.float))
        PhaseData = img_phase.get_data().astype(np.float)
        PhaseData = PhaseData/1000.0  # Because GE gave me scaled angles for the phase

        #PhaseData = np.angle(img_phase.get_data().astype(np.float))
        voxel_sizes = img_phase.header.get_zooms()

        filter_size = np.round(0.15*img_phase.shape[0])

        ### We now assume there is no empty spacing between slices
       
        unwarped_file = 'swi_phase_unwarped.nii.gz'
        if not os.path.exists(unwarped_file): 
            print 'Phase unwrapping'
            PhaseData_unwrap = self.LaplacianUnwrap(PhaseData, voxel_sizes, 8)
            nib.save(nib.Nifti1Image(PhaseData_unwrap, img_phase.get_affine()), unwarped_file)
        else:
            PhaseData_unwrap = nib.load(unwarped_file).get_data().astype(np.float)
        
        #PhaseData_unwrap = PhaseData
        hemofilt_file = 'swi_phase_filt.nii.gz'
        if not os.path.exists(hemofilt_file): 
            print 'Homodyne filtering'
            ## Don't remember why I do that.. :|
            PhaseData = MagnData * np.exp(1j * PhaseData_unwrap)
            PhaseData = fftshift(fftn(PhaseData))
            ##
            Corrected, orig, smooth = self.Homodyne2D(PhaseData, filter_size)
            PhaseData = -1*np.angle(Corrected)
            nib.save(nib.Nifti1Image(PhaseData, img_phase.get_affine()), hemofilt_file)
            
            nib.save(nib.Nifti1Image(np.real(orig), img_phase.get_affine()), 'test_orig.nii.gz')
            nib.save(nib.Nifti1Image(np.real(smooth), img_phase.get_affine()), 'test_smooth.nii.gz')
        else:
            PhaseData = nib.load(hemofilt_file).get_data().astype(np.float) 

        #Filter the phase data()
        swi_mask = 'swi_ss_mask.nii.gz'
        swi_ss_file = 'swi_ss.nii.gz'
        swi_phase_ss = 'swi_phase_ss.nii.gz'
        if not os.path.exists(swi_phase_ss) or not os.path.exists(swi_ss_file): 
            os.system("bet " + self._input_mag + ' ' + swi_ss_file + ' -A2 ' + hemofilt_file + ' -R -m -f 0.3')
            os.system("3dcalc -overwrite -a " + hemofilt_file + " -b " + swi_mask + " -expr '(step(b)*a)' -prefix " + swi_phase_ss)
        
        #print "Start denoising"
        #os.system("dipy_nlmeans.py -std 30 -mask " + swi_mask + " -o swi_ss " + swi_ss_file)
        #print "end"

        MagnData = nib.load(swi_ss_file).get_data().astype(np.float)        
        PhaseData = nib.load(swi_phase_ss).get_data().astype(np.float)
        MagnData_file = 'swi_factor2_final.nii.gz'
        PhaseMask_file = 'swi_phase_mask.nii.gz'
        if not os.path.exists(MagnData_file): 
            print 'PADRE processing'
            MagnData4times, MagnData, PhaseMask = self.PADRE(MagnData,-1*PhaseData,1.6,1.2,0); 
            nib.save(nib.Nifti1Image(MagnData, img_mag.get_affine()), MagnData_file)
            nib.save(nib.Nifti1Image(MagnData4times, img_mag.get_affine()), 'swi_factor4_final.nii.gz')
            nib.save(nib.Nifti1Image(PhaseMask, img_mag.get_affine()), PhaseMask_file)
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input_mag", type=str,
                        help="The input SWI magnetude file name (relative path)")

    parser.add_argument("input_phase", type=str,
                        help="The input SWI magnetude file name (relative path)")

    exec_shell = SWIPhaseBoost(from_cmd=True)
    exec_shell.parse_arg(parser.parse_args())
    exec_shell.run()
