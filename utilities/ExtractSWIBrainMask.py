#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extract the brain mask from a SWI image.

Usage:
  ExtractSWIBrainMask mag <mag_path> <inverted> <brain> <brain_mask>
  <low_thresh> [options]
  ExtractSWIBrainMask phase <mag_path> <phase_path> <brain> <brain_mask> [
  options]
"""

import argparse
from subprocess import call

import nibabel as nib
import numpy as np
from scipy.ndimage.morphology import binary_fill_holes


class ExtractSWIBrainMask:
    def __init__(self, mode=None, mag_path=None, inverted_path=None,
                 brain=None, brain_mask=None, phase_path=None,
                 low_threshold=80, from_cmd=False):

        self._mode = mode
        self._mag_path = mag_path
        self._inverted_path = inverted_path
        self._brain = brain
        self._brain_mask = brain_mask
        self._phase_path = phase_path
        self._low_thresh = low_threshold

        if not from_cmd:
            self.valid_arg()

    def parse_arg(self, args):
        self._mode = args.mode
        self._mag_path = args.mag_path
        self._inverted_path = args.inverted_path
        self._brain = args.brain
        self._brain_mask = args.brain_mask
        self._phase_path = args.phase_path
        self._low_thresh = args.low_thresh

        self.valid_arg()

    def valid_arg(self):

        mode_type = ['phase', 'mag']
        if self._mode not in mode_type:
            raise Exception('Mode needs to be : mag or phase')

        if not self._mag_path:
            raise Exception('You need to provide an input file name for the '
                            'magnitude image.')

        if not self._brain:
            raise Exception('You need to provide an output file name for the '
                            'skullstrip image.')

        if not self._brain_mask:
            raise Exception('You need to provide an output file name for the'
                            'brain mask image.')

        if self._mode == 'mag':

            if not self._inverted_path or not self._low_thresh:
                raise Exception("For 'mag' mode, you need to provide :"
                                "--inverted_path and --low_thresh.")

        elif self._mode == 'phase':

            if not self._phase_path:
                raise Exception("For 'phase' mode, you need to provide :"
                                "--phase_path.")

    def run(self):

        if self._mode == 'mag':

            bet_mag(self._mag_path, self._inverted_path, self._brain,
                    self._brain_mask, self._low_thresh)

        elif self._mode == 'phase':

            bet_phase(self._phase_path, self._mag_path, self._brain,
                      self._brain_mask)


def bet_mag(mag_path, invert_path, brain_path, brain_mask_path, low_thresh):
    img = nib.load(mag_path)
    data = img.get_data()

    # Those lines are equivalent to this AFNI call :
    # 3dcalc -expr "step(a-$lowval) * ($maxval-a)"
    # Colors are inverted based on the threshold value and the maximum
    # value of the data
    norm_data = data - low_thresh
    invert_data = np.zeros(data.shape)
    invert_data[norm_data > 0] = np.max(data) - data[norm_data > 0]

    nib.save(nib.Nifti1Image(invert_data, img.get_affine(),
                             img.get_header()), invert_path)

    call(['3dSkullStrip',
          '-input', invert_path,
          '-orig_vol',
          '-prefix', brain_path])

    brain = nib.load(brain_path)
    brain_data = brain.get_data()
    mask_data = np.zeros(brain.shape)
    mask_data[brain_data > 0] = 1
    mask_data = binary_fill_holes(mask_data)

    header = brain.get_header()
    header.set_data_dtype('uint8')
    nib.save(nib.Nifti1Image(mask_data.astype(np.uint8),
                             brain.get_affine(), header), brain_mask_path)


def bet_phase(phase_path, mag_path, brain_path, brain_mask_path):
    # Generate MASK from phase image.
    call(['3dSkullStrip',
          '-input', phase_path,
          '-orig_vol',
          '-mask_vol',
          '-use_edge',
          '-touchup',
          '-prefix', brain_mask_path])

    mask = nib.load(brain_mask_path)
    mask_data = mask.get_data()

    # AFNI maximum value is the brain mask, other values are representing the
    # CSF, the skin, etc.
    mask_data[mask_data < np.max(mask_data)] = 0

    # Re-save as the brain mask with only the max value.
    header = mask.get_header()
    header.set_data_dtype('uint8')
    nib.save(nib.Nifti1Image(mask_data.astype(np.uint8),
                             mask.get_affine(), header), brain_mask_path)

    mag = nib.load(mag_path)
    brain_data = mag.get_data() * mask_data
    nib.save(nib.Nifti1Image(brain_data, mag.get_affine(), mag.get_header()),
             brain_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("mode", type=str, choices=['mag', 'phase'],
                        help="The mode (mag|phase). "
                             " If mag => mag_path, --inverted_path "
                             " and --low_thresh are required."
                             " If phase => mag_path and --phase_path"
                             " are required.")

    parser.add_argument("mag_path", type=str,
                        help="The input file name to the magnitude image (SWI)"
                             " (relative path)")

    parser.add_argument("brain", type=str,
                        help="The output file name, which is skullstripped. "
                             "(relative path)")

    parser.add_argument("brain_mask", type=str,
                        help="The output file name, which is the mask of the"
                             "skullstripped (relative path).")

    # Mode == mag
    parser.add_argument("-i", "--inverted_path", type=str,
                        help="The output file name to the inverted magnitude"
                             " image (SWI with white vessels (relative path).")

    parser.add_argument("-l", "--low_thresh", type=int, default=80,
                        help="The lower threshold for the inversion of the "
                             "magnitude. [default: 80]")

    # Mode == phase
    parser.add_argument("-p", "--phase_path", type=str,
                        help="The input file name to the phase image (SWI) "
                             "(relative path)")

    exec_shell = ExtractSWIBrainMask(from_cmd=True)
    exec_shell.parse_arg(parser.parse_args())
    exec_shell.run()
