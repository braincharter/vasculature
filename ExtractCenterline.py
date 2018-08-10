#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extract the medial centerline from a binary image.

Usage:
  ExtractCenterline <input> <output> [options]
"""

import argparse

import nibabel as nib
import numpy as np
import skimage.morphology

class ExtractCenterline:
    def __init__(self, input_filename=None,
                 output_filename=None,
                 from_cmd=False, ):

        self._input = input_filename
        self._output = output_filename

        if not from_cmd:
            self.valid_arg()

    def parse_arg(self, args):
        self._input = args.input
        self._output = args.output

    def valid_arg(self):

        if not self._input:
            raise Exception('You need to provide a file name for the input.')

        if not self._output:
            raise Exception('You need to provide a file name for the output.')

    def run(self):

        img = nib.load(self._input)
        data = img.get_data().astype(np.bool)

        skel = skimage.morphology.skeletonize_3d(data)        
        
        nib.save(nib.Nifti1Image(skel, img.get_affine(),
                                 img.get_header()), self._output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str,
                        help="The input file name (relative path)")

    parser.add_argument("output", type=str,
                        help="The output file name (relative path)")

    exec_shell = ExtractCenterline(from_cmd=True)
    exec_shell.parse_arg(parser.parse_args())
    exec_shell.run()
