#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Invert the color of the image base on the maximum in the mask and only change
color in this mask.

Usage:
  InvertContrast <input> <mask> <output> [options]
"""

import argparse

import nibabel as nib
import numpy as np


class InvertContrast:
    def __init__(self, input_filename=None,
                 mask_filename=None,
                 output_filename=None,
                 from_cmd=False, ):

        self._input = input_filename
        self._mask = mask_filename
        self._output = output_filename

        if not from_cmd:
            self.valid_arg()

    def parse_arg(self, args):
        self._input = args.input
        self._mask = args.mask
        self._output = args.output

    def valid_arg(self):

        if not self._input:
            raise Exception('You need to provide a file name for the input.')

        if not self._mask:
            raise Exception('You need to provide a file name for the mask.')

        if not self._output:
            raise Exception('You need to provide a file name for the output.')

    def run(self):

        img = nib.load(self._input)
        mask_data = nib.load(self._mask).get_data().astype(np.bool)

        data = img.get_data()
        new_arr = np.ma.array(data, mask=~mask_data)

        result = np.max(new_arr) - data
        inverse = result * mask_data

        nib.save(nib.Nifti1Image(inverse, img.get_affine(),
                                 img.get_header()), self._output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str,
                        help="The input file name (relative path)")

    parser.add_argument("mask", type=str,
                        help="The input file name for the mask of the brain "
                             "(relative path)")

    parser.add_argument("output", type=str,
                        help="The output file name (relative path)")

    exec_shell = InvertContrast(from_cmd=True)
    exec_shell.parse_arg(parser.parse_args())
    exec_shell.run()
