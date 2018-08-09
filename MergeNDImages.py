#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Merge X Nd images into a single image. Final = N+1_D = NdxX

Usage:
  MergeNDImages <output> <images>... [options]
"""

import argparse

import nibabel as nib
import numpy as np


class MergeNDimages:
    def __init__(self, images_filenames=None, output_filename=None,
                 from_cmd=False):

        self._images = images_filenames
        self._output = output_filename

        if not from_cmd:
            self.valid_arg()

    def parse_arg(self, args):
        self._images = args.images
        self._output = args.output

    def valid_arg(self):

        if not self._images:
            raise Exception('You need to provide a list[] of images as input.')

        if type(self._images) is not list:
            raise Exception('You need to provide a list[] of images as input.')

        if not self._output:
            raise Exception('You need to provide a file name for the output.')

    def run(self):

        if len(self._images) < 2:
            raise Exception("Provide at least 2 images to merge them.")

        first_img = nib.load(self._images[0])
        first_data = first_img.get_data()
        merge = first_data[..., np.newaxis]

        for image in self._images[1:]:
            data = nib.load(image).get_data()

            if first_data.shape != data.shape:
                raise Exception("Your images must have the same shape.")

            merge = np.concatenate((merge, data[..., np.newaxis]),
                                   axis=len(first_data.shape))

        # No needs to modify the header. Nibabel will append the new
        # dimension to this header, and the rest of the header still the same.
        nib.save(nib.Nifti1Image(merge, first_img.get_affine(),
                                 first_img.get_header()), self._output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("images", type=str, nargs='+',
                        help="A list of inputs file name to merge."
                             " (relative path)")

    parser.add_argument("output", type=str,
                        help="The output file name (relative path)")

    exec_shell = MergeNDimages(from_cmd=True)
    exec_shell.parse_arg(parser.parse_args())
    exec_shell.run()
