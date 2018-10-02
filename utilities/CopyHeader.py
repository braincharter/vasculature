#! /usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'biza2301'

import argparse

import nibabel as nib


class CopyHeader:
    def __init__(self, master_filename=None, input_filename=None,
                 output_filename=None, from_cmd=False, ):

        self._master = master_filename
        self._input = input_filename
        self._output = output_filename

        if not from_cmd:
            self.valid_arg()

    def parse_arg(self, args):
        self._master = args.master
        self._input = args.input
        self._output = args.output

    def valid_arg(self):

        if not self._input:
            raise Exception('You need to provide an input file name.')

        if not self._master:
            raise Exception('You need to provide an input master file name.')

        if not self._output:
            raise Exception('You need to provide an output file name.')

    def run(self):
        master = nib.load(self._master)
        input_data = nib.load(self._input)

        nib.save(nib.Nifti1Image(input_data.get_data(), master.get_affine(),
                                 master.get_header()), self._output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("master", type=str,
                        help="The master input file name where the header will"
                             "be taken (relative path)")

    parser.add_argument("input", type=str,
                        help="The input file name where the data array will"
                             "be taken(relative path)")

    parser.add_argument("output", type=str,
                        help="The output file name where the array is equal "
                             "to the input and header to the master. "
                             "(relative path)")

    exec_shell = CopyHeader(from_cmd=True)
    exec_shell.parse_arg(parser.parse_args())
    exec_shell.run()
