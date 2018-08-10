#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Compute the vessel enhancing diffusion (VED)

Based on paper :

Frangi, AF, Niessen, WJ, Vincken, KL, & Viergever, MA (1998). Multiscale
Vessel Enhancement Filtering. In Wells, WM, Colchester, A, & Delp, S,
Editors, MICCAI '98 Medical Image Computing and Computer-Assisted
Intervention, Lecture Notes in Computer Science, pages 130-137, Springer
Verlag, 1998.

Manniesing, R, Viergever, MA, & Niessen, WJ (2006). Vessel Enhancing
Diffusion: A Scale Space Representation of Vessel Structures. Medical Image
Analysis, 10(6), 815-825.
"""

import argparse
import os, sys
from subprocess import call


class ComputeVED:
    def __init__(self, input_filename=None, output_filename=None,
                 sigma_min=0.3, sigma_max=6.0, number_scale=10,
                 dark_blood=False, alpha=0.5, beta=1.0, c=0.00001,
                 number_iterations=1, sensitivity=5.0, wstrength=25.0,
                 epsilon=0.1,
                 generate_iteration_files=False,
                 generate_scale=False,
                 generate_hessian=False,
                 out_folder=None,
                 frangi_only=False,
                 scale_object=False,
                 from_cmd=False):

        self._input = input_filename
        self._output = output_filename

        self._sigma_min = sigma_min
        self._sigma_max = sigma_max
        self._number_scale = number_scale

        self._dark_blood = dark_blood
        self._alpha = alpha
        self._beta = beta
        self._c = c

        self._number_iterations = number_iterations
        self._sensitivity = sensitivity
        self._wstrength = wstrength
        self._epsilon = epsilon

        self._generate_iteration_files = generate_iteration_files
        self._generate_scale = generate_scale
        self._generate_hessian = generate_hessian
        self._out_folder = out_folder

        self._frangi_only = frangi_only
        self._scale_object = scale_object

        if not from_cmd:
            self.valid_arg()

    def parse_arg(self, args):
        self._input = args.input
        self._output = args.output

        self._sigma_min = args.sigma_min
        self._sigma_max = args.sigma_max
        self._number_scale = args.number_scale

        self._dark_blood = args.dark_blood
        self._alpha = args.alpha
        self._beta = args.beta
        self._c = args.c

        self._number_iterations = args.number_iterations
        self._sensitivity = args.sensitivity
        self._wstrength = args.wstrength
        self._epsilon = args.epsilon

        self._generate_iteration_files = args.generate_iteration_files
        self._generate_scale = args.generate_scale
        self._generate_hessian = args.generate_hessian
        self._out_folder = args.out_folder

        self._frangi_only = args.frangi_only
        self._scale_object = args.scale_object

    def valid_arg(self):

        if not self._input:
            raise Exception('You need to provide a file name for the input')

        if not self._output:
            raise Exception('You need to provide a file name for the output')

    def run(self):

        if self._generate_iteration_files or self._generate_scale \
                or self._generate_hessian:

            if not self._out_folder:
                raise Exception(
                    "Since you wants to generate extra files, you need "
                    "to provide an output folder with: out_folder.")

        files_created_ved = []
        shorts = []
        
        prefix = str.split(self._output, ".")[0] + "_"

        if self._generate_iteration_files or self._frangi_only:
            vesselness_path = os.path.join(self._out_folder, 
                                           (prefix + 'frangi_only_vesselness'
                                            '_measure.nii.gz'))
            files_created_ved.append(vesselness_path)
            shorts.append('frangi_only_vesselness_measure.nii.gz')

        if self._generate_scale:
            best_scale_path = os.path.join(self._out_folder, 
                                           (prefix + 'ved_generated_'
                                            'best_scale.nii.gz'))
            files_created_ved.append(best_scale_path)
            shorts.append('ved_generated_best_scale.nii.gz')

        if self._generate_hessian:
            best_hessian_path = os.path.join(self._out_folder, 
                                             (prefix + 'ved_generated_'
                                              'best_Hessian.nii.gz'))
            files_created_ved.append(best_hessian_path)
            shorts.append('ved_generated_best_Hessian.nii.gz')

        print(self._generate_iteration_files)

        # range to save last iteration.
        if self._generate_iteration_files:
            iteration_path = [os.path.join(self._out_folder, 
                                           (prefix + 'ved_iteration'
                                            '_{}.nii.gz'.format(i)))
                              for i in range(1, self._number_iterations + 1)]
            files_created_ved.extend(iteration_path)

            shorts_path = [os.path.join('ved_iteration_{}.nii.gz'.format(i))
                              for i in range(1, self._number_iterations + 1)]
            shorts.extend(shorts_path)

        if files_created_ved:
            if not os.path.exists(self._out_folder):
                os.mkdir(self._out_folder)

        kwargs = {}

        print("shorts : ", shorts)
        print("outputs : ", files_created_ved)

        # All parameters need to be pass as string or None if it is a flag.
        kwargs['--input'] = self._input
        kwargs['--output'] = self._output

        # Multi-scale parameters.
        kwargs['--sigmaMin'] = str(self._sigma_min)
        kwargs['--sigmaMax'] = str(self._sigma_max)
        kwargs['--numberOfScale'] = str(self._number_scale)

        # Frangi parameters.
        if self._dark_blood:
            kwargs['--darkBlood'] = None

        kwargs['--alpha'] = str(self._alpha)
        kwargs['--beta'] = str(self._beta)
        kwargs['--c'] = str(self._c)

        # VED parameters.
        kwargs['--numberOfIteration'] = str(self._number_iterations)
        kwargs['--sensitivity'] = str(self._sensitivity)
        kwargs['--wStrength'] = str(self._wstrength)
        kwargs['--epsilon'] = str(self._epsilon)

        # Flags
        if self._frangi_only:
            kwargs['--frangiOnly'] = None
        if self._scale_object:
            kwargs['--scaleObject'] = None
        if self._generate_scale:
            kwargs['--generateScale'] = None
        if self._generate_hessian:
            kwargs['--generateHessian'] = None
        if self._generate_iteration_files:
            kwargs['--generateIterationFiles'] = None

        cmd_string = [sys.path[0] + '/itkVEDMain']
        
        for k in kwargs:
            cmd_string.append(k)
            if kwargs[k]:
                cmd_string.append(kwargs[k])

        print(cmd_string)
        ret_code = call(cmd_string)

        print("VED binary execution is finished with code: " + str(ret_code))

        if files_created_ved:
            print("TEST")
            print("shorts are: " + str(shorts))
            print("files_created_ved are: " + str(files_created_ved))
            for short_path, long_path in zip(shorts, files_created_ved):
                os.rename(short_path, long_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str,
                        help="The input file name (relative path)")
    parser.add_argument("output", type=str,
                        help="The output file name (relative path)")

    # Multi-scale parameters.
    parser.add_argument("-m", "--sigma_min", type=float, default=0.3,
                        help="The minimum sigma used in the "
                             "multi-scale analysis. Usually the "
                             "smallest resolution of the "
                             "acquisition. [default: 0.3]")

    parser.add_argument("-M", "--sigma_max", type=float, default=6.0,
                        help="The maximum sigma used in the "
                             "multi-scale analysis. Usually the "
                             "biggest expected size of a vessel. "
                             "[default: 6.0]")

    parser.add_argument("-n", "--number_scale", type=int, default=10,
                        help="The number of scales created between"
                             " the sigma min/max in the multi-scale "
                             "analysis. [default: 10]")

    # Frangi parameters.
    parser.add_argument("-d", "--dark_blood", action="store_true",
                        help="Flag to extract black blood vessel.")

    parser.add_argument("-a", "--alpha", type=float, default=0.5,
                        help="The alpha parameter used in Frangi "
                             "vesselness equation to limit blob-like "
                             "structure. [default: 0.5]")

    parser.add_argument("-b", "--beta", type=float, default=1.0,
                        help="The beta parameter used in Frangi "
                             "vesselness equation to limit "
                             "plate-like structure. [default: 1]")

    parser.add_argument("-c", "--c", type=float, default=0.00001,
                        help="The c parameter used in Frangi "
                             "vesselness equation to limit the "
                             "background comparison. "
                             "[default: 0.00001]")

    # VED parameters.
    parser.add_argument("-t", "--number_iterations", type=int,
                        default=1,
                        help="The number of diffusion iterations. "
                             "[default: 1]")

    parser.add_argument("-s", "--sensitivity", type=float,
                        default=5.0,
                        help="The sensitivity used in VED param. "
                             "[default: 5.0]")

    parser.add_argument("-w", "--wstrength", type=float, default=25.0,
                        help="The weight strength used in VED "
                             "equation. [default: 25.0]")

    parser.add_argument("-e", "--epsilon", type=float, default=0.1,
                        help="The epsilon used in VED equation. "
                             "[default: 0.1]")

    # Flags
    parser.add_argument("-f", "--frangi_only", action="store_true",
                        help="Flag to stop the pipeline after Frangi "
                             "equation (No post-processing with "
                             "VED).")

    parser.add_argument("-O", "--scale_object", action="store_true",
                        help="Flag to rescale vesselness measure "
                             "based on eigen amplitude.")

    parser.add_argument("-S", "--generate_scale", action="store_true",
                        help="Flag to generate an output for the "
                             "best sigma value per voxel.")

    parser.add_argument("-H", "--generate_hessian",
                        action="store_true",
                        help="Flag to generate an output with the "
                             "hessian matrix per voxel.")

    parser.add_argument("-I", "--generate_iteration_files",
                        action="store_true",
                        help="Flag to generate output iteration "
                             "files and vesselness.")

    parser.add_argument("-D", "--out_folder", type=str,
                        help="he output folder for all the optional "
                             "generated files. This is required if "
                             "any of those flags are set : "
                             "[generate_scale, generate_hessian, "
                             "generate_iteration_files, frangi_only]")

    exec_shell = ComputeVED(from_cmd=True)
    exec_shell.parse_arg(parser.parse_args())
    exec_shell.run()
