#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extract the magnitude and velocity images from a Dicom file. Also apply the
slope and intercept to the data to be in float.

Usage:
  ConvertFourDFlowImage <input> <mag> <velocity> [options]
"""

from collections import namedtuple
import argparse

import nibabel as nib
import numpy as np
from pydicom.filereader import read_file

image_info = namedtuple('ImageInfo', ['intercept', 'slope', 'start_index'])


class ConvertFourDFlowImage:
    def __init__(self, input_filename=None, mag=None, velocity=None,
                 from_cmd=False):
        self._input = input_filename
        self._mag = mag
        self._velocity = velocity

        if not from_cmd:
            self.valid_arg()

    def parse_arg(self, args):
        self._input = args.input
        self._mag = args.mag
        self._velocity = args.velocity

    def valid_arg(self):
        if not self._input:
            raise Exception('You need to provide a file name for the input.')

        if not self._mag:
            raise Exception('You need to provide an output file name for the '
                            'magnitude image.')

        if not self._velocity:
            raise Exception('You need to provide an output file name for the '
                            'velocity image.')

    def run(self):
        data = read_file(self._input)

        if data.Manufacturer != 'Philips Medical Systems':
            raise Exception("This code can only read Philips Dicom, because "
                            "tags may vary between manufacturer.")

        list_seq = data.PerFrameFunctionalGroupsSequence
        pixel_spacing = [float(i) for i in
                         list_seq[0].PixelMeasuresSequence[0].PixelSpacing]

        # The 4th spacing is the delta time from TriggerTime. To obtain
        # the delta between two times: delta = seq[1].TriggerTime -
        # seq[0].TriggerTime, but since seq[0] is always 0. We can take
        # TriggerTime from seq[1], it is already equals to the delta.
        image_spacing = (pixel_spacing[0],
                         pixel_spacing[1],
                         float(data.SpacingBetweenSlices),
                         float(list_seq[1][0x2005, 0x0140f][0].TriggerTime))

        dimension = data[0x2001, 0x1018].value
        phase = data[0x2001, 0x1017].value
        seq_length = dimension * phase

        _mag_name = ['ORIGINAL', 'PRIMARY', 'M_FFE', 'M', 'FFE']
        _vel_name = ['ORIGINAL', 'PRIMARY', 'VELOCITY MAP', 'P', 'PCA']

        image_info_mag = obtain_sequence_info(list_seq, seq_length,
                                              _mag_name)
        image_info_vel = obtain_sequence_info(list_seq, seq_length,
                                              _vel_name)

        # Data in dicom is stored like (dimension * phase, row, column)
        dicom_shape = (dimension, phase, data.Rows, data.Columns)

        create_image(data, dicom_shape, seq_length, image_spacing,
                     image_info_mag, self._mag)
        create_image(data, dicom_shape, seq_length, image_spacing,
                     image_info_vel, self._velocity)


def obtain_sequence_info(list_seq, seq_length, seq_name):
    """
    Return a namedtuple with the slope, intercept and starting indexing of
    the image.

    :param list_seq:
        The list of sequences of the Dicom
    :param seq_length:
        The length of a sequence
    :param seq_name:
        The name of the sequence in the dicom.
    :return:
        A namedtuple which contains the image information.
    """

    for x in range(0, len(list_seq), seq_length):

        if list_seq[x][0x2005, 0x0140f][0].ImageType == seq_name:
            pixel = list_seq[x].PixelValueTransformationSequence[0]

            return image_info(pixel.RescaleIntercept, pixel.RescaleSlope, x)


def create_image(data, dicom_shape, seq_length, spacing, img_info, output):
    """
    Create a nifti image with his associates header and affine matrix for
    magnitude (0) and for velocity (1)

    :param data:
        The dicom object which contains the data array and the dicom header.
    :param dicom_shape:
        The shape of the data to reconstruct (dimension, phase, row, column)
    :param seq_length:
        The length of a single sequence (normally dimension * phase).
    :param spacing:
        Data spacing for dimension: x, y, z, t.
    :param img_info:
        Namedtuple with the image information.
    :param output:
        The output name for the image
    """
    top_bound = img_info.start_index + seq_length

    arr_data = np.squeeze(
        data.pixel_array[img_info.start_index:top_bound, :, :])

    arr_data = np.reshape(arr_data, dicom_shape)
    arr_data = arr_data * img_info.slope + img_info.intercept

    # Re-orient the data to be in the right orientation (same as Philips
    # MatLab code)
    arr_data = np.rot90(np.flipud(np.transpose(arr_data, (2, 3, 0, 1))))

    affine = create_affine_matrix(data, spacing)

    # Need to pass arr_data.shape, because it now ordered as:
    # (row, column, dimension, phase)
    hdr = create_header(img_info.slope, img_info.intercept, arr_data.shape,
                        spacing)

    nib.save(nib.Nifti1Image(arr_data, affine, hdr), output)


def create_affine_matrix(data, spacing):
    """
    This needs to be called to fill Nifti affine matrix while saving.

    To create the matrix, the code is based on :

    http://stackoverflow.com/questions/21759013/dicom-affine-matrix-transformation-from-image-space-to-patient-space-in-matlab

    :param data:
        The dicom header with the information to create the affine matrix.
    :param spacing:
        Data spacing for dimension: x, y, z, t.
    :return:
        An affine matrix for the nifti.
    """

    list_seq = data.PerFrameFunctionalGroupsSequence

    cos_xy = np.array([float(i) for i in
                       list_seq[0].PlaneOrientationSequence[0]
                      .ImageOrientationPatient])

    cos_x = cos_xy[0:3]
    cos_y = cos_xy[3:6]

    patient_pos_in_first_seq = np.array([float(i) for i in list_seq[0]
                                        .PlanePositionSequence[0]
                                        .ImagePositionPatient])

    dimension = data[0x2001, 0x1018].value  # Can't access private tag

    # Deducing an orthogonal direction for Z.
    if dimension == 1:
        cos_z = np.cross(cos_x, cos_y)
    else:
        patient_pos_in_last_seq = np.array([float(i) for i in
                                            list_seq[len(list_seq) - 1]
                                           .PlanePositionSequence[0]
                                           .ImagePositionPatient])

        cos_z = (patient_pos_in_first_seq -
                 patient_pos_in_last_seq) / (dimension * spacing[2])

    # Concatenate the R matrix with the translation T.
    top = np.transpose(np.vstack([cos_x * spacing[0],
                                  cos_y * spacing[1],
                                  cos_z * spacing[2],
                                  patient_pos_in_first_seq]))

    affine = np.vstack([top, np.array([[0.0, 0.0, 0.0, 1.0]])])

    # Overwrite affine diagonal to be pure resolution and not approximate.
    affine[0, 0] = -spacing[0]
    affine[1, 1] = spacing[1]
    affine[2, 2] = spacing[2]

    return affine


def create_header(slope, inter, array_shape, spacing):
    """
    :param slope:
        The rescale slope to get array in float
    :param inter:
        The rescale intercept to get array in float
    :param array_shape:
        The shape of the original array (row, columns, dimension, phase)
    :param spacing:
        Data spacing for dimension: x, y, z, t.
    :return:
        Header of the nifti file.
    """
    hdr = nib.Nifti1Header()

    hdr['scl_slope'] = slope
    hdr['scl_inter'] = inter

    hdr.set_data_shape(array_shape)

    # Divide by 1000.0 to have it in second.
    hdr.set_zooms((spacing[0], spacing[1], spacing[2], spacing[3] / 1000.0))

    return hdr


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str,
                        help="The input file name (relative path)")
    parser.add_argument("mag", type=str,
                        help="The output file name for the magnitude image"
                             "(relative path)")

    parser.add_argument("velocity", type=str,
                        help="The output file name  for the velocity image"
                             "(relative path)")

    exec_shell = ConvertFourDFlowImage(from_cmd=True)
    exec_shell.parse_arg(parser.parse_args())
    exec_shell.run()
