#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compute the distance between 2 ROIs, which is restricted by the binary mask

Usage:
  DistanceBetweenROI <input> <output> <start_roi> <end_roi> [options]

Options:
    -m --max_distance <m>   Maximum search distance in mm. [default: 15]
    -d --distance_type <d>  Type of distance (addition | euclidean)
                            [default: euclidean]
"""
import argparse

import nibabel as nib
import numpy as np
from scipy import ndimage
from scipy.spatial import distance


class DistanceBetweenROI:
    def __init__(self, input_filename=None, output_filename=None,
                 start_roi=None, end_roi=None, max_distance=15,
                 distance_type='euclidean', from_cmd=False):
        self._input = input_filename
        self._output = output_filename
        self._start_roi = start_roi
        self._end_roi = end_roi
        self._max_distance = max_distance
        self._distance_type = distance_type

        if not from_cmd:
            self.valid_arg()

    def parse_arg(self, args):
        self._input = args.input_filename
        self._output = args.output_filename
        self._start_roi = args.start_roi
        self._end_roi = args.end_roi
        self._max_distance = args.max_distance
        self._distance_type = args.distance_type

    def valid_arg(self):

        if not self._input:
            raise Exception('You need to provide a file name for the input.')

        if not self._output:
            raise Exception('You need to provide a file name for the output.')

        if not self._start_roi:
            raise Exception('You need to provide a file name for the input'
                            'start ROI.')

        if not self._end_roi:
            raise Exception('You need to provide a file name for the input'
                            'end ROI.')

    def run(self):
        img = nib.load(self._input)

        pix_dim = np.array(img.get_header().get_zooms())

        if not np.allclose(pix_dim, pix_dim[::-1]):
            raise Exception('Your data must be isotropic to extract distance.')

        data = img.get_data()
        start_data = nib.load(self._start_roi).get_data()
        end_data = nib.load(self._end_roi).get_data()

        start_seed = self.get_seed(data, start_data, pix_dim[0])
        end_seed = self.get_seed(data, end_data, pix_dim[0])

        distance_between_roi = self.flood_fill_distance(data, end_data,
                                                        start_seed,
                                                        end_seed)
        distance_between_roi *= pix_dim[0]

        print 'Maximum inside vessel length between ROIs: {} ' \
              '(approximate mm).'\
              .format(np.max(end_data * distance_between_roi))

        nib.save(nib.Nifti1Image(distance_between_roi, img.get_affine(),
                                 img.get_header()), self._output)

    def get_seed(self, mask, roi_data, resolution):
        """
        :param mask:
            A binary mask
        :param roi_data:
            A binary ROI mask indicating where the distance process starts or
            ends.
        :param resolution:
            Resolution of the original mask data.
        :return:
            Return the coordinate of the seed for the current Roi.
        """

        intersect = np.clip(mask + roi_data - 1, 0, 1).astype(np.int16)
        inter_label = ndimage.label(intersect, np.ones((3, 3, 3)))[0]

        # Label 0 = background, 1 to N = intersection. We limit to only 1 label
        if len(np.unique(inter_label)) > 2:
            raise Exception('There are multiples intersection regions between '
                            'binary mask and ROIs. Please redefine your ROIs.')

        if len(np.unique(inter_label)) == 1:
            seed = self.search_nearest_point_from_roi(mask, roi_data,
                                                      resolution)

            print 'Your ROI does not intersect the binary mask, but nearest mask ' \
                  'point from ROI is {}.'.format(seed)
        else:
            seed = ndimage.measurements.center_of_mass(intersect)

        seed = [int(round(x)) for x in seed]

        if not mask[tuple(seed)]:
            raise Exception('The value of one of the seed is not included in '
                            'the binary mask.')

        return seed

    def search_nearest_point_from_roi(self, mask, roi, resolution):
        """
        If there is no intersection between the binary mask and a ROI, this
        will obtain the coordinate inside the binary mask which is the
        nearest from the ROI.

        :param mask:
            A binary mask
        :param roi:
            A ROI volume where coordinates will be paired with the
            binary mask to get the nearest coordinate from this volume.
        :param resolution:
            Resolution of the original mask data.
        :return:
            Return the nearest coordinate in binary mask from the ROI.
        """

        # A potential optimization could be to reduce the mask volume into a
        # bounding box of size 'max_distance' around the ROI. Which limits the
        # number of coordinates used for mask_vec. But cdist is fast and this
        # optimization won't be profitable.

        mask_vec = np.argwhere(mask)
        roi_vec = np.argwhere(roi)

        dist_table = distance.cdist(mask_vec, roi_vec, 'euclidean')

        if np.min(dist_table) * resolution > self._max_distance:
            raise Exception('Your ROI is too far from the binary mask. Please '
                            'redefines it or increase the maximum distance '
                            'parameter (max_distance: {:.2f}, current distance'
                            ' from ROI: {:.2f}.'
                            .format(self._max_distance,
                                    np.min(dist_table) * resolution))

        # Unravel_index[0] = mask_vec and [1] = roi_vec
        return mask_vec[np.unravel_index(np.argmin(dist_table),
                                         dist_table.shape)[0]]

    def flood_fill_distance(self, mask, roi_end, start_seed, end_seed):
        """
        Use a flood fill method to run through the binary mask and increment
        the distance value from the start seed.

        :param mask:
            A binary mask
        :param roi_end
            A binary ROI mask indicating where the processing ends.
        :param start_seed:
            Coordinate in N-D of the start seed (must be a list).
        :param end_seed:
            Coordinate in N-D of the end seed (must be a list).
        :return:
            The distance mask with value in voxel.
        """

        result = np.zeros(mask.shape)
        connectivity = create_connectivity_list(1)

        # Add zero as distance value for the initial seed. Tuple = (x, y, z, D)
        start_seed.append(0)
        stack = set()
        stack.add(tuple(start_seed))

        pass_inside_end = False
        while stack:
            col, row, dim, dist = stack.pop()

            for e in [tuple(e) for e in np.asarray((col, row, dim)) -
                      connectivity if is_inside(mask.shape, tuple(e))]:

                # Check if the current point is in the vessel mask and also if
                # it is not processed in the distance map.
                if mask[e] and result[e] < 1:

                    if self._distance_type == 'addition':
                        result[e] = dist + 1
                        # Append distance to tuple of coordinates.
                        e = e + (dist + 1,)
                    else:
                        orientation = np.asarray(e) - np.asarray(
                            (col, row, dim))
                        d_increase = dist + np.sqrt(np.sum(orientation ** 2))
                        result[e] = d_increase
                        # Append distance to tuple of coordinates.
                        e = e + (d_increase,)

                    # Check if the end seed has been reached
                    if e[0:3] == tuple(end_seed):
                        return result

                    if roi_end[e[0:3]]:
                        pass_inside_end = True

                    stack.add(e)

        if not pass_inside_end:
            raise Exception(
                'Flood fill did not reach the end ROI. The ROI might '
                'be disconnected from the start ROI.')
        else:
            print 'Flood fill have reached a point in the end ROI, but not its' \
                  ' center of mass.'

        return result


def is_inside(shape, coord):
    for shape_, coord_ in zip(shape, coord):
        if coord_ < 0 or coord_ >= shape_:
            return False

    return True


def create_connectivity_list(size):
    return [(i, j, k)
            for i in range(-size, size + 1)
            for j in range(-size, size + 1)
            for k in range(-size, size + 1)
            if (not (i == j == k == 0))]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str,
                        help="The input file name (relative path)")
    parser.add_argument("output", type=str,
                        help="The output file name (relative path)")

    parser.add_argument("start_roi", type=str,
                        help="The start ROI file name (relative path)")
    parser.add_argument("end_roi", type=str,
                        help="The end ROI file name (relative path)")

    parser.add_argument("-m", "--max_distance", type=float, default=15.0,
                        help="Maximum search distance in mm. [default: 15]")

    parser.add_argument("-d", "--distance_type", type=str,
                        choices=['addition', 'euclidean'],
                        default='euclidean',
                        help="Type of distance (addition | euclidean) "
                             "[default: euclidean]")

    exec_shell = DistanceBetweenROI(from_cmd=True)
    exec_shell.parse_arg(parser.parse_args())
    exec_shell.run()
