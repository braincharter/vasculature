#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extract the diameter from a binary image

Based on paper:

New algorithms for Euclidean distance transformation on an n-dimensional
digitized picture with applications," T. Saito and J. Toriwaki, Pattern
Recognition 27 (1994) 1551-1565.

Reverse Distance Transformation and Skeletons Based upon the Euclidean
Metric For n-Dimensional Binary Pictures, T. Saito and J. Toriwaki, IEICE
Trans. Inf. & Syst., Vol E77-D, No. 9, Sept. 1994

Exact medial axis with euclidean distance, E. Remy and E. Thiel, Image and
Vision Computing 23 (2005) 167-175.

A new method for the model-independent assessment of thickness in
three-dimensional images" T. Hildebrand and P. Rüesgsegger, J. of Microscopy,
185 (1996) 67-75.

http://www.optinav.com/Local_Thickness.htm

Usage:
  ExtractDiameter <input> <output> [options]

Options:
  -t --threshold <t>      Threshold value to create a binary mask
                          [default: 0]
  -d --distance <d>       Flag to save the distance map
  -r --ridge <r>          Flag to save the ridge distance map
  -u --thickness <u>      Flag to save the unclean thickness image
"""

import argparse
import math

import nibabel as nib
import numpy as np
from scipy import ndimage


class ExtractDiameter:
    def __init__(self, input_filename=None, output_filename=None,
                 threshold=0, distance_output=None, ridge_output=None,
                 thickness_output=None, from_cmd=False):

        self._input = input_filename
        self._output = output_filename
        self._threshold = threshold
        self._distance_output = distance_output
        self._ridge_output = ridge_output
        self._thickness_output = thickness_output

        if not from_cmd:
            self.valid_arg()

    def parse_arg(self, args):
        self._input = args.input
        self._output = args.output
        self._threshold = args.threshold
        self._distance = args.distance_output
        self._ridge = args.ridge_output
        self._thickness = args.thickness_output

    def valid_arg(self):

        if not self._input:
            raise Exception('You need to provide a file name for the input.')

        if not self._output:
            raise Exception('You need to provide a file name for the output.')

    def run(self):
        img = nib.load(self._input)

        pix_dim = np.array(img.get_header().get_zooms())
        if not np.allclose(pix_dim, pix_dim[::-1]):
            raise Exception('Your data must be isotropic to extract diameter. Now it is pix_dim=' + str(pix_dim) )

        data = img.get_data()
        data[data > self._threshold] = 1

        edt = compute_distance_map(data)

        if self._distance_output:
            nib.save(nib.Nifti1Image(edt, img.get_affine(),
                                     img.get_header()), self._distance_output)

        ridge_distance = compute_ridge_distance(edt)

        if self._ridge_output:
            nib.save(nib.Nifti1Image(ridge_distance, img.get_affine(),
                                     img.get_header()), self._ridge_output)

        thickness = compute_thickness(ridge_distance)

        if self._thickness_output:
            nib.save(nib.Nifti1Image(thickness * data, img.get_affine(),
                                     img.get_header()), self._thickness_output)

        diameter_clean = compute_clean_thickness(data, thickness, pix_dim[0])

        nib.save(nib.Nifti1Image(diameter_clean, img.get_affine(),
                                 img.get_header()), self._output)


def compute_distance_map(data):
    """
    To extract local thickness, the algorithm starts by creating the distance
    map from a binary image (segmentation). And here, we use the scipy EDT
    function (similar to the one used in the Fiji plugin).

    The code is based on this article:

    New algorithms for Euclidean distance transformation on an n-dimensional
    digitized picture with applications," T. Saito and J. Toriwaki, Pattern
    Recognition 27 (1994) 1551-1565.

    :param data
        A binary mask.
    :return:
        A volume with distance value in each voxel of the binary mask.

    """

    return ndimage.distance_transform_edt(data)


def compute_ridge_distance(edt):
    """
    The next step is to remove redundant points from the distance map to reduce
    the number of voxels to be tested for sphere fitting.

    The code is based on those articles:

    Reverse Distance Transformation and Skeletons Based upon the Euclidean
    Metric For n-Dimensional Binary Pictures, T. Saito and J. Toriwaki,
    IEICE Trans. Inf. & Syst., Vol E77-D, No. 9, Sept. 1994

    Exact medial axis with euclidean distance, E. Remy and E. Thiel,
    Image and Vision Computing 23 (2005) 167-175.

    :param edt
        A volume with distance value in each voxel of the binary mask.
    :return:
        A volume where redundant voxels from input has been removed.
    """

    # Get the rounded integer value of the square max value
    r_sq_max = int(round(np.max(edt) ** 2)) + 1

    # Get square distance for a sphere as an integer for each possible radius.
    dist_sq_values = np.round(np.unique(edt) ** 2).astype(np.int16)

    # Create an index array which associates potential distances with the max.
    dist_sq_index = np.zeros((r_sq_max,), dtype=np.int16)
    ind_ds = 0
    for x in xrange(r_sq_max):
        if x in dist_sq_values:
            dist_sq_index[x] = ind_ds
            ind_ds += 1

    # Create LUT based on existing distance values in the original mask.
    r_sq_lut = create_lookup_table(dist_sq_values)

    ridge_distance = np.zeros(edt.shape)
    connectivity = create_connectivity_list(1)

    for x, y, z in np.argwhere(edt):
        dist_value = edt[x, y, z]
        redundant_point = False

        # Get the associate index to search in LUT
        sk0_sq_ind = dist_sq_index[int(round(dist_value ** 2))]

        # Loop on each neighbor to check if the current voxel
        # is a redundant voxel in the distance map.
        for n in [tuple(e) for e in np.asarray((x, y, z)) - connectivity if
                  is_inside(edt.shape, tuple(e))]:

            # Get the distance of the neighbor in number of voxel (1, 2 or 3),
            # Also reduced by 1, to have LUT index between 0 and 2 (python
            # indexing).
            num_comp = np.sum(np.abs(np.asarray(n) - (x, y, z))) - 1

            # Check if current point (x,y,z) is in the neighbor sphere at
            # position n.
            if int(round(edt[n] ** 2)) >= r_sq_lut[num_comp][sk0_sq_ind]:
                redundant_point = True
                break

        if not redundant_point:
            ridge_distance[x, y, z] = dist_value

    return ridge_distance


def compute_thickness(ridge_distance):
    """
    The next step is to compute the local thickness by fitting a sphere in each
    points that remains in the ridge distance map.

    The code is based on this article:

    A new method for the model-independent assessment of thickness in
    three-dimensional images" T. Hildebrand and P. Rüesgsegger,
     J. of Microscopy, 185 (1996) 67-75.


    :param ridge_distance
        A volume where redundant voxels from input has been removed.
    :return:
        A volume with diameter value in each voxel (Bigger than binary mask).
    """

    w, h, d = ridge_distance.shape
    thickness = np.copy(ridge_distance)

    for i, j, k in np.argwhere(ridge_distance):
        r_squared = int(round(ridge_distance[i, j, k] ** 2))
        r_int = int(math.ceil(ridge_distance[i, j, k]))

        i_start = 0 if i - r_int < 0 else i - r_int
        i_stop = w if i + r_int >= w else i + r_int + 1

        j_start = 0 if j - r_int < 0 else j - r_int
        j_stop = h if j + r_int >= h else j + r_int + 1

        k_start = 0 if k - r_int < 0 else k - r_int
        k_stop = d if k + r_int >= d else k + r_int + 1

        for k1 in range(k_start, k_stop):
            r1_squaredk = (k1 - k) ** 2
            for j1 in range(j_start, j_stop):
                r1_squaredjk = r1_squaredk + (j1 - j) ** 2

                # Check if we are still within the original sphere
                if r1_squaredjk > r_squared:
                    continue

                for i1 in range(i_start, i_stop):
                    r1_squaredijk = r1_squaredjk + (i1 - i) ** 2

                    # Check if we are still within the original sphere
                    if r1_squaredijk <= r_squared:
                        # If the current sphere includes this neighbor and
                        # radius is bigger, modify neighbor value.
                        thickness[i1, j1, k1] = max(r_squared,
                                                    thickness[i1, j1, k1])

    # This function does not return the masked version of the thickness,
    # because it needs to be cleaned with clean thickness function. But when
    # it is written on the disk, the output is masked using the original
    # volume. (This is done to be identical to BoneJ plugin output).
    return 2.0 * np.sqrt(thickness)


def compute_clean_thickness(mask, thickness, image_resolution):
    """
    Next step is to clean the border of local thickness images.

    The code is based on this article:

    Clean thickness code from Local thickness OptiNav.
    http://www.optinav.com/Local_Thickness.htm

    :param mask
        A binary mask.
    :param thickness
        A volume with diameter value in each voxel (Bigger than binary mask).
    :param image_resolution
        Resolution of the original mask data.
    :return:
        A volume with approximate diameter value in mm, same size as the
        original mask.

    """
    connectivity = create_connectivity_list(1)

    # Create a flag array (0 = background, -1 = border, Thickness = inside)
    flag_array = np.copy(thickness)
    for x, y, z in np.argwhere(thickness):
        for n in [tuple(e) for e in np.asarray((x, y, z)) - connectivity if
                  is_inside(thickness.shape, tuple(e))]:

            # Set this point to border flag if one of neighbors is background.
            if thickness[n] == 0:
                flag_array[x, y, z] = -1
                break

    result_array = np.copy(thickness)
    for x, y, z in np.argwhere(flag_array == -1):
        positive_sum = 0
        nb_point = 0
        for n in [tuple(e) for e in np.asarray((x, y, z)) - connectivity if
                  is_inside(thickness.shape, tuple(e))]:
            if flag_array[n] > 0:
                positive_sum += thickness[n]
                nb_point += 1

        if nb_point > 0:
            result_array[x, y, z] = positive_sum / nb_point

    # Return the clean thickness masked with the original volume and also
    # multiply with the resolution of the input to obtain value as mm.
    return result_array * mask * image_resolution


def create_lookup_table(dist_sq_values):
    """
    For each offset from the origin (dx, dy, dz), and each radius-squared, r_sq
    find the smallest radius-squared, r1_sq, such that a ball of radius r1
    centered at (dx, dy, dz) includes a ball of r_sq centered at the origin.
    These balls refer to a 3D Integer grid. The set of (dx, dy, dz) points
    considered is a cube center at the origin.

    The code is based on this article:

    Exact medial axis with euclidean distance, E. Remy and E. Thiel,

    :param dist_sq_values
        A list of all possible distances values.
    :return:
        A look-up table to compare distance for redundant points.
    """
    # Create a Look-up table based on direction v (all directions are in 3D)
    # and the possible radius from the structure.
    side = scan_cube(1, 0, 0, dist_sq_values)
    diagonal = scan_cube(1, 1, 0, dist_sq_values)
    corner = scan_cube(1, 1, 1, dist_sq_values)

    return side, diagonal, corner


def scan_cube(dx, dy, dz, dist_sq_values):
    """
    For a list of r_sq, find the smallest r1_sq values such that a "ball" of
    radius r1 centered at (dx, dy, dz) includes the origin.

    The code is based on this article:

    Exact medial axis with euclidean distance, E. Remy and E. Thiel,

    A shape point p is the centre of a maximal disk if there is no other
    shape point q such that the ball BK1(q, DT[q]) entirely covers the
    ball BK1(p, DT[p]). (see figure 2 in paper to understand)

    :param dx
        Direction in x.
    :param dy
        Direction in y.
    :param dz
        Direction in z.
    :param dist_sq_values
        A list of all possible distances values
    :return
        A look-up table in a certain direction.
    """
    r1_radii = np.zeros((len(dist_sq_values),), dtype=np.int16)

    # Obtain the direction vector v.
    vx_abs = -int(math.fabs(dx))
    vy_abs = -int(math.fabs(dy))
    vz_abs = -int(math.fabs(dz))

    # For each possible square radius of the structure
    for t in xrange(len(dist_sq_values)):

        sq_r = dist_sq_values[t]
        r = 1 + int(math.sqrt(sq_r))
        max_sq_r = 0

        # Iterate on each value of the selected radius
        for k in xrange(r + 1):
            # Get offset in K
            sq_k = k ** 2
            # Radius k in direction v_z
            sq_rad_k = (k - vz_abs) ** 2

            for j in xrange(r + 1):
                # Get offset in K and J
                sq_kj = sq_k + j ** 2

                # Check if offset are still within the structure
                if sq_kj > sq_r:
                    continue

                # Radius I (offset between r and kj) in direction v_x
                rad_i = int(math.sqrt(sq_r - sq_kj)) - vx_abs

                # Radius K + Radius J in direction v_y + squared Radius I
                sq_r1 = sq_rad_k + (j - vy_abs) ** 2 + rad_i ** 2

                # Get maximal radius r1 into the template for current t
                if sq_r1 > max_sq_r:
                    max_sq_r = sq_r1

        r1_radii[t] = max_sq_r

    return r1_radii


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

    parser.add_argument("-t", "--threshold", type=int, default=0,
                        help="Threshold value to create a binary mask")

    parser.add_argument("-d", "--distance_output", type=str,
                        help="Flag to save the distance map")
    parser.add_argument("-r", "--ridge_output", type=str,
                        help="Flag to save the ridge distance map")
    parser.add_argument("-u", "--thickness_output", type=str,
                        help="Flag to save the unclean thickness image")

    exec_shell = ExtractDiameter(from_cmd=True)
    exec_shell.parse_arg(parser.parse_args())
    exec_shell.run()
