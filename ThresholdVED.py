#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Threshold the VED result (take-off X% of the value)

Usage:
  ThresholdVED <input> <output> --mode=<m> [--threshold=<t>]
  [--process_per_slice] [options]

  ThresholdVED <input> <output> otsu [--process_per_slice] [options]
  ThresholdVED <input> <output> yen [--process_per_slice] [options]
  ThresholdVED <input> <output> triangle [--process_per_slice] [options]
  ThresholdVED <input> <output> renyi [--process_per_slice] [options]
  ThresholdVED <input> <output> shanbhag [--process_per_slice] [options]

  ThresholdVED <input> <output> max [--threshold=<t>]
  [--process_per_slice] [options]
  ThresholdVED <input> <output> median [--threshold=<t>]
  [--process_per_slice] [options]

Options:
  --mode <m>
        Flag to select between different threshold methods [default: triangle]
  --process_per_slice
        Flag to process image slice per slice [default: False]
  --threshold <t>
        A thresholding percentage (%) which is used only if mode = max or
        median
"""

import argparse
import math

import nibabel as nib
import numpy as np
from scipy.ndimage.morphology import binary_fill_holes
from skimage import filters
from skimage.exposure import histogram


class ThresholdVED:
    def __init__(self, input_filename=None, output_filename=None, mode=None,
                 process_per_slice=False, threshold=0.20, from_cmd=False):

        self._input = input_filename
        self._output = output_filename
        self._mode = mode
        self._process_per_slice = process_per_slice
        self._threshold = threshold

        if not from_cmd:
            self.valid_arg()

    def parse_arg(self, args):
        self._input = args.input
        self._output = args.output
        self._mode = args.mode
        self._process_per_slice = args.process_per_slice
        self._threshold = args.threshold

    def valid_arg(self):

        mode_type = ['max', 'median', 'triangle', 'otsu', 'yen', 'renyi',
                     'shanbhag']

        if self._mode not in mode_type:
            raise Exception('You need to select a mode between those :'
                            'max, median, triangle, otsu, yen, renyi and '
                            'shanbhag.')

        if not self._input:
            raise Exception('You need to provide a file name for the input.')

        if not self._output:
            raise Exception('You need to provide a file name for the output.')

    def run(self):

        if self._mode == 'max' or self._mode == 'median':

            if not self._threshold:
                raise Exception("You need to provide a % threshold when using"
                                "mode == (max|median).")

        img = nib.load(self._input)
        scale_data = self.get_scaled_data(img.get_data())

        mask_data = np.zeros(scale_data.shape)
        mask_data[scale_data > 0] = 1
        mask_data = binary_fill_holes(mask_data)


        nib.save(nib.Nifti1Image(mask_data, img.get_affine(),
                                 img.get_header()), self._output)

    def get_scaled_data(self, data):
        """
        :param data:
            The image volume (3D) with vesselness measure in each voxel.
        :return:
            Return the thresholded image.
        """

        if self._mode == "renyi" or self._mode == "shanbhag":
            # Those mode need the data to be scaled between 0 and 255.
            data *= 255.0 / data.max()

        if self._process_per_slice:
            result = np.zeros(data.shape)
            for slices in range(data.shape[2]):
                result[:, :, slices] = self.apply_method(data[:, :, slices])

            return result

        return self.apply_method(data)

    def apply_method(self, data):
        """
        Apply different methods based on user selection.

        :param data:
            The image volume (2-3D) with vesselness measure in each voxel.
        :return:
            Return the threshold value.
        """
        val = 0

        # Threshold value based on histogram.
        if self._mode == 'otsu':  # Good with TOF
            val = filters.threshold_otsu(data)

        if self._mode == 'triangle':  # Good with TOF (best on TOF)
            val = threshold_triangle(data)

        if self._mode == "renyi":  # Good for SWI/TOF (best on swi)
            val = threshold_renyi_entropy(data)

        if self._mode == 'yen':  # Good with SWI
            val = filters.threshold_yen(data)

        if self._mode == 'shanbhag':  # Good for SWI
            val = threshold_shanbhag(data)

        # Threshold value data-driven based on a %, Good on full volume,
        # don't use those with per slice processing.
        if self._mode == 'max':
            # Good percentage to used with MAX mode.
            # For TOF : 0.15 to 0.25.
            # For SWI : 0.80 to 0.90.
            val = (np.max(data) * self._threshold)

        if self._mode == 'median':
            # Good percentage to used with MEDIAN mode.
            # For TOF : 0.30 to 0.40
            # For SWI : 0.10 to 0.15.
            val = (np.median(data) * self._threshold) + np.median(data)

        print "Threshold value : {}".format(val)

        return np.clip(data - val, a_min=0, a_max=1)


def threshold_triangle(image, nbins=256):
    """
    This function has been pasted from skimage 0.13_dev on :
    https://github.com/scikit-image/scikit-image/blob/master/skimage/filters/thresholding.py

    Return threshold value based on the triangle algorithm.
    Parameters
    ----------
    image : (N, M[, ..., P]) ndarray
        Grayscale input image.
    nbins : int, optional
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.
    Returns
    -------
    threshold : float
        Upper threshold value. All pixels with an intensity higher than
        this value are assumed to be foreground.
    References
    ----------
    .. [1] Zack, G. W., Rogers, W. E. and Latt, S. A., 1977,
       Automatic Measurement of Sister Chromatid Exchange Frequency,
       Journal of Histochemistry and Cytochemistry 25 (7), pp. 741-753
       DOI:10.1177/25.7.70454
    .. [2] ImageJ AutoThresholder code,
       http://fiji.sc/wiki/index.php/Auto_Threshold
    Examples
    --------
    from skimage.data import camera
    image = camera()
    thresh = threshold_triangle(image)
    binary = image > thresh
    """
    # nbins is ignored for integer arrays
    # so, we recalculate the effective nbins.
    hist, bin_centers = histogram(image.ravel(), nbins)
    nbins = len(hist)

    # Find peak, lowest and highest gray levels.
    arg_peak_height = np.argmax(hist)
    peak_height = hist[arg_peak_height]
    arg_low_level, arg_high_level = np.where(hist > 0)[0][[0, -1]]

    # Flip is True if left tail is shorter.
    flip = arg_peak_height - arg_low_level < arg_high_level - arg_peak_height
    if flip:
        hist = hist[::-1]
        arg_low_level = nbins - arg_high_level - 1
        arg_peak_height = nbins - arg_peak_height - 1

    # If flip == True, arg_high_level becomes incorrect
    # but we don't need it anymore.
    del arg_high_level

    # Set up the coordinate system.
    width = arg_peak_height - arg_low_level
    x1 = np.arange(width)
    y1 = hist[x1 + arg_low_level]

    # Normalize.
    norm = np.sqrt(peak_height ** 2 + width ** 2)
    peak_height /= norm
    width /= norm

    # Maximize the length.
    # The ImageJ implementation includes an additional constant when calculating
    # the length, but here we omit it as it does not affect the location of the
    # minimum.
    length = peak_height * x1 - width * y1
    arg_level = np.argmax(length) + arg_low_level

    if flip:
        arg_level = nbins - arg_level - 1

    return bin_centers[arg_level]


def threshold_renyi_entropy(image, nbins=256):
    """
    Code translated from:
    https://github.com/fiji/Auto_Threshold/blob/master/src/main/java/fiji/threshold/Auto_Threshold.java

    Source of this code:
    http://imagej.net/Auto_Threshold#RenyiEntropy

    THRESHOLD SELECTION USING RENYI'S ENTROPY,
    Sahoo P., Wilkins C. and Yeager J.
    http://www.sciencedirect.com/science/article/pii/S0031320396000659

    Supplementary information:
    https://en.wikipedia.org/wiki/R%C3%A9nyi_entropy

    :param image:
        The image array (in 3D or a single slice 2D).
    :param nbins:
        Number of bins used in the histogram.
    :return:
        Return the threshold value to use.
    """
    # Get P1, the probability to be a background voxel if gray-scale
    # value is into the bin X and P2 is the probability to be a
    # classified as the object for the similar bin (so P2 = 1 - P1).
    # Also this method set : self.hist, self.norm_hist and
    # self.entropy_range.
    hist, norm_hist, entropy_range, p1, p2 = init_threshold(image, nbins)

    # Define the t_star with all 3 Renyi equation, based on P1 and P2
    # distribution (P1 = Background, P2 = the object we want to preserve)
    t_star = [None] * 3
    t_star[0] = renyi_entropy(hist, norm_hist, entropy_range, p1, p2, 0.5)
    t_star[1] = shannon_entropy(hist, norm_hist, entropy_range, p1, p2)
    t_star[2] = renyi_entropy(hist, norm_hist, entropy_range, p1, p2, 2.0)
    t_star = sorted(t_star)

    # Obtain the beta weighting for the threshold value equation.
    beta = renyi_beta(t_star)

    # Determine the optimal threshold value.
    omega = p1[t_star[2]] - p1[t_star[0]]
    optimal_thresh = t_star[0] * (p1[t_star[0]] + 0.25 * omega * beta[0])
    optimal_thresh += 0.25 * t_star[1] * omega * beta[1]
    optimal_thresh += t_star[2] * (p2[t_star[2]] + 0.25 * omega * beta[2])

    return optimal_thresh


def init_threshold(image, nbins=256):
    """
    Initialize the Renyi Entropy, by creating normalized histogram and
    the bin range for the Entropy loops.

    This function fills the member variable : hist, norm_hist and
    entropy_range.

    :param image:
        The image array (in 3D or a single slice 2D).
    :param nbins:
        The number of bins to be divided in the histogram (max 256).
    :return:
        Return the probability to be backgroud => P1 and also the
        probability to be an object structure => P2.
    """
    hist, _ = histogram(image.ravel(), nbins)

    hist = hist.astype(np.float32)
    norm_hist = hist / np.sum(hist)

    # P1 is the probability to be background at a certain bin.
    # The probability of P1 needs to sum at 1.
    p1 = np.cumsum(norm_hist)

    # P2 is the probability for being an object for the first bin,
    # so this is the inverse probability of the background (P1)
    p2 = 1.0 - p1

    # nonzero returns a tuple of coord, so we access the first element,
    # and after we need the first bin, so again the first element.
    first_bin = np.nonzero(p1)[0][0]
    # Here, we need the last bin. So we access the last coordinate.
    last_bin = np.nonzero(p2)[0][-1]

    entropy_range = range(first_bin, last_bin + 1)

    return hist, norm_hist, entropy_range, p1, p2


def shannon_entropy(hist, norm_hist, entropy_range, p1, p2):
    """
    Since Alpha is always 1.0, This function is the Shannon entropy.
    https://en.wiktionary.org/wiki/Shannon_entropy

    Also, this is called instead of renyi_entropy with alpha = 1.0,
    because its computed only on non-zero histogram bins.

    H = -Sum_i_n(p_i * log(p_i))

    Where N is the number of class and P_i = Probability to be from the
    class I [Pr(X == I)].

    :param p1:
        Array with the probability to be in class 1 for all bins.
    :param p2:
        Array with the probability to be in class 2 for all bins.
    :return:
        Return the bin (t_star) which maximize this equation
    """

    thresh_val = 0
    max_entropy = 0
    nbins = len(hist)
    non_zero_el = np.nonzero(hist)

    for bin_it in entropy_range:
        back_id = np.intersect1d(non_zero_el, range(bin_it + 1))
        obj_id = np.intersect1d(non_zero_el, range(bin_it + 1, nbins))

        ent_back = np.sum(-norm_hist[back_id] *
                          np.log(norm_hist[back_id] / p1[bin_it]))
        ent_back /= p1[bin_it]

        if len(obj_id):
            ent_obj = np.sum(-norm_hist[obj_id] *
                             np.log(norm_hist[obj_id] / p2[bin_it]))
            ent_obj /= p2[bin_it]

        # Here, its the sum of P_i_n (P1 + P2).
        if max_entropy < ent_back + ent_obj:
            max_entropy = ent_back + ent_obj
            thresh_val = bin_it

    return thresh_val


def renyi_entropy(hist, norm_hist, entropy_range, p1, p2, alpha):
    """
    H_alpha(x) =     1
                 _________ * log(Sum_i_n((p_i/p_class)^alpha))

                 1 - alpha


    Where N is the number of class and P_i = Probability to be from the
    class I [Pr(X == I)].


    In case of the first t_star (0.5), this equals to do a sqrt (that
    what they done in the Fiji Java function), but the ^alpha still the
    best usage to keep the function general.

    :param p1:
        Array with the probability to be in class 1 for all bins.
    :param p2:
        Array with the probability to be in class 2 for all bins.
    :param alpha:
        The alpha parameter for the Renyi Entropy equation.
    :return:
        Return the bin (t_star) which maximize this equation
    """
    thresh_val = 0
    max_entropy = 0.0
    nbins = len(hist)

    for bin_it in entropy_range:
        # Define bins 0 to bin_it as background and check entropy of this
        # group. Use the inverse bin_it + 1 to nBins as object.
        ent_back = np.sum(norm_hist[0:(bin_it + 1)] ** alpha)
        ent_back /= p1[bin_it] ** alpha

        ent_obj = np.sum(norm_hist[(bin_it + 1):nbins] ** alpha)
        ent_obj /= p2[bin_it] ** alpha

        # Here, its the sum of the log, in the Renyi Equation define in
        # the function doc. log(x) + log(y) = log(x*y)
        if (ent_back * ent_obj) > 0.0:
            tot_ent = 1.0 / (1.0 - alpha) * math.log(ent_back * ent_obj)

        if tot_ent > max_entropy:
            max_entropy = tot_ent
            thresh_val = bin_it

    return thresh_val


def renyi_beta(t_star):
    """
    See this article for information about Beta:
    http://www.sciencedirect.com/science/article/pii/S0031320396000659

    This function defines the weighting to use on the optimal threshold
    definition function. Those weighting are define to emphasis on t_star
    (bin) who are nearby, this helps to converge to a threshold value. So
    gaussian (1,2,1) or a left/right bell (3,1,0)/(0,1,3).

    :param t_star:
        A list with the sorted t_star.
    :return:
        A list with the beta value.
    """
    if abs(t_star[0] - t_star[1] <= 5):
        if abs(t_star[1] - t_star[2] <= 5):
            return [1, 2, 1]
        else:
            return [0, 1, 3]
    else:
        if abs(t_star[1] - t_star[2] <= 5):
            return [3, 1, 0]
        else:
            return [1, 2, 1]


def threshold_shanbhag(image, nbins=256):
    """
    This method is based on :

    https://github.com/fiji/Auto_Threshold/blob/master/src/main/java/fiji/threshold/Auto_Threshold.java
    http://imagej.net/Auto_Threshold#RenyiEntropy

    :param image:
        The image array (in 3D or a single slice 2D).
    :param nbins:
        Number of bins used in the histogram.
    :return:
        Return the threshold value to use.
    """
    hist, norm_hist, entropy_range, p1, p2 = init_threshold(image, nbins)

    thresh_val = -1.0
    min_entropy = np.max(image)

    for bin_it in entropy_range:
        term = 0.5 / p1[bin_it]
        ent_back = term * np.sum(-norm_hist[1:(bin_it + 1)] *
                                 np.log(1.0 - term * p1[0:bin_it]))

        term2 = 0.5 / p2[bin_it]
        ent_obj = term2 * np.sum(-norm_hist[(bin_it + 1):nbins] *
                                 np.log(1.0 - term2 * p2[(bin_it + 1):nbins]))

        if abs(ent_back - ent_obj) < min_entropy:
            min_entropy = abs(ent_back - ent_obj)
            thresh_val = bin_it

    return thresh_val


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str,
                        help="The input file name (relative path)")
    parser.add_argument("output", type=str,
                        help="The output file name (relative path)")

    parser.add_argument("mode", type=str, default='triangle',
                        choices=['max', 'median', 'triangle', 'otsu', 'yen',
                                 'renyi', 'shanbhag'],
                        help="The thresholding method [default: triangle]")

    parser.add_argument("-s", "--process_per_slice", action="store_true",
                        help="Flag to process image slice per slice "
                             "[default: False]")

    parser.add_argument("-t", "--threshold", type=float, default=0.20,
                        help="A thresholding percentage which is used "
                             "only if mode = max or median [default: 0.20]")

exec_shell = ThresholdVED(from_cmd=True)
exec_shell.parse_arg(parser.parse_args())
exec_shell.run()
