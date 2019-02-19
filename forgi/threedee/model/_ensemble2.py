"""
This code is still experimental
"""

from __future__ import absolute_import, unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip)


from collections import MutableSequence, Sequence
import sys
import math
import numpy as np
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.descriptors as ftmd
import scipy.stats
import matplotlib.pyplot as plt
import warnings
import pymbar

all = ["autocorrelate_data", "Ensemble", "is_stationary_adf"]
################## Functions for performing multiple analysis steps on data #######################


def overall_analysis(y, measure_name="???"):
    import pymbar
    #skip = int(len(y)/20000)
    # if skip > 1:
    #    y_sk = y[::skip]
    #print("Is the distribution stationary?")
    #adf_stationary =  is_stationary_adf(y_sk)
    # if adf_stationary:
    #    print("Augmented Dickey-Fuller Test passed")
    # else:
    #    print("Augmented Dickey Fuller Test failed")
    #burn_in, production, integr_corr_time = split_and_bin_timeseries(y)
    # need_split=False
    # if integr_corr_time>0.8*len(y):
    #    print("Whole process is considered burn-in.")
    #    need_split=True
    # elif not is_stationary_adf(production):
    #    need_split=True
    #    print("Augmented Dickey Fuller: Production is not stationary! ")
    # if not need_split:
    #    plot(y, measure_name)
    #    sys.exit()
    #splitlength = int(len(y)/500)
    #if splitlength<2000: splitlength = 2000
    #print("Using splitlength of {}".format(splitlength))

    # Change point detection. See Section 3 of
    # https://github.com/amanahuja/change-detection-tutorial

    fig, ax = plt.subplots(3)
    ax[0].plot(np.arange(len(y)), y)
    # Calculate whole-data statistic inefficiency
    """parts = [(y, 0)]
    newparts=[]
    final = {}
    while True:
        for part in parts:
            y1, y2, s, o = split_timeseries_into_two(part[0], part[1])
            if not o:
                final[s[2]] = y2
            else:
                if len(y1)>100:
                  newparts.append((y1, s[0]))
                if len(y2)>100:
                  newparts.append((y2, s[2]))
        if not newparts:
            break
        parts=newparts
        newparts=[]
    fig, ax = plt.subplots(2)
    ax[0].plot(np.arange(len(y)), y)
    for start in sorted(final.keys()):
        y=final[start]
        ax[1].plot(np.arange(start, start+len(y)), y)
    plt.show()
    sys.exit()"""
    inefficiencies = []
    ineff = pymbar.timeseries.statisticalInefficiency(y)
    inefficiencies.append(ineff)
    window_lengths = [len(y)]

    def two_halfs(x, start, end):
        mid = math.ceil((start + end) / 2)
        return pymbar.timeseries.statisticalInefficiency(x[start:mid]), pymbar.timeseries.statisticalInefficiency(x[mid:end]), mid
    indices = [0, len(y)]
    new_indices = set(indices)
    for split in range(1, 30):
        curr_ineff = []
        plt_x = []
        plt_y = []
        for i in range(len(indices) - 1):
            try:
                ineff1, ineff2, mid = two_halfs(y, indices[i], indices[i + 1])
                plt_x += [indices[i], mid, mid, indices[i + 1]]
                plt_y += [ineff1, ineff1, ineff2, ineff2]
            except pymbar.utils.ParameterError:
                break
            curr_ineff += [ineff1, ineff2]
            new_indices.add(mid)

        ax[2].plot(plt_x, plt_y, label="Split {}".format(split))
        if curr_ineff:
            inefficiencies.append(np.mean(np.array(curr_ineff)))
            window_lengths.append(mid - indices[i])
        indices = sorted(new_indices)
    ax[1].plot(window_lengths, inefficiencies, 'o-')
    plt.show()
    ie = np.array(inefficiencies)
    delta_ie = ie[1:] - ie[:-1]
    dd_ie = delta_ie[1:] - delta_ie[:-1]
    for i, dd in enumerate(dd_ie):
        if i == 0:
            continue
        if dd < 0:  # dd_ie[1]/50:
            break
    fig, ax = plt.subplots(3)
    ax[0].plot(np.arange(len(ie)), ie, "o-")
    ax[0].plot([i, i + 1, i + 2], [ie[i], ie[i + 1], ie[i + 2]], "ro")
    ax[1].plot(np.arange(len(delta_ie)), delta_ie, "o-")
    ax[1].plot([0, (len(delta_ie))], [0, 0])
    ax[1].plot([i, i + 1], [delta_ie[i], delta_ie[i + 1]], "ro")
    ax[2].plot(np.arange(len(dd_ie)), dd_ie, "o-")
    ax[2].plot([0, (len(dd_ie))], [0, 0])
    ax[2].plot([i], [dd_ie[i]], "ro")
    window_length = len(y) // (2**(i + 1))
    print("Window length", window_length)
    ineffs = []
    plt_x = []
    plt_y = []
    for start in range(0, len(y) - window_length // 2, window_length):
        try:
            ie = pymbar.timeseries.statisticalInefficiency(
                y[start:start + window_length])
        except pymbar.utils.ParameterError:
            ie = window_length
        ineffs.append(ie)
        plt_x += [start, min(len(y), start + window_length)]
        plt_y += [ie, ie]
    fig, ax = plt.subplots(2)
    ax[0].plot(np.arange(len(y)), y)
    ax[1].plot(plt_x, plt_y, "o-")
    plt.show()
    sys.exit()

    fig, ax = plt.subplots(2)

    def two_halfs(x, start, end):
        mid = math.ceil((start + end) / 2)
        return pymbar.timeseries.statisticalInefficiency(x[start:mid]), pymbar.timeseries.statisticalInefficiency(x[mid:end]), mid
    indices = [0, len(y)]
    new_indices = set(indices)
    for split in range(1, 7):
        plt_x = []
        plt_y = []
        for i in range(len(indices) - 1):
            neff1, neff2, mid = two_halfs(y, indices[i], indices[i + 1])
            plt_x += [indices[i], mid, mid, indices[i + 1]]
            plt_y += [neff1, neff1, neff2, neff2]
            new_indices.add(mid)
        ax[1].plot(plt_x, plt_y, label="Split {}".format(split))
        indices = sorted(new_indices)
    plt_x = [0, len(y) // 3, len(y) // 3, 2 * len(y) //
             3, 2 * len(y) // 3, len(y)]
    plt_y = [pymbar.timeseries.statisticalInefficiency(y[0: len(y) // 3]), pymbar.timeseries.statisticalInefficiency(y[0: len(y) // 3]),
             pymbar.timeseries.statisticalInefficiency(
                 y[len(y) // 3: 2 * len(y) // 3]), pymbar.timeseries.statisticalInefficiency(y[len(y) // 3: 2 * len(y) // 3]),
             pymbar.timeseries.statisticalInefficiency(y[2 * len(y) // 3:]), pymbar.timeseries.statisticalInefficiency(y[2 * len(y) // 3:])]
    ax[1].plot(plt_x, plt_y, "o-")
    ax[1].legend()
    plt_x = [0, len(y) // 3, len(y) // 3, 2 * len(y) //
             3, 2 * len(y) // 3, len(y)]
    plt_y = [pymbar.timeseries.statisticalInefficiency(y[0: len(y) // 3]), pymbar.timeseries.statisticalInefficiency(y[0: len(y) // 3]),
             pymbar.timeseries.statisticalInefficiency(
                 y[len(y) // 3: 2 * len(y) // 3]), pymbar.timeseries.statisticalInefficiency(y[len(y) // 3: 2 * len(y) // 3]),
             pymbar.timeseries.statisticalInefficiency(y[2 * len(y) // 3:]), pymbar.timeseries.statisticalInefficiency(y[2 * len(y) // 3:])]
    ax[1].plot(plt_x, plt_y, "o-")
    ax[1].legend()
    for start in range(0, len(y), 700):
        for end in range(start, len(y),  700):
            neff = [(end - start) /
                    pymbar.timeseries.statisticalInefficiency(y[start:end])]
            ax[1].plot([start, end], [neff, neff])
    plt.show()
    sys.exit()
    # for y in (x[::-1], np.concatenate((x[30000:60000],x[60000:],x[:30000])), x, np.concatenate((x[:60000],x[:30000][::-1])), np.concatenate((x[60000:],x[30000:60000],x[60000:])), np.concatenate((x[30000:60000],x[60000:],x[30000:60000]))):
    if True:
        for step in [1000, 2000]:
            neffs = []
            neff_back = []
            plt_x = []
            print("Step", step)
            for start in range(0, len(y) - step, step):
                neffs += [(start) /
                          pymbar.timeseries.statisticalInefficiency(y[:start])]
                neff_back += [(len(y) - start) /
                              pymbar.timeseries.statisticalInefficiency(y[start:])]
                plt_x += [start]
                if (start % 500) == 0:
                    print(start, "*", end="")
                    sys.stdout.flush()  # , flush=True)

            ax[1].plot(plt_x, neffs, "o-", label="Forward {}".format(step))
            ax[1].plot(plt_x, neff_back, label="Backward {}".format(step))
            ax[2].plot(plt_x, np.array(neffs) +
                       np.array(neff_back), label=step)
        ax[1].legend()
    plt.show()
    sys.exit()

    # Use this as window size
    num_blocks = int(len(y) / ineff)
    plt_x = []
    plt_y = []
    for start in range(0, len(y) - ineff, ineff):
        plt_x += [start, start + ineff]
        plt_y += [pymbar.timeseries.statisticalInefficiency(
            y[start: start + ineff])] * 2
    ax.plot(plt_x, plt_y)
    plt.show()
    sys.exit()

    x = y

    for y in (x[::-1], np.concatenate((x[30000:60000], x[60000:], x[:30000])), x):
        indices, data = optimal_neff_step_detection(y)
        fig, ax = plt.subplots(2)
        old_i = 0
        print(indices)
        for i in indices:
            ax[0].plot(np.arange(old_i, i), y[old_i:i],
                       label="{}-{}".format(old_i, i))
            old_i = i
        if indices:
            ax[0].plot(np.arange(i, len(y)), y[i:], label="Rest")
        # ax[0].legend()
        old_len = 0
        for z in data:
            ax[1].plot(np.arange(old_len, old_len + len(z)), z)
            old_len += len(z)

    plt.show()
    for y in (x[::-1], np.concatenate((x[30000:60000], x[60000:], x[:30000])), x):
        y = y[::25]

        # Choose window size long enough to have a constant mean
        oldmean = y[:50].mean()
        for windowsize in range(100, int(len(y)) / 10, 50):
            newmean = y[:windowsize].mean()
            print(newmean, oldmean, abs(newmean - oldmean), newmean * 0.01)
            if abs(newmean - oldmean) < newmean * 0.01:
                break
            oldmean = newmean

        # Use sliding windows and optimize the number of effective samples from the beginning.

        print("windowsize", windowsize)
        z = []
        acorr = []
        neffs = []
        for start in range(2 * windowsize, len(y) - windowsize, windowsize // 2):
            #print("stat-ineff window", pymbar.timeseries.statisticalInefficiency(y[start:start+windowsize]))
            neffs.append((start + windowsize) /
                         pymbar.timeseries.statisticalInefficiency(y[:start + windowsize]))
            #mean_sofar = y[:start+windowsize].mean()
            #mean_window = y[start:start+windowsize].mean()
            #std_sofar = y[:start+windowsize].std()
            #z_score =  (mean_window - mean_sofar)/(std_sofar/math.sqrt(windowsize))
            # z.append(z_score)
            if (start % 5000) == (2 * windowsize % 5000):
                # print(start, z_score, mean_sofar, mean_window, std_sofar, math.sqrt(windowsize))
                print("*", end="")
            # if abs(z_score)>0.05:
            #   break
        print("start is", start)
        fig, ax = plt.subplots(2, 2)
        ax[0, 0].plot(np.arange(len(y)), y)
        #ax[0,0].plot(np.arange(start,len(y)), y[start:])
        #ax[0,1].plot(np.arange(len(z)), z)
        #ax[0,1].set_title("Z score")
        #ax[1,0].plot(np.arange(len(sie)), sie)
        # ax[1,0].set_title("inefficiency")
        ax[1, 0].plot(np.arange(len(neffs)) *
                      windowsize + 2 * windowsize, neffs)
        ax[1, 0].set_title("Neff so far")
        acorr1 = autocorrelate_data(y[:1200])
        acorr2 = autocorrelate_data(y[:1200 + windowsize])
        ax[0, 1].plot(np.arange(len(acorr1)), acorr1, label="Until Step")
        ax[0, 1].set_title("Acorr 1")
        ax[0, 1].plot(np.arange(len(acorr2)), acorr2, label="After step")
        ax[0, 1].plot([0, len(acorr2)], [0, 0])
        ax[0, 1].legend()
        mean1 = np.mean(y[:1200])
        unbiased1 = y[:1200] - mean1
        ynorm1 = np.sum(unbiased1**2)
        ax[1, 1].plot(np.arange(1200), unbiased1 / ynorm1, label="Until Step")
        mean2 = np.mean(y[:1200 + windowsize])
        unbiased2 = y[:1200 + windowsize] - mean2
        ynorm2 = np.sum(unbiased2**2)
        ax[1, 1].set_title("Normalized for calculation of autocorrelation")
        ax[1, 1].plot(np.arange(1200 + windowsize),
                      unbiased2 / ynorm2, label="After Step")
        ax[1, 1].plot([0, 1200], [(unbiased2 / ynorm2)[:1200].mean(), (unbiased2 /
                                                                       ynorm2)[:1200].mean()], label="Mean for first part, after step")
        ax[1, 1].plot([0, 1200 + windowsize], [0, 0])
        ax[1, 1].legend()
    plt.show()
    sys.exit()

    fig, ax = plt.subplots(2, 3)
    for splitlength in [100, 200, 500, 1000, 2000, 5000]:
        means = np.zeros(int(len(y) / 100))
        stds = np.zeros(int(len(y) / 100))
        for i, start in enumerate(range(0, len(y), splitlength)):
            if len(y) - start > 200:  # Enough samples
                means[i] = y[start:start + splitlength].mean()
                stds[i] = y[start:start + splitlength].std()
        means = means[:i]
        stds = stds[:i]

        ax[0, 0].plot(np.arange(len(means)), means)
        ax[1, 0].plot(np.arange(len(stds)), stds)
    mom = np.array([means[:e].mean() for e in range(1, len(means))])
    mos = np.array([stds[:e].mean() for e in range(1, len(stds))])
    ax[0, 1].plot(np.arange(len(mom)), mom)
    ax[1, 1].plot(np.arange(len(mos)), mos)
    dmom = mom[1:] - mom[:-1]
    dmos = mos[1:] - mos[:-1]
    ax[0, 2].plot(np.arange(len(dmom)), dmom)
    ax[1, 2].plot(np.arange(len(dmos)), dmos)
    plt.show()
    for i, start in enumerate(range(0, len(y), splitlength)):
        if abs(y[:start].mean() - y[start:start + splitlength].mean()) > 2 * y[:start].std():
            break

    fig, ax = plt.subplots()
    ax.plot(np.arange(start), y[:start])
    ax.plot(np.arange(start, len(y)), y[start:])
    plt.show()


def plot(y, measure_name="", ax=None):
    """
    Plot several graphs for this timeseries.
    """
    corr = autocorrelate_data(y)
    burn_in, production, integr_corr_time = split_and_bin_timeseries(y)

    # The MEASURE
    if ax is None:
        fig, ax = plt.subplots(3, 2)

    ax[0, 0].set_xlabel("step")
    ax[0, 0].set_ylabel(measure_name)
    ax[0, 0].set_title(measure_name)
    if integr_corr_time > 0:
        ax[0, 0].plot(np.arange(len(y[:integr_corr_time])),
                      y[:integr_corr_time], label="burn-in", color="red")
        ax[0, 0].plot([integr_corr_time],
                      y[integr_corr_time], "ro", markersize=5)
    if integr_corr_time < len(y):
        ax[0, 0].plot(np.arange(len(y[:integr_corr_time]), len(y)),
                      y[integr_corr_time:], label="production", color="blue")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ax[0, 0].legend()

    # THE AUTOCORRELATION #Uses code from https://github.com/choderalab/pymbar/blob/master/pymbar/timeseries.py under GPL license
    corr_chodera = np.zeros(len(corr))
    A_n = y[integr_corr_time:]  # TODO: add skip
    mu_A = A_n.mean()
    dA_n = A_n.astype(np.float64) - mu_A
    sigma2_AA = (dA_n * dA_n).mean()
    t = 1
    increment = 1
    N = A_n.size
    while (t < N - 1):
        C = np.sum(dA_n[0:(N - t)] * dA_n[t:N] + dA_n[0:(N - t)]
                   * dA_n[t:N]) / (2.0 * float(N - t) * sigma2_AA)
        corr_chodera[t] = C
        if (C <= 0.0) and (t > 3):
            break
        t += increment
        increment += 1

    ax[1, 0].plot(np.arange(len(corr)), corr,
                  label="whole dataset", color="violet")
    ax[1, 0].plot(np.arange(len(corr_chodera)), corr_chodera,
                  color="lightblue", label="pymbar (production)")

    ax[1, 0].set_xlabel("lag")
    ax[1, 0].set_ylabel("autocorrelation")
    ax[1, 0].plot([0, len(corr)], [0, 0])
    ax[1, 0].legend()

    # BLOCK AVERAGE
    if len(burn_in) > 0:
        ax[2, 0].plot(np.arange(len(burn_in)), burn_in,
                      "o-", label="burn-in", color="red")
    if len(production) > 0:
        ax[2, 0].plot(np.arange(len(production)) + len(burn_in),
                      production, "o-", label="production", color="blue")
    ax[2, 0].set_xlabel("effective step")
    ax[2, 0].set_ylabel("Block average of " + measure_name)
    ax[2, 0].set_title("Block average for {} independent samples".format(
        len(production) + len(burn_in)))
    ax[2, 0].legend()

    # THE DEVIATION OF THE MEASURE
    dy = y[1:] - y[:-1]
    ax[0, 1].plot(np.arange(len(dy)), dy)
    ax[0, 1].set_title("Change of {}".format(measure_name))
    ax[0, 1].set_xlabel("step")
    ax[0, 1].set_ylabel("Delta {}".format(measure_name))

    # Autocorr of deviation #Uses code from https://github.com/choderalab/pymbar/blob/master/pymbar/timeseries.py under GPL license
    corrD = autocorrelate_data(dy)
    corr_chodera = np.zeros(len(corr))
    burn_in, production, integr_corr_time = split_and_bin_timeseries(dy)
    A_n = dy[integr_corr_time:]  # TODO: add skip
    mu_A = A_n.mean()
    dA_n = A_n.astype(np.float64) - mu_A
    sigma2_AA = (dA_n * dA_n).mean()
    t = 1
    increment = 1
    N = A_n.size
    while (t < N - 1):
        C = np.sum(dA_n[0:(N - t)] * dA_n[t:N] + dA_n[0:(N - t)]
                   * dA_n[t:N]) / (2.0 * float(N - t) * sigma2_AA)
        corr_chodera[t] = C
        if (C <= 0.0) and (t > 3):
            break
        t += increment
        increment += 1

    ax[1, 1].plot(np.arange(len(corrD)), corrD)
    ax[1, 1].plot(np.arange(len(corr_chodera)),
                  corr_chodera, label="pymbar (production)")
    ax[1, 1].set_xlabel("lag")
    ax[1, 1].set_ylabel("autocorrelation")
    ax[1, 1].plot([0, len(corrD)], [0, 0])
    ax[1, 1].legend()
    # BLOCK AVERAGE
    if len(burn_in) > 0:
        ax[2, 1].plot(np.arange(len(burn_in)), burn_in, "o-", label="burn-in")
    if len(production) > 0:
        ax[2, 1].plot(np.arange(len(production)) + len(burn_in),
                      production, "o-", label="production")
    ax[2, 1].set_xlabel("effective step")
    ax[2, 1].set_ylabel("Block average of Delta " + measure_name)
    ax[2, 1].set_title("Block average for {} independent samples".format(
        len(production) + len(burn_in)))
    ax[2, 1].legend()

    plt.show()


def split_and_bin_timeseries(y, skip=None):
    """
    Split the time series into an equilibration and equilibrium part (if applicable) and
    bin the equilibrium part in statistically independent samples.

    :param y: A time series as numpy 1D array
    :param skip: INT>0 or None.
                 If None, decide wether or not to skip datapoints in y for faster processing
                 depending on the length of y
                 If INT i: Only use every ith entry in y (use 1 for disabeling skip)
    """
    import pymbar
    if pymbar.version.version.split(".")[0] != '3':
        warnings.warn("pymbar version {} found. This code is tested with version 3 of pymbar, "
                      "available at github.")
    if skip is None:
        skip = int(len(y) / 20000)
    if skip > 1:
        x = y[::skip]
    else:
        x = y
        skip = 1

    t, g, neff = pymbar.timeseries.detectEquilibration(x, fast=True, nskip=1)

    production = blocked_average_reverse(
        y[t * skip:], math.ceil(g * skip))[::-1]
    burn_in = blocked_average_reverse(y[:t * skip], math.ceil(g * skip))[::-1]
    return burn_in, production, t * skip


def split_timeseries_into_two(y, overall_start):
    ie = []
    ineff = pymbar.timeseries.statisticalInefficiency(y)
    ie.append(ineff)

    def two_halfs(x, start, end):
        mid = math.ceil((start + end) / 2)
        return pymbar.timeseries.statisticalInefficiency(x[start:mid]), pymbar.timeseries.statisticalInefficiency(x[mid:end]), mid
    indices = [0, len(y)]
    new_indices = set(indices)
    for split in range(1, 30):
        curr_ineff = []
        for i in range(len(indices) - 1):
            try:
                ineff1, ineff2, mid = two_halfs(y, indices[i], indices[i + 1])
            except pymbar.utils.ParameterError:
                break
            curr_ineff += [ineff1, ineff2]
            new_indices.add(mid)
        if curr_ineff:
            ie.append(np.mean(np.array(curr_ineff)))
        indices = sorted(new_indices)
    ie = np.array(ie)
    delta_ie = ie[1:] - ie[:-1]
    dd_ie = delta_ie[1:] - delta_ie[:-1]
    for i, dd in enumerate(dd_ie):
        if i == 0:
            continue
        if dd < abs(dd_ie[0]) / 500:  # Almost zero
            break
    fig, ax = plt.subplots(4)
    ax[0].plot(np.arange(len(ie)), ie, "o-")
    ax[0].plot([i, i + 1, i + 2], [ie[i], ie[i + 1], ie[i + 2]], "ro")
    ax[1].plot(np.arange(len(delta_ie)), delta_ie, "o-")
    ax[1].plot([0, (len(delta_ie))], [0, 0])
    ax[1].plot([i, i + 1], [delta_ie[i], delta_ie[i + 1]], "ro")
    ax[2].plot(np.arange(len(dd_ie)), dd_ie, "o-")
    ax[2].plot([0, (len(dd_ie))], [0, 0])
    ax[2].plot([i], [dd_ie[i]], "ro")
    ax[3].plot(np.arange(len(y)), y)
    plt.show()
    window_length = len(y) // (2**(i + 1))
    print("Window length", window_length)
    ineffs = []
    plt_x = []
    plt_y = []
    for start in range(0, len(y) - window_length // 3, window_length // 2):
        try:
            ie = pymbar.timeseries.statisticalInefficiency(
                y[start:start + window_length])
        except pymbar.utils.ParameterError:
            ie = window_length
        ineffs.append(ie)
        plt_x += [start, min(len(y), start + window_length)]
        plt_y += [ie, ie]
    max_i = np.argmax(np.array(ineffs))
    splitpoint = window_length // 2 * max_i + window_length // 2

    fig, ax = plt.subplots(2)
    ax[0].plot(np.arange(len(y)), y)
    ax[1].plot(plt_x, plt_y, "o-")
    ax[0].plot([splitpoint], [y[splitpoint]], "ro")
    ax[1].plot([splitpoint], [ineffs[max_i]], "ro")
    plt.show()
    if splitpoint < 100:
        y1 = []
        start1 = 0
        end1 = 0
    else:
        _, _, start1 = split_and_bin_timeseries(y[:splitpoint], 25)
        print("First")
        if start1 < splitpoint / 4:
            y1 = y[start1:splitpoint]
        else:
            y1 = y[:splitpoint]
            start1 = 0
        end1 = splitpoint
    if len(y) - splitpoint < 100:
        y2 = []
        start2 = len(y)
    else:
        _, _, start2 = split_and_bin_timeseries(y[splitpoint:], 25)
        print("Second")
        if start2 < (len(y) - splitpoint) / 4:
            y2 = y[splitpoint + start2:]
            start2 = splitpoint + start2
        else:
            y2 = y[splitpoint:]
            start2 = splitpoint
    n_eff1 = len(y1) / pymbar.timeseries.statisticalInefficiency(y1)
    n_eff2 = len(y2) / pymbar.timeseries.statisticalInefficiency(y2)
    n_eff_total = len(y) / pymbar.timeseries.statisticalInefficiency(y)
    if n_eff_total >= n_eff1 + n_eff2:
        print("No split", n_eff_total, ">=", n_eff1 + n_eff2)
        return [], y, [overall_start, overall_start, overall_start, overall_start + len(y)], False
    else:
        print("Split", n_eff_total, "<", n_eff1 + n_eff2)
        return y1, y2, [overall_start + start1, overall_start + end1, overall_start + start2, overall_start + len(y)], True


########################## Functions working with timeseries as 1D arrays #########################
"""def optimal_neff_step_detection(x):
    '''
    We use the idea from J.Chodera, JTCT 2016, to optimize the effective number of independent
    samples for step detection.

    The assumption is that the integrated autocorrelation of a timeseries increases,
    if a step appears in a previousely stationary timeseries.

    Note that the autocorrelation function uses the mean and the variance of the sample.
    If several subsequent values are far from the mean, the autocorrelation increases.
    If a step is included in the data, this slightly shifts the mean, so everything before the
    step starts to contribute to the autocorrelation. Further more, the new window has a high autocorrelation.
    '''
    import pymbar
    # Choose window size long enough to have a constant mean
    oldmean=x[:50].mean()
    for windowsize in range(100,len(x)//10, 50):
        newmean = x[:windowsize].mean()
        print(oldmean, newmean, abs(newmean-oldmean), newmean*0.01)
        if abs(newmean-oldmean)<newmean*0.01:
            break
        oldmean=newmean
    else:
        return [], [x]
    #Find the first splitpoint from the left using sliding windows
    oldNeff = 0
    neffs = []
    for start in range(2*windowsize, len(x)-windowsize, windowsize//4):
        ineff = pymbar.timeseries.statisticalInefficiency(x[:start+windowsize])
        newNeff = (start+windowsize)/ineff
        neffs.append(newNeff)
        if newNeff<oldNeff:
            break;
        oldNeff = newNeff
    else:
        ineff = pymbar.timeseries.statisticalInefficiency(x)
        return [], [blocked_average_reverse(x, math.ceil(ineff))[::-1]] #No step-points found.
    #Optional refinement
    fig, ax = plt.subplots(2)
    ax[0].plot(np.arange(start), x[:start])
    ax[0].plot(np.arange(start,start+windowsize), x[start:start+windowsize])
    ax[0].plot(np.arange(start+windowsize, len(x)), x[start+windowsize:])
    ax[1].plot(np.arange(len(neffs)), neffs)
    plt.show()
    #Split into 2 And process the second part
    ineff = pymbar.timeseries.statisticalInefficiency(x[:start])
    y = blocked_average_reverse(x[:start], math.ceil(ineff))[::-1]
    rest_indices, rest_data = optimal_neff_step_detection(x[start:])
    return [start]+[ start+p for p in rest_indices], [y]+rest_data
"""


def autocorrelate_data(y, mean=None):
    """
    Calculate the autocorrelation for the time series y for different lags.

    :param y: A 1D np.array
    :param mean: None (use mean of y) or a FLOAT (Use this value as "mean" for normalization)
    :returns: A 1D np.array
    """
    if mean is None:
        mean = np.mean(y)
    yunbiased = y - mean
    ynorm = np.sum(yunbiased**2)
    corr = np.correlate(yunbiased, yunbiased, mode="same") / ynorm
    return corr[len(y) / 2:]


def is_stationary_adf(y):
    """
    Whether a time series y is stationary or not. Uses the Augmented Dickey-Fuller Test.

    :param y: A 1D np.array
    :returns: BOOLEAN
    """
    import statsmodels.tsa.stattools as smtools
    adf = smtools.adfuller(y)
    return adf[0] < adf[4]['5%']


def blocked_average_reverse(x, bin_width):
    """
    The blocked average of x in blocks of length tau, starting with the last block and
    ignoring the first datapoints of x that do not fully fit into blocks.

    See Yang et al, Journal of Chemical Physics 2004, doi:10.1063/1.1638996

    :param x: A numpy 1D array representing a timeseries
    :param bin_width: INT. How many samples should be averaged per block

    :returns: Numpy 1D array.
    """
    length, = x.shape  # x must be 1 dimensional
    num_blocks = int(length / bin_width)

    try:
        return np.mean(x[::-1][:num_blocks * bin_width].reshape((num_blocks, -1)), axis=1)
    except ValueError:
        if num_blocks == 0:
            return np.array([])
        else:
            raise

############################### Ensembles of CG-Objects ###########################################


def ensemble_from_filenames(filenames):
    cgs = []
    for fn in filenames:
        cgs.append(ftmc.CoarseGrainRNA(fn))
    return Ensemble(cgs)


class EnsembleBase(Sequence):
    """
    Baseclass for Ensemble and Ensemble View
    """

    AVAILABLE_DESCRIPTORS = ["rog", "anisotropy", "asphericity"]

    def _get_descriptor(self, descriptor, domain=None):
        """
        :param ensemble: an Ensemble or EnsembleView object
        :param descriptor: A STRING. One of AVAILABLE_DESCRIPTORS
        :param domain: An iterable of cg element names or None (whole cg)
        :returns: A np.array
        """
        if descriptor not in self.AVAILABLE_DESCRIPTORS:
            raise ValueError("Descriptor {} not available.".format(descriptor))
        if descriptor == "rog":
            if domain:
                return np.array([ftmd.radius_of_gyration(cg.get_poss_for_domain(domain, "vres")) for cg in self])
            else:
                return np.array([cg.radius_of_gyration() for cg in self])
        elif descriptor == "anisotropy":
            if domain:
                return np.array([ftmd.anisotropy(cg.get_poss_for_domain(domain, "vres")) for cg in self])
            else:
                return np.array([ftmd.anisotropy(cg.get_ordered_stem_poss()) for cg in self])
        elif descriptor == "asphericity":
            if domain:
                return np.array([ftmd.asphericity(cg.get_poss_for_domain(domain, "vres")) for cg in self])
            else:
                return np.array([ftmd.asphericity(cg.get_ordered_stem_poss()) for cg in self])

    def autocorrelation(self, descriptor="rog", domain=None, mean=None):
        """
        Return the normalized autocorrelation as a 1D array for the given measure along the
        sequence of structures in this ensemble.

        Uses Code from http://stackoverflow.com/a/17090200/5069869

        :param ensemble: A ensemble or EnsembleView object
        :param descriptor: A STRING. One of AVAILABLE_DESCRIPTORS
        :param domain: An iterable of cg element names or None (whole cg)
        :param mean: A FLOAT or None.
                      If this is None, all datapoints will be shifted by their mean
                      before calculating the autocorrelation.
                      If this is numeric, datapoints will be shifted by this value instead.

        """
        y = self.get_descriptor(descriptor, domain)
        return autocorrelate_data(y, mean)


class Ensemble(EnsembleBase):
    def __init__(self, cgs):
        """
        An Ensemble is a sequence of Coarse grained RNAs, all of which must correspond
                    to the same RNA 2D structure.

        :param cgs: An ordered iterable of coarse grain RNAs.
        """
        self._cgs = cgs
        # Cached Data
        self._descriptors = {}

    # Methods for accessing the stored cg structures
    def __getitem__(self, i):
        if isinstance(i, int):
            return self._cgs[i]
        if isinstance(i, slice):
            if i.step is not None:
                raise IndexError("Steps are not supperted when using slice notation "
                                 "on Ensemble objects")
            try:
                return EnsembleView(self, i.start, i.stop)
            except ValueError as e:
                raise IndexError(e)
        raise IndexError("Unsupported index {}" / format(i))

    def __len__(self):
        return len(self._cgs)

    def update(self):
        """
        Clear all cached values.

        Has to be called, whenever a CoarseGrainRNA that is
        part of the ensemble is modified in place.
        """
        self._descriptors = {}

    def get_descriptor(self, descriptor, domain=None):
        """
        Get a numpy 1D array of one descriptor, e.g. the Radius of Gyration,
        for all RNAs in this ensemble.

        :param descriptor: STRING. One of AVAILABLE_DESCRIPTORS
        :param domain: None (whole RNA) or a collection of cg. elements.
                       If domain is present, calculate the descriptor only for these cg-elements.
        """
        if domain is None:
            if descriptor not in self._descriptors:
                self._descriptors[descriptor] = self._get_descriptor(
                    descriptor, domain)
            return self._descriptors[descriptor]
        else:  # No caching with domains.
            return self._get_descriptor(descriptor, domain)


class EnsembleView(EnsembleBase):
    def __init__(self, ens, start, end):
        """
        A View on a subsequence of an Ensemble.
        Implements the same interface like an ensemble.
        """
        if start is None:
            start = 0
        if end is None:
            end = len(ens)
        if not isinstance(start, int):
            raise ValueError("Start of EnsembleView must be integer.")
        if not isinstance(end, int):
            raise ValueError("End of EnsembleView must be integer.")
        if start < 0 or start > end or end < 0:
            raise ValueError(
                "Start and End of slice EnsembleView must be positive and end>start!.")
        if end > len(ens):
            end = len(ens)
        self._ensemble = ens
        self._start = start
        self._end = end

    def __len__(self):
        return self._end - self._start

    def __getitem__(self, i):
        if i > len(self):
            raise IndexError(i)
        return self._ensemble[self._start + i]

    def get_descriptor(self, descriptor, domain=None):
        """
        See Ensemble.get_descriptor.

        Gets the descriptor only for part of the ensemble that is in this view.
        """
        if domain is None:
            if descriptor not in self._ensemble._descriptors:
                return self._get_descriptor(descriptor, domain)
            else:
                return self._ensemble._descriptors[descriptor][self._start:self._end]
        else:  # No caching with domains.
            return self._get_descriptor(descriptor, domain)


from math import sqrt
import multiprocessing as mp
import numpy as np
from six.moves import range
from six.moves import zip


def t_scan(L, window=1000, num_workers=-1):
    """
    Computes t statistic for i to i+window points versus i-window to i
    points for each point i in input array. Uses multiple processes to
    do this calculation asynchronously. Array is decomposed into window
    number of frames, each consisting of points spaced at window
    intervals. This optimizes the calculation, as the drone function
    need only compute the mean and variance for each set once.

    Parameters
    ----------
    L : numpy array
        1 dimensional array that represents time series of datapoints
    window : int / float
        Number of points that comprise the windows of data that are
        compared
    num_workers : int
        Number of worker processes for multithreaded t_stat computation
        Defult value uses num_cpu - 1 workers


    Returns
    -------
    t_stat : numpy array
        Array which holds t statistic values for each point. The first
        and last (window) points are replaced with zero, since the t
        statistic calculation cannot be performed in that case.

    """
    size = L.size
    window = int(window)
    frames = list(range(window))
    n_cols = (size // window) - 1

    t_stat = np.zeros((window, n_cols))

    if num_workers == 1:
        results = [_t_scan_drone(L, n_cols, frame, window) for frame in frames]
    else:
        if num_workers == -1:
            num_workers = mp.cpu_count() - 1
        pool = mp.Pool(processes=num_workers)
        results = [pool.apply_async(_t_scan_drone, args=(
            L, n_cols, frame, window)) for frame in frames]
        results = [r.get() for r in results]
        pool.close()

    for index, row in results:
        t_stat[index] = row

    t_stat = np.concatenate((
        np.zeros(window),
        t_stat.transpose().ravel(order='C'),
        np.zeros(size % window)
    ))

    return t_stat


def _t_scan_drone(L, n_cols, frame, window=1e3):
    """
    Drone function for t_scan. Not Intended to be called manually.
    Computes t_scan for the designated frame, and returns result as
    array along with an integer tag for proper placement in the
    aggregate array
    """
    size = L.size
    window = int(window)
    root_n = sqrt(window)

    output = np.zeros(n_cols)
    b = L[frame:window + frame]
    b_mean = b.mean()
    b_var = b.var()
    for i in range(window + frame, size - window, window):
        a = L[i:i + window]
        a_mean = a.mean()
        a_var = a.var()
        output[i // window - 1] = root_n * \
            (a_mean - b_mean) / sqrt(a_var + b_var)
        b_mean, b_var = a_mean, a_var

    return frame, output


def mz_fwt(x, n=2):
    """
    Computes the multiscale product of the Mallat-Zhong discrete forward
    wavelet transform up to and including scale n for the input data x.
    If n is even, the spikes in the signal will be positive. If n is odd
    the spikes will match the polarity of the step (positive for steps
    up, negative for steps down).

    This function is essentially a direct translation of the MATLAB code
    provided by Sadler and Swami in section A.4 of the following:
    http://www.dtic.mil/dtic/tr/fulltext/u2/a351960.pdf

    Parameters
    ----------
    x : numpy array
        1 dimensional array that represents time series of data points
    n : int
        Highest scale to multiply to


    Returns
    -------
    prod : numpy array
        The multiscale product for x

    """
    N_pnts = x.size
    lambda_j = [1.5, 1.12, 1.03, 1.01][0:n]
    if n > 4:
        lambda_j += [1.0] * (n - 4)

    H = np.array([0.125, 0.375, 0.375, 0.125])
    G = np.array([2.0, -2.0])

    Gn = [2]
    Hn = [3]
    for j in range(1, n):
        q = 2**(j - 1)
        Gn.append(q + 1)
        Hn.append(3 * q + 1)

    S = np.concatenate((x[::-1], x))
    S = np.concatenate((S, x[::-1]))
    prod = np.ones(N_pnts)
    for j in range(n):
        n_zeros = 2**j - 1
        Gz = _insert_zeros(G, n_zeros)
        Hz = _insert_zeros(H, n_zeros)
        current = (1.0 / lambda_j[j]) * np.convolve(S, Gz)
        current = current[N_pnts + Gn[j]:2 * N_pnts + Gn[j]]
        prod *= current
        if j == n - 1:
            break
        S_new = np.convolve(S, Hz)
        S_new = S_new[N_pnts + Hn[j]:2 * N_pnts + Hn[j]]
        S = np.concatenate((S_new[::-1], S_new))
        S = np.concatenate((S, S_new[::-1]))
    return prod


def _insert_zeros(x, n):
    """
    Helper function for mz_fwt. Splits input array and adds n zeros
    between values.
    """
    newlen = (n + 1) * x.size
    out = np.zeros(newlen)
    indices = list(range(0, newlen - n, n + 1))
    out[indices] = x
    return out


def find_steps(array, threshold):
    """
    Finds local maxima by segmenting array based on positions at which
    the threshold value is crossed. Note that this thresholding is
    applied after the absolute value of the array is taken. Thus,
    the distinction between upward and downward steps is lost. However,
    get_step_sizes can be used to determine directionality after the
    fact.

    Parameters
    ----------
    array : numpy array
        1 dimensional array that represents time series of data points
    threshold : int / float
        Threshold value that defines a step


    Returns
    -------
    steps : list
        List of indices of the detected steps

    """
    steps = []
    array = np.abs(array)
    above_points = np.where(array > threshold, 1, 0)
    ap_dif = np.diff(above_points)
    cross_ups = np.where(ap_dif == 1)[0]
    cross_dns = np.where(ap_dif == -1)[0]
    for upi, dni in zip(cross_ups, cross_dns):
        steps.append(np.argmax(array[upi:dni]) + upi)
    return steps


def get_step_sizes(array, indices, window=1000):
    """
    Calculates step size for each index within the supplied list. Step
    size is determined by averaging over a range of points (specified
    by the window parameter) before and after the index of step
    occurrence. The directionality of the step is reflected by the sign
    of the step size (i.e. a positive value indicates an upward step,
    and a negative value indicates a downward step). The combined
    standard deviation of both measurements (as a measure of uncertainty
    in step calculation) is also provided.

    Parameters
    ----------
    array : numpy array
        1 dimensional array that represents time series of data points
    indices : list
        List of indices of the detected steps (as provided by
        find_steps, for example)
    window : int, optional
        Number of points to average over to determine baseline levels
        before and after step.


    Returns
    -------
    step_sizes : list
        List of the calculated sizes of each step
    step_error : list

    """
    step_sizes = []
    step_error = []
    indices = sorted(indices)
    last = len(indices) - 1
    for i, index in enumerate(indices):
        if i == 0:
            q = min(window, indices[i + 1] - index)
        elif i == last:
            q = min(window, index - indices[i - 1])
        else:
            q = min(window, index - indices[i - 1], indices[i + 1] - index)
        a = array[index:index + q]
        b = array[index - q:index]
        step_sizes.append(a.mean() - b.mean())
        step_error.append(sqrt(a.var() + b.var()))
    return step_sizes, step_error


def damped_randomwalk(start, length, mean=0, std=1, damp=0.9, transition=1000):
    x = np.zeros(length)
    x[0] = start
    for i in range(1, length):
        if i < transition:
            x[i] = x[i - 1] * damp + \
                (np.random.randn() * std + mean * (i / transition))
        else:
            x[i] = x[i - 1] * damp + (np.random.randn() * std + mean)
    return x


def rand_or_stay(start, length, mean=0, std=1, stayrate=0.1):
    x = np.zeros(length)
    x[0] = start
    for i in range(1, length):
        if np.random.rand() < stayrate:
            x[i] = x[i - 1]
        else:
            x[i] = np.random.randn() * std + mean
    return x


if __name__ == "__main__":
    import pymbar
    """x1 = damped_randomwalk(0,20000)
    x2 = rand_or_stay(x1[-1], 20000, stayrate = 0.8)
    x3 = damped_randomwalk(x2[-1], 20000, 3)
    x4 = damped_randomwalk(x3[-1], 20000, -1, 2)
    x5 = rand_or_stay(x4[-1], 20000, 2, 5, 0.02)
    x = np.concatenate([x1,x2,x3,x4,x5])
    overall_analysis(x)
    sys.exit()
    print("Testing pymbar")
    x1 = np.random.randn(20000)
    x2 = np.zeros(20000)
    x2[0] = x1[-1]
    for i in range(1,20000):
        x2[i]=x2[i-1]*0.99+np.random.randn()
    x3 = np.zeros(20000)
    x3[0] = x2[-1]
    for i in range(1,20000):
        x3[i]=x3[i-1]*0.79+np.random.randn()*2
    x4 = np.zeros(20000)
    x4[0] = x3[-1]
    for i in range(1,20000):
        x4[i]=x4[i-1]*0.95+np.random.randn()*0.5
    x5 = np.zeros(20000)
    x5[0] = x4[-1]
    for i in range(1,20000):
        x5[i]=x5[i-1]/i+np.random.randn()+1.1*(1-1/i)
    x = np.concatenate((x1,x2,x3,x4,x5))
    fig, ax = plt.subplots(2)
    t = t_scan(x)#overall_analysis(rogs, "ROG")
    ax[0].plot(np.arange(len(x)), x)
    ax[1].plot(np.arange(len(t)), t)
    overall_analysis(x)

    fig, ax = plt.subplots(2)
    ax[0].plot(np.arange(len(x1)), x1, label = "rand, {}".format(len(x1)/pymbar.timeseries.statisticalInefficiency(x1) ))
    #ax[0].plot(np.arange(len(x2)), x2, label = "rand with higher std, {}".format(len(x2)/pymbar.timeseries.statisticalInefficiency(x2) ))
    ax[0].plot(np.arange(len(x3)), x3, label = "damped random walk, {}".format(len(x3)/pymbar.timeseries.statisticalInefficiency(x3) ))
    ax[0].legend()
    ax[1].plot(np.arange(len(x1)-1), x1[1:]-x1[:-1])
    #ax[1].plot(np.arange(len(x2)-1), x2[1:]-x2[:-1])
    ax[1].plot(np.arange(len(x3)-1), x3[1:]-x3[:-1])
    plt.show()
    #ax[1].plot(np.arange(len(t)), t)
    plt.show()
    sys.exit()"""
    if len(sys.argv) == 2:
        rogs = []
        energies = []
        with open(sys.argv[1]) as f:
            rog_field = None
            energy_field = None
            for line in f:
                line = line.strip()
                if line[0] == "#":
                    continue
                fields = line.split("\t")
                if rog_field is None:  # The header line
                    for i, field in enumerate(fields):
                        if field == "Sampling_Energy":
                            energy_field = i
                        if field == "ROG":
                            rog_field = i
                            break
                    else:
                        print(
                            "Data seems not to contain any ROG field. Exiting", file=sys.stderr)
                        sys.exit()
                else:
                    e = fields[energy_field]
                    energies.append(float(e))
                    f = fields[rog_field].split()
                    try:
                        # Not an ROG Field? Expected FLOAT-SPACE-"A" format
                        assert f[1] == "A"
                    except:
                        print(f)
                        raise
                    rogs.append(float(f[0]))
        rogs = np.array(rogs)
        t = t_scan(rogs)
        overall_analysis(rogs, "ROG")
        fig, ax = plt.subplots(2)
        ax[0].plot(np.arange(len(rogs)), rogs)
        ax[1].plot(np.arange(len(t)), t)
        plt.show()
    else:
        ens = ensemble_from_filenames(sys.argv[1:])
        rogs = ens.get_descriptor("rog")
        overall_analysis(rogs, "ROG")
        rogs = ens.get_descriptor("anisotropy")
        overall_analysis(rogs, "Anisotropy")
        rogs = ens.get_descriptor("asphericity")
        overall_analysis(rogs, "Asphericity")
