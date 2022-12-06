import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.pyplot as plt
import os


def _1gaussian(x, amplitude, centre, sigma, baseline=0):
    """
    A single gaussian function

    Parameters
    ----------
    x: array
        x axis data
    amplitude1: float
        amplitudelitude of the function
    centre1: float
        centretre of the function (mean)
    sigma1: float
        width of the function (standard deviation)

    Returns
    -------
    function: numpy array
        y values for the function
    """
    return (
        amplitude
        * (1 / (sigma * (np.sqrt(2 * np.pi))))
        * (np.exp(-((x - centre) ** 2) / ((2 * sigma) ** 2)))
        + baseline
    )


def _2gaussian(x, amplitude1, centre1, sigma1, amplitude2, centre2, sigma2, baseline=0):
    """
    A double gaussian function

    Parameters
    ----------
    x: array
        x axis data
    amplituden: float
        amplitudelitude of a component gaussian function
    centren: float
        centretre of a component gaussian function (mean)
    sigman: float
        width of a component gaussian function (standard deviation)

    Returns
    -------
    function: numpy array
        y values for the function
    """
    return _1gaussian(x, amplitude1, centre1, sigma1, baseline) + _1gaussian(
        x, amplitude2, centre2, sigma2
    )


def _3gaussian(
    x,
    amplitude1,
    centre1,
    sigma1,
    amplitude2,
    centre2,
    sigma2,
    amplitude3,
    centre3,
    sigma3,
    baseline=0,
):
    return _1gaussian(x, amplitude1, centre1, sigma1, baseline) + _2gaussian(
        x, amplitude2, centre2, sigma2, amplitude3, centre3, sigma3
    )


def _4gaussian(
    x,
    amplitude1,
    centre1,
    sigma1,
    amplitude2,
    centre2,
    sigma2,
    amplitude3,
    centre3,
    sigma3,
    amplitude4,
    centre4,
    sigma4,
    baseline=0,
):
    return _1gaussian(x, amplitude1, centre1, sigma1, baseline) + _3gaussian(
        x,
        amplitude2,
        centre2,
        sigma2,
        amplitude3,
        centre3,
        sigma3,
        amplitude4,
        centre4,
        sigma4,
    )


def _5gaussian(
    x,
    amplitude1,
    centre1,
    sigma1,
    amplitude2,
    centre2,
    sigma2,
    amplitude3,
    centre3,
    sigma3,
    amplitude4,
    centre4,
    sigma4,
    amplitude5,
    centre5,
    sigma5,
    baseline=0,
):
    return _1gaussian(x, amplitude1, centre1, sigma1, baseline) + _4gaussian(
        x,
        amplitude2,
        centre2,
        sigma2,
        amplitude3,
        centre3,
        sigma3,
        amplitude4,
        centre4,
        sigma4,
        amplitude5,
        centre5,
        sigma5,
    )


def _6gaussian(
    x,
    amplitude1,
    centre1,
    sigma1,
    amplitude2,
    centre2,
    sigma2,
    amplitude3,
    centre3,
    sigma3,
    amplitude4,
    centre4,
    sigma4,
    amplitude5,
    centre5,
    sigma5,
    amplitude6,
    centre6,
    sigma6,
    baseline=0,
):
    return _1gaussian(x, amplitude1, centre1, sigma1, baseline) + _5gaussian(
        x,
        amplitude2,
        centre2,
        sigma2,
        amplitude3,
        centre3,
        sigma3,
        amplitude4,
        centre4,
        sigma4,
        amplitude5,
        centre5,
        sigma5,
        amplitude6,
        centre6,
        sigma6,
        baseline=0,
    )


def _7gaussian(
    x,
    amplitude1,
    centre1,
    sigma1,
    amplitude2,
    centre2,
    sigma2,
    amplitude3,
    centre3,
    sigma3,
    amplitude4,
    centre4,
    sigma4,
    amplitude5,
    centre5,
    sigma5,
    amplitude6,
    centre6,
    sigma6,
    amplitude7,
    centre7,
    sigma7,
    baseline=0,
):
    return _1gaussian(x, amplitude1, centre1, sigma1, baseline) + _6gaussian(
        x,
        amplitude2,
        centre2,
        sigma2,
        amplitude3,
        centre3,
        sigma3,
        amplitude4,
        centre4,
        sigma4,
        amplitude5,
        centre5,
        sigma5,
        amplitude6,
        centre6,
        sigma6,
        amplitude7,
        centre7,
        sigma7,
    )


def _8gaussian(
    x,
    amplitude1,
    centre1,
    sigma1,
    amplitude2,
    centre2,
    sigma2,
    amplitude3,
    centre3,
    sigma3,
    amplitude4,
    centre4,
    sigma4,
    amplitude5,
    centre5,
    sigma5,
    amplitude6,
    centre6,
    sigma6,
    amplitude7,
    centre7,
    sigma7,
    amplitude8,
    centre8,
    sigma8,
    baseline=0,
):
    return _1gaussian(x, amplitude1, centre1, sigma1, baseline) + _7gaussian(
        x,
        amplitude2,
        centre2,
        sigma2,
        amplitude3,
        centre3,
        sigma3,
        amplitude4,
        centre4,
        sigma4,
        amplitude5,
        centre5,
        sigma5,
        amplitude6,
        centre6,
        sigma6,
        amplitude7,
        centre7,
        sigma7,
        amplitude8,
        centre8,
        sigma8,
    )


def _9gaussian(
    x,
    amplitude1,
    centre1,
    sigma1,
    amplitude2,
    centre2,
    sigma2,
    amplitude3,
    centre3,
    sigma3,
    amplitude4,
    centre4,
    sigma4,
    amplitude5,
    centre5,
    sigma5,
    amplitude6,
    centre6,
    sigma6,
    amplitude7,
    centre7,
    sigma7,
    amplitude8,
    centre8,
    sigma8,
    amplitude9,
    centre9,
    sigma9,
    baseline=0,
):
    return _1gaussian(x, amplitude1, centre1, sigma1, baseline) + _8gaussian(
        x,
        amplitude2,
        centre2,
        sigma2,
        amplitude3,
        centre3,
        sigma3,
        amplitude4,
        centre4,
        sigma4,
        amplitude5,
        centre5,
        sigma5,
        amplitude6,
        centre6,
        sigma6,
        amplitude7,
        centre7,
        sigma7,
        amplitude8,
        centre8,
        sigma8,
        amplitude9,
        centre9,
        sigma9,
    )


def _10gaussian(
    x,
    amplitude1,
    centre1,
    sigma1,
    amplitude2,
    centre2,
    sigma2,
    amplitude3,
    centre3,
    sigma3,
    amplitude4,
    centre4,
    sigma4,
    amplitude5,
    centre5,
    sigma5,
    amplitude6,
    centre6,
    sigma6,
    amplitude7,
    centre7,
    sigma7,
    amplitude8,
    centre8,
    sigma8,
    amplitude9,
    centre9,
    sigma9,
    amplitude10,
    centre10,
    sigma10,
    baseline=0,
):
    return _1gaussian(x, amplitude1, centre1, sigma1, baseline) + _9gaussian(
        x,
        amplitude2,
        centre2,
        sigma2,
        amplitude3,
        centre3,
        sigma3,
        amplitude4,
        centre4,
        sigma4,
        amplitude5,
        centre5,
        sigma5,
        amplitude6,
        centre6,
        sigma6,
        amplitude7,
        centre7,
        sigma7,
        amplitude8,
        centre8,
        sigma8,
        amplitude9,
        centre9,
        sigma9,
        amplitude10,
        centre10,
        sigma10,
    )


def fit_gaussian_peaks(
    time,
    signal,
    initial_guess,
    boundaries,
    num_peaks,
):
    """
    TODO: This kind of function could be useful, but a better adapted function
        for peak deconvolution should be written. The function could take
        similar concepts to this one, but with a different concept for its
        implementation.

    Fitting sums of gaussian peaks to data using supplied peak indices.

    Parameters
    ----------
    time: array
        time values
    sig: array
        signal to be fit to
    peaks: list of peak positions
        list peak positions in the time

    initial_guess: list
        Initial guess for the peak amplitudelitude, position and width
        e.g. see _1gaussian() function arguments.

    lowerbounds: list
        Lower bounds for the peak amplitudelitude, position and width
        e.g. see _1gaussian() function arguments.
    upperbounds: list
        Upper bounds for the peak amplitudelitude, position and width
        e.g. see _1gaussian() function arguments.

    Returns
    -------
    popt: ndarray
        list of fitted values [[amplitudelitude, centretre, width],]
    pcov: array
        correlation matrix
    """

    if num_peaks == 1:
        popt, pcov = curve_fit(
            _1gaussian, time, signal, p0=initial_guess, bounds=boundaries
        )
    elif num_peaks == 2:
        popt, pcov = curve_fit(
            _2gaussian, time, signal, p0=initial_guess, bounds=boundaries
        )
    elif num_peaks == 3:
        popt, pcov = curve_fit(
            _3gaussian, time, signal, p0=initial_guess, bounds=boundaries
        )  # ,method='trf')
    elif num_peaks == 4:
        popt, pcov = curve_fit(
            _4gaussian, time, signal, p0=initial_guess, bounds=boundaries
        )
    elif num_peaks == 5:
        popt, pcov = curve_fit(
            _5gaussian, time, signal, p0=initial_guess, bounds=boundaries
        )
    elif num_peaks == 6:
        popt, pcov = curve_fit(
            _6gaussian, time, signal, p0=initial_guess, bounds=boundaries
        )
    elif num_peaks == 7:
        popt, pcov = curve_fit(
            _7gaussian, time, signal, p0=initial_guess, bounds=boundaries
        )
    elif num_peaks == 8:
        popt, pcov = curve_fit(
            _8gaussian, time, signal, p0=initial_guess, bounds=boundaries
        )
    elif num_peaks == 9:
        popt, pcov = curve_fit(
            _9gaussian, time, signal, p0=initial_guess, bounds=boundaries
        )
    elif num_peaks == 10:
        popt, pcov = curve_fit(
            _10gaussian, time, signal, p0=initial_guess, bounds=boundaries
        )
    else:
        print("Error: number of peaks to large")
        popt, pcov = [0, 0, 0], 0

    return popt, pcov


def deconvolute_region(chromatogram, region, peak_diff, num_peaks=1):

    """
    TODO: Combine the ideas in this function with fit_gaussian_peaks()

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
        Chromatogram
    region: list
        region of chromatogram under operation [lower bound, upper bound]

    Returns
    -------
    popt: ndarray
        list of fitted values [[amplitude, centre, width],]
    pcov: array
        correlation matrix
    """

    upper = region[1]
    lower = region[0]

    inds = np.where((chromatogram.time > lower) & (chromatogram.time < upper))[0]

    time = chromatogram.time[inds]
    signal = chromatogram.signal[inds]

    signal = signal - np.average(signal[-5:-1])

    for p in range(0, len(peaks)):  # extend the boundary
        initial_guess[0] = np.amax(signal)
        initial_guess[1] = peaks[p]
        lowerbounds[1] = time[0]
        upperbounds[1] = time[-1]
        guess.extend(initial_guess)
        lbds.extend(lowerbounds)
        ubds.extend(upperbounds)

    # peak_list = np.array([*chromatogram.peaks])
    # peak_inds = np.where((peak_list > lower)&(peak_list < upper))[0]
    #
    # peaks = peak_list[peak_inds]
    #
    # while len(peaks) < num_peaks:
    #    peaks = np.append(peaks, np.average(peaks))
    #
    # if len(peaks) > num_peaks:
    #    peaks = peaks[:num_peaks]

    return fit_gaussian_peaks(time, signal, peaks)


def deconvolute_peak(
    chromatogram,
    peak_folder,
    indices,
    info_dict,
    plotting=True,
):
    """
    TODO: this function is quite similar in scope to deconvolute_region().
    Refactor with the other two deconvolution macros.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
        Chromatogram
    region: list
        region of chromatogram under operation [lower bound, upper bound]
    """
    from ChromProcess import Classes

    time = chromatogram.time[indices[0] : indices[-1]]
    signal = chromatogram.signal[indices[0] : indices[-1]]
    # baseline = np.interp(time, [time[0], time[-1]], [min(signal), min(signal)])
    signal = signal  # - baseline
    # lower_bounds = []
    # upper_bounds = []
    # for p in range(0,num_peaks): # extend the boundary
    #    standard_lower_bounds[0] = min(signal)*0.1
    #    standard_lower_bounds[1] = time[0]
    #    standard_upper_bounds[0] = max(signal)*3
    #    standard_upper_bounds[1] = time[-1]
    #    lower_bounds.extend(standard_lower_bounds)
    #    upper_bounds.extend(standard_upper_bounds)
    #
    if "peaks" in info_dict:
        initial_guess = []
        boundaries = [[],[]]
        for peak_dict in info_dict["peaks"].values():
            if "amplitude" in peak_dict:
                boundaries[0].append(peak_dict["amplitude"][0])
                boundaries[1].append(peak_dict["amplitude"][2])
                initial_guess.append(peak_dict["amplitude"][1])
            else:
                boundaries[0].append(info_dict["standard_amplitude"][0])
                boundaries[1].append(info_dict["standard_amplitude"][2])
                initial_guess.append(info_dict["standard_amplitude"][1])

            boundaries[0].append(peak_dict["centre"][0])
            boundaries[1].append(peak_dict["centre"][2])
            initial_guess.append(peak_dict["centre"][1])
            if "sigma" in peak_dict:
                boundaries[0].append(peak_dict["sigma"][0])
                boundaries[1].append(peak_dict["sigma"][2])
                initial_guess.append(peak_dict["sigma"][1])
            else:
                boundaries[0].append(info_dict["standard_sigma"][0])
                boundaries[1].append(info_dict["standard_sigma"][2])
                initial_guess.append(info_dict["standard_sigma"][1])
            #add values for the automatic baseline adaptation
        initial_guess.append(0)
        boundaries[0].append(0)
        boundaries[1].append(min(signal)+0.001)

    else: #this code was used for the old file system
        if "lower_fit_boundaries" in info_dict:
            boundaries = [
                [*info_dict["lower_fit_boundaries"], 0],
                [*info_dict["upper_fit_boundaries"], min(signal) + 0.001],
            ]
        else:
            lower_bounds = [0, time[0], 0]
            upper_bounds = [1e7, time[-1], 0.03]
            boundaries = [
                info_dict["number_of_peaks"] * lower_bounds,
                info_dict["number_of_peaks"] * upper_bounds,
            ]

        if not "initial_guess" in info_dict:
            info_dict["initial_guess"] = []
            for n in range(1, info_dict["number_of_peaks"] + 1):
                info_dict["initial_guess"].extend(
                    [
                        max(signal),
                        time[0]
                        + ((time[-1] - time[0]) / (info_dict["number_of_peaks"] + 2)) * n,
                        0.009,
                    ]
                )
            initial_guess = [*info_dict["initial_guess"], 0]

    popt, pcov = fit_gaussian_peaks(
        time,
        signal,
        initial_guess,
        boundaries,
        info_dict["number_of_peaks"],
    )
    peaks = []
    baseline = np.ones(len(signal)) * popt[-1]
    fig, ax = plt.subplots(1, 1)
    ax.plot(time, signal)
    final_curve = np.zeros(len(signal))

    ax.plot(time, baseline)
    for p in range(0, info_dict["number_of_peaks"]):
        amp = popt[0 + p * 3]
        centre = popt[1 + p * 3]
        sigma = popt[2 + p * 3]
        gauss = _1gaussian(time, amp, centre, sigma)
        pk_idx = np.argmin(abs(time - centre))
        start_idx = np.argmin(abs(time - centre + 3 * sigma))
        end_idx = np.argmin(abs(time - centre - 3 * sigma))
        integral = np.trapz(gauss, x=time)
        retention_time = time[pk_idx]
        start = time[start_idx]
        end = time[end_idx]
        peaks.append(
            Classes.Peak(
                retention_time,
                start,
                end,
                integral=integral,
                indices=[],
                height=max(gauss),
            )
        )
        final_curve += gauss
        ax.plot(time, gauss + baseline)

    if plotting == True:
        ax.plot(time, final_curve + baseline)
        fig.set_size_inches(18.5, 10.5)
        plt.tight_layout()
        plt.savefig(f"{peak_folder}\\{chromatogram.filename}.png")
        plt.close()
    else:
        plt.close()

    mse = sum((signal - final_curve) ** 2) / len(signal)
    return popt, pcov, mse, peaks
