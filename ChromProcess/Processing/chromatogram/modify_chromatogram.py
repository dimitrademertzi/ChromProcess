import numpy as np
from ChromProcess.Processing.chromatogram import find_peaks


def add_peaks_to_chromatogram(peaks, chromatogram):
    """
    Add peaks to a chromatogram (modifies the chromatogram in place).

    Parameters
    ----------
    peaks: list of Peak objects

    chromatogram: Chromatogram object

    Returns
    ------ None
    """

    for peak in peaks:
        rt = peak.retention_time
        indices = np.where(
            (chromatogram.time >= peak.start) & (chromatogram.time <= peak.end)
        )[0]
        peak.indices = indices
        chromatogram.peaks[rt] = peak


def integrate_chromatogram_peaks(chromatogram, baseline_subtract=False):
    """
    Integrate all of the peaks in a chromatogram (modifies them in place).

    Parameters
    ----------
    chromatogram: Classes.Chromatogram object
        Chromatogram containing peaks.
    baseline_subtract: bool
        Whether to perform a local baseline subtraction on the peak.

    Returns
    ------
    None
    """

    for p in chromatogram.peaks:
        chromatogram.peaks[p].get_integral(
            chromatogram, baseline_subtract=baseline_subtract
        )


def internal_standard_integral(chromatogram, is_start, is_end):
    """
    Finds and adds internal standard information into a chromatogram.

    Parameters
    ----------
    series: Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.
    is_start: float
        Start of the region of the chromatogram in which the internal standard
        is found.
    is_end: float
        End of the region of the chromatogram in which the internal standard
        is found.

    Returns
    ------
    None
    """

    peaks = find_peaks.find_peaks_in_region(
        chromatogram, is_start, is_end, threshold=0.1
    )

    peak = peaks[0]

    peak.indices = np.where(
        (chromatogram.time >= peak.start) & (chromatogram.time <= peak.end)
    )[0]

    peak.get_integral(chromatogram)

    chromatogram.internal_standard = peak


def internal_standard_integral_look_ahead(chromatogram, is_start, is_end):
    """
    Finds and adds internal standard information into a chromatogram.

    Parameters
    ----------
    series: Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.
    is_start: float
        Start of the region of the chromatogram in which the internal standard
        is found.
    is_end: float
        End of the region of the chromatogram in which the internal standard
        is found.

    Returns
    ------
    None
    """
    from ChromProcess.Utils.peak_finding import find_peaks_scipy
    from ChromProcess.Utils import indices_from_boundary
    from ChromProcess import Classes

    inds = indices_from_boundary(chromatogram.time, is_start, is_end)
    time = chromatogram.time[inds]
    signal = chromatogram.signal[inds]
    picked_peaks = find_peaks_scipy(
        signal,
        threshold=0.01,
        min_dist=5,
        max_inten=1e100,
        prominence=0.1,
        wlen=1001,
        look_ahead=25,
        smooth_window=25,
    )

    peaks = []
    for x in range(0, len(picked_peaks["Peak_indices"])):
        pk_idx = picked_peaks["Peak_indices"][x]
        start_idx = picked_peaks["Peak_start_indices"][x]
        end_idx = picked_peaks["Peak_end_indices"][x]

        retention_time = time[pk_idx]
        start = time[start_idx]
        end = time[end_idx]
        peaks.append(
            Classes.Peak(retention_time, start, end, indices=[], height=signal[pk_idx])
        )

    peak = peaks[0]

    peak.indices = np.where(
        (chromatogram.time >= peak.start) & (chromatogram.time <= peak.end)
    )[0]

    peak.get_integral(chromatogram)

    chromatogram.internal_standard = peak
