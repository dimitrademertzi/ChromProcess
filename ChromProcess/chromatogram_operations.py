def pickPeaks(time, signal, threshold = 0.1):
    '''
    Finds peaks in signal and returns their times and indices

    Parameters
    ----------
    time, signal: 1D arrays
        Signal from which peaks will be picked and its corresponding time axis.
    threshold: float
        threshold for peak picking as a proportion of the maximum value of the
        signal.

    Returns
    -------
    retention_times: list
        list of peak retention times (floats)

    peak_windows:
        Time windows for peaks
    '''
    import numpy as np
    from ChromProcess import processing_functions as p_f
    
    peaks = p_f.Peak_finder(signal, thres = threshold, min_dist = 0.1,
                            max_inten = 1e100, min_inten = -1e100)

    retention_times = []
    peak_windows = []

    for p in range(0,len(peaks["Peak_indices"])):

        peak_idx = peaks["Peak_indices"][p]
        start_idx  = peaks["Peak_start_indices"][p]
        end_idx = peaks["Peak_end_indices"][p]

        rt = time[peak_idx]
        peak_inds = np.where( (time >= time[start_idx]) &
                              (time <= time[end_idx]) )[0]

        if len(peak_inds) > 0:
            retention_times.append(rt)
            peak_windows.append(time[peak_inds])


    return retention_times, peak_windows

def find_peaks_scipy(time, intensity, threshold=0.1, min_dist=1, max_inten = 1e100, prominence = 0.7, wlen = 101, peak_window = 12):
    '''Peak detection for the GC, meant to replace both pickPeaks and Peak_finder.
    Rather than using the first derivative this function is based on scipy.signal.find_peaks, 
    this allows us to specify prominence and find peaks in the more noisy GC spectra.

    Parameters
    ----------
    time, signal: 1D arrays
        Signal from which peaks will be picked and its corresponding time axis.
    thres : float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    min_dist : int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.
    max_inten: float
        peaks will not be detected above this threshold.
    min_inten: float
        peaks will not be detected below this threshold.
    Prominence: number or ndarray or sequence, optional
        Required prominence of peaks. Either a number, None, an array matching 
        x or a 2-element sequence of the former. The first element is always 
        interpreted as the minimal and the second, if supplied, as the maximal required prominence.
    wlen: int, optional
        Used for calculation of the peaks prominences, 
        thus it is only used if one of the arguments prominence or width is given. 
        See argument wlen in peak_prominences for a full description of its effects.
    peak_window: int
        Number of spaces the function will search through to find the start or end of a peak

    Returns
    -------
    ndarray
        Array containing the indexes of the peaks that were detected
    '''
    
    import numpy as np
    from ChromProcess import processing_functions as p_f
    from scipy.signal import find_peaks

    smooth_intensity = p_f.savitzky_golay(intensity, 7, 3, deriv=0, rate=1) #smooth out the function before starting the picking
    peaks, properties = find_peaks(smooth_intensity,  
                                    distance = min_dist, 
                                    height = [max(smooth_intensity)*threshold,max_inten], 
                                    prominence = prominence, 
                                    wlen = wlen)

    peak_starts = []
    for n in range(0,len(peaks)):
        cursor = peaks[n]-2
        if cursor < peak_window:
            print(f'peak window overlapping with peak start at RT {time[cursor]}')
            while smooth_intensity[cursor] > 0 and not cursor-1 > 0:
                cursor -= 1
        else:
            while smooth_intensity[cursor] > np.min(smooth_intensity[cursor-peak_window:cursor-1]) and not cursor-1 < peak_window:
                cursor -= 1
        
        peak_starts.append(cursor)

    peak_ends = []

    for n in range(0,len(peaks)):
        cursor = peaks[n]+2
        if cursor+peak_window > len(smooth_intensity):
            print(f'peak window overlapping with peak end at RT {time[cursor]}')
            while smooth_intensity[cursor] > 0 and not cursor + 1 < len(smooth_intensity):
                cursor += 1
        else:
            while smooth_intensity[cursor] > np.min(smooth_intensity[cursor+1:cursor+peak_window]) and not cursor+1+peak_window > len(smooth_intensity):
                cursor += 1
        
        peak_ends.append(cursor)

    retention_times = []
    peak_windows = []

    for p in range(0,len(peaks)):

        peak_idx = peaks[p]
        start_idx  = peak_starts[p]
        end_idx = peak_ends[p]

        rt = time[peak_idx]
        peak_inds = np.where( (time >= time[start_idx]) &
                              (time <= time[end_idx]) )[0]

        if len(peak_inds) > 0:
            retention_times.append(rt)
            peak_windows.append(time[peak_inds])

    return retention_times, peak_windows

def getpeaksfromfile(chromatogram, Peakfile, peak_window = 12):
    '''
    Get the peak retention times from a PeakCollections file rather than directly from a chromatogram.
    Peak boundaries will be generated. 

    chromatogram: Chromatogram object
        Chromatogram object to which the peaks should be added
    Peakfile: string
        Location of peak collection csv
    peak_window: int 
        Number of spaces the function will search through to find the start or end of a peak

    Returns
    -------
    ndarray
        Array containing the indexes of the peaks that were detected
    '''

    import numpy as np
    import pandas as pd
    from ChromProcess import processing_functions as p_f
    from ChromProcess import chromatogram_operations as chrom_ops

    df = pd.read_csv(Peakfile,index_col=False)
    peak_retention_times = df.iloc[3:,0].values #This function is very specific to the file type
    peaks = np.empty(0,dtype='int64') #the functions written to find the start and end of the peak relies on indexes, so we need to convert the retention times to their corresponding indexes.
    for p_rt in peak_retention_times:
        idx = np.where(chromatogram.time == float(p_rt)) #This covers any exact matches
        if len(idx[0]) == 0: #no exact match found, this should mean that a value was manually inserted, the function will now search the closest time value.
            search_idx = np.searchsorted(chromatogram.time, float(p_rt), "left")
            max_signal = np.max(chromatogram.signal[search_idx-4:search_idx+4]) #search the 9 nearest value for the highest peak, this gives a bit of play for time input
            idx = np.where(chromatogram.signal == max_signal)
        peaks = np.append(peaks, int(idx[0]))
    smooth_intensity = p_f.savitzky_golay(chromatogram.signal, 7, 3, deriv=0, rate=1) #smooth out the function before starting the picking
    



    peak_starts = []
    for n in range(0,len(peaks)):
        cursor = peaks[n]-2
        if cursor < peak_window:
            print(f'peak window overlapping with peak start at RT {chromatogram.time[cursor]}')
            while smooth_intensity[cursor] > 0 and not cursor-1 > 0:
                cursor -= 1
        else:
            while smooth_intensity[cursor] > np.min(smooth_intensity[cursor-peak_window:cursor-1]) and not cursor-1 < peak_window:
                cursor -= 1
        
        peak_starts.append(cursor)

    peak_ends = []

    for n in range(0,len(peaks)):
        cursor = peaks[n]+2
        if cursor+peak_window > len(smooth_intensity):
            print(f'peak window overlapping with peak end at RT {chromatogram.time[cursor]}')
            while smooth_intensity[cursor] > 0 and not cursor + 1 < len(smooth_intensity):
                cursor += 1
        else:
            while smooth_intensity[cursor] > np.min(smooth_intensity[cursor+1:cursor+peak_window]) and not cursor+1+peak_window > len(smooth_intensity):
                cursor += 1
        
        peak_ends.append(cursor)

    retention_times = []
    peak_windows = []

    for p in range(0,len(peaks)):

        peak_idx = peaks[p]
        start_idx  = peak_starts[p]
        end_idx = peak_ends[p]

        rt = chromatogram.time[peak_idx]
        peak_inds = np.where( (chromatogram.time >= chromatogram.time[start_idx]) &
                              (chromatogram.time <= chromatogram.time[end_idx]) )[0]

        if len(peak_inds) > 0:
            retention_times.append(rt)
            peak_windows.append(chromatogram.time[peak_inds])

    peak_features = []
    for rt, times in zip(retention_times, peak_windows):
        feature = [times[0], rt, times[-1]]
        peak_features.append(feature)

    chrom_ops.addPeaksToChromatogram(peak_features, chromatogram)

   

def regionPeakPick(chromatogram, region, threshold = 0.1, spectra_type = 'GC-MS', min_dist=1, max_inten = 1e100, prominence = 0.7, wlen = 11, peak_window = 12):
    '''
    Pick peaks in a chromatogram within the region given.

    Parameters
    ----------
    chromatogram: Chromatogram object

    region: list
        [lower bound, upper bound]
    threshold: float
        Peaks below this fraction of the highest 
        intensity of the chromatogram will not be picked.

    Returns
    -------
    None
    '''
    import numpy as np
    from ChromProcess import chromatogram_operations as chrom_ops

    low_limit = region[0]
    high_limit = region[1]
    
    if low_limit > chromatogram.time.max():
        # there is nothing to find outside the chromatogram
        pass
    elif high_limit < chromatogram.time.min():
        # there is nothing to find outside the chromatogram
        pass
    else:
        # get the indices of the chromatogram region
        inds = np.where( (chromatogram.time > low_limit) &
                        (chromatogram.time < high_limit) )[0]
        # get peaks in the chromatogram region 
        if spectra_type == 'GC':
            rts, peak_times = chrom_ops.GCpickPeaks(
                                chromatogram.time[inds],
                                chromatogram.signal[inds],
                                threshold = threshold, 
                                min_dist=1, 
                                max_inten = 1e100, 
                                prominence = 0.5, 
                                wlen = 1001, 
                                peak_window = 12
                                )
        else:
            rts, peak_times = chrom_ops.pickPeaks(
                                            chromatogram.time[inds],
                                            chromatogram.signal[inds],
                                            threshold = threshold
                                            )
        peak_features = []
        for rt, times in zip(rts,peak_times):
            feature = [times[0], rt, times[-1]]
            peak_features.append(feature)

        chrom_ops.addPeaksToChromatogram(peak_features, chromatogram)
       
def getPeakFromBounds(chrom, start, end):
    '''
    Create a peak using the boundaries 
    defined within chrom

    chrom: Chromatogram object
    start: float
        Start of the peak
    end: float
        End of the peak
    
    peak: Peak object
        peak created from the bounds
    '''
    import numpy as np
    from ChromProcess import Classes

    if start > end:
        print(f'peak start ({start}) > peak end, ({end}) returning None')
        return None

    inds = np.where(
            (chrom.time > start) & 
            (chrom.time < end)
            )[0]

    timeseg = chrom.time[inds]
    sigseg  = chrom.signal[inds]

    peak_idx = np.argmax(sigseg)
    retention_time = timeseg[peak_idx] 
 
    peak = Classes.Peak(retention_time, inds)

    return peak

def addPeaksToChromatogram(peak_times, chromatogram):
    '''
    peak_times: list
        [start, peak, end] in units of the chromatogram's
        retention time axis.

    chromatogram: Chromatogram object
    '''
    import numpy as np
    from ChromProcess import Classes

    # iterate through indices found and use them to
    # create peak objects which are added to the 
    # chromatogram
    for p in peak_times:
        # find indices of the peaks upper and lower
        # bounds in relationship to the entire 
        # chromatogram.
        idx = np.where((chromatogram.time >= p[0]) &
                        (chromatogram.time <= p[-1]))[0]

        peak = Classes.Peak(p[1], idx)
        chromatogram.peaks[p[1]] = peak

def integratePeaks(chromatogram, baseline_subtract = False):
    '''
    Parameters
    ----------

    chromatogram: Chromatogram object

    Return
    ------
    None
    '''
    # put integral information into the peaks.
    for p in chromatogram.peaks:
        chromatogram.peaks[p].get_integral(chromatogram, baseline_subtract = baseline_subtract)

def mz_background_subtraction(chromatogram, threshold = 500):
    '''
    Gets the ion chromatograms of the analysis and reconstitutes the total ion
    chromatogram ommitting mass signals which do not exceed a threshold.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object.
        Chromatogram to be processed.
    threshold: float
        Ion chromatograms which do not exceed this threshold will be removed
        before the total ion chromatogram is reconsituted.

    Returns: 1D numpy array.
        original signal if not mass spectra information present.
        processed signal if ms info present.
    '''
    import numpy as np

    if len(chromatogram.mass_intensity) == 0:
        return chromatogram.signal
    else:
        inds = chromatogram.mass_intensity < threshold
        mass_intens = np.copy(chromatogram.mass_intensity)
        mass_intens[inds] = 0.0

        new_chromatogram = np.zeros(len(chromatogram.time))
        for s in range(0,len(chromatogram.point_counts)):

            start = chromatogram.scan_indices[s]
            end = start + chromatogram.point_counts[s]

            inten = mass_intens[start:end]

            new_chromatogram[s] = np.sum(inten)

        return new_chromatogram

def get_mz_background(chromatogram, threshold = 500):
    '''
    Gets the ion chromatograms of the analysis and reconstitutes the total ion
    chromatogram ommitting mass signals which do not exceed a threshold.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object.
        Chromatogram to be processed.
    threshold: float
        Ion chromatograms which do not exceed this threshold will be removed
        before the total ion chromatogram is reconsituted.

    Returns: None
        Modifies the chromatogram in-place.
    '''
    import numpy as np

    if len(chromatogram.mass_intensity) == 0:
        print("MS_intensity_threshold_chromatogram: no mass spectra information in chromatogram")
        return chromatogram.signal
    else:
        inds = chromatogram.mass_intensity > threshold
        mass_intens = np.copy(chromatogram.mass_intensity)
        mass_intens[inds] = 0.0

        new_chromatogram = np.zeros(len(chromatogram.time))
        for s in range(0,len(chromatogram.point_counts)):

            start = chromatogram.scan_indices[s]
            end = start + chromatogram.point_counts[s]

            inten = mass_intens[start:end]

            new_chromatogram[s] = np.sum(inten)

        return new_chromatogram

def internalRefIntegral(chromatogram, internal_ref_region):
    '''
    Get integrals of internal references for a Chromatogram_Series
    Parameters
    ----------
    series: Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.
    '''
    from ChromProcess import Classes
    from ChromProcess import chromatogram_operations as chrom_ops

    if len(chromatogram.peaks) == 0:
        chrom_ops.regionPeakPick(chromatogram, internal_ref_region, threshold = 0.1)
    else:
        pass

    for p in chromatogram.peaks:
        if p > internal_ref_region[0] and p < internal_ref_region[1]:
            rt = chromatogram.peaks[p].retention_time
            indices = chromatogram.peaks[p].indices
            chromatogram.internal_reference = Classes.Peak(rt, indices)
            chromatogram.internal_reference.get_integral(chromatogram)

def getIonChromatogramsFromRegion(chromatogram, lower, upper, threshold = 0.1):
    '''
    Get the ion chromatograms for a region of a chromatogram. Requires mass
    spectra information to be present in the chromatogram.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
        Chromatogram containing information
    lower: float
        Lower bound for retention time region in chromatogram
    upper: float
        Upper bound for retention time region in chromatogram

    threshold: float
        Threshold for mass spectra extraction relative to the maximum signal
        in the region.

    Returns
    -------

    ion_dict: dict
        dictionary of ion chromatograms
        {m/z: intensities over time}
    '''

    import numpy as np

    if not chromatogram.mass_spectra:
        return None
    else:
        ion_dict = {}
        inds = np.where((chromatogram.time > lower)&(chromatogram.time < upper))[0]

        time = chromatogram.time[inds]
        scan_inds = chromatogram.scan_indices[inds]
        p_counts = chromatogram.point_counts[inds]

        for s in range(0,len(time)):

            inten = chromatogram.mass_intensity[scan_inds[s]:scan_inds[s]+p_counts[s]]
            masses = chromatogram.mass_values[scan_inds[s]:scan_inds[s]+p_counts[s]]

            if len(inten) > 0:

                filt_inds = np.where(inten > threshold*np.amax(inten))[0]
                inten = inten[filt_inds]
                masses = masses[filt_inds]

                round = np.round(masses, 2)

                for m in range(0,len(round)):
                    if round[m] in ion_dict:
                        ion_dict[round[m]][s] = inten[m]
                    else:
                        ion_dict[round[m]] = np.zeros(len(time))
                        ion_dict[round[m]][s] = inten[m]
            else:
                pass

    return ion_dict

def RegionSVD(chromatogram, region, ic_threshold = 0.1,
                          components = 5):

    '''
    Reduce a GC-MS chromatogram region dimensionality using singular value
    decomposition. The chromatogram must contain ion chromatogram information.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
        Chromatogram to work on.
    region: list, tuple or array
        Two-element list of floats defining the lower and upper bounds of the
        region.

    Returns
    -------

    U: numpy array
        Unitary matrix ith left singular vectors

    S: numpy array
        array of singular values

    Vh: numpy array
        Unitary matrix with right singular vectors


    '''
    import numpy as np
    from ChromProcess import chromatogram_operations as chrom_ops
    from ChromProcess import simple_functions as simp_func

    ic_dict = chrom_ops.getIonChromatogramsFromRegion(chromatogram,
                                                      region[0], region[1],
                                                      threshold = ic_threshold)
    if ic_dict == None:
        return None
    else:
        # convert dict into nump arrays
        ic_stack = np.array(list(ic_dict.values()))
        trans_IC = ic_stack.T


        U, S, Vh = simp_func.runSVD(trans_IC)

        # inverse transform
        #IT = np.dot(U*S, Vh)

        return U[:,:components], S[:components], Vh[:components]
