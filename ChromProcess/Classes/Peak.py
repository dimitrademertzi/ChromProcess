
class Peak:
    def __init__(self,retention_time, indices):
        '''
        Creates a Peak object using a retention time
        and the indices of the places in the data's
        parent chromatogram from which the time and
        signal of the peak can be obtained.

        Parameters
        ----------
        retention_time: float
        indices: list
        '''

        self.retention_time = retention_time

        self.indices = indices

        self.integral = 0.0

        self.height = 0.0

        self.ion_chromatograms = {}

        self.ion_integrals = {}

        self.mass_spectrum = []

        self.deconvolution = []

    def peak_table_entry_text(self, chromatogram):
        '''
        For writing the peak's row in a peak table.
        Returns
        -------
        entry: str
            Row for the peak table output
        '''

        st_ind = 0
        end_ind = 0
        peak_start = 0.0
        peak_end = 0.0
        if len(self.indices) > 0:
            st_ind = self.indices[0]
            end_ind = self.indices[-1]
            peak_start = chromatogram.time[st_ind]
            peak_end = chromatogram.time[end_ind]

        rtn_time = self.retention_time
        integral = self.integral

        entry = f"{rtn_time},{integral},{peak_start},{peak_end}\n"

        return entry

    def get_integral(self, chromatogram, baseline_subtract = False):
        '''
        Get the integral of the peak using a chromatogram. Note that an
        arbitray chromatogram can be passed to this method, meaning it is not
        secure. The baseline substraction substracts a baseliner interpolated
        linearly between the start and the end of the peak.

        Parameters
        ----------
        chromatogram: ChromProcess Chromatogram object
        baseline_subtract: bool
        '''

        import numpy as np

        time = chromatogram.time[self.indices]
        signal = chromatogram.signal[self.indices]

        if baseline_subtract:
            time_bound = [time[0], time[-1]]
            signal_bound = [signal[0], signal[-1]]
            linterp = np.interp(time, time_bound, signal_bound)
            self.integral = ( np.trapz(signal - linterp, x = time) )
        else:
            self.integral = ( np.trapz(signal, x = time) )

        return self.integral

    def get_height(self, chromatogram):
        '''
        Get the height of the peak.

        Parameters
        ----------
        chromatogram: ChromProcess Chromatogram object
        '''

        self.height = chromatogram.signal[self.indices]

    def get_mass_spectrum(self, chromatogram):
        '''
        Get the mass spectrum at the apex of the peak.

        Parameters
        ----------
        chromatogram: ChromProcess Chromatogram object
        '''
        import numpy as np

        if len(chromatogram.scan_indices) == 0:
            # no mass spectral information in the chromatogram
            pass
        else:
            time = chromatogram.time

            ind = np.where(time == self.retention_time)[0]

            start = chromatogram.scan_indices[ind][0]
            end = start + chromatogram.point_counts[ind][0]

            self.mass_spectrum = [
                                np.round(chromatogram.mass_values[start:end],2),
                                chromatogram.mass_intensity[start:end]
                                ]


