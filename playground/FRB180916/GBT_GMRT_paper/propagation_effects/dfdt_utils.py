"""Utilities for linear drift rate measurements."""

import numpy as np
import scipy.signal

def delay_from_dm(dm, freq_emitted):
    """Return the delay in seconds caused by dispersion in the
    interstellar medium.
    Parameters
    ----------
    dm : float
        Dispersion measure,  in pc cm-3.
    freq_emitted : float
        Observing frequency, in MHz.
    Returns
    -------
    float
        Dispersive delay, in seconds.
    """
    if type(freq_emitted) == type(0.0):
        if freq_emitted > 0.0:
            return dm / (0.000241 * freq_emitted * freq_emitted)
        else:
            return 0.0
    else:
        return np.where(
            freq_emitted > 0.0, dm / (0.000241 * freq_emitted * freq_emitted),
            0.0
        )


def boxcar_kernel(width):
    """Returns the boxcar kernel of given width normalized by
    sqrt(width) for S/N reasons.
    Parameters
    ----------
    width : int
        Width of the boxcar.
    Returns
    -------
    boxcar : array_like
        Boxcar of width `width` normalized by sqrt(width).
    """
    width = int(round(width, 0))
    return np.ones(width, dtype="float32") / np.sqrt(width)


def find_burst(ts, width_factor=4, min_width=1, max_width=128):
    """Find burst peak and width using boxcar convolution.
    Parameters
    ----------
    ts : array_like
        Time-series.
    width_factor : int, optional
        Windowing factor for on and off-pulse determination.
    min_width : int, optional
        Minimum width to search from, in number of time samples.
        1 by default.
    max_width : int, optional
        Maximum width to search up to, in number of time samples.
        128 by default.
    Returns
    -------
    peak : int
        Index of the peak of the burst in the time-series.
    width : int
        Width of the burst in number of samples.
    snr : float
        S/N of the burst.
    """
    min_width = int(min_width)
    max_width = int(max_width)

    # do not search widths bigger than timeseries
    widths = list(range(min_width,
                        min(max_width + 1, int((len(ts) - 50) // 6))))

    # envelope finding
    snrs = np.empty_like(widths, dtype=float)
    peaks = np.empty_like(widths, dtype=int)

    # borders for on and off-pulse determination
    outer = 3 * width_factor // 2
    inner = width_factor // 2

    for i in range(len(widths)):
        convolved = scipy.signal.convolve(ts, boxcar_kernel(widths[i]))
        peaks[i] = np.nanargmax(convolved)
        # peak should not be on the edge of time-series
        if (peaks[i] > 0.9 * ts.shape[0]) or (peaks[i] < 0.1 * ts.shape[0]):
            snrs[i] = np.nan
        else:
            # get RMS for S/N weighting, as in PRESTO's single_pulse_search.py
            baseline = np.concatenate(
                [
                    convolved[0 : max(0, peaks[i] - 3 * widths[i])],
                    convolved[peaks[i] + 3 * widths[i] :],
                ]
            )

            # cutoff of at least 50 samples is a bit arbitrary, but seems
            # reasonable
            if baseline.shape[0] > 50:
                rms = np.std(baseline)
            else:
                rms = np.nan

            snrs[i] = convolved[peaks[i]] / rms

    best_idx = np.nanargmax(snrs)

    return peaks[best_idx], widths[best_idx], snrs[best_idx]


def add_text(ax, text, color="black"):
    """Add text in top left corner of plot."""
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    xtext = 0.03 * (xlim[1] - xlim[0]) + xlim[0]
    ytext = 0.95 * (ylim[1] - ylim[0]) + ylim[0]

    ax.text(xtext, ytext, text, ha="left", va="top", color=color)


def shift_channels(intensity, bins):
    """Shift each channel to the left by the corresponding value in
    bins.
    Parameters
    ----------
    bins : array_like
        Array containing the number of bins to shift each channel
        by.
    Note
    ----
    *** Shifting happens in-place! ***
    """
    nchan, nsamp = intensity.shape
    
    nchan = len(bins)

    pad = np.ones_like(bins)

    for ii in range(nchan):
        intensity[ii, :] == np.roll(intensity[ii, :], -bins[ii], axis=0)

def dedisperse(intensity, center_frequencies, dt_s, dm=0.0,
               reference_frequency=600.0):
    """Shift channels according to the delays predicted by the given
    dispersion measure.
    Parameters
    ----------
    intensity : array_like
        Intensity array of shape (nchan, nsamp).
    center_frequencies : array_like
        Center frequencies of nchan channels, in MHz.
    dt_s : float
        Sampling time, in seconds.
    dm : float
        Dispersion measure, in pc cm-3.
    reference_frequency : float
        Reference frequency for DM calculation, in MHz. 600 MHz by
        default.
    Note
    ----
    *** Dedispersion happens in place! ***
    """
    dm = float(dm)
    reference_frequency = float(reference_frequency)

    # calculate the integer shifts for the full DM and
    # then subtract the existing integer shifts. Reduces rounding error.
    ref_delay = delay_from_dm(dm, reference_frequency)
    delays = delay_from_dm(dm, center_frequencies)

    # relative delay
    rel_delays = delays - ref_delay
    rel_bindelays = np.round(rel_delays / dt_s).astype("int")

    shift_channels(intensity, rel_bindelays)
    
