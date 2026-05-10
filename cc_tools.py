import numpy as np
import torch
from torch import fft

def nextpow2(i):
    n = 1
    while(n < i): n *= 2
    return n

def torch_xcorr(signal_1, signal_2):
    if len(signal_1.shape)<2 | len(signal_2.shape)<2:
        print('input dimension must be ntrace*npts !')
        return 0
    else:
        signal_length = signal_1.shape[-1]
        x_cor_sig_length = signal_length*2 - 1
        fast_length = nextpow2(x_cor_sig_length)

        # The last signal_ndim axes will be transformed
        fft_1 = fft.rfft(signal_1, fast_length, dim=-1)
        fft_2 = fft.rfft(signal_2, fast_length, dim=-1)

        # Take the complex conjugate of one of the spectrums. 
        # Which one you choose depends on domain specific conventions
        fft_multiplied = torch.conj(fft_1) * fft_2

        # back to time domain.
        prelim_correlation = fft.irfft(fft_multiplied, dim=-1)

        # Shift the signal to make it look like a proper crosscorrelation,
        # and transform the output to be purely real
        final_result = torch.roll(prelim_correlation, 
                                  fast_length//2, dims=-1)[:,  fast_length//2-x_cor_sig_length//2:fast_length//2-x_cor_sig_length//2+x_cor_sig_length]
        
        return final_result

def computeCC(data, dt, max_lag, isource=0, ch_buffer_in=-1):
    """Function to compute CC between single DAS channel as virtual source"""
    nch = data.shape[0]
    ch_buffer = data.shape[0] if ch_buffer_in == -1 else ch_buffer_in
    nchunk = int(nch/ch_buffer+1.5) if ch_buffer < nch else 1
    max_lag_npts = int(max_lag/dt)
    cc = np.zeros((nch,max_lag_npts*2+1))
    npts = data.shape[1]
    dataDASsource = torch.from_numpy(data[isource,:npts])
    
    # Processing channel chunks for CC computation
    first_ch = 0
    for ichunk in range(nchunk):
        second_ch = min(first_ch+ch_buffer,nch)
        dataDAS = torch.from_numpy(data[first_ch:second_ch,:].copy()) # DAS data

        cc_all = torch_xcorr(dataDASsource.reshape(1, npts).repeat(dataDAS.shape[0], 1), dataDAS)/npts
        cc[first_ch:second_ch,:] = cc_all[:, npts-max_lag_npts:npts+max_lag_npts+1].numpy()

        first_ch += ch_buffer
    return cc

from scipy.signal import butter, filtfilt, detrend, convolve
from scipy.signal.windows import tukey

def running_absolute_mean(trace, nwin):
    '''
    reference: refer to noisepy package, but need to be improved faster on GPU
    :param trace: 1d array shape: npts
    :param nwin: # of points in moving window
    :return: smoothed data
    '''
    npts = len(trace)
    tmp = np.zeros(npts + 2 * nwin)
    tmp[nwin:-nwin] = np.abs(trace)
    tmp[:nwin] = tmp[nwin]
    tmp[-nwin:] = tmp[-nwin - 1]
    return np.nan_to_num(trace/convolve(tmp, np.ones(nwin) / nwin, mode='same')[nwin: -nwin])

def temporal_normalization(data, fs, window_time):
    '''
    running absolute mean normalization or one-bit, depending on window_time
    :param data: shape: nch * npts
    :param fs: sampling frequency
    :param window_time: running window length, in seconds. recommended: half the longest period
    :return: normalized data
    '''
    if window_time == 0: # 1-bit
        return np.sign(data)
    else:
        nwin = int(fs * window_time)
        nch = data.shape[0]
        for i in range(nch):
            data[i,:] = running_absolute_mean(data[i,:], nwin)
        return data
    
def spectral_whitening(rfftdata, df, window_freq, f1, f2):
    '''
    phase-only or running absolute mean spectral whitening, depending on window_freq
    :param rfftdata: shape: nch * npts, !!torch.tensor!!
    :param df: frequency interval
    :param window_freq: running window length, in Hz.
    :return: whitened spectra
    '''
    idxf1 = int(np.floor(f1 / df))
    idxf2 = int(np.ceil(f2 / df))
    rfftdata_angle = torch.angle(rfftdata)

    if window_freq == 0: # phase-only
        rfftdata = torch.exp(1j * rfftdata_angle)

    else: # running absolute mean
        nwin = int(window_freq / df)
        nch = rfftdata.shape[0]
        for i in range(nch):
            rfftdata[i,:] = torch.from_numpy(running_absolute_mean(rfftdata[i,:].cpu().numpy(), nwin))

    rfftdata[:, :idxf1] = torch.cos(torch.linspace(np.pi / 2, np.pi, idxf1, device=rfftdata.device)) ** 2 * rfftdata[:, :idxf1]
    rfftdata[:, idxf2:] = torch.cos(torch.linspace(np.pi, np.pi / 2, rfftdata.shape[-1] - idxf2, device=rfftdata.device)) ** 2 * rfftdata[:, idxf2:]

    return rfftdata


import torch
from torch import fft

def nextpow2(i):
    n = 1
    while(n < i): n *= 2
    return n

def torch_xcorr(signal_1, signal_2, whitening_params=None):
    if len(signal_1.shape)<2 | len(signal_2.shape)<2:
        print('input dimension must be ntrace*npts !')
        return 0
    else:
        signal_length = signal_1.shape[-1]
        x_cor_sig_length = signal_length*2 - 1
        fast_length = nextpow2(x_cor_sig_length)

        # The last signal_ndim axes will be transformed
        fft_1 = fft.rfft(signal_1, fast_length, dim=-1)
        fft_2 = fft.rfft(signal_2, fast_length, dim=-1)
        
        if(whitening_params is not None):
            fs, window_freq, f1, f2 = whitening_params
            df = fs/fast_length
            fft_1 = spectral_whitening(fft_1, df, window_freq, f1, f2)
            fft_2 = spectral_whitening(fft_2, df, window_freq, f1, f2)

        # Take the complex conjugate of one of the spectrums. 
        # Which one you choose depends on domain specific conventions
        fft_multiplied = torch.conj(fft_1) * fft_2

        # back to time domain.
        prelim_correlation = fft.irfft(fft_multiplied, dim=-1)

        # Shift the signal to make it look like a proper crosscorrelation,
        # and transform the output to be purely real
        final_result = torch.roll(prelim_correlation, 
                                  fast_length//2, dims=-1)[:,  fast_length//2-x_cor_sig_length//2:fast_length//2-x_cor_sig_length//2+x_cor_sig_length]
        
        return final_result
    
def computeCC(data, dt, max_lag, isource=0, ch_buffer_in=-1, whitening_params=None):
    """Function to compute CC between single DAS channel as virtual source"""
    nch = data.shape[0]
    ch_buffer = data.shape[0] if ch_buffer_in == -1 else ch_buffer_in
    nchunk = int(nch/ch_buffer+1.5) if ch_buffer < nch else 1
    max_lag_npts = int(max_lag/dt)
    cc = np.zeros((nch,max_lag_npts*2+1))
    npts = data.shape[1]
    dataDASsource = torch.from_numpy(data[isource,:npts])
    
    # Processing channel chunks for CC computation
    first_ch = 0
    for ichunk in range(nchunk):
        second_ch = min(first_ch+ch_buffer,nch)
        dataDAS = torch.from_numpy(data[first_ch:second_ch,:].copy()) # DAS data
        if dataDAS.shape[0] == 0:
            break
        cc_all = torch_xcorr(dataDASsource.reshape(1, npts).repeat(dataDAS.shape[0], 1), dataDAS, whitening_params)/npts
        #cc[first_ch:second_ch,:] = cc_all[:, npts-max_lag_npts:npts+max_lag_npts+1].numpy()
        cc[first_ch:second_ch,:] = cc_all[:, (npts-1)-max_lag_npts:(npts-1)+max_lag_npts+1].numpy()

        first_ch += ch_buffer
    return cc
