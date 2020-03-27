# haskell-fft

Straight-forward implementation of 1-D Fast Fourier Transform (FFT), its inverse, and Welch Power Spectral Density (PSD) estimation in Haskell. Implements the Cooley-Tukey FFT algorithm for real and complex timeseries data. Variants of the following numpy functions are implemented:
- numpy.fft.fft
- numpy.fft.ifft
- numpy.fft.rfft
- numpy.fft.irfft
- numpy.fft.fftfreq
- scipy.signal.hann
- scipy.signal.welch

# Requirements

[Haskell Stack](https://docs.haskellstack.org/en/stable/README/)

# Installation

```
stack install
stack test
```

# Example Usage

FFT (and inverse) of real-valued timeseries:
```
import FFT

fourier_components = rfft [1..10]
frequencies = fftfreq 10 1
timeseries = irfft fourier_components
```

FFT (and inverse) of complex-valued timeseries:
```
import FFT
import Data.Complex

fourier_components = fft [mkPolar i (pi/i) | i <- [1..10]]
frequencies = fftfreq 10 1
timeseries = ifft fourier_components
```

Welch Power Spectral Density (PSD) estimation:
```
import Signal
(frequencies, psd) = welch [1..100] 10 2
```
