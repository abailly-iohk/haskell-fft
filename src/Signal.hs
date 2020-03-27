module Signal where

import Data.Complex as C
import FFT
import Utils

-- Hann window function.
-- scipy.signal.hann
-- https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.hann.html
hann :: (Floating a) => Int -> [a]
hann n = [0.5 - 0.5 * cos ((2 * pi * i) / (n' - 1)) | i <- map fromIntegral [0..n-1]]
  where
    n' = fromIntegral n

-- Welch Power Spectral Density (PSD) Estimation.
-- scipy.signal.welch
-- https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.welch.html
--
-- xs: Timeseries data
-- windowSize: Length of each window to compute FFT.
-- overlap: Number of points to overlap between windows.
-- Returns a tuple of sample frequencies and the power spectral density of xs
--
-- Unlike the scipy implementation, returns PSD corresponding to all frequencies
-- (including those of the negative frequencies, which are simply duplicates of
-- the positive frequencies for real data).
welch :: (RealFloat a) => [a] -> Int -> Int -> ([a], [a])
welch xs windowSize overlap = (fftfreq windowSize 1.0, welchPsd)
  where
    welchPsd = map (/ norm) $ sumFrequencies $ magnitudeSquared <$> rfft <$> applyHann <$> windows
    norm = hannWindowNorm * (fromIntegral . length $ windows)
    -- sum powers over all windows correponding to the frequency bins
    sumFrequencies = foldr (zipWith (+)) (repeat 0)
    -- compute square magnitude of frequencies computed by FFT for a given window
    magnitudeSquared fs = (^2) . C.magnitude <$> fs
    -- multiply by Hann windowing function for a given data window
    applyHann window = zipWith (*) hannWindow window
    -- create windows of size exactly `windowSize` with the specified overlap
    windows = filter (\w -> length w == windowSize) [slice start windowSize xs | start <- [0,stride..n-1]]
    slice start len = drop start . take (start + len)
    n = length xs
    stride = windowSize - overlap
    -- Hann window function of length `windowSize`
    hannWindow = hann windowSize
    hannWindowNorm = sum $ (^2) <$> hannWindow
