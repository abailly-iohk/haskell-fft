module Signal where

import Data.Complex as C
import Data.Vector (Vector, slice)
import qualified Data.Vector as Vector
import FFT
import Utils

-- Hann window function.
-- scipy.signal.hann
-- https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.hann.html
hann :: Floating a => Int -> Vector a
hann n = Vector.fromList [0.5 - 0.5 * cos ((2 * pi * i) / (n' - 1)) | i <- map fromIntegral [0 .. n - 1]]
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
welch :: RealFloat a => Vector a -> Int -> Int -> (Vector a, Vector a)
welch xs windowSize overlap = (fftfreq windowSize 1.0, welchPsd)
 where
  welchPsd = fmap (/ norm) $ sumFrequencies $ magnitudeSquared <$> rftHann
  rftHann = rfft <$> hannWindows
  hannWindows = applyHann <$> windows
  norm = hannWindowNorm * (fromIntegral . length $ windows)
  -- sum powers over all windows correponding to the frequency bins
  sumFrequencies ys = Vector.foldr (Vector.zipWith (+)) (Vector.replicate (length ys) 0) ys
  -- compute square magnitude of frequencies computed by FFT for a given window
  magnitudeSquared fs = (^ 2) . C.magnitude <$> fs
  -- multiply by Hann windowing function for a given data window
  applyHann window = Vector.zipWith (*) hannWindow window
  -- create windows of size exactly `windowSize` with the specified overlap
  windows = Vector.fromList $ filter (\w -> length w == windowSize) [slice start windowSize xs | start <- [0, stride .. n - 1]]
  n = length xs
  stride = windowSize - overlap
  -- Hann window function of length `windowSize`
  hannWindow = hann windowSize
  hannWindowNorm = sum $ (^ 2) <$> hannWindow
