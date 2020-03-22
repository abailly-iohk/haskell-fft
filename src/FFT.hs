module FFT where

import Control.Applicative

import Data.Matrix as M
import Data.Complex as C

import Utils

-- 1-D Discrete Fourier Transform (DFT), the naive way.
dft :: (RealFloat a) => [Complex a] -> [Complex a]
dft xs =
  let n = length xs
      theta i j = 2 * pi * (fromIntegral i) * (fromIntegral j) / (fromIntegral n)
      w = M.matrix n n (\(i,j) -> C.mkPolar 1 $ theta (i-1) (j-1))
      xs' = M.fromList n 1 xs
  in M.toList $ vdot w xs'

-- Inverse transform for `dft`.
idft :: (RealFloat a) => [Complex a] -> [Complex a]
idft fs =
  let n = length fs
      theta i j = -2 * pi * (fromIntegral i) * (fromIntegral j) / (fromIntegral n)
      norm = fromIntegral n
      winv = M.matrix n n (\(i,j) -> C.mkPolar (1/norm) $ theta (i-1) (j-1))
      fs' = M.fromList n 1 fs
  in M.toList $ vdot winv fs'

-- 1-D DFT for real-valued vector. Unlike the numpy equivalent, returns all
-- frequency components including the negative frequency ones (which are simply
-- complex conjugates of the positive frequency components).
rdft :: (RealFloat a) => [a] -> [Complex a]
rdft xs = dft $ toComplex <$> xs
  where
    toComplex magnitude = C.mkPolar magnitude 0

-- Inverse transform for `rdft`.
irdft :: (RealFloat a) => [Complex a] -> [a]
irdft = map C.realPart . idft

-- numpy.fft.fft
-- 1-D Fast Fourier Transform (FFT).
fft :: (RealFloat a) => [Complex a] -> [Complex a]
fft xs
  | n == 0 = []
  | n == 1 = xs
  | even n = positiveFreqs ++ negativeFreqs
  | otherwise = dft xs
  where
    n = length xs
    positiveFreqs = zipWith3 (\evenX oddX w -> evenX + w * oddX) evenFreqs oddFreqs weights
    negativeFreqs = zipWith3 (\evenX oddX w -> evenX - w * oddX) evenFreqs oddFreqs weights
    evenFreqs = fft $ [xi | (i,xi) <- zip [0..] xs, even i]
    oddFreqs = fft $ [xi | (i,xi) <- zip [0..] xs, odd i]
    weights = [C.conjugate . C.mkPolar 1 $ theta k | k <- [0 .. (n `div` 2)-1]]
    theta k = 2 * pi * (fromIntegral k) / (fromIntegral n)

-- numpy.fft.ifft
-- Inverse transform for `fft`.
ifft :: (RealFloat a) => [Complex a] -> [Complex a]
ifft fs = normalize . ifft' $ fs
  where
    normalize xs = map (/ (fromIntegral (length xs))) xs
    ifft' fs
      | n == 0 = []
      | n == 1 = fs
      | even n = positiveFreqs ++ negativeFreqs
      | otherwise =  unnormalize . idft $ fs -- undo normalization done by idft do avoid double normalization
      where
        n = length fs
        positiveFreqs = zipWith3 (\evenX oddX w -> evenX + w * oddX) evenFreqs oddFreqs weights
        negativeFreqs = zipWith3 (\evenX oddX w -> evenX - w * oddX) evenFreqs oddFreqs weights
        evenFreqs = ifft' $ [xi | (i,xi) <- zip [0..] fs, even i]
        oddFreqs = ifft' $ [xi | (i,xi) <- zip [0..] fs, odd i]
        weights = [C.mkPolar 1 $ theta k | k <- [0 .. (n `div` 2)-1]]
        theta k = 2 * pi * (fromIntegral k) / (fromIntegral n)
        unnormalize = map (* (fromIntegral n))

-- numpy.fft.rfft
-- 1-D FFT for real-valued vector. Unlike the numpy equivalent, returns all
-- frequency components including the negative frequency ones (which are simply
-- complex conjugates of the positive frequency components).
rfft :: (RealFloat a) => [a] -> [Complex a]
rfft xs = fft $ toComplex <$> xs
  where
    toComplex magnitude = C.mkPolar magnitude 0

-- numpy.fft.irfft
-- Inverse transform for `rfft`.
irfft :: (RealFloat a) => [Complex a] -> [a]
irfft = map C.realPart . ifft

-- numpy.fft.fftfreq
-- Computes frequency bins.
-- n: Window length
-- d: Sample spacing (inverse of sampling rate)
fftfreq :: (Fractional a) => Int -> a -> [a]
fftfreq n d = map ((/ norm) . fromIntegral) $ freqs n
  where
    freqs n
      | n == 0 = []
      | otherwise = [0..n_2-1] ++ [i - n | i <- [n_2..n-1]]
    n_2 = if even n then n `div` 2 else (n + 1) `div` 2
    norm = (fromIntegral n) * d
