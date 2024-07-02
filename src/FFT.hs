{-# LANGUAGE OverloadedLists #-}

module FFT where

import Control.Applicative

import Data.Complex as C
import Data.Matrix as M

import Data.Vector (Vector, ifilter)
import qualified Data.Vector as Vector
import Utils

-- 1-D Discrete Fourier Transform (DFT), the naive way.
dft :: RealFloat a => Vector (Complex a) -> Vector (Complex a)
dft xs =
  let n = length xs
      theta i j = 2 * pi * fromIntegral i * fromIntegral j / fromIntegral n
      w = M.matrix n n (\(i, j) -> C.mkPolar 1 $ theta (i - 1) (j - 1))
      xs' = M.colVector xs
   in M.getMatrixAsVector $ vdot w xs'

-- Inverse transform for `dft`.
idft :: RealFloat a => Vector (Complex a) -> Vector (Complex a)
idft fs =
  let n = length fs
      theta i j = (-2) * pi * fromIntegral i * fromIntegral j / fromIntegral n
      norm = fromIntegral n
      winv = M.matrix n n (\(i, j) -> C.mkPolar (1 / norm) $ theta (i - 1) (j - 1))
      fs' = M.colVector fs
   in M.getMatrixAsVector $ vdot winv fs'

-- 1-D DFT for real-valued vector. Unlike the numpy equivalent, returns all
-- frequency components including the negative frequency ones (which are simply
-- complex conjugates of the positive frequency components).
rdft :: RealFloat a => Vector a -> Vector (Complex a)
rdft xs = dft $ toComplex <$> xs
 where
  toComplex magnitude = C.mkPolar magnitude 0

-- Inverse transform for `rdft`.
irdft :: RealFloat a => Vector (Complex a) -> Vector a
irdft = fmap C.realPart . idft

-- numpy.fft.fft
-- 1-D Fast Fourier Transform (FFT).
fft :: RealFloat a => Vector (Complex a) -> Vector (Complex a)
fft xs
  | n == 0 = mempty
  | n == 1 = xs
  | even n = positiveFreqs <> negativeFreqs
  | otherwise = dft xs
 where
  n = length xs
  positiveFreqs = Vector.zipWith3 (\evenX oddX w -> evenX + w * oddX) evenFreqs oddFreqs weights
  negativeFreqs = Vector.zipWith3 (\evenX oddX w -> evenX - w * oddX) evenFreqs oddFreqs weights
  evenFreqs = fft $ ifilter (\i _ -> even i) xs
  oddFreqs = fft $ ifilter (\i _ -> odd i) xs
  weights = Vector.fromList [C.conjugate . C.mkPolar 1 $ theta k | k <- [0 .. n `div` 2 - 1]]
  theta k = 2 * pi * fromIntegral k / fromIntegral n

-- numpy.fft.ifft
-- Inverse transform for `fft`.
ifft :: RealFloat a => Vector (Complex a) -> Vector (Complex a)
ifft fs = normalize . ifft' $ fs
 where
  normalize xs = fmap (/ fromIntegral (length xs)) xs
  ifft' fs
    | n == 0 = mempty
    | n == 1 = fs
    | even n = positiveFreqs <> negativeFreqs
    | otherwise = unnormalize . idft $ fs -- undo normalization done by idft do avoid double normalization
   where
    n = length fs
    positiveFreqs = Vector.zipWith3 (\evenX oddX w -> evenX + w * oddX) evenFreqs oddFreqs weights
    negativeFreqs = Vector.zipWith3 (\evenX oddX w -> evenX - w * oddX) evenFreqs oddFreqs weights
    evenFreqs = ifft' $ ifilter (\i _ -> even i) fs
    oddFreqs = ifft' $ ifilter (\i _ -> odd i) fs
    weights = Vector.fromList [C.mkPolar 1 $ theta k | k <- [0 .. n `div` 2 - 1]]
    theta k = 2 * pi * fromIntegral k / fromIntegral n
    unnormalize = fmap (* fromIntegral n)

-- numpy.fft.rfft
-- 1-D FFT for real-valued vector. Unlike the numpy equivalent, returns all
-- frequency components including the negative frequency ones (which are simply
-- complex conjugates of the positive frequency components).
rfft :: RealFloat a => Vector a -> Vector (Complex a)
rfft xs = fft $ toComplex <$> xs
 where
  toComplex magnitude = C.mkPolar magnitude 0

-- numpy.fft.irfft
-- Inverse transform for `rfft`.
irfft :: RealFloat a => Vector (Complex a) -> Vector a
irfft = fmap C.realPart . ifft

-- numpy.fft.fftfreq
-- Computes frequency bins.
-- n: Window length
-- d: Sample spacing (inverse of sampling rate)
fftfreq :: Fractional a => Int -> a -> Vector a
fftfreq n d = fmap ((/ norm) . fromIntegral) $ freqs n
 where
  freqs n
    | n == 0 = mempty
    | otherwise = Vector.fromList $ [0 .. n_2 - 1] <> [i - n | i <- [n_2 .. n - 1]]
  n_2 = if even n then n `div` 2 else (n + 1) `div` 2
  norm = fromIntegral n * d
