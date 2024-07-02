{-# LANGUAGE FlexibleInstances #-}

module FFTSpec where

import Data.Complex as C (Complex, magnitude)
import Data.Vector (Vector)
import qualified Data.Vector as Vector
import FFT (dft, fft, idft, ifft, irdft, irfft, rdft, rfft)
import System.Exit (exitFailure, exitSuccess)
import Test.Hspec (Spec)
import Test.Hspec.QuickCheck (prop)
import Test.QuickCheck (Arbitrary (arbitrary), quickCheckAll)

spec :: Spec
spec = do
  prop "can roundtrip discrete fourier transform" prop_dftInverse
  prop "can roundtrip real discrete fourier transform" prop_rdftInverse
  prop "can roundtrip real fast fourier transform" prop_rfftInverse
  prop "can roundtrip fast fourier transform" prop_fftInverse
  prop "dft is close to fft" prop_fft
  prop "real dft is close to real fft" prop_rfft

-- Helpers

allClose :: (Ord a, Fractional a) => Vector a -> Vector a -> Bool
allClose xs ys = all (\diff -> diff < 1e-4) diffs
 where
  diffs = Vector.zipWith (\x y -> abs $ x - y) xs ys

allCloseComplex :: RealFloat a => Vector (Complex a) -> Vector (Complex a) -> Bool
allCloseComplex xs ys = all (\diff -> diff < 1e-4) diffs
 where
  diffs = Vector.zipWith (\x y -> C.magnitude $ x - y) xs ys

instance Arbitrary a => Arbitrary (Vector a) where
  arbitrary = Vector.fromList <$> arbitrary

-- Properties

prop_dftInverse :: Vector (Complex Double) -> Bool
prop_dftInverse xs = allCloseComplex xs $ idft (dft xs)

prop_rdftInverse :: Vector Double -> Bool
prop_rdftInverse xs = allClose xs $ irdft (rdft xs)

prop_fftInverse :: Vector (Complex Double) -> Bool
prop_fftInverse xs = allCloseComplex xs $ ifft (fft xs)

prop_rfftInverse :: Vector Double -> Bool
prop_rfftInverse xs = allClose xs $ irfft (rfft xs)

prop_fft :: Vector (Complex Double) -> Bool
prop_fft xs = allCloseComplex (fft xs) (dft xs)

prop_rfft :: Vector Double -> Bool
prop_rfft xs = allCloseComplex (rfft xs) (rdft xs)
