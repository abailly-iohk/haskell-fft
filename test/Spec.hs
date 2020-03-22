{-# LANGUAGE TemplateHaskell #-}

import Test.QuickCheck
import System.Exit
import Data.Complex as C
import FFT

-- Helpers

allClose :: (Ord a, Fractional a) => [a] -> [a] -> Bool
allClose xs ys = and $ fmap (\diff -> diff < 1e-4) diffs
  where
    diffs = [abs $ x - y | (x,y) <- zip xs ys]

allCloseComplex :: (RealFloat a) => [Complex a] -> [Complex a] -> Bool
allCloseComplex xs ys = and $ fmap (\diff -> diff < 1e-4) diffs
  where
    diffs = [C.magnitude $ x - y | (x,y) <- zip xs ys]

-- Properties

prop_dftInverse :: [Complex Double] -> Bool
prop_dftInverse xs = allCloseComplex xs $ idft (dft xs)

prop_rdftInverse :: [Double] -> Bool
prop_rdftInverse xs = allClose xs $ irdft (rdft xs)

prop_fftInverse :: [Complex Double] -> Bool
prop_fftInverse xs = allCloseComplex xs $ ifft (fft xs)

prop_rfftInverse :: [Double] -> Bool
prop_rfftInverse xs = allClose xs $ irfft (rfft xs)

prop_fft :: [Complex Double] -> Bool
prop_fft xs = allCloseComplex (fft xs) (dft xs)

prop_rfft :: [Double] -> Bool
prop_rfft xs = allCloseComplex (rfft xs) (rdft xs)

-- https://hackage.haskell.org/package/QuickCheck-2.13.2/docs/Test-QuickCheck-All.html
return []
runTests :: IO Bool
runTests = $quickCheckAll

main :: IO ()
main = do
  success <- and <$> sequence [runTests]
  if success then exitSuccess else exitFailure
