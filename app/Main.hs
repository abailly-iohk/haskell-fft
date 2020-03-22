module Main where

import FFT

main :: IO ()
main = do
  putStrLn "Frequency components for signal [1..4]:"
  print $ rfft [1..4]
  putStrLn "Corresponding sample frequencies:"
  print $ fftfreq 4 1
