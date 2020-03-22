module Utils where

import Data.Matrix
import qualified Data.Complex as C

type ComplexMatrix a = Matrix (C.Complex a)

-- Complex conjugate of a matrix
conjugate :: (Num a) => ComplexMatrix a -> ComplexMatrix a
conjugate = fmap C.conjugate

-- Similar to numpy.vdot. Takes complex conjugate of the first matrix
-- and computes matrix product.
vdot :: (RealFloat a) => ComplexMatrix a -> ComplexMatrix a -> ComplexMatrix a
vdot m n = (conjugate m) * n
