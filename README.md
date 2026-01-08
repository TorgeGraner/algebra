# Algebra
Calculates the jordan normal form of a given trigonalizable Matrix with coefficients in a ring of characteristic zero using exact arithmetic. 

The ring must be implemented as a class implementing the operators ```+, -, *, /``` and ```%```, aswell as a constructor from an integer. Note that division might throw an exception if the divisor does not divide the divident.

Since some functions, such as matrix inversion, require multiplicative inverses, when using a Ring R, results might have to be multiplied by a scalar. Take for example the following integer matrix, which is invertible over the rationals but not over the integers:

```math
A:=\begin{bmatrix}
  -2 & 22  & 1 & -16 & 0\\\
  1  & 13  & 0 & -10 & 1\\\
  2  & -18 & 0 & 12  & 0\\\
  0  & -9  & 0 & 6   & 0\\\
  0  & 11  & 0 & -8  & 0
\end{bmatrix}\quad\text{with}\quad A^{-1}=\begin{bmatrix}
  0 & 0 & -1 & 2 & 0\\\
  0 & 0 & 0 & -\frac{4}{3} & -1\\\
  1 & 0 & 1 & -2 & -2\\\
  0 & 0 & 0 & -\frac{11}{6} & \frac{3}{2}\\\
  0 & 1 &-\frac{1}{2} & 0 &-2
\end{bmatrix}\in\mathbb{Q}^{5\times5}\quad\text{and}\quad 6A^{-1}=\begin{bmatrix}
  0 & 0 & 3  & -6  & 0\\\
  0 & 0 & 0  & -8  & -6\\\
  6 & 0 & 6  & -12 & -12\\\
  0 & 0 & 0  & -11 & -9\\\
  0 & 6 & -3 & 0   &-12
\end{bmatrix}\in\mathbb{Z}^{5\times5}
```

Using the scaled inverse therefore allows us to stay in the ring of integers, even though later results might be scaled accordingly. This problem does not occur when using a field.

## Algorithm
The algorithm calculates the characteristic polynomial of the given matrix and searches (very naively) for a factorization. For each linear factor, it determines a base for the corresponding generalized Eigenspaces and combines them to a complete base. If the characteristic polynomial does not break down completely, return a partial base.


## Notes
Note that all datastructures are defined in the .hpp, since template classes cannot be seperated into .hpp and .cpp. The problem comes from the compiler not being able to infer the same template parameter for .hpp and .cpp. This problem could be fixed by always including the corresponding .cpp together with the .hpp.

Even though the native type int can be used as a ring, under given circumstances this can produce problems. For example the statement ```Matrix<Fractionalize<int>> A = 0``` would not compile, since the compiler would get confused with the current type. Use the wrapper class Integer instead. 

The algorithm does not work well with floating point arithmetic, even with a wrapper class, since rounding error can still influence the rank of a matrix, leading to faulty logic producing errors.

Some array is allocated somewhere that is never deallocated.
