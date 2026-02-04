# Algebra
Calculates the jordan normal form of a given trigonalizable Matrix with coefficients in any euclidean ring using exact arithmetic. 

Currently decomposes a given integer matrix over ℤ, ℚ and the finite field $𝔽_p$, aswell as a wrapper class of double, modeling $ℝ$.

## Build
In project directory:
* Create debug/release subdirectory (e.g. ```mkdir release``` or ```mkdir debug```)

In corresponding subdirectory:
* Invoke ```cmake .. -DCMAKE_BUILD_TYPE = Release``` or ```cmake .. -DCMAKE_BUILD_TYPE = Debug```
* Then ```make```

## Execution
Enter a matrix in form of a string, e.g. this could be either

```25 -16 30 -44 -12 13 -7 18 -26 -6 -18 12 -21 36 12 -9 6 -12 21 6 11 -8 15 -22 -3```

or 

```
25 -16 30 -44 -12
13 -7 18 -26 -6
-18 12 -21 36 12
-9 6 -12 21 6
11 -8 15 -22 -3
```

as long the number of integers entered is a perfect square.

## Algorithm
The algorithm first calculates the characteristic polynomial of the given matrix. If the characteristic of $R$ is zero, this is done using Faddeev-LeVerrier, else (as for $𝔽_p$) the naive way as the determinant of the polynomial-valued matrix $$(X\cdot I-A)$$ also works. Factorization is done by brute force (bad).

For each eigenvalue $\lambda\in R$ of $A$, the kernels of the endomorphisms $\psi^k:=(A-\lambda\cdot I)^k$ are calculated. From top to bottom, for each base vector ${v\in\ker\psi^k/\ker\psi^{k-1}}$ (that is not already in the result) the corresponding jordan chain

$$v,\psi(v),\dots,\psi^{k-1}(v)$$

is added to the result. The matrix containing all given jordan chains form a jordan base $S$, such that

$$S^{-1}\cdot A\cdot S$$

is in jordan normal form (with ones above the diagonal). If the characteristic polynomial does not break down completely, return a partial base.

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
\end{bmatrix}\in ℚ^{5\times5}\quad\text{and}\quad 6A^{-1}=\begin{bmatrix}
  0 & 0 & 3  & -6  & 0\\\
  0 & 0 & 0  & -8  & -6\\\
  6 & 0 & 6  & -12 & -12\\\
  0 & 0 & 0  & -11 & -9\\\
  0 & 6 & -3 & 0   &-12
\end{bmatrix}\in ℤ^{5\times5}
```

Using the scaled inverse therefore stays in the ring of integers, even though later results might be scaled accordingly. For obvious reasons this problem does not occur when working in a field.

## Notes
Works with any euclidean ring realized as class, correctly implementing the operators ```+, -, *, /``` and ```%``` (such that ```%``` is the actual rest of division not something like ```fmod```), aswell as ```+=, -=, *=, /=, %=``` a constructor from an integer and a static function ```int characteristic()``` returning the characteristic. Division in a ring might throw an exception if the divisor does not divide the dividend.

Note that this programs is primarily used to get familiar with C++ template and memory management. I am aware that direct (de-)allocation of memory and use of raw pointers is not incentivized.

All datastructures are defined in the .hpp, since template classes cannot be seperated into .hpp and .cpp. The problem comes from the compiler not being able to infer the same template parameter for .hpp and .cpp. This problem could be fixed by always including the corresponding .cpp together with the .hpp.

## Limitations
Does not scale well due to extreme growth of coefficients as is usual with exact algorithms (at least with integer and rational coefficients).

Even though the native type int can be used as a template parameter, under given circumstances (when using nested types as template parameters) this can produce unwanted problems. Use the wrapper class Integer instead (actually wrapping the native type int64_t).

## Root calculus
Defines numbers only by its minimal polynomial e.g. $i$ as the root of $X^2+1$. If 

$$P(X)=0\quad\text{with}\quad P(X)=\sum_{i=0}^na_iX^i\quad\text{then}\quad X^n=-\frac{1}{a_n}\sum_{i=0}^{n-1}a_iX^i.$$

This can be abused to find a polynomial, whose root is e.g. $X+Y$. Since $K[X,Y]/(P,Q)$ is a finite dimensional vector space, there must be a linear dependence between the powers $(X+Y)^k$. Finding this dependence gives us the coefficients of the minimal (see below) polynomial of $X+Y$. To generate the powers of $1/X$, factor the greatest possible power of $X$ from $P$, subtract the constant term and divide by $X$, i.e.

$$\sum_{k=0}^na_kX^k=0\quad\Leftrightarrow\quad X^{\ell}\sum_{k=\ell+1}^na_kX^{k-\ell}+a_{\ell}X^{\ell}=0\quad\Leftrightarrow\quad-\sum_{k=\ell+1}^na_kX^{k-\ell-1}=\frac{a_{\ell}}{X}$$

Since there is no algebraic distinction between the roots of a polynomial, this result is the shared minimal polynomial of ALL possible sums of roots of $P$ and $Q$, which can be undesired. Especially, the given polynomial does not have to be irreducible. 

Let e.g. $P(X)=X^2+1$ and $Q(Y)=Y^2+1$, such that $X$ and $Y$ are both $\pm i$. Since $(\pm i)\cdot(\pm i)=\pm1$ (as $i\cdot(-i)$ is permitted), the calculated polynomial to $XY$ is $Z^2-1=(Z-1)(Z+1)$.
