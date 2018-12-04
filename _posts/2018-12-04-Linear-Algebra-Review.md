---
layout: post
title: Linear Algebra Review
---
The materials here are mostly directly transcribed from Standford CS229 notes. Writing them in $\LaTeX$ to get myself familiar with the formulas and for future reference.

[TOC]

## Matrix Multiplication

Given two matrices $A \in \mathbb{R}^{m \times n}$ and $B \in \mathbb{R}^{n \times p}$, their product is a matrix.
$$
\begin{align}
                   C  &= AB \in \mathbb{R}^{m \times p} \\
\text{where}\space C_{ij} &= \sum_{k=1}^n A_{ik}B_{kj}
\end{align}
$$

### Vector-Vector Products

Given two vectors $x,y \in \mathbb{R}^n$, their ***inner product*** or ***dot product*** is a real number (scalar).
$$
x^Ty \in \mathbb{R} = \begin{bmatrix}x_1 & x_2 & \cdots & x_n\end{bmatrix}
                      \begin{bmatrix}y_1 \\ y_2 \\ \vdots \\ y_n\end{bmatrix}
                    = \sum_{i=1}^nx_iy_i = y^Tx
$$

Given two vector $x \in \mathbb{R}^m, y \in \mathbb{R}^n$, their ***outer product*** is a matrix.
$$
xy^T \in \mathbb{R}^{m \times n}
= \begin{bmatrix}x_1 \\ x_2 \\ \vdots \\ x_m\end{bmatrix}
  \begin{bmatrix}y_1 & y_2 & \cdots & y_n\end{bmatrix}
= \begin{bmatrix}x_1y_1 & x_1y_2 & \cdots & x_1y_n \\
                 x_2y_1 & x_2y_2 & \cdots & x_2y_n \\
                 \vdots & \vdots & \ddots & x_my_n \\
                 x_my_1 & x_my_2 & \cdots & x_my_n
  \end{bmatrix}
$$

A useful example of outer product is to represent a matrix $A \in \mathbb{R}^{m \times n}$ whose columns are all equal to some vector $x \in \mathbb{R}^m$ by using $n$-dimensional all-one vector $\mathbf{1} \in \mathbb{R}^n$ compactly as,
$$
A = \begin{bmatrix}\vert & \vert & & \vert \\
                   x & x & \cdots & x \\
                   \vert & \vert & & \vert
    \end{bmatrix}
  = \begin{bmatrix}x_1 & x_1 & \cdots & x_1 \\
                   x_2 & x_2 & \cdots & x_2 \\
                   \vdots & \vdots & \ddots & \vdots \\
                   x_m & x_m & \cdots & x_m
    \end{bmatrix}
  = \begin{bmatrix}x_1 \\ x_2 \\ \vdots \\ x_m\end{bmatrix}
    \begin{bmatrix}1 & 1 & \cdots & 1\end{bmatrix}
  = x\mathbf{1}^T
$$

### Matrix-Vector Products

Given a matrix $A \in \mathbb{R}^{m \times n}$ and a vector $x \in \mathbb{R}^{n}$, their product is a vector $y=Ax \in \mathbb{R}^m$.

If we write $A$ by rows,
$$
y=Ax=\begin{bmatrix}- & a_1^T & - \\
                    - & a_2^T & - \\
                      & \vdots & \\
                    - & a_m^T & -
     \end{bmatrix}x
 =\begin{bmatrix}a_1^Tx \\
                 a_2^Tx \\
                 \vdots \\
                 a_m^Tx
  \end{bmatrix}
$$

the $i$th entry of $y$ is the inner product of the $i$th row of $A$ and $x$ ($y_i=a_i^Tx$).

If we write $A$ by columns,
$$
y=Ax=\begin{bmatrix}\vert & \vert &        & \vert \\
                    a_1   & a_2   & \cdots & a_n \\
                    \vert & \vert &        & \vert
     \end{bmatrix}
     \begin{bmatrix}x_1 \\ x_2 \\ \vdots \\ x_n\end{bmatrix}
 =\begin{bmatrix}\vert \\ a_1 \\ \vert \end{bmatrix}x_1 +
  \begin{bmatrix}\vert \\ a_2 \\ \vert \end{bmatrix}x_2 + \dots +
  \begin{bmatrix}\vert \\ a_n \\ \vert \end{bmatrix}x_n
$$
we can see $y$ is a ***linear combination*** of the columns of $A$, with the coefficients given by the entries of $x$.

We can also multiply the matrix on the left by a row vector. Given $A \in \mathbb{R}^{m \times n}$, $x \in \mathbb{R}^{m}$ and $y \in \mathbb{R}^n$,
$$
y^T=x^TA=x^T
    \begin{bmatrix}\vert & \vert &        & \vert \\
                   a_1   & a_2   & \cdots & a_n \\
                   \vert & \vert &        & \vert
    \end{bmatrix}
   =\begin{bmatrix}x^Ta_1 & x^Ta_2 & \cdots & x^Ta_n 
    \end{bmatrix}
$$
Writing $A$ by rows,
$$
\begin{align}
y^T=x^TA
  &=\begin{bmatrix}x_1 & x_2 & \cdots & x_m
    \end{bmatrix}
    \begin{bmatrix}- & a_1^T & - \\
                    - & a_2^T & - \\
                      & \vdots & \\
                    - & a_m^T & -
    \end{bmatrix} \\
  &=x_1\begin{bmatrix}- & a_1^T & -\end{bmatrix} +
    x_2\begin{bmatrix}- & a_2^T & -\end{bmatrix} + \dots +
    x_m\begin{bmatrix}- & a_m^T & -\end{bmatrix}
\end{align}
$$
we can see $y^T$ is a linear combination of the rows of $A$, with the coefficients given by the entries of $x$.

### Matrix-Matrix Products

We can view matrix-matrix multiplication in four different ways.

First, the product $C=AB$ can be viewed as a set of vector-vector products. The $(i, j)$th entry of $C$ is equal to the inner product of the $i$th row of $A$ and the $j$th column of $B$.
$$
C=AB
 =\begin{bmatrix}- & a_1^T & - \\
                 - & a_2^T & - \\
                   & \vdots & \\
                 - & a_m^T & -
  \end{bmatrix}
  \begin{bmatrix}\vert & \vert &        & \vert \\
                 b_1   & b_2   & \cdots & b_p \\
                 \vert & \vert &        & \vert
  \end{bmatrix}
 =\begin{bmatrix}a_1^Tb_1 & a_1^Tb_2 & \cdots & a_1^Tb_p \\
                 a_2^Tb_1 & a_2^Tb_2 & \cdots & a_2^Tb_p \\
                 \vdots   & \vdots   & \ddots & \vdots \\
                 a_m^Tb_1 & a_m^Tb_2 & \cdots & a_m^Tb_p 
  \end{bmatrix}
$$
Alternatively, we can view the same multiplication as the sum of all outer products of the $i$th column of $A$ and $i$th row of $B$ (i.e. the sum of $n$ matrices of dimension $m \times p$).
$$
C=AB
 =\begin{bmatrix}\vert & \vert &        & \vert \\
                 a_1   & a_2   & \cdots & a_n \\
                 \vert & \vert &        & \vert
  \end{bmatrix}
  \begin{bmatrix}- & b_1    & - \\
                 - & b_2    & - \\
                   & \vdots & \\
                 - & b_n    & -
  \end{bmatrix}
 =a_1b_1 + a_2b_2 + \dots + a_nb_n
 =\sum_{i=1}^{n} a_ib_i
$$
Second, the matrix-matrix multiplication can also be viewed as a set of matrix-vector products.

If we represent $B$ by columns,
$$
C=AB=A\begin{bmatrix}\vert & \vert &        & \vert \\
                     b_1   & b_2   & \cdots & b_p \\
                     \vert & \vert &        & \vert
      \end{bmatrix}
    =\begin{bmatrix}\vert & \vert &        & \vert \\
                     Ab_1 & Ab_2   & \cdots & Ab_p \\
                     \vert & \vert &        & \vert
     \end{bmatrix}
$$

If we represent $A$ by rows,
$$
C=AB=\begin{bmatrix}- & a_1^T  & - \\
                    - & a_2^T  & - \\
                      & \vdots & \\
                    - & a_m^T  & -
     \end{bmatrix}B
    =\begin{bmatrix}- & a_1^TB  & - \\
                    - & a_2^TB  & - \\
                      & \vdots  & \\
                    - & a_m^TB  & -
     \end{bmatrix}
$$

Properties of matrix multiplication:

* Matrix multiplication is associative: $(AB)C=A(BC)$
* Matrix multiplication is distributive: $A(B+C)=AB+AC$
* Matrix multiplication is, in general, *not* commutative: $AB \ne BA$

## Operations and Properties

### The Identity Matrix and Diagonal Matrices

The ***identity matrix***, denoted $I \in \mathbb{R}^{n \times n}$, is a square matrix with ones on the diagonal and zeros everywhere else.
$$
I_{ij}=\begin{cases}1 & i=j \\
                    0 & i \ne j
       \end{cases}
$$
It has the property that for all $A \in \mathbb{R}^{m \times n}$,
$$
AI = A = IA
$$
A ***diagonal matrix*** is a matrix where all non-diagonal elements are 0, denoted $D=\text{diag}(d_1,d_2,\dots,d_n)$, with
$$
D_{ij}=\begin{cases}d_i & i=j \\
                    0   & i \ne j
       \end{cases}
$$
Clearly, $I=\text{diag}(1,1,\dots\,1)$.

### The Transpose

Given a matrix $A \in \mathbb{R}^{m \times n}$, its ***transpose*** is $A^T \in \mathbb{R}^{n \times m}$ where $(A^T)_{ij}=A_{ji}$ (i.e. the result from "flipping" the rows and columns).

Properties of transpose:

* $(A^T)^T=A$
* $(AB)^T=B^TA^T$
* $(A+B)^T=A^T+B^T$

### Symmetric Matrices

A square matrix $A \in \mathbb{R}^{n \times n}$ is ***symmetric*** if $A=A^T$. It is ***anti-symmetric*** if $A=-A^T$.

It is easy to show the matrix $A+A^T$ is symmetric and $A-A^T$ is anti-symmetric. So any square matrix $A \in \mathbb{R}^{n \times n}$ can be represented as the sum of a symmetric matrix and an anti-symmetric matrix.
$$
A =\frac{1}{2}(A+A^T) + \frac{1}{2}(A-A^T)
$$
The set of all symmetric matrices of size $n$ is commonly denoted as $\mathbb{S}^n$, so that $A \in \mathbb{S}^n$ means $A$ is a symmetric $n \times n$ matrix.

### The Trace

The ***trace*** of a square matrix $A \in \mathbb{R}^{n \times n}$, denoted $\text{tr}(A)$ (or just $\text{tr}A$), is the sum of diagonal elements in the matrix.
$$
\text{tr}A=\sum_{i=1}^{n}A_{ii}
$$
Properties of trace:

* For $A \in \mathbb{R}^{n \times n}$, $\text{tr}A=\text{tr}A^T$
* For $A,B \in \mathbb{R}^{n \times n}$, $\text{tr}(A+B)=\text{tr}A+\text{tr}B $
* For $A \in \mathbb{R}^{n \times n}$, $t \in \mathbb{R}$, $\text{tr}(tA)=t\space\text{tr}A$
* For $A, B$ such that $AB$ is square, $\text{tr}AB=\text{tr}BA$
* For $A, B, C$ such that $ABC$ is square, $\text{tr}ABC=\text{tr}BCA=\text{tr}CAB$, and so on for the product of more matrices.

### Norms

A ***norm*** of a vector $\|x\|$ is informally a measure of the "length" of the vector.

More formally, a norm is any function $f:\mathbb{R}^n \rightarrow \mathbb{R}$ that satisfies four properties:

1. For all $x \in \mathbb{R}^n$, $f(x) \ge 0$ (non-negativity).
2. $f(x)=0$ if and only if $x=0$ (definiteness).
3. For all $x \in \mathbb{R}^n, t \in \mathbb{R}$, $f(tx)=|t|f(x)$ (homogeneity).
4. For all $x,y \in \mathbb{R}^n$, $f(x+y) \le f(x)+f(y)$ (triangle inequality).

The family of $\ell_p$ norms parameterized by a real number $p \ge 1$, is defined as,
$$
\|x\|_p = \left( \sum_{i=1}^n |x_i|^p \right)^{1/p}
$$

The $\ell_1$ norm,
$$
\|x\|_1 = \sum_{i=1}^n |x_i|
$$
The $\ell_2$ norm (Euclidean norm),
$$
\|x\|_2 = \sqrt{\sum_{i=1}^n x_i^2}
$$
The $\ell_\infty$ norm,
$$
\|x\|_\infty = \text{max}_i|x_i|
$$
Norms can also be defined for matrices, such as the Frobenius norm,
$$
\|A\|_F = \sqrt{ \sum_{i=1}^{m}\sum_{j=1}^{n}A_{ij}^2 }
        = \sqrt{\text{tr}(A^TA)}
$$

### Linear Independence and Rank

A set of vectors $\{x_1, x_2, \dots, x_n\} \subset \mathbb{R}^m$ is said to be ***(linearly) independent*** if no vector can be represented as a linear combination of the remaining vectors.

The ***column rank*** of a matrix is the size of the largest subset of columns of $A$ that constitute a linearly independent set. The ***row rank*** is the largest number of rows of $A$ that constitute a linearly independent set. For any matrix $A \in \mathbb{R}^{m \times n}$, the column rank is equal to the row rank, so they are collectively referred to as the ***rank*** of $A$, denoted as $\text{rank}(A)$.

Properties of the rank:

* For $A \in \mathbb{R}^{m \times n}$, $\text{rank}(A) \le \text{min}(m,n)$. If $\text{rank}(A) = \text{min}(m,n)$ then $A$ is said to be ***full rank***.
* For $A \in \mathbb{R}^{m \times n}$, $\text{rank}(A)=\text{rank}(A^T)$
* For $A \in \mathbb{R}^{m \times n}, B \in \mathbb{R}^{n \times p}$, $\text{rank}(AB) \le \text{min}(\text{rank}(A),\text{rank}(B))$
* For $A, B \in \mathbb{R}^{m \times n}$, $\text{rank}(A+B) \le \text{rank}(A)+\text{rank}(B)$

### The Inverse

The ***inverse*** of a square matrix $A \in \mathbb{R}^{n \times n}$, denoted by $A^{-1}$, is the unique matrix such that
$$
A^{-1}A = I = AA^{-1}
$$
Not all matrices have inverses. If $A^{-1}$ exists, we say that $A$ is ***invertible*** or ***non-singular***.

In order for a square matrix $A$ to have an inverse $A^{−1}$, then $A$ must be full rank.

Properties of the inverse:

* $(A^{-1})^{-1}=A$
* $(AB)^{-1}=B^{-1}A^{-1}$
* $(A^{-1})^T=(A^T)^{-1}$. For this reason this matrix is often denoted $A^{−T}$.

### Orthogonal Matrices

Two vectors $x, y \in \mathbb{R}^n$ are ***orthogonal*** if $x^Ty=0$. A vector $x \in \mathbb{R}^n$ is ***normalized*** if $\|x\|_2=1$.

A square matrix $U \in \mathbb{R}^{n \times n}$ is ***orthogonal*** if all its columns are orthogonal to each other and are normalized (the columns are then referred to as being ***orthonormal***).

It follows immediately from the definition of orthogonality and normality that
$$
U^TU = I = UU^T
$$
In other words, the inverse of an orthogonal matrix is its transpose ($U^{-1}=U^T$).

A nice property of orthogonal matrices is that operating on a vector with an orthogonal matrix will not change its Euclidean norm, i.e.,
$$
\|Ux\|_2 = \|x\|_2
$$

### Range and Nullspace of a Matrix

The ***span*** of a set of vectors $\{x_1, x_2, \dots, x_n\}$ is the set of all vectors that can be expressed as a linear combination of $\{x_1, x_2, \dots, x_n\}$.
$$
\text{span}(\{x_1, \dots, x_n\}) = \left\{v:v=\sum_{i=1}^{n}\alpha_ix_i,\space\alpha_i \in \mathbb{R}\right\}
$$

The ***projection*** of a vector $y \in \mathbb{R}^m$ onto the span of $\{x_1, \dots, x_n\}$ is the vector $v \in \text{span}(\{x_1, \dots, x_n\})$, such that $v$ is as close as possible to $y$, as measured by the Euclidean norm $\|v-y\|_2$.
$$
\text{Proj}(y;\{x_1,\dots,x_n\})
    =\text{argmin}_{v \in \text{span}(\{x_1,\dots,x_n\})} \|y-v\|_2
$$
The ***range*** (also called the columnspace) of a matrix $A \in \mathbb{R}^{m \times n}$, denoted $\mathcal{R}(A)$, is the span of the columns of $A$.
$$
\mathcal{R}(A)=\{ v \in \mathbb{R}^m : v=Ax, x \in \mathbb{R}^n \}
$$
With some assumptions (namely that $A$ is full rank and $n \lt m$), the projection of a vector $y \in \mathbb{R}^m$ onto the range of $A$ is given by,
$$
\text{Proj}(y;A)=\text{argmin}_{v \in \mathcal{R}(A)} \|v-y\|_2
                =A(A^TA)^{-1}A^Ty
$$
This equation is almost the same formula as the least squares problem. In fact we are minimizing the same objective (except for the squaring of the norm which does not affect the optimal point).

When $A$ contains only a single column, $a \in \mathbb{R}^m$, this gives the special case for a projection of a vector onto a line:
$$
\text{Proj}(y;a)=\frac{aa^T}{a^Ta} y
$$
The ***nullspace*** of a matrix $A \in \mathbb{R}^{m \times n}$, denoted $\mathcal{N}(A)$, is the set of all vectors that equal 0 when multiplied by $A$.
$$
\mathcal{N}(A)=\{ x \in \mathbb{R}^n : Ax=0 \}
$$
It turns out that
$$
\{ w:w=u+v,u \in \mathcal{R}(A^T),v \in \mathcal{N}(A) \}=\mathbb{R}^n
\text{ and } \mathcal{R}(A^T) \cap \mathcal{N}(A)=\{\mathbf{0}\}
$$
In other words, $\mathcal{R}(A^T)$ and $\mathcal{N}(A)$ are disjoint subsets that together span the entire space of $\R^n$. Sets of this type are called ***orthogonal complements*** and denoted by $\mathcal{R}(A^T)=\mathcal{N}(A)^\bot$.

### The Determinant

The ***determinant*** of a square matrix $A \in \R^{n \times n}$ is a function $\det : \R^{n \times n} \rightarrow \R$, and is denoted $|A|$ or $\det{A}$.

We start by looking at the geometric interpretation of the determinant. Consider the set of points $S \subset \R^n$ formed by taking all possible linear combinations of the row vectors $a_1, \dots, a_n \in \R^n$ of $A$, where the coefficients of the linear combination are all between 0 and 1.
$$
S=\{ v \in \R^n : v=\sum_{i=1}^{n}\alpha_i a_i \text{ where } 0 \le \alpha_i \le 1, i=1,\dots,n \}
$$
It turns out the absolute value of the determinant of $A$ is a measure of the "volume" of the set $S$.

For two-dimensional matrices, $S$ generally has the shape of a *parallelogram*. In three dimensions, $S$ corresponds to the object *parallelepiped* (a three-dimensional box with skewed sides, such that every face has the shape of a parallelogram). In even higher dimensions, $S$ is an object known as an $n$-dimensional *parallelotope*.

Algebraically, the determinant satisfies the following three properties:

1. The determinant of the identity is 1, $|I| = 1$. (Geometrically, the volume of a unit hypercube is 1).
2. Given a matrix $A \in \R^{n×n}$, if we multiply a single row in $A$ by a scalar $t \in \R$, then the determinant of the new matrix is $t|A|$. (Geometrically, multiplying one of the sides of the set $S$ by a factor $t$ causes the volume to increase by a factor $t$.)
3. If we exchange any two rows $a_i^T$ and $a_j^T$ of $A$, then the determinant of the new matrix is $-|A|$.

Such a function satisfying the above conditions exists and is unique.

Properties of the determinant:

* For $A \in \R^{n \times n}$, $|A|=|A^T|$.
* For $A, B \in \R^{n \times n}$, $|AB|=|A||B|$.
* For $A \in \R^{n \times n}$, $|A|=0$ if and only if $A$ is singular (i.e. non-invertible). (If $A$ is singular, then it does not have full rank, and hence its columns are linearly dependent. In this case, the set $S$ corresponds to a "flat sheet" within the $n$-dimensional space and hence has zero volume.)
* For $A \in \R^{n \times n}​$ and $A​$ non-singular, $|A^{-1}|=\dfrac{1}{|A|}​$.

First we define, for $A \in \R^{n \times n}, A_{\setminus i,\setminus j} \in \R^{(n-1) \times (n-1)}$ to be the *matrix* that results from deleting the $i$th row and $j$th column from $A$.

The general (recursive) formula for the determinant is
$$
\begin{align}
|A|&=\sum_{i=1}^n (-1)^{i+j} a_{ij}|A_{\setminus i,\setminus j}| & (\text{for any j} \in 1,\dots,n) \\
   &=\sum_{j=1}^n (-1)^{i+j} a_{ij}|A_{\setminus i,\setminus j}| & (\text{for any i} \in 1,\dots,n)
\end{align}
$$
with the initial case that $|A|=a_{11}$ for $A \in \R^{1 \times 1}$. The expanded formula for $A \in \R^{n \times n}$ will have $n!$ terms.
$$
\begin{align}
|\begin{bmatrix}a_{11}\end{bmatrix}| &= a_{11} \\
\left|\begin{bmatrix}a_{11} & a_{12} \\
                     a_{21} & a_{22}
\end{bmatrix}\right| &= a_{11}a_{22}-a_{12}a_{21} \\
\left|\begin{bmatrix}a_{11} & a_{12} & a_{13} \\
                     a_{21} & a_{22} & a_{23} \\
                     a_{31} & a_{32} & a_{33}
\end{bmatrix}\right| &= \begin{split}
                        a_{11}&a_{22}a_{33}+a_{12}a_{23}a_{31}+a_{13}a_{21}a_{32} \\
                        &-a_{11}a_{23}a_{32}-a_{12}a_{21}a_{33}-a_{13}a_{22}a_{31}
                        \end{split}
\end{align}
$$
The ***classical adjoint*** (often just called the adjoint) of a matrix $A \in \R^{n \times n}$, denoted $\text{adj}(A)$, is defined as,
$$
\begin{align}
&\text{adj}(A) \in \R^{n \times n}, & (\text{adj}(A))_{ij}=(-1)^{i+j}|A_{\setminus j,\setminus i}|
\end{align}
$$
(note the switch in the indices $A_{\setminus j,\setminus i}$). It can be shown that for any non-singular $A \in \R^{n \times n}$,
$$
A^{-1}=\frac{1}{|A|}\text{adj}(A)
$$

### Quadratic Forms and Positive Semidefinite Matrices

Given a square matrix $A \in \R^{n \times n}$ and a vector $x \in \R^n$, the scalar value $x^TAx$ is called a ***quadratic form***.
$$
x^TAx = \sum_{i=1}^n x_i(Ax)_i = \sum_{i=1}^n x_i \left(\sum_{j=1}^n A_{ij}x_j\right) = \sum_{i=1}^n\sum_{j=1}^n A_{ij}x_ix_j
$$
Note that,
$$
x^TAx = (x^TAx)^T = x^TA^Tx = x^T\left(\frac{1}{2}A+\frac{1}{2}A^T\right)x
$$
where the first equality uses the fact that the transpose of a scalar is itself, and the last equality follows from the fact that we are just averaging two equal quantities. We can see that only the symmetric part of $A$ contributes to the quadratic form. For this reason, we often implicitly assume that the matrices appearing in a quadratic form are symmetric.

* A symmetric matrix $A \in \mathbb{S}^n$ is ***positive definite*** (PD) if $x^TAx>0 \ \forall x \in \R^n, x \ne \mathbf{0}$. This is usually denoted $A \succ 0$ (or just $A>0$), and the set of all positive definite matrices is often denoted $\mathbb{S}^n_{++}$.
* A symmetric matrix $A \in \mathbb{S}^n$ is ***positive semidefinite*** (PSD) if $x^TAx \ge 0 \ \forall x \in \R^n$. This is usually denoted $A \succeq 0$ (or just $A \ge 0$), and the set of all positive semidefinite matrices is often denoted $\mathbb{S}^n_{+}$.
* A symmetric matrix $A \in \mathbb{S}^n$ is ***negative definite*** (ND) if $x^TAx<0 \ \forall x \in \R^n, x \ne \mathbf{0}$. This is usually denoted $A \prec 0$ (or just $A<0$).
* A symmetric matrix $A \in \mathbb{S}^n$ is ***negative semidefinite*** (NSD) if $x^TAx \le 0 \ \forall x \in \R^n$. This is usually denoted $A \preceq 0$ (or just $A \le 0$).
* A symmetric matrix $A \in \mathbb{S}^n$ is ***indefinite***, if it is neither positive semidefinite nor negative semidefinite.

It should be obvious that if $A$ is PD, then $-A$ is ND and vice versa. If $A$ is PSD, then $-A$ is NSD and vice versa. If $A$ is indefinite, then so is $-A$.

One important property of PD and ND matrices is that they are always full rank, and hence, invertible.

Finally, there is one type of PD matrix that worth mentioning. Given any matrix $A \in \R^{m \times n}$ (not necessarily symmetric or even square), the matrix $G=A^TA$ (sometimes called a ***Gram matrix***) is always positive semidefinite. Further, if $m \ge n$ (and assume $A$ is full rank), then $G=A^TA$ is positive definite.

### Eigenvalues and Eigenvectors

Given a square matrix $A \in \R^{n \times n}$, we say that $\lambda \in \C$ is an ***eigenvalue*** of $A$ and $x \in C^n$ is the corresponding ***eigenvector*** if
$$
\begin{align}
&Ax = \lambda x, &x \ne 0
\end{align}
$$
As for any eigenvector $x \in \C^n$ and scalar $t \in \C$, $A(cx)=cAx=c\lambda x=\lambda(cx)$, so $cx$ is also an eigenvector. So when we talk about "the" eigenvector, we usually refer to the normalized eigenvector (i.e. length 1).

Rewriting the equation, we get $(\lambda I - A)x=0,\space x \ne 0$. This has a non-zero solution if and only if $(\lambda I-A)$ has a non-empty nullspace, which means $(\lambda I-A)$ is singular, i.e.
$$
|(\lambda I-A)|=0
$$
This would expand to an $n$-degree polynomial in $\lambda$ and we can solve for the $n$ (possibly complex) roots for the eigenvalues. Then the corresponding eigenvectors can be found by solving the linear equation $(\lambda_i I - A)x=0$. In practice, we will not compute the eigenvalues and eigenvectors in this way as more efficient method exists.

Properties of eigenvalues and eigenvectors:

* The trace of $A$ is equal to the sum of its eigenvalues,

$$
\text{tr}A = \sum_{i=1}^n \lambda_i
$$

* The determinant of $A$ is equal to the product of its eigenvalues,

$$
|A| = \prod_{i=1}^n \lambda_i
$$

* The rank of $A$ is equal to the number of non-zero eigenvalues of $A$.

* If $A$ is non-singular, then $1/\lambda_i$ is an eigenvalue of $A^{-1}$ with associated eigenvector $x_i$, i.e. $A^{-1}x_i=(1/\lambda_i)x_i$.

* The eigenvalues of a diagonal matrix $D=\text{diag}(d_1, \dots, d_n)$ are just the diagonal entries $d_1, \dots, d_n$.

We can write all the eigenvector equations simultaneously as
$$
AX = X\Lambda
$$
where the columns of $X \in \R^{n \times n}$ are the eigenvectors of $A$ and $\Lambda$ is a diagonal matrix whose entries are the eigenvalues of $A$, i.e.
$$
\begin{align}
&X \in \R^{n \times n}=\begin{bmatrix}
                      |   & |   &        & | \\
                      x_1 & x_2 & \cdots & x_n \\
                      |   & |   &        & |
                      \end{bmatrix},
                      &\Lambda = \text{diag}(\lambda_1,\dots,\lambda_n)
\end{align}
$$
If the eigenvectors of $A$ are linearly independent, then $X$ will be invertible, so $A=X \Lambda X^{-1}$. A matrix that can be written in this form is called ***diagonalizable***. 

### Eigenvalues and Eigenvectors of Symmetric Matrices

For a symmetric matrix $A \in \mathbb{S}^n$,

* all the eigenvalues of $A$ are real.
* the eigenvectors of $A$ are orthonormal (i.e. the matrix $X$ defined above is an orthogonal matrix).

So $A \in \mathbb{S}^n=U \Lambda U^{-1}=U \Lambda U^T$.

The quadratic form is then,
$$
x^TAx = x^TU \Lambda U^Tx = y^T \Lambda y = \sum_{i=1}^n \lambda_iy_i^2
$$
where $y=U^Tx$ (and since $U$ is full rank, so any vector $y \in \R^n$ can be represented in this form).

The sign of the expression depends entirely on the $\lambda_i$'s. If all $\lambda_i>0$, then the matrix is PD; If all $\lambda_i \ge 0$, then the matrix is PSD. Likewise, If all $\lambda_i<0$ or $\lambda_i \le 0$, then the matrix is ND or NSD. Finally, if $A$ has both positive and negative eigenvalues, then it is indefinite.

Eigenvalues and eigenvectors are often used in maximizing or minimizing some function of a matrix. Given a matrix $A \in \mathbb{S}^n$ and assuming its eigenvalues are ordered as $\lambda_1 \ge \lambda_2 \ge \dots \ge \lambda_n$.

Consider the maximization problem,
$$
\text{max}_{x \in \R^n}x^TAx \quad\text{subject to } \Vert x\Vert_2^2=1
$$

the optimal $x$ is $x_1$ (the eigenvector corresponding to $\lambda_1$), and the maximal value of the form is $\lambda_1$.

Similarly, for the minimization problem,
$$
\text{min}_{x \in \R^n}x^TAx \quad\text{subject to } \Vert x\Vert_2^2=1
$$
the optimal solution is $x_n$ (the eigenvector corresponding to $\lambda_n$), and the minimal value is $\lambda_n$.

## Matrix Calculus

### The Gradient

Suppose that $f:\R^{m \times n} \rightarrow \R$ is a function that takes as input a matrix $A$ of size $m \times n$ and returns a real value. The the ***gradient*** of $f$ (with respect to $A \in \R^{m \times n}$) is the matrix of partial derivatives, defined as,
$$
\nabla_Af(A) \in \R^{m \times n} =
\begin{bmatrix}
\frac{\partial f(A)}{\partial A_{11}} & \frac{\partial f(A)}{\partial A_{12}} & \cdots & \frac{\partial f(A)}{\partial A_{1n}} \\
\frac{\partial f(A)}{\partial A_{21}} & \frac{\partial f(A)}{\partial A_{22}} & \cdots & \frac{\partial f(A)}{\partial A_{2n}} \\
\vdots & \vdots & \ddots & \vdots \\
\frac{\partial f(A)}{\partial A_{m1}} & \frac{\partial f(A)}{\partial A_{m2}} & \cdots & \frac{\partial f(A)}{\partial A_{mn}}
\end{bmatrix}
$$
i.e., an $m \times n$ matrix with
$$
(\nabla_Af(A))_{ij} = \frac{\partial f(A)}{\partial A_{ij}}
$$
Note that the size of $\nabla_Af(A)$ is always the same as the size of $A$.

It is very important to remember that the gradient of a function is *only* defined if the function is real-valued. For example, we *cannot* take the gradient of $f(x)=Ax$ with respect to $x$, since this quantity is vector-valued.

* $\nabla_x(f(x)+g(x)) = \nabla_x f(x) + \nabla_x g(x)$
* For $t \in \R$, $\nabla_x(t\space f(x)) = t \nabla_x f(x)$

We need to explicitly define the variable which we are differentiating with respect to. In the case we are differentiating the function $f$ w.r.t. its argument $z$ and evaluating the gradient at point $Ax$, we denote as $\nabla_z f(Ax)$ (often just written as $\nabla f(Ax)$). In the other case we are differentiating the composite function $g(x)=f(Ax)$ w.r.t. $x$ directly, we denote as $\nabla_x f(Ax)$.

### The Hessian

Suppose that $f:\R^n \rightarrow \R$ is a function that takes a vector in $\R^n$ and returns a real number. The the ***Hessian*** matrix with respect to $x$, written $\nabla_x^2 f(x)$ or simply as $H$ is the $n \times n$ matrix of partial derivatives,
$$
\nabla_x^2 f(x) \in \R^{n \times n} =
\begin{bmatrix}
\frac{\partial^2 f(x)}{\partial x_1^2} & \frac{\partial^2 f(x)}{\partial x_1 \partial x_2} & \cdots & \frac{\partial^2 f(x)}{\partial x_1 \partial x_n} \\
\frac{\partial^2 f(x)}{\partial x_2 \partial x_1} & \frac{\partial^2 f(x)}{\partial x_2^2} & \cdots & \frac{\partial^2 f(x)}{\partial x_2 \partial x_n} \\
\vdots & \vdots & \ddots & \vdots \\
\frac{\partial^2 f(x)}{\partial x_n \partial x_1} & \frac{\partial^2 f(x)}{\partial x_n \partial x_2} & \cdots & \frac{\partial^2 f(x)}{\partial x_n^2}
\end{bmatrix}
$$
i.e., an $n \times n$ matrix with
$$
(\nabla_x^2 f(x))_{ij} = \frac{\partial^2 f(x)}{\partial x_i \partial x_j}
$$
The Hessian is always symmetric, since
$$
\frac{\partial^2 f(x)}{\partial x_i \partial x_j} =
\frac{\partial^2 f(x)}{\partial x_j \partial x_i}
$$
Note that similar to the gradient, the Hessian is defined *only* when $f(x)$ is real-valued.

While we can take the gradient w.r.t. a matrix $A \in \R^n$, we will only consider taking the Hessian w.r.t. a vector $x \in \R^n$.

### Gradients and Hessians of Quadratic and Linear Functions

* $\nabla_x b^Tx = b$
* $\nabla_x x^TAx = 2Ax$ (if $A$ symmetric)
* $\nabla_x^2 x^TAx = 2A$ (if $A$ symmetric)

### Least Squares

Given a matrix $A \in \R^{m \times n}$ (assume $A$ is full rank) and a vector $b \in \R^m$ such that $b \notin \mathcal{R}(A)$. In this situation we will not be able to find a vector $x \in \R^n$ such that $Ax=b$. Instead we want to find a vector $x$ such that $Ax$ is as close as possible to $b$, as measured by the square of the Euclidean norm.
$$
\begin{align}
\Vert Ax-b\Vert_2^2 &= (Ax-b)^T(Ax-b) \\
                    &= x^TA^TAx-2b^TAx+b^Tb
\end{align}
$$
Taking the gradient w.r.t $x$ and setting it to zero.
$$
\begin{align}
\nabla_x(x^TA^TAx-2b^TAx+b^Tb) &= 0 \\
                  2A^TAx-2A^Tb &= 0 \\
                         A^TAx &= A^Tb \\
                             x &= (A^TA)^{-1}A^Tb
\end{align}
$$

### Gradients of the Determinant

Given $A \in \R^{n \times n}$, we want to find the gradient of the determinant w.r.t. matrix $A$, i.e. $\nabla_A |A|$.

Recall that
$$
|A|=\sum_{i=1}^n (-1)^{i+j}A_{ij}|A_{\setminus i,\setminus j}| \quad\text{(for any } j \in 1,\dots,n)
$$
so
$$
\frac{\partial}{\partial A_{kl}}|A| =
\frac{\partial}{\partial A_{kl}} \sum_{i=1}^n (-1)^{i+j}A_{ij}|A_{\setminus i,\setminus j}| =
(-1)^{k+l} |A_{\setminus k,\setminus l}|
=(\text{adj}(A))_{lk}
$$
Therefore
$$
\nabla_A|A| = (\text{adj}(A))^T = (|A|A^{-1})^T = |A|A^{-T}
$$
Next consider the function $f: \mathbb{S}^n_{++} \rightarrow \R$, $f(A)=\log|A|$. Note we have to restrict the domain of $f$ to be the positive definite matrices, which ensures $|A|>0$ and hence $\log|A| \in \R$.

Using the chain rule, we can see
$$
\frac{\partial\log|A|}{\partial A_{ij}} =
\frac{\partial\log|A|}{\partial|A|} \frac{\partial|A|}{\partial A_{ij}} =
\frac{1}{|A|} \frac{\partial|A|}{\partial A_{ij}}
$$
so
$$
\begin{align}
\nabla_A\log|A| &= \frac{1}{|A|}\nabla_A|A| \\
                & = A^{-T} = A^{-1} \text{ (transpose dropped as }A\text{ is symmetric)}
\end{align}
$$

### Eigenvalues as Optimzation

Consider the following equality constrained optimization problem again.
$$
\text{max}_{x \in \R^n} x^TAx \quad\text{subject to }\Vert x\Vert_2^2=1 \text{, for }A\in \mathbb{S}^n
$$
A standard way of solving optimization problems with equality constraints is by forming the objective function ***Lagrangian*** (where $\lambda$ is the Lagrange multiplier).
$$
\mathcal{L}(x,\lambda) = x^TAx-\lambda x^Tx
$$
It can be established that for $x^*$ to be optimal, the gradient of the Lagrangian has to be zero at $x^*$ (not the only condition, but it is required).
$$
\begin{align}
\nabla_x\mathcal{L}(x,\lambda) &= \nabla_x(x^TAx-\lambda x^Tx) \\
                               &= 2Ax-2\lambda x = 0 \\
                               &\Rightarrow Ax=\lambda x \quad\text{(the eigenvalues linear equation)}
\end{align}
$$
This shows that the only points which can possibly maximize (or minimize) $x^TAx$ given the constraint $x^Tx=1$ are the eigenvectors of $A$.
