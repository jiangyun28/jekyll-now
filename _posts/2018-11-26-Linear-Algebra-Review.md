---
layout: post
title: Linear Algebra Review
---
# Linear Algebra Review

My own notes using $\LaTeX$ for reference.

[TOC]

## Matrix Multiplication

Given two matrices $A \in \mathbb{R}^{m \times n}$ and $B \in \mathbb{R}^{n \times p}$, their product is a matrix
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
* For $A \in \R^{n \times n}$, $|A|=0$ if and only if $A$ is singular (i.e. non-invertible).
* For $A \in \R^{n \times n}​$ and $A​$ non-singular, $|A^{-1}|=\dfrac{1}{|A|}​$.

ff