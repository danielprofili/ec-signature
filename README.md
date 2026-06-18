# ec-signature

## Introduction

A set of Maple functions which implement the algorithm described in
[1]; namely: given a pair of real symmetric matrices with symbolic
entries and a desired arrangement of their eigenvalues, compute a
quantifier-free condition on the entries of the matrices so that their
eigenvalues are arranged in the given way.

This implementation uses the theory of generalized real root counting
and the signature of matrices. See [1] for more details.

## Usage
The main function is `ECsig`.

```
# In : F, G, symbolic real symmetric matrices.
# Out: Csig, Asig such that c = EC(F,G) iff c = Csig . sig(Asig)
```

## Example
Let
```
  F = [ 1  -2  ]       G =  [ 1        p^2+p  p + 1 ]
      [ -2  1  ] ,          [ p^2 + p  -p     p - 1 ]
                            [ p + 1    p - 1  p     ]
```

where `p` is a parameter. Suppose we would like to find a
quantifier-free condition
on `p` so that the eigenvalue configuration of `F` and `G` is
`(1,1)`; i.e., the eigenvalues of `F` (red) and `G` (blue) are
arranged as follows:

![EC of 011](http://danielprofili.github.io/resources/ec11.png)

Make sure that Maple's working directory is the same directory as the
`ECSignature.mpl` script. In a Maple worksheet, enter the following.

```
restart:
with(LinearAlgebra):
read(`ECSignature.mpl`);
F := Matrix([[1,-2],[-2,1]]);
G := Matrix([[1, p^2 + p, p + 1], [p^2 + p, -p, p-1], [-p + 1, p - 1, p]]);
Csig, Asig := ECSig(F,G);
```

![Maple output](http://danielprofili.github.io/resources/sig-output.png)

The outputs are a numeric matrix `Csig` and a column vector `Asig` of
polynomials in `p`, such that if `c` is any valid eigenvalue
configuration, then 

```
c = Csig . sig(Asig),
```

where `sig(Asig)` means to take the signature (difference between
number of positive roots and negative roots) of each entry in `Asig`.

## Reference
1. Hong, Hoon, **Daniel Profili**, and J. Rafael Sendra. ["Conditions
   for eigenvalue configurations of two real symmetric matrices
   (signature approach)."](https://arxiv.org/pdf/2401.00866).
    SIAM Journal on Matrix Analysis and Applications (pending). (2026).

