# e2c2: Edwards Elliptic Curve Cryptography

[![Build Status](https://travis-ci.org/genos/e2c2.svg?branch=master)](https://travis-ci.org/genos/e2c2)

The group of rational points on an elliptic curve over a finite fields has
proven very useful in cryptography since Koblitz and Miller first suggested its
use in the 1980s. Due to the lack of subexponential algorithms for the Discrete
Logarithm Problem in this group, elliptic curve cryptography enjoys a level of
security comparable to other ElGamal-type systems with much smaller key sizes.

There is some room for improvement, however. Typically the rules for the group
operation on an elliptic curve involve a number of special cases:

* What if one point is the point at infinity?
* What if the two points are the same?
* What if they're inverses of each other?

In each of these cases, the need to handle an exception to the typical
geometric understanding can lead to implementations giving off more information
than intended---leaking information through a "side-channel."

In 2007, Dr. Harold Edwards put forth a new form of elliptic curve; despite his
paper not focusing on cryptography, his new normal form has very desirable
cryptographic properties: the addition law is *unified* and *complete*. In
other words, Edwards curves do not leak as much side-channel information as
curves in typical Weierstrass (or other) forms.  Moreover, in many cases the
addition laws involve less operations, making for faster computations.  While
this is not the case over binary fields, the benefits of the law's completeness
make the loss of speed seem negligible. In fact, some authors argue that with
specialized hardware the speed difference can be greatly reduced, while the
completeness of the binary Edwards curve group law actually makes it faster
than Weierstrass implementations that must constantly check for special cases.
Add to this the reduced code complexity, and binary Edwards curves look much
more promising from an implementation point of view.

This code library consists the second version of a C++ proof-of-concept
implementation of Edwards Curves (over both binary fields and fields of odd
prime characteristic) and both affine and projective points over them, built
using Victor Shoup's NTL.


## Acknowledgements
Mike Blackmon was a great help in getting this code working.
