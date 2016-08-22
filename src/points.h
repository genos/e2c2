/**
 * @file points.h
 * @brief Affine points on Edwards Curves
 *
 * This file contains the interface and implementation (since this base "class"
 * is really a template) of affine and projective points on Edwards Curves.
 */

#ifndef _POINTS_H
#define _POINTS_H

#include <iostream>         // Readable output
#include <NTL/ZZ.h>         // Arbitrarily large integers
#include <NTL/ZZ_pE.h>      // Field elements from @f$ \mathbf{F}_{p^n} @f$
#include <NTL/GF2E.h>       // Field elements from @f$ \mathbf{F}_{2^n} @f$
#include "utilities.h"      // Utilities header for e2c2 project
#include "curves.h"         // e2c2 curve interface and implementation
#include "mol.h"            // Birational map and utilities from MOL paper
#include "utilities.h"      // e2c2 utilities


/// Namespace for our library
namespace e2c2 {
//------- Helper Function -------//

/**
 * @p counterBTTiming(const NTL::ZZ& k, const NTL::ZZ& m)
 * @brief This function counteracts the timing attack outlined in Brumley &
 * Tuveri's paper
 */
inline auto counterBTTiming(const NTL::ZZ& k, const NTL::ZZ& m) -> NTL::ZZ {
    if (NTL::NumBits(k + m) == NTL::NumBits(m)) {
        return k + 2 * m;
    } else
    { return k + m; }
}


//------- Class Skeletons -------//

/**
 * @class Affine
 * @brief Class skeleton for affine points on Edwards Curves
 *
 * Collects relevant information and functionality for affine points on
 * Edwards Curves
 */
template <class Elt, class Curve>
class Affine {
public:
    /// x-coordinate
    Elt x;

    /// y-coordinate
    Elt y;

    /// curve to which this point belongs
    Curve curve;

    /// Default constructor
    Affine() : x(), y(), curve() {}

    /// Destructor
    ~Affine() {}

    /// Constructor given all information
    Affine(const Elt& x, const Elt& y, const Curve& curve) :
        x(x), y(y), curve(curve) {
        if (!curveEquation(curve, x, y)) {
            throw InvalidParametersException();
        }
    }

    /// Constructor if given just a curve (left to specific classes)
    Affine(const Curve& curve) : x(), y(), curve(curve) {
        *this = aff_id(curve);
    }

    /// Checking whether we have the neutral element
    auto isID() const -> bool {
        return *this == Affine(this->curve);
    }

    /// Two affine points are equal if all relevant info is the same...
    auto operator==(const Affine& that) const -> bool {
        return (x == that.x) && (y == that.y) && (curve == that.curve);
    }

    /// ...and are not equal otherwise
    auto operator !=(const Affine& that) const -> bool {
        return !(*this == that);
    }

    /// Assignment by addition; left to fleshed-out classes
    auto operator+=(const Affine& that) -> Affine& {
        if (this->curve != that.curve) {
            throw DifferentCurvesException();
        } else {
            *this = aff_add(*this, that);
            return *this;
        }
    }

    /// Negation of a point; left to full class
    auto operator-() const -> Affine {
        return aff_neg(*this);
    }

    /// Assignment by subtraction; makes use of += and -
    auto operator-=(const Affine& that) -> Affine& {
        return *this += -that;
    }

    /// Addition via +=
    auto operator+(const Affine& that) const -> Affine {
        return Affine(*this) += that;
    }

    /// Subtraction via -=
    auto operator-(const Affine& that) const -> Affine {
        return Affine(*this) -= that;
    }

    /// Point doubling; left to full class
    auto pointDouble() const -> Affine {
        return aff_double(*this);
    }

    /// Montgomery Ladder for scalar multiplication
    auto montgomery(const NTL::ZZ& k) const -> Affine {
        // Work with positive scalars
        if (k < 0) {
            return (-(*this)).montgomery(-k);
        }

        // Counteract Brumley & Tuveri's timing attack
        NTL::ZZ kk = counterBTTiming(k, curve.cardinality());

        Affine aTmp(curve);  // identity element
        Affine bTmp(*this);
        for (auto i = NumBits(kk) - 1; i >= 0; i--) {
            if (bit(k, i)) {
                aTmp += bTmp;
                bTmp = bTmp.pointDouble();
            } else {
                bTmp += aTmp;
                aTmp = aTmp.pointDouble();
            }
        }
        return aTmp;
    }

    /// Assignment by scalar multiplication (using Montgomery Ladder)
    auto operator*=(const NTL::ZZ& k) -> Affine& {
        return *this = montgomery(k);
    }

    /// Assignment by scalar multiplication (using Montgomery Ladder)
    template <class N>
    auto operator*=(const N& k) -> Affine& {
        return *this = montgomery(NTL::to_ZZ(k));
    }

    /// Scalar multiplication (using Montgomery Ladder), via *= (e.g. k *
    /// point)
    friend auto operator*(const NTL::ZZ& k, const Affine& point) ->
    Affine {
        return Affine(point) *= k;
    }

    /// Scalar multiplication (using Montgomery Ladder), via *= (e.g. k *
    /// point)
    template <class N>
    friend auto operator*(const N& k, const Affine& point) -> Affine {
        return Affine(point) *= NTL::to_ZZ(k);
    }

    /// Output
    friend auto operator<<(std::ostream& out, const Affine&
                           point) -> std::ostream& {
        return (out << "(" << point.x << ", " << point.y << ")");
    }
};


/**
 * @class Projective
 * @brief Class skeleton for projective points on Edwards Curves
 *
 * Collects relevant information and functionality for projective points on
 * Edwards Curves
 */
template <class Elt, class Curve>
class Projective {
public:
    /// x-coordinate
    Elt x;

    /// y-coordinate
    Elt y;

    /// z-coordinate
    Elt z;

    /// curve to which this point belongs
    Curve curve;

    /// Default constructor
    Projective() : x(), y(), z(), curve() {}

    /// Destructor
    ~Projective() {}

    /// Constructor given all information
    Projective(const Elt& x, const Elt& y, const Elt& z,
               const Curve& curve) : x(x), y(y), z(z), curve(curve) {
        if (!curveEquation(curve, x / z, y / z)) {
            throw InvalidParametersException();
        }
    }

    /// Constructor if given just a curve (left to specific class)
    Projective(const Curve& curve) : x(), y(), z(), curve(curve) {
        *this = proj_id(curve);
    }

    /// Constructor from an affine point of the same type
    explicit Projective(const Affine<Elt, Curve>& a) : x(a.x), y(a.y), z(),
        curve(a.curve) {
        z = 1;
        if (!curveEquation(curve, x / z, y / z)) {
            throw InvalidParametersException();
        }
    }

    /// Equivalence class representative: z = 1
    auto equivalenceClassRep() const -> Projective {
        Elt one;
        one = 1;
        return Projective(x / z, y / z, one, curve);
    }

    /// Checking whether we have the neutral element
    auto isID() const -> bool {
        return *this == Projective(this->curve);
    }

    /// Two projective points are equal iff all relevant info is the
    /// same...
    auto operator==(const Projective& that) const -> bool {
        return (x / z == that.x / that.z) && (y / z == that.y / that.z) &&
               (curve == that.curve);
    }

    /// ...and are not equal otherwise
    auto operator !=(const Projective& that) const -> bool {
        return !(*this == that);
    }

    /// Assignment by addition; left to fleshed-out class specificities
    auto operator+=(const Projective& that) -> Projective& {
        if (this->curve != that.curve) {
            throw DifferentCurvesException();
        }

        *this = proj_add(*this, that);
        return *this;
    }

    /// Negation of a point; left to class specificities
    auto operator-() const -> Projective {
        return proj_neg(*this);
    }

    /// Assignment by subtraction; makes use of += and -
    auto operator-=(const Projective& that) -> Projective& {
        return *this += -that;
    }

    /// Addition via +=
    auto operator+(const Projective& that) const -> Projective {
        return Projective(*this) += that;
    }

    /// Subtraction via -=
    auto operator-(const Projective& that) const -> Projective {
        return Projective(*this) -= that;
    }

    /// Point doubling; left to class specifics
    auto pointDouble() const -> Projective {
        return proj_double(*this);
    }

    /// Montgomery Ladder for scalar multiplication
    auto montgomery(const NTL::ZZ& k) const -> Projective {
        // Work with positive scalars
        if (k < 0) {
            return (-(*this)).montgomery(-k);
        }

        // Counteract Brumley & Tuveri's timing attack
        NTL::ZZ kk = counterBTTiming(k, curve.cardinality());

        Projective aTmp(curve);  // identity element
        Projective bTmp(*this);
        for (auto i = NumBits(kk) - 1; i >= 0; i--) {
            if (bit(k, i)) {
                aTmp += bTmp;
                bTmp = bTmp.pointDouble();
            } else {
                bTmp += aTmp;
                aTmp = aTmp.pointDouble();
            }
        }
        return aTmp;
    }

    /// Assignment by scalar multiplication (using Montgomery Ladder)
    auto operator*=(const NTL::ZZ& k) -> Projective& {
        return *this = montgomery(k);
    }

    /// Assignment by scalar multiplication (using Montgomery Ladder)
    template <class N>
    auto operator*=(const N& k) -> Projective& {
        return *this = montgomery(NTL::to_ZZ(k));
    }

    /// Scalar multiplication (using Montgomery Ladder), via *= (e.g. k *
    /// point)
    friend auto operator*(const NTL::ZZ& k,
                          const Projective& point) -> Projective {
        return Projective(point) *= k;
    }

    /// Scalar multiplication (using Montgomery Ladder), via *= (e.g. k *
    /// point)
    template <class N>
    friend auto operator*(const N& k, const Projective& point) ->
    Projective {
        return Projective(point) *= NTL::to_ZZ(k);
    }


    /// Output
    friend auto operator<<(std::ostream& out,
                           const Projective& point) -> std::ostream& {
        Projective tmp = point.equivalenceClassRep();
        return (out << "(" << tmp.x << " : " << tmp.y << " : " << tmp.z <<
        ")");
    }

};


//------- Affine Specialization Functions -------//

/**
 * @var OddAff
 * @brief Affine points on an odd curve
 *
 * This typedef builds the four functions necessary to flesh out our Affine
 * template for odd affine points.
 */
using OddAff = Affine<NTL::ZZ_pE, OddCurve>;

/**
 * @p aff_id(const OddCurve& curve)
 * @brief Affine neutral element on odd curve.
 */
inline auto aff_id(const OddCurve& curve) -> OddAff {
    return OddAff(NTL::ZZ_pE::zero(), curve.c, curve);
}

/**
 * @p aff_add(const OddAff& a1, const OddAff& a2)
 * @brief Affine addition for points on an odd curve.
 */
inline auto aff_add(const OddAff& a1, const OddAff& a2) -> OddAff {
    NTL::ZZ_pE w, num_x, num_y, den_x, den_y;
    w = a1.curve.d * a1.x * a2.x * a1.y * a2.y;
    num_x = a1.x * a2.y + a1.y * a2.x;
    num_y = a1.y * a2.y - a1.x * a2.x;
    den_x = a1.curve.c * (1 + w);
    den_y = a1.curve.c * (1 - w);

    return OddAff(num_x / den_x,
    num_y / den_y,
    a1.curve);
}

/**
 * @p aff_neg(const OddAff& a)
 * @brief Negation of an affine point on an odd curve.
 */
inline auto aff_neg(const OddAff& a) -> OddAff {
    return OddAff(-a.x, a.y, a.curve);
}

/**
 * @p aff_double(const OddAff& a)
 * @brief Affine point doubling on odd curve.
 */
inline OddAff aff_double(const OddAff& a) {
    auto xx = NTL::sqr(a.x);
    auto yy = NTL::sqr(a.y);
    auto num_x = 2 * a.x * a.y * a.curve.c;
    auto num_y = (yy - xx) * a.curve.c;
    auto den_x = xx + yy;
    auto den_y = 2 * NTL::sqr(a.curve.c) - (xx + yy);

    return OddAff(num_x / den_x, num_y / den_y, a.curve);
}


/**
 * @var BinaryAff
 * @brief Affine points on a binary curve
 *
 * This typedef builds the four functions necessary to flesh out our Afine
 * template for binary affine points.
 */
using BinaryAff = Affine<NTL::GF2E, BinaryCurve>;

/**
 * @p aff_id(const BinaryCurve& curve)
 * @brief Affine neutral element on binary curve.
 */
inline auto aff_id(const BinaryCurve& curve) -> BinaryAff {
    return BinaryAff(NTL::GF2E::zero(), NTL::GF2E::zero(), curve);
}

/**
 * @p aff_add(const BinaryAff& a1, const BinaryAff& a2)
 * @brief Affine addition for points on a binary curve.
 */
inline auto aff_add(const BinaryAff& a1, const BinaryAff& a2) -> BinaryAff {
    auto w1 = a1.x + a1.y;
    auto w2 = a2.x + a2.y;
    auto a = NTL::sqr(a1.x) + a1.x;
    auto b = NTL::sqr(a1.y) + a1.y;
    auto c = a1.curve.d * w1 * w2;
    auto d = a2.x * a2.y;

    return BinaryAff(
        a1.y + (c + a1.curve.c * (w1 + a2.x) + a * (d + a2.x)) /
        (a1.curve.c + a * w2),
        a1.x + (c + a1.curve.c * (w1 + a2.y) + b * (d + a2.y)) /
        (a1.curve.c + b * w2),
        a1.curve);
}

/**
 * @p aff_neg(const BinaryAff& a)
 * @brief Negation of an affine point on a binary curve.
 */
inline auto aff_neg(const BinaryAff& a) -> BinaryAff {
    return BinaryAff(a.y, a.x, a.curve);
}

/**
 * @p aff_double(const BinaryAff& a)
 * @brief Affine point doubling on binary curve.
 */
inline auto aff_double(const BinaryAff& a) -> BinaryAff {
    auto aa = NTL::sqr(a.x);
    auto b = NTL::sqr(aa);
    auto c = NTL::sqr(a.y);
    auto d = NTL::sqr(c);
    auto f = a.curve.c;
    auto g = (a.curve.d / a.curve.c) * (b + d);
    auto j = aa + c;
    auto k = g + a.curve.d * j;
    auto z = f + j + g;

    return BinaryAff((k + aa + d) / z,
    (k + c + b) / z,
    a.curve);
}

/**
 * @p birMapAff(const NTL::GF2E& u, const NTL::GF2E& v,
 * const BinaryCurve& curve)
 * @brief Birational Map from Weierstrass curve to Binary Edwards curve.
 */
inline auto birMapAff(const NTL::GF2E& u, const NTL::GF2E& v,
                      const NTL::GF2E& a2, const BinaryCurve& curve) -> BinaryAff {
    NTL::GF2E x, y;
    mol_bm_aff(x, y, u, v, NTL::GF2E::degree(), curve.c, curve.d, a2);
    return BinaryAff(x, y, curve);
}


/**
 * @var TwistedAff
 * @brief Affine points on a twisted curve
 *
 * This typedef builds the four functions necessary to flesh out our Affine
 * template for twisted affine points.
 */
using TwistedAff = Affine<NTL::ZZ_pE, TwistedCurve>;

/**
 * @p aff_id(const TwistedCurve& curve)
 * @brief Affine neutral element on twisted curve.
 */
inline auto aff_id(const TwistedCurve& curve) -> TwistedAff {
    return TwistedAff(NTL::ZZ_pE::zero(), NTL::to_ZZ_pE(1), curve);
}

/**
 * @p aff_add(const TwistedAff& a1, const TwistedAff& a2)
 * @brief Affine addition for points on a twisted curve.
 */
inline auto aff_add(const TwistedAff& a1, const TwistedAff& a2) ->
TwistedAff {
    auto w = a1.curve.d * a1.x * a2.x * a1.y * a2.y;
    auto num_x = a1.x * a2.y + a1.y * a2.x;
    auto num_y = a1.y * a2.y - a1.curve.d * a1.x * a2.x;
    auto den_x = 1 + w;
    auto den_y = 1 - w;

    return TwistedAff(num_x / den_x, num_y / den_y, a1.curve);
}

/**
 * @p aff_neg(const TwistedAff& a)
 * @brief Negation of an affine point on a twisted curve.
 */
inline auto aff_neg(const TwistedAff& a) -> TwistedAff {
    return TwistedAff(-a.x, a.y, a.curve);
}

/**
 * @p aff_double(const TwistedAff& a)
 * @brief Affine point doubling on a twisted curve.
 */
inline auto aff_double(const TwistedAff& a) -> TwistedAff {
    auto b = NTL::sqr(a.x + a.y);
    auto c = NTL::sqr(a.x);
    auto d = NTL::sqr(a.y);
    auto e = a.curve.c * c;
    auto f = e + d;
    auto j = f - 2;
    auto z = f * j;

    return TwistedAff(((b - c - d) * j) / z,
    (f * (e - d)) / z,
    a.curve);
}


//------- Projective Specialization Functions -------//

/**
 * @var OddProj
 * @brief Projective points on an odd curve
 *
 * This typedef builds the four functions necessary to flesh out our
 * Projective template for odd affine points.
 */
using OddProj = Projective<NTL::ZZ_pE, OddCurve>;

/**
 * @p proj_id(const OddCurve& curve)
 * @brief Projective neutral element on odd curve.
 */
inline auto proj_id(const OddCurve& curve) -> OddProj {
    return OddProj(NTL::ZZ_pE::zero(),
    curve.c,
    NTL::to_ZZ_pE(1),
    curve);
}

/**
 * @p proj_add(const OddProj& p1, const OddProj& p2)
 * @brief Projective addition for points on an odd curve.
 */
inline auto proj_add(const OddProj& p1, const OddProj& p2) -> OddProj {
    // From Bernstein & Lange, "Faster Addition and Doubling on Elliptic
    // Curves"
    auto a = p1.z * p2.z;
    auto b = NTL::sqr(a);
    auto c = p1.x * p2.x;
    auto d = p1.y * p2.y;
    auto e = p1.curve.d * c * d;
    auto f = b - e;
    auto g = b + e;
    return OddProj(a * f * ((p1.x + p1.y) * (p2.x + p2.y) - c - d),
    a * g * (d - c),
    p1.curve.c * f * g,
    p1.curve);
}

/**
 * @p proj_neg(const OddProj& p)
 * @brief Negation of an projective point on an odd curve.
 */
inline auto proj_neg(const OddProj& p) -> OddProj {
    return OddProj(-p.x, p.y, p.z, p.curve);
}

/**
 * @p proj_double(const OddProj& p)
 * @brief Projective point doubling on odd curve.
 */
inline auto proj_double(const OddProj& p) -> OddProj {
    /// From Bernstein & Lange, "Faster Addition and Doubling on Elliptic
    /// Curves"
    auto b = NTL::sqr(p.x + p.y);
    auto c = NTL::sqr(p.x);
    auto d = NTL::sqr(p.y);
    auto e = c + d;
    auto h = NTL::sqr(p.curve.c * p.z);
    auto j = e - 2 * h;
    return OddProj(p.curve.c * (b - e) * j,
    p.curve.c * e * (c - d),
    e * j,
    p.curve);
}


/**
 * @var BinaryProj
 * @brief Projective points on an binary curve
 *
 * This typedef builds the four functions necessary to flesh out our
 * Projective template for binary affine points.
 */
using BinaryProj = Projective<NTL::GF2E, BinaryCurve>;

/**
 * @p proj_id(const BinaryCurve& curve)
 * @brief Projective neutral element on binary curve.
 */
inline auto proj_id(const BinaryCurve& curve) -> BinaryProj {
    return BinaryProj(NTL::GF2E::zero(),
    NTL::GF2E::zero(),
    NTL::to_GF2E(1),
    curve);
}

/**
 * @p proj_add(const BinaryProj& p1, const BinaryProj& p2)
 * @brief Projective addition for points on a binary curve.
 */
inline auto proj_add(const BinaryProj& p1, const BinaryProj& p2) ->
BinaryProj {
    /// from Bernstein, Lange, and Farashahi, "Binary Edwards Curves"
    auto w1 = p1.x + p1.y;
    auto w2 = p2.x + p2.y;
    auto a = p1.x * (p1.x + p1.z);
    auto b = p1.y * (p1.y + p1.z);
    auto c = p1.z * p2.z;
    auto d = w2 * p2.z;
    auto e = p1.curve.c * NTL::sqr(c);
    auto h = (p1.curve.c * p2.z + p1.curve.d * w2) * w1 * c;
    auto i = p1.curve.c * c * p1.z;
    auto u = e + a * d;
    auto v = e + b * d;
    auto s = u * v;
    return BinaryProj(
        s * p1.y + (h + p2.x * (i + a * (p2.y + p2.z))) * v * p1.z,
        s * p1.x + (h + p2.y * (i + b * (p2.x + p2.z))) * u * p1.z,
        s * p1.z,
        p1.curve);
}

/**
 * @p proj_neg(const BinaryProj& p)
 * @brief Negation of an projective point on a binary curve.
 */
inline auto proj_neg(const BinaryProj& p) -> BinaryProj {
    return BinaryProj(p.y, p.x, p.z, p.curve);
}

/**
 * @p proj_double(const BinaryProj& p)
 * @brief Projective point doubling on binary curve.
 */
inline auto proj_double(const BinaryProj& p) -> BinaryProj {
    /// from Bernstein, Lange, and Farashahi, "Binary Edwards Curves"
    auto a = NTL::sqr(p.x);
    auto b = NTL::sqr(a);
    auto c = NTL::sqr(p.y);
    auto d = NTL::sqr(c);
    auto e = NTL::sqr(p.z);
    auto f = p.curve.c * NTL::sqr(e);
    auto g = (p.curve.d / p.curve.c) * (b + d);
    auto h = a * e;
    auto i = c * e;
    auto j = h + i;
    auto k = g + p.curve.d * j;
    return BinaryProj(k + h + d, k + i + b, f + j + g, p.curve);
}

/**
 * @p birMapProj(const NTL::GF2E& u, const NTL::GF2E& v,
 * const BinaryCurve& curve)
 * @brief Birational Map from Weierstrass curve to Binary Edwards curve.
 */
inline auto birMapProj(const NTL::GF2E& u, const NTL::GF2E& v,
                       const NTL::GF2E& a2, const BinaryCurve& curve) -> BinaryProj {
    NTL::GF2E x, y, z;
    mol_bm_proj(x, y, z, u, v, NTL::GF2E::degree(), curve.c, curve.d,  a2);
    return BinaryProj(x, y, z, curve);
}


/**
 * @var TwistedProj
 * @brief Projective points on an twisted curve
 *
 * This typedef builds the four functions necessary to flesh out our
 * Projective* template for twisted affine points.
 */
using TwistedProj = Projective<NTL::ZZ_pE, TwistedCurve>;

/**
 * @p proj_id(const TwistedCurve& curve)
 * @brief Projective neutral element on twisted curve.
 */
inline auto proj_id(const TwistedCurve& curve) -> TwistedProj {
    return TwistedProj(NTL::ZZ_pE::zero(),
    NTL::to_ZZ_pE(1),
    NTL::to_ZZ_pE(1),
    curve);
}

/**
 * @p proj_add(const TwistedProj& p1, const TwistedProj& p2)
 * @brief Projective addition for points on a twisted curve.
 */
inline auto proj_add(const TwistedProj& p1, const TwistedProj& p2) ->
TwistedProj {
    /// From Bernstein, Birkner, Joye, Lange, Peters, "Twisted Edwards
    /// Curves"
    auto a = p1.z * p2.z;
    auto b = NTL::sqr(a);
    auto c = p1.x * p2.x;
    auto d = p1.y * p2.y;
    auto e = p1.curve.d * c * d;
    auto f = b - e;
    auto g = b + e;
    return TwistedProj(a * f * ((p1.x + p1.y) * (p2.x + p2.y) - c - d),
    a * g * (d - p1.curve.c * c),
    f * g,
    p1.curve);
}

/**
 * @p proj_neg(const TwistedProj& p)
 * @brief Negation of an projective point on a twisted curve.
 */
inline auto proj_neg(const TwistedProj& p) -> TwistedProj {
    return TwistedProj(-p.x, p.y, p.z, p.curve);
}

/**
 * @p proj_double(const TwistedProj& p)
 * @brief Projective point doubling on twisted curve.
 */
inline auto proj_double(const TwistedProj& p) -> TwistedProj {
    /// From Bernstein, Birkner, Joye, Lange, Peters, "Twisted Edwards
    /// Curves"
    auto b = NTL::sqr(p.x + p.y);
    auto c = NTL::sqr(p.x);
    auto d = NTL::sqr(p.y);
    auto e = p.curve.c * c;
    auto f = e + d;
    auto h = NTL::sqr(p.z);
    auto j = f - 2 * h;
    return TwistedProj((b - c - d) * j, f * (e - d), f * j, p.curve);
}
}



#endif  // _POINTS_H
