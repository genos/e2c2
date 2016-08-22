/**
 * @file curves.h
 * @brief Edwards Curves over finite fields of prime characteristic
 * @author Graham Enos
 *
 * This file contains the interface and implementation (since this base "class"
 * is really a template) of Edwards Curves over finite fields of prime
 * characteristic, viz @f$ \mathbf{F}_{p^n} @f$ in C++
 */

#ifndef _CURVES_H
#define _CURVES_H

#include <iostream>         // Readable output
#include <NTL/ZZ.h>         // Arbitrarily large integers
#include <NTL/ZZ_pE.h>      // Field elements from @f$ \mathbf{F}_{p^n} @f$
#include <NTL/GF2E.h>       // Field elements from @f$ \mathbf{F}_{2^n} @f$
#include "utilities.h"      // Utilities header for e2c2 project


/// Namespace for our library
namespace e2c2 {

//------- Class Skeletons -------//

/**
 * @p CurveID
 * @brief quick typesafe identifiers for curves
 *
 * @details Used in curve definitions, made an enum class for typesafety
 * and speed; used to differentiate curves that are defined similarly
 */
enum class CurveID {
    Odd,
    Binary,
    Twisted
};

/**
 * @p Curve
 * @brief Base "class" (really class template) for Edwards Curves
 *
 * @tparam Elt          field element
 * @tparam ID           id of curve
 *
 * @details Collects all relevant information about the curve and provides
 * basic functionality
 */
template <class Elt, CurveID ID>
class Curve {
public:
    /// Parameter #1 of curve; @f$ c @f$ in papers on odd curves,
    /// @f$ d_1 @f$ in binary papers, @f$ a @f$ in twisted papers,
    Elt c;

    /// Parameter #2 of curve; @f$ d @f$ in papers on odd and twisted
    /// curves, @f$ d_2 @f$ in binary papers
    Elt d;

    /// Cardinality of curve
    NTL::ZZ m;

    /// Default constructor
    Curve() : c(), d(), m() {}

    /// Destructor
    ~Curve() {}

    /// Constructor given full tuple as input
    Curve(const Elt& c, const Elt& d, const NTL::ZZ& m) :
        c(c), d(d), m(m) {
        if (!parametersValid(*this)) {
            throw InvalidParametersException();
        }
    }

    /// Equality test
    auto operator==(const Curve& that) const -> bool {
        return (c == that.c) && (d == that.d) && (m == that.m);
    }

    /// Inequality test
    auto operator!=(const Curve& that) const -> bool {
        return not(*this == that);
    }

    /// Number of rational points on curve over field
    auto cardinality() const -> NTL::ZZ {
        return m;
    }

    /// Output of curve's information
    friend auto operator<<(std::ostream& out, const Curve& curve) ->
    std::ostream& {
        out << getName(curve) << std::endl
        << "Curve parameter #1: "  << curve.c << std::endl
        << "Curve parameter #2: "  << curve.d << std::endl
        << "Cardinality: " << curve.m;
        return out;
    }
};


//------- Edwards Curves -------//

/**
 * @var OddCurve
 * @brief Edwards Curve over finite field with odd characteristic
 *
 */
using OddCurve = Curve<NTL::ZZ_pE, CurveID::Odd>;


/**
 * @p getName(const OddCurve&)
 * @brief output of name
 */
inline auto getName(const OddCurve&) -> std::string {
    return std::string("Edwards Curve over odd field");
}

/**
 * @p parametersValid(const OddCurve& curve)
 * @brief Parameter validation
 */
inline auto parametersValid(const OddCurve& curve) -> bool {
    return (curve.c != 0 && NTL::sqr(curve.c) + curve.c != curve.d);
}

/**
 * @p curveEquation(const OddCurve& curve, const NTL::ZZ_pE& x,
 *      const NTL::ZZ_pE& y)
 * @brief Validation of point coordinates
 */
inline auto curveEquation(const OddCurve& curve, const NTL::ZZ_pE& x,
                          const NTL::ZZ_pE& y) -> bool {
    auto xx = NTL::sqr(x), yy = NTL::sqr(y);
    return (xx + yy == NTL::sqr(curve.c) * (1 + curve.d * xx * yy));
}


/**
 * @var BinaryCurve
 * @brief Edwards Curve over finite field of characteristic two
 */
using BinaryCurve = Curve<NTL::GF2E, CurveID::Binary>;

/**
 * @p getName(const BinaryCurve&)
 * @brief output of name
 */
inline auto getName(const BinaryCurve&) -> std::string {
    return "Binary Edwards Curve";
}

/**
 * @p parametersValid(const BinaryCurve& curve)
 * @brief Parameter validation
 */
inline auto parametersValid(const BinaryCurve& curve) -> bool {
    return (curve.c * curve.d * (1 - curve.d * NTL::power(curve.c, 4)) !=
    0);
}

/**
 * @p curveEquation(const BinaryCurve& curve, const NTL::GF2E& x,
 *      const NTL::GF2E& y)
 * @brief Validation of point coordinates
 */
inline auto curveEquation(const BinaryCurve& curve, const NTL::GF2E& x,
                          const NTL::GF2E& y) -> bool {
    auto xx = NTL::sqr(x), yy = NTL::sqr(y);
    return (curve.c * (x + y) + curve.d * (xx + yy)
    == (x + xx) * (y + yy));
}


/**
 * @var TwistedCurve
 * @brief Twisted Edwards Curve over a non-binary field
 */
using TwistedCurve = Curve<NTL::ZZ_pE, CurveID::Twisted>;

/**
 * @p getName(const TwistedCurve&)
 * @brief output of name
 */
inline auto getName(const TwistedCurve&) -> std::string {
    return "Twisted Edwards Curve";
}

/**
 * @p parametersValid(const TwistedCurve& curve)
 * @brief Parameter validation
 */
inline auto parametersValid(const TwistedCurve& curve) -> bool {
    return (curve.c != curve.d) || ((curve.c == 0) || (curve.d == 0));
}

/**
 * @p curveEquation(const TwistedCurve& curve, const NTL::ZZ_pE& x,
 *      const NTL::ZZ_pE& y)
 * @brief Validation of point coordinates
 */
inline auto curveEquation(const TwistedCurve& curve, const NTL::ZZ_pE& x,
                          const NTL::ZZ_pE& y) -> bool {
    auto xx = NTL::sqr(x), yy = NTL::sqr(y);
    return (curve.c * xx + yy == 1 + curve.d * xx * yy);
}
}

#endif  // _CURVES_H
