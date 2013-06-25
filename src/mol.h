/**
 * @file mol.h
 * @brief Birational Map and utilities from MOL Paper
 * @author Graham Enos
 *
 * This file contains the interface to the C++ implementation of the birational
 * map and assorted utility functions given in Moloney, O'Mahony, & Laurent's
 * paper, available here: http://eprint.iacr.org/2010/208
 */

#ifndef _MOL_H
#define _MOL_H

#include <iostream>         // Readable output
#include <NTL/ZZ.h>         // Arbitrarily large integers
#include <NTL/ZZ_pE.h>      // Field elements from @f$ \mathbf{F}_{p^n} @f$
#include <NTL/GF2E.h>       // Field elements from @f$ \mathbf{F}_{2^n} @f$
#include "curves.h"         // Edwards curves (and variations) 
#include "utilities.h"      // Utilities header for e2c2 project


/// Namespace for our library
namespace e2c2 {

    /// Per MOL paper, \f$\sqrt{\alpha} = \alpha^{2^{m-1}}\f$
    inline const NTL::GF2E gf2m_sqrt(const NTL::GF2E& alpha, const long& m) {
        auto s = alpha;
        for (auto i = 0L; i < m - 1; ++i) {
            sqr(s, s);
        }
        return s;
    }


    /// Half-Trace function
    inline const NTL::GF2E half_trace(const NTL::GF2E& alpha, const long m) {
        auto ht = NTL::GF2E::zero();
        auto e = NTL::to_ZZ(1);
        for (auto i = 0L; i <= (m - 1) / 2; ++i) {
            ht += power(alpha, e);
            e <<= 2;
        }
        return ht;
    }


    /// MOL's Algorithm 1 to compute d1
    inline const NTL::GF2E mol_alg_1(const long& n, const NTL::GF2E& a2,
                            const NTL::GF2E& a6) {
        auto t = trace(a2), r = trace(a6);
        auto a6_2 = gf2m_sqrt(a6, n);
        auto a6_4 = gf2m_sqrt(a6_2, n);
        auto x = NTL::GF2E::zero(), d1 = NTL::GF2E::zero();
        set_parameter(x, "[0 1]");
        auto w = x + trace(x);
        if (t == 0 && r == 1) {
            set(d1);
        } else {
            if (t == 1 && r == 0) {
                d1 = a6_4;
            } else {
                if (t == 1 && r == 1 && a6 != 1) {
                    if (trace(inv(a6 + 1)) == 1) {
                        d1 = a6_2 + a6_4;
                    } else {
                        d1 = a6_4 + 1;
                    }
                } else {
                    if (t == 1 && a6 == 1) {
                        if (trace(inv(w)) == 1) {
                            d1 = w;
                        } else {
                            if (trace(inv(w + 1)) == 1) {
                                d1 = inv(w + 1);
                            } else {
                                d1 = inv(w + 1) + 1;
                            }
                        }
                    } else {
                        if (t == 0 && r == 0) {
                            if (trace(inv(a6 + 1)) == 0) {
                                d1 = a6_4 + 1;
                            } else {
                                auto i = 1;
                                auto s = a6_2;
                                while (trace(NTL::power(a6,
                                                NTL::power_long(2, i) + 1))
                                        == 0) {
                                    s *= s;
                                    ++i;
                                }
                                d1 = inv(s + 1);
                            }
                        }
                    }
                }
            }
        }
        return d1;
    }


    /// Construct a Binary Edwards Curve from Weierstrass Parameters, degree of
    /// field extension, and supplied cardinality of curve
    inline BinaryCurve from_weierstrass(const long n, const NTL::ZZ& m,
            const NTL::GF2E& a2, const NTL::GF2E& a6) {
        auto c = mol_alg_1(n, a2, a6);
        auto d = NTL::sqr(c) + c + gf2m_sqrt(a6, n) / NTL::sqr(c);
        return BinaryCurve(c, d, m);
    }


    /// MOL's birational map from Weierstrass curve to Affine Binary Edwards
    inline void mol_bm_aff(NTL::GF2E& x, NTL::GF2E& y, const NTL::GF2E& u,
                  const NTL::GF2E& v, const long m, const NTL::GF2E& d1,
                  const NTL::GF2E& d2, const NTL::GF2E& a2) {
        auto b = half_trace(sqr(d1) + d2 + a2, m), tmp = sqr(d1) + d1 + d2;
        auto z = sqr(u) + d1 * u + sqr(d1) * tmp;
        x = d1 * (b * u + v + (sqr(d1) + d1) * tmp);
        y = (x + d1 * u);
        x /= z;
        y /= z;
    }


    /// MOL's birational map from Weierstrass curve to Projective Binary
    /// Edwards
    inline void mol_bm_proj(NTL::GF2E& x, NTL::GF2E& y, NTL::GF2E& z,
            const NTL::GF2E& u, const NTL::GF2E& v, const long m,
            const NTL::GF2E& d1, const NTL::GF2E& d2, const NTL::GF2E& a2) {
        auto b = half_trace(sqr(d1) + d2 + a2, m), tmp = sqr(d1) + d1 + d2;
        x = d1 * (b * u + v + (sqr(d1) + d1) * tmp);
        y = (x + d1 * u);
        z = sqr(u) + d1 * u + sqr(d1) * tmp;
    }
}
#endif // _MOL_H
