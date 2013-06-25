/**
 * @file curves_test.cc
 * @brief A quick test of curve functionality
 * @author Graham Enos
 *
 * This file gives a demonstration of current curve functionality in e2c2.
 */

#include <cstdlib>          // User input
#include <iostream>         // Readable output
#include "e2c2.h"

using namespace std;
using namespace NTL;
using namespace e2c2;


/**
 * @p odd_test()
 * @brief Test of odd curve functionality
 */
void odd_test() {
    ZZ_p::init(power2_ZZ(255) - 19);  /// Sets F_p
    ZZ_pE::init(ZZ_pX(1, 1));  /// Sets F_(p^n) = F_(p^1)

    /// Bernstein's "Curve25519"
    auto c = to_ZZ_pE(1), d = to_ZZ_pE(121665) / to_ZZ_pE(121666);
    auto m = 8 * (power2_ZZ(252) +
                  to_ZZ("27742317777372353535851937790883648493"));
    OddCurve o(c, d, m);
    cout << o << endl;
}


/**
 * @p binary_test()
 * @brief Test of binary curve functionality
 */
void binary_test() {
    /// Our irred. polynomial is x^163 + x^7 + x^6 + x^3 + 1, per FIPS 186-3
    GF2E::init(GF2X(163, 1) + GF2X(7, 1) + GF2X(6, 1) + GF2X(3, 1) +
               GF2X(0, 1));
    GF2X::HexOutput = true;  /// more compact output 
    auto n = 163;
    auto m = to_ZZ("5846006549323611672814742442876390689256843201587");
    auto a2 = to_GF2E(1), a6 = GF2E::zero();
    /// a6 = b in Fips 186-3 language
    set_parameter(a6, "20a601907b8c953ca1481eb10512f78744a3205fd", true);
    auto b_163 = from_weierstrass(n, m, a2, a6);

    cout << b_163 << endl;
    cout << "Here are the (c, d) parameters again: " << endl;
    cout << "(" << b_163.c << ", " << b_163.d << ")" << endl;

    /// Purposefully raising an InvalidParametersException
    cout << endl << endl << "Hold on tight; I'm going to try to"
        << " make a binary curve with c = d = 0..." << endl;
    BinaryCurve b(GF2E::zero(), GF2E::zero(), ZZ::zero());
}


/**
 * @p twisted_test()
 * @brief Test of twisted curve functionality
 */
void twisted_test() {
    ZZ_p::init(power2_ZZ(255) - 19);  /// Sets F_p
    ZZ_pE::init(ZZ_pX(1, 1));  /// Sets F_(p^n) = F_(p^1)

    /// Bernstein's "Curve25519," in twisted form
    auto c = to_ZZ_pE(121666), d = to_ZZ_pE(121665);
    auto m = 8 * (power2_ZZ(252) +
                  to_ZZ("27742317777372353535851937790883648493"));
    TwistedCurve t(c, d, m);
    cout << t << endl;
}


/**
 * @p main(int argc, char *argv[])
 * @brief Runs appropriate curve test
 */
int main(int argc, char *argv[]) {
    try {
        if (argc == 1)
            binary_test();
        else {
            switch(atoi(argv[1])) {
            case 0:
                odd_test();
                break;
            case 1:
                binary_test();
                break;
            case 2:
                twisted_test();
                break;
            default:
                cout << "Please select a type of curve." << endl;
            }
        }
    } catch (InvalidParametersException& e) {
        cout << e.what() << endl;
    }

    return 0;
}
