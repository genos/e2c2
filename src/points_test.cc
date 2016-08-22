/**
 * @file points_test.cc
 * @brief A quick test of point functionality
 * @author Graham Enos
 *
 * This file gives a demonstration of current point functionality in e2c2.
 */

#include <cstdlib>          // User input
#include <iostream>         // Readable output
#include "e2c2.h"

using namespace std;
using namespace NTL;
using namespace e2c2;


/**
 * @p odd_test()
 * @brief Test of odd point functionality
 */
auto odd_test() -> void {
    ZZ_p::init(power2_ZZ(255) - 19);  /// Sets F_p
    ZZ_pE::init(ZZ_pX(1, 1));  /// Sets F_(p^n) = F_(p^1)

    /// Bernstein's "Curve25519"
    auto c = to_ZZ_pE(1), d = to_ZZ_pE(121665) / to_ZZ_pE(121666);
    auto m = 8 * (power2_ZZ(252) +
    to_ZZ("27742317777372353535851937790883648493"));
    OddCurve o(c, d, m);
    cout << o << endl;

    OddAff id(o);
    cout << "id = " << id << endl;
    cout << "17 * id = " << 17 * id << endl;

    auto x = to_ZZ_pE(1), y = ZZ_pE::zero(), z = to_ZZ_pE(1);
    OddProj point1(x, y, z, o), point2(point1);
    for (auto i = 0; i < 4; ++i) {
        cout << "ProjectivePoint 2 = " << point2 << endl;
        cout << "ProjPoint2 + ProjPoint1 = " << point2 + point1 << endl;
        point2 += point1;
    }

    cout << (point1 + OddProj(OddAff(x, y, o))) << endl;

    cout << -point1 << endl;

    cout << "3 * point2 = " << 3 * point2
    << endl;

}


/**
 * @p binary_test()
 * @brief Test of binary point functionality
 */
auto binary_test() -> void {
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

    BinaryAff id(b_163);
    cout << "id = " << id << endl;
    cout << "17 * id = " << 17 * id << endl;

    GF2E x = to_GF2E(1);
    GF2E y = to_GF2E(1);
    GF2E z = to_GF2E(1);
    BinaryProj point1(x, y, z, b_163);
    BinaryProj point2(point1);
    for (auto i = 0; i < 4; ++i) {
        cout << "ProjectivePoint 2 = " << point2 << endl;
        cout << "ProjPoint2 + ProjPoint1 = " << point2 + point1 << endl;
        point2 += point1;
    }

    cout << (point1 + BinaryProj(BinaryAff(x, y, b_163))) << endl;

    cout << -point1 << endl;

    cout << "3 * point2 = " << 3 * point2
    << endl;
    /// Purposefully raising an InvalidParametersException
    cout << endl << endl << "Hold on tight; I'm going to try to"
    << " make a binary point with (x : y : z) = (1 : 0 : 1)..." << endl;
    BinaryProj(to_GF2E(1), GF2E::zero(), to_GF2E(1), b_163);
}


/**
 * @p twisted_test()
 * @brief Test of twisted point functionality
 */
auto twisted_test() -> void {
    ZZ_p::init(power2_ZZ(255) - 19);  /// Sets F_p
    ZZ_pE::init(ZZ_pX(1, 1));  /// Sets F_(p^n) = F_(p^1)

    /// Bernstein's "Curve25519," in twisted form
    auto c = to_ZZ_pE(121666), d = to_ZZ_pE(121665);
    auto m = 8 * (power2_ZZ(252) +
    to_ZZ("27742317777372353535851937790883648493"));
    TwistedCurve t(c, d, m);
    cout << t << endl;

    TwistedAff id(t);
    cout << "id = " << id << endl;
    cout << "17 * id = " << 17 * id << endl;

    ZZ_pE x = ZZ_pE::zero();
    ZZ_pE y = to_ZZ_pE(-1);
    ZZ_pE z = to_ZZ_pE(1);
    TwistedProj point1(x, y, z, t);
    TwistedProj point2(point1);
    for (auto i = 0; i < 4; ++i) {
        cout << "ProjectivePoint 2 = " << point2 << endl;
        cout << "ProjPoint2 + ProjPoint1 = " << point2 + point1 << endl;
        point2 += point1;
    }

    cout << (point1 + TwistedProj(TwistedAff(x, y, t))) << endl;

    cout << -point1 << endl;

    cout << "3 * point2 = " << 3 * point2
    << endl;
    /// Purposefully raising an InvalidParametersException
    cout << endl << endl << "Hold on tight; I'm going to try to"
    << " make a binary point with (x : y : z) = (1 : 0 : 1)..." << endl;
    TwistedProj(to_ZZ_pE(1), ZZ_pE::zero(), to_ZZ_pE(1), t);
}


/**
 * @p main(int argc, char *argv[])
 * @brief Runs appropriate point test
 */
auto main(int argc, char *argv[]) -> int {
    try {
        if (argc == 1) {
            binary_test();
        } else {
            switch (atoi(argv[1])) {
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
