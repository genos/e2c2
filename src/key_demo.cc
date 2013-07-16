/**
 * @file key_demo.cc
 * @brief Demo of DH Key Exchange with projective points
 * @author Graham Enos
 *
 * This is a sample file that makes use of the classes and methods in
 * b_projective.cc and b_edwards.cc to conduct a DH key exchange demonstration.
 */

#include <iostream>
#include <sstream>
#include "e2c2.h"

using namespace std;
using namespace NTL;
using namespace e2c2;

/**
 * @brief A quick key exchange demo
 *
 * This function shows how to use our projective points on a Binary Edwards
 * Curve to conduct a Diffie-Hellman key exchange.
 */
auto main(int argc, char *argv[]) -> int
{
    /// Our irred. polynomial is x^163 + x^7 + x^6 + x^3 + 1, per FIPS 186-3
    GF2E::init(GF2X(163, 1) + GF2X(7, 1) + GF2X(6, 1) + GF2X(3, 1) +
               GF2X(0, 1));
    GF2X::HexOutput = true;  /// more compact output 
    auto n = 163;
    GF2E a2 = to_GF2E(1), a6;
    /// a6 = b in Fips 186-3 language
    set_parameter(a6,
                 string("20a601907b8c953ca1481eb10512f78744a3205fd"),
                 true);
    auto m = to_ZZ("5846006549323611672814742442876390689256843201587");
    auto c = mol_alg_1(n, a2, a6);
    auto d = NTL::sqr(c) + c + gf2m_sqrt(a6, n) / NTL::sqr(c);
    auto b_163 = from_weierstrass(n, m, a2, a6);
    
    /// Weierstrass params; need to be changed to Edwards Curve
    GF2E g_x, g_y, x, y, z;
    set_parameter(g_x, "3f0eba16286a2d57ea0991168d4994637e8343e36", true);
    set_parameter(g_y, "0d51fbc6c71a0094fa2cdd545b11c5c0c797324f1", true);
    mol_bm_proj(x, y, z, g_x, g_y, n, c, d, a2);
    BinaryProj P(x, y, z, b_163);
    BinaryProj id(b_163);

    /// Key Exchange
    cout << "Alice and Bob wish to communicate in private, " <<
        "so they need a shared secret key." << endl << endl << endl;
    /// Seeding pseudorandom generator
    SetSeed(argc <= 1 ? to_ZZ(1729) : to_ZZ(argv[1]));
    /// Alice's first steps
    auto a = RandomLen_ZZ(NumBits(m));
    cout << "Alice first picks a secret random number a = " << a <<
        " with the same number of bits as m (the size of our group)..." <<
        endl;
    auto aP = a * P;
    cout << "...then sends Bob the point pA = a * P (the generator) = "
        << aP << endl << endl;
    /// Bob's first steps
    auto b = RandomLen_ZZ(NumBits(m));
    cout << "Bob also picks a secret random number b = " << b << endl;
    auto bP = b * P;
    cout << "Then he sends Alice the point bP = b * P = " << bP << endl
        << endl;
    /// Alice reconstructs private key
    auto key_a = a * bP;
    cout << "Then Alice takes a and multiplies it by bP to get key_a = " <<
        key_a << endl;
    /// Bob does the same
    auto key_b = b * aP;
    cout << "Similarly, Bob calculates key_b = b * aP = " << key_b << endl;
    /// Quick check
    cout << "Are these in fact the same key? " << endl;
    if (key_a == key_b) {
        cout << "Yes!" << endl <<
            "So now they share a secret key, and can communicate securely."
             << endl << endl << endl;
    } else {
        cout << "NO...uh oh, my code is wrong somewhere..." << endl;
    }
    return 0;
}
