/**
 * @file ecoh_echo.cc
 * @brief Elliptic Curve only key derivation function
 * @author Graham Enos
 *
 */
#include <algorithm>
#include <iostream>
#include <sstream>
#include "e2c2.h"

using namespace std;
using namespace NTL;
using namespace e2c2;

//------- Constants --------//
/// Number of blocks
const size_t NUM_BLOCKS = 32;
/// Length of each block
const size_t B_LEN = 192;
/// Bit length of integer representation
const size_t I_LEN = 64;


//------- Utilities -------//
/// Helps omega() decide which coordinate to use
enum class Coord { X, Y };

auto ZZtoHex(const ZZ& z) -> string {
    /// The hex representation of z
    static const size_t out_len = (B_LEN + I_LEN) >> 4;
    static const string digits = "0123456789abcdef";
    stringstream ss;
    if (z <= 0) {
        ss << "0x";
        while (ss.str().length() < out_len + 2) {
            ss << "0";
        }
        return ss.str();
    } else {
        auto zz = z;
        while (zz > 0) {
            ss << digits.at(to_long(zz & 0xf));
            zz  >>= 4;
        }
        while (ss.str().length() < out_len) {
            ss << "0";
        }
        ss << "x0";
        auto s = ss.str();
        reverse(s.begin(), s.end());
        return s;
    }
}

auto StoZZ(const string& s) -> ZZ {
    /// Treat a string of characters as a bunch of bytes (ASCII), then use
    /// those bytes to build a ZZ
    size_t i = 0;
    auto z = ZZ::zero();
    for (auto c : s) {
        ++i;
        z <<= 8;
        z |= to_ZZ(c);
    }
    /// Extra NOPS (to protect against side-channel/timing attacks)
    while (i < (B_LEN >> 3)) {
        ++i;
    }
    return z;
}

auto GF2EtoZZ(const GF2E& x) -> ZZ {
    /// Use the bit representation of x to create a ZZ
    GF2X::HexOutput = false;
    auto z = ZZ::zero();
    stringstream ss;
    ss << x;
    auto s = ss.str();
    /// Only leave behind 0s and 1s (remove spaces and brackets)
    for (auto c : string("[ ]")) {
        auto i = remove(s.begin(), s.end(), c);
        while (i != s.end()) {
            s.erase(i);
            i = remove(s.begin(), s.end(), c);
        }
    }
    /// Set bits where there are 1s
    for (auto c : s) {
        z <<= 1;
        z |= to_ZZ(c == '1') & 0x1 ;
    }
    return z;
}

auto ZZtoGF2E(const ZZ& z) -> GF2E {
    /// Use the bit representation of z to create a GF2E
    auto x = GF2E::zero();
    stringstream ss;
    ss << "[";
    for (auto i = 0; i < NumBits(z); ++i) {
        ss << bit(z, i) << " ";
    }
    ss << "]";
    ss >> x;
    return x;
}


//------- ECOH's Echo Helpers -------//
auto block(const size_t i) -> ZZ {
    /// The ith block O_i is (random)^{B_LEN - I_LEN} || {I_LEN bit repr of i},
    /// similar to the original ECOH submission
    auto O = (RandomBits_ZZ(B_LEN - I_LEN)) << I_LEN;
    O |= to_ZZ(i) & ((to_ZZ(1) << I_LEN) - 1);
    return O;
}

auto omega(const BinaryAff& Q, const BinaryAff& G, const Coord z) -> ZZ {
    /// Output, based on coordinate, of bit-length B_LEN
    if (z == Coord::X) {
        return GF2EtoZZ((Q + (GF2EtoZZ(Q.x) / 2) * G).x) / 2 %
        power2_ZZ(B_LEN);
    } else {  /// z == Coord::Y
        return GF2EtoZZ((Q + (GF2EtoZZ(Q.y) / 2) * G).y) / 2 %
        power2_ZZ(B_LEN);
    }
}

auto pi(const GF2E& z, const GF2E& a2, const GF2E& a6,
        const BinaryCurve& E) -> BinaryAff {
    /// Uses Icart's f from "How to Hash into Elliptic Curves" and Moloney,
    /// O'Mahony, & Laurent's birational map from "Efficient Implementation of
    /// Elliptic Curve Point Operations Using Binary Edwards Curves"
    const auto alpha = a2 + z + sqr(z);
    const auto u = power(power(alpha, 4) + power(alpha, 3) + a6,
    (2 * GF2E::cardinality() - 1) / 3) + alpha;
    const auto v = z * u + sqr(alpha);
    return birMapAff(u, v, a2, E);
}


//------- ECOH's Echo -------//
auto ecoh_echo(const ZZ& password, const ZZ& salt, const GF2E& a2,
               const GF2E& a6, const BinaryCurve& E, const BinaryAff& G) -> ZZ {
    /// Our PBKDF

    /// Helpers defined as lambda expressions
    auto phi = [&G](const BinaryAff & P) {
        return omega(P, G, Coord::X);
    };
    auto psi = [&G](const BinaryAff & P) {
        return omega(P, G, Coord::Y) & 1;
    };

    /// Memory efficient version
    SetSeed(salt ^ password);
    auto O = block(0);
    auto X_1_tmp = to_ZZ(NUM_BLOCKS);
    auto X_2_tmp = to_ZZ(NUM_BLOCKS);
    auto P = pi(ZZtoGF2E(O), a2, a6, E);
    auto Q = P;
    for (size_t i = 1; i < NUM_BLOCKS; ++i) {
        O = block(i);
        P = pi(ZZtoGF2E(O ^ i ^ phi(P)), a2, a6, E);
        X_1_tmp ^= (psi(P) << i);
        X_2_tmp ^= O;
        Q += P;
    }
    auto X1 = pi(ZZtoGF2E(X_1_tmp), a2, a6, E);
    auto X2 = pi(ZZtoGF2E(X_1_tmp ^ X_2_tmp), a2, a6, E);
    Q += X1 + X2;
    return phi(Q);
}


//------- An Example -------//
auto main(int argc, char *argv[]) -> int {
    //------- Error Checking -------//
    if (argc != 3) {
        cerr << "usage: " << argv[0] << " salt password" << endl;
        return EXIT_FAILURE;
    }

    //------- Setup -------//
    /// Our irred. polynomial is x^163 + x^7 + x^6 + x^3 + 1, per FIPS 186-3
    GF2E::init(GF2X(163, 1) + GF2X(7, 1) + GF2X(6, 1) + GF2X(3, 1) +
    GF2X(0, 1));
    GF2X::HexOutput = true;  /// more compact output
    auto a2 = to_GF2E(1), a6 = GF2E::zero();
    /// a6 = b in Fips 186-3 language
    set_parameter(a6, "20a601907b8c953ca1481eb10512f78744a3205fd", true);
    auto E = from_weierstrass(163,
    to_ZZ("5846006549323611672814742442876390689256843201587"),
    a2,
    a6);
    auto u = GF2E::zero(), v = GF2E::zero();
    set_parameter(u, "3f0eba16286a2d57ea0991168d4994637e8343e36", true);
    set_parameter(v, "0d51fbc6c71a0094fa2cdd545b11c5c0c797324f1", true);
    auto x = GF2E::zero(), y = GF2E::zero();
    BinaryAff G = birMapAff(u, v, a2, E);

    //------- Run ECOH's ECHO -------//
    auto salt = to_ZZ(argv[1]);
    auto password = StoZZ(argv[2]);
    cout << ZZtoHex(ecoh_echo(password, salt, a2, a6, E, G)) << endl;

    return EXIT_SUCCESS;
}
