/*
 * This file is part of the BSGS distribution (https://github.com/JeanLucPons/BSGS).
 * Copyright (c) 2020 Jean Luc PONS.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Optimized by NoMachine
 */

#include "SECP256K1.h"
#include <string.h>

Secp256K1::Secp256K1() {}

void Secp256K1::Init() {
    // Prime for the finite field
    Int P;
    P.SetBase16("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F");
    Int::SetupField(&P);

    // Generator point (G) and order (n)
    G.x.SetBase16("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798");
    G.y.SetBase16("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8");
    G.z.SetInt32(1);
    order.SetBase16("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141");

    Int::InitK1(&order);

    // Compute Generator table (pre-computed multiples of G for faster scalar multiplication)
    Point N(G);
    for (int i = 0; i < 32; i++) {
        GTable[i * 256] = N;
        N = DoubleDirect(N); // Double the current point to get the next multiple of G
        for (int j = 1; j < 255; j++) {
            GTable[i * 256 + j] = N;
            N = AddDirect(N, GTable[i * 256]); // Add G to the current point
        }
        GTable[i * 256 + 255] = N; 
    }
}

Secp256K1::~Secp256K1() {}

// Direct addition of two points on the elliptic curve
Point Secp256K1::AddDirect(Point &p1, Point &p2) {
    Int _s, _p, dy, dx, dxInv;
    Point r;
    r.z.SetInt32(1);

    // dy = y2 - y1 (difference in y-coordinates)
    dy.ModSub(&p2.y, &p1.y);

    // dx = x2 - x1 (difference in x-coordinates)
    dx.ModSub(&p2.x, &p1.x);

    // dxInv = 1 / dx (modular inverse of dx)
    dxInv = dx;
    dxInv.ModInv();

    // _s = dy / dx (slope of the line through p1 and p2)
    _s.ModMulK1(&dy, &dxInv);

    // _p = _s^2 (square of the slope)
    _p.ModSquareK1(&_s);

    // r.x = _s^2 - x1 - x2 (x-coordinate of the resulting point)
    r.x.ModSub(&_p, &p1.x);
    r.x.ModSub(&p2.x);

    // r.y = _s(x1 - r.x) - y1 (y-coordinate of the resulting point)
    r.y.ModSub(&p2.x, &r.x);
    r.y.ModMulK1(&_s);
    r.y.ModSub(&p2.y);

    return r;
}

// Optimized mixed Jacobian-affine addition of two points
Point Secp256K1::Add2(Point &p1, Point &p2) {
    Int u, v, vs2, vs3, us2, a, tmp;
    Point r;

    // U = Y2*Z1 - Y1 (difference in y-coordinates adjusted for Z1)
    u.ModMulK1(&p2.y, &p1.z);
    u.ModSub(&p1.y);

    // V = X2*Z1 - X1 (difference in x-coordinates adjusted for Z1)
    v.ModMulK1(&p2.x, &p1.z);
    v.ModSub(&p1.x);

    // V² and V³ (used in subsequent calculations)
    vs2.ModSquareK1(&v);
    vs3.ModMulK1(&vs2, &v);

    // U² (square of the difference in y-coordinates)
    us2.ModSquareK1(&u);

    // A = U²*Z1 - V³ - 2*V²*X1 (intermediate value for x-coordinate calculation)
    a.ModMulK1(&us2, &p1.z);
    tmp.ModMulK1(&vs2, &p1.x);
    tmp.ModAdd(&tmp);
    a.ModSub(&vs3);
    a.ModSub(&tmp);

    // X3 = V*A (final x-coordinate of the resulting point)
    r.x.ModMulK1(&v, &a);

    // Y3 = (V²*X1 - A)*U - V³*Y1 (final y-coordinate of the resulting point)
    tmp.ModMulK1(&vs2, &p1.x);
    tmp.ModSub(&a);
    tmp.ModMulK1(&u);
    r.y.ModMulK1(&vs3, &p1.y);
    r.y.ModNeg();
    r.y.ModAdd(&tmp);

    // Z3 = V³*Z1 (final z-coordinate of the resulting point)
    r.z.ModMulK1(&vs3, &p1.z);

    return r;
}

// General addition of two points on the elliptic curve
Point Secp256K1::Add(Point &p1, Point &p2) {
    Int u, v, u1, u2, v1, v2, vs2, vs3, us2, w, a, tmp;
    Point r;

    // u1 = Y2*Z1, u2 = Y1*Z2 (adjusted y-coordinates for both points)
    u1.ModMulK1(&p2.y, &p1.z);
    u2.ModMulK1(&p1.y, &p2.z);
    u.ModSub(&u1, &u2);
    
    // v1 = X2*Z1, v2 = X1*Z2 (adjusted x-coordinates for both points)
    v1.ModMulK1(&p2.x, &p1.z);
    v2.ModMulK1(&p1.x, &p2.z);
    v.ModSub(&v1, &v2);
    
    // w = Z1*Z2 (product of z-coordinates)
    w.ModMulK1(&p1.z, &p2.z);

    // us2 = U², vs2 = V² (squares of differences in coordinates)
    us2.ModSquareK1(&u);
    vs2.ModSquareK1(&v);

    // vs3 = V³ (cube of the difference in x-coordinates)
    vs3.ModMulK1(&vs2, &v);
    
    // A = U²*Z1 - V³ - 2*V²*X1 (intermediate value for x-coordinate calculation)
    a.ModMulK1(&us2, &w);
    tmp.ModMulK1(&vs2, &v2);
    tmp.ModAdd(&tmp);
    a.ModSub(&vs3);
    a.ModSub(&tmp);

    // X3 = V*A (final x-coordinate of the resulting point)
    r.x.ModMulK1(&v, &a);

    // Y3 = (V²*X1 - A)*U - V³*Y1 (final y-coordinate of the resulting point)
    tmp.ModMulK1(&vs3, &u2);
    r.y.ModSub(&vs2, &a);
    r.y.ModMulK1(&r.y, &u);
    r.y.ModSub(&tmp);

    // Z3 = V³*Z1 (final z-coordinate of the resulting point)
    r.z.ModMulK1(&vs3, &w);

    return r;
}

// Doubling a point on the elliptic curve using direct affine coordinates
Point Secp256K1::DoubleDirect(Point &p) {
    Int s, x_sq, inv_2y, tmp;
    Point r;
    r.z.SetInt32(1);

    // s = 3x² (slope formula for doubling)
    x_sq.ModSquareK1(&p.x);
    s.ModAdd(&x_sq, &x_sq);
    s.ModAdd(&x_sq);
    
    // s = (3x²)/(2y) (compute the slope)
    inv_2y.ModAdd(&p.y, &p.y);
    inv_2y.ModInv();
    s.ModMulK1(&s, &inv_2y);

    // r.x = s² - 2x (x-coordinate of the resulting point)
    r.x.ModSquareK1(&s);
    tmp.ModAdd(&p.x, &p.x);
    r.x.ModSub(&tmp);

    // r.y = s*(x - r.x) - y (y-coordinate of the resulting point)
    tmp.ModSub(&p.x, &r.x);
    r.y.ModMulK1(&s, &tmp);
    r.y.ModSub(&p.y);

    return r;
}

// Doubling a point on the elliptic curve using Jacobian coordinates
Point Secp256K1::Double(Point &p) {
    Int z2, x2, w, s, s2, b, h, tmp, _8b;
    Point r;

    // z2 = Z², x2 = X² (squares of z and x coordinates)
    z2.ModSquareK1(&p.z);
    x2.ModSquareK1(&p.x);
    
    // w = 3x² + a*z⁴ (slope formula for doubling, a=0 in secp256k1)
    w.ModAdd(&x2, &x2);
    w.ModAdd(&x2);
    w.ModAdd(&z2);

    // s = Y*Z, b = X*Y*s (intermediate values for coordinate calculation)
    s.ModMulK1(&p.y, &p.z);
    b.ModMulK1(&p.x, &p.y);
    b.ModMulK1(&s);
    
    // h = w² - 8b (intermediate value for x-coordinate calculation)
    h.ModSquareK1(&w);
    _8b.ModAdd(&b, &b);
    _8b.ModDouble();
    _8b.ModDouble();
    h.ModSub(&_8b);

    // X3 = 2*h*s (final x-coordinate of the resulting point)
    r.x.ModMulK1(&h, &s);
    r.x.ModAdd(&r.x);

    // s2 = S², tmp = 8*y²*s² (intermediate values for y-coordinate calculation)
    s2.ModSquareK1(&s);
    tmp.ModSquareK1(&p.y);
    tmp.ModMulK1(&s2);
    tmp.ModDouble();
    tmp.ModDouble();
    tmp.ModDouble();

    // Y3 = (4*b - h)*w - 8*y²*s² (final y-coordinate of the resulting point)
    r.y.ModAdd(&b, &b);
    r.y.ModAdd(&r.y, &r.y);
    r.y.ModSub(&h);
    r.y.ModMulK1(&w);
    r.y.ModSub(&tmp);

    // Z3 = 8*S³ (final z-coordinate of the resulting point)
    r.z.ModMulK1(&s2, &s);
    r.z.ModDouble();
    r.z.ModDouble();
    r.z.ModDouble();

    return r;
}

// Compute the y-coordinate corresponding to a given x-coordinate
Int Secp256K1::GetY(Int x, bool isEven) {
    Int y, x2;
    x2.ModSquareK1(&x); // x²
    y.ModMulK1(&x2, &x); // x³
    y.ModAdd(7); // x³ + 7
    y.ModSqrt(); // sqrt(x³ + 7)

    // Adjust y to ensure it matches the parity (even or odd)
    if (y.IsEven() ^ isEven) {
        y.ModNeg();
    }
    return y;
}

// Verify if a point lies on the elliptic curve
bool Secp256K1::EC(Point &p) {
    Int lhs, rhs, x2;
    x2.ModSquareK1(&p.x); // x²
    rhs.ModMulK1(&x2, &p.x); // x³
    rhs.ModAdd(7); // x³ + 7
    lhs.ModSquareK1(&p.y); // y²
    lhs.ModSub(&rhs); // y² - (x³ + 7)
    return lhs.IsZero(); // Check if the result is zero
}

// Compute the public key from a private key
Point Secp256K1::ComputePublicKey(Int *privKey) {
    Point Q;
    Q.Clear();
    
    int i = 0;
    // Skip leading zero bytes in the private key
    while (i < 32 && privKey->GetByte(i) == 0) {
        i++;
    }

    if (i < 32) {
        uint8_t b = privKey->GetByte(i);
        if (b > 0) {
            Q = GTable[256 * i + (b - 1)]; // Use precomputed multiples of G
        }
        i++;
    }

    for (; i < 32; i++) {
        uint8_t b = privKey->GetByte(i);
        if (b > 0) {
            Q = Add2(Q, GTable[256 * i + (b - 1)]); // Add precomputed multiples of G
        }
    }

    Q.Reduce(); // Normalize the resulting point
    return Q;
}
