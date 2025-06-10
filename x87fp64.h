//=========================================================
//  x87fp64.h
//
//  64-bit floating-point support
//=========================================================
//
// BSD 3-Clause License
//
// Copyright (c) 2025, Aaron Giles
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#pragma once

#ifndef X87FP64_H
#define X87FP64_H

#include "x87common.h"
#include "x87fp80.h"


//===========================================================================
//
// x87::fp64_t
//
// Class for handling x87 floating-point values using 64-bit floats.
// This is generally "good enough" accuracy for video games and many other
// scenarios, and is much faster than manually doing the 80-bit math.
//
//===========================================================================

namespace x87
{

struct fp64_t
{
public:
    //
    // construction/destruction
    //
    explicit fp64_t() { }
    fp64_t(fp64_t const &v64) { m_value.d = v64.m_value.d; }
    explicit fp64_t(fp80_t const &v80) { m_value.d = v80.as_double(); }
    explicit fp64_t(uint64_t man, uint16_t se) { m_value.d = fp80_t(man, se).as_double(); }
    fp64_t(double _val) { m_value.d = _val; }
    explicit fp64_t(float _val) { m_value.d = double(_val); }
    explicit fp64_t(int64_t _val) { m_value.d = double(_val); }
    explicit fp64_t(int32_t _val) { m_value.d = double(_val); }
    explicit fp64_t(int16_t _val) { m_value.d = double(_val); }

    //
    // operators
    //
    fp64_t &operator=(fp64_t const &src) { m_value.d = src.m_value.d; return *this; }
    fp64_t &operator+=(fp64_t const &rhs) { m_value.d += rhs.m_value.d; return *this; }
    fp64_t &operator-=(fp64_t const &rhs) { m_value.d -= rhs.m_value.d; return *this; }
    fp64_t &operator*=(fp64_t const &rhs) { m_value.d *= rhs.m_value.d; return *this; }
    fp64_t &operator/=(fp64_t const &rhs) { m_value.d /= rhs.m_value.d; return *this; }

    //
    // comparison operations
    //
    bool operator==(fp64_t const &rhs) const { return (this->as_double() == rhs.as_double()); }
    bool operator!=(fp64_t const &rhs) const { return (this->as_double() != rhs.as_double()); }
    bool operator<(fp64_t const &rhs) const { return (this->as_double() < rhs.as_double()); }
    bool operator<=(fp64_t const &rhs) const { return (this->as_double() <= rhs.as_double()); }
    bool operator>(fp64_t const &rhs) const { return (this->as_double() > rhs.as_double()); }
    bool operator>=(fp64_t const &rhs) const { return (this->as_double() >= rhs.as_double()); }

    //
    // pieces
    //
    uint64_t mantissa() const { return m_value.i & FP64_MANTISSA_MASK; }
    int32_t exponent() const { return ((m_value.i & FP64_EXPONENT_MASK) >> FP64_EXPONENT_SHIFT) - FP64_EXPONENT_BIAS; }
    uint8_t sign() const { return m_value.i >> FP64_SIGN_SHIFT; }
    uint64_t as_fpbits64() const { return m_value.i; }
    uint32_t as_fpbits32() const { float_int32_t u = { float(m_value.d) }; return u.i; }

    //
    // conversion
    //
    int16_t as_int16() const { return int16_t(m_value.d); }
    int16_t as_int16(x87cw_t round) const { fpround_t r(round); return int16_t(m_value.d); }
    int32_t as_int32() const { return int32_t(m_value.d); }
    int32_t as_int32(x87cw_t round) const { fpround_t r(round); return int32_t(m_value.d); }
    int64_t as_int64() const { return int64_t(m_value.d); }
    int64_t as_int64(x87cw_t round) const { fpround_t r(round); return int64_t(m_value.d); }
    float as_float() const { return float(m_value.d); }
    float as_float(x87cw_t round) const { fpround_t r(round); return float(m_value.d); }
    double as_double(x87cw_t round = X87CW_ROUNDING_NEAREST) const { return m_value.d; }
    fp80_t as_fp80() const { return fp80_t(m_value.d); }

    //
    // misc
    //
    bool isnormal() const { return ((((m_value.i >> FP64_EXPONENT_SHIFT) + 1) & 0x7fe) != 0 || this->iszero()); }
    bool isminexp() const { return ((m_value.i & FP64_EXPONENT_MASK) == 0); }
    bool ismaxexp() const { return ((m_value.i & FP64_EXPONENT_MASK) == FP64_EXPONENT_MASK); }
    bool isnan() const { return ((m_value.i & FP64_ABS_MASK) > 0x7ff0000000000000ull); }
    bool isqnan() const { return ((m_value.i & FP64_ABS_MASK) >= 0x7ff8000000000000ull); }
    bool issnan() const { return this->isnan() && !isqnan(); }
    bool isinf() const { return ((m_value.i & FP64_ABS_MASK) == 0x7ff0000000000000ull); }
    bool ispinf() const { return (m_value.i == 0x7ff0000000000000ull); }
    bool isninf() const { return (m_value.i == 0xfff0000000000000ull); }
    bool iszero() const { return ((m_value.i & FP64_ABS_MASK) == 0); }
    bool isdenorm() const { return ((m_value.i & FP64_EXPONENT_MASK) == 0 && !this->iszero()); }
    fp64_t &copysign(fp64_t const &src) { m_value.i = (m_value.i & FP64_ABS_MASK) | (src.m_value.i & FP64_SIGN_MASK); return *this; }

    //
    // static constants
    //
    static fp64_t const_zero()  { static int64_double_t const c = { 0x0000000000000000 }; return fp64_t(c.d); }
    static fp64_t const_nzero() { static int64_double_t const c = { 0x8000000000000000 }; return fp64_t(c.d); }
    static fp64_t const_one()   { static int64_double_t const c = { 0x3ff0000000000000 }; return fp64_t(c.d); }
    static fp64_t const_none()  { static int64_double_t const c = { 0xbff0000000000000 }; return fp64_t(c.d); }
    static fp64_t const_2t()    { static int64_double_t const c = { 0x400a934f0979a371 }; return fp64_t(c.d); }
    static fp64_t const_2e()    { static int64_double_t const c = { 0x3ff71547652b82fe }; return fp64_t(c.d); }
    static fp64_t const_pi()    { static int64_double_t const c = { 0x400921fb54442d18 }; return fp64_t(c.d); }
    static fp64_t const_lg2()   { static int64_double_t const c = { 0x3fd34413509f79ff }; return fp64_t(c.d); }
    static fp64_t const_ln2()   { static int64_double_t const c = { 0x3fe62e42fefa39ef }; return fp64_t(c.d); }
    static fp64_t const_snan()  { static int64_double_t const c = { 0x7ff0000000000001 }; return fp64_t(c.d); }
    static fp64_t const_qnan()  { static int64_double_t const c = { 0x7ff8000000000000 }; return fp64_t(c.d); }
    static fp64_t const_pinf()  { static int64_double_t const c = { 0x7ff0000000000000 }; return fp64_t(c.d); }
    static fp64_t const_ninf()  { static int64_double_t const c = { 0xfff0000000000000 }; return fp64_t(c.d); }
    static fp64_t const_indef() { static int64_double_t const c = { 0xfff8000000000000 }; return fp64_t(c.d); }

    //
    // static unary ops
    //
    static fp64_t abs(fp64_t const &src) { return fp64_t(std::abs(src.as_double())); }
    static fp64_t chs(fp64_t const &src) { return fp64_t(-src.as_double()); }
    static fp64_t sqrt(fp64_t const &src) { return fp64_t(std::sqrt(src.as_double())); }
    static fp64_t floor(fp64_t const &src) { return fp64_t(std::floor(src.as_double())); }
    static fp64_t ceil(fp64_t const &src) { return fp64_t(std::ceil(src.as_double())); }

    //
    // static transcendental ops
    //
    static fp64_t ldexp(fp64_t const &a, int32_t factor);

    //
    // x87 ops
    //
    static uint16_t x87_fxtract(fp64_t const &src, fp64_t &dst1, fp64_t &dst2);
    static uint16_t x87_fscale(fp64_t const &src1, fp64_t const &src2, fp64_t &dst);
    static uint16_t x87_fprem(fp64_t const &src1, fp64_t const &src2, fp64_t &dst);
    static uint16_t x87_fprem1(fp64_t const &src1, fp64_t const &src2, fp64_t &dst);
    static uint16_t x87_f2xm1(fp64_t const &src, fp64_t &dst);
    static uint16_t x87_fyl2x(fp64_t const &src1, fp64_t const &src2, fp64_t &dst);
    static uint16_t x87_fyl2xp1(fp64_t const &src1, fp64_t const &src2, fp64_t &dst);
    static uint16_t x87_fsin(fp64_t const &src, fp64_t &dst);
    static uint16_t x87_fcos(fp64_t const &src, fp64_t &dst);
    static uint16_t x87_fsincos(fp64_t const &src, fp64_t &dst1, fp64_t &dst2);
    static uint16_t x87_fptan(fp64_t const &src, fp64_t &dst1, fp64_t &dst2);
    static uint16_t x87_fpatan(fp64_t const &src1, fp64_t const &src2, fp64_t &dst);

    //
    // static misc ops
    //
    static fp64_t make_qnan(fp64_t const &src) { x87_assert(src.isnan()); fp64_t result(src); result.m_value.i |= 0x0008000000000000ull; return result; }
    static fp64_t from_fpbits32(uint32_t bits) { int32_float_t u = { bits }; return fp64_t(u.d); }
    static fp64_t from_fpbits64(uint64_t bits) { fp64_t result; result.m_value.i = bits; return result; }
    static bool samesign(fp64_t const &src1, fp64_t const &src2) { return (((src1.m_value.i ^ src2.m_value.i) & FP64_SIGN_MASK) == 0); }

protected:
    //
    // internal state
    //
    double_int64_t m_value;
};

//
// core math operations
//
inline fp64_t operator+(fp64_t const &a, fp64_t const &b) { return fp64_t(a.as_double() + b.as_double()); }
inline fp64_t operator-(fp64_t const &a, fp64_t const &b) { return fp64_t(a.as_double() - b.as_double()); }
inline fp64_t operator*(fp64_t const &a, fp64_t const &b) { return fp64_t(a.as_double() * b.as_double()); }
inline fp64_t operator/(fp64_t const &a, fp64_t const &b) { return fp64_t(a.as_double() / b.as_double()); }

//
// transcendental ops
//
inline fp64_t fp64_t::ldexp(fp64_t const &a, int32_t factor)
{
    int32_t exp = a.exponent() + factor + FP64_EXPONENT_BIAS;
    if (exp > 0 && exp < FP64_EXPONENT_MAX_BIASED)
        return fp64_t::from_fpbits64(a.m_value.i + (int64_t(factor) << FP64_EXPONENT_SHIFT));
    else if (exp <= 0)
        return const_zero();
    else
        return const_pinf();
}

}

#endif
