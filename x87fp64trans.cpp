//=========================================================
//  x87fp64trans.cpp
//
//  64-bit floating-point support (transcendental funcs)
//
//  The routines in this code are based on existing work from several sources:
//
//  * fxtract/fscale/f2xm1 implementations are by Aaron Giles
//  * fprem/fprem1 implementation was derived from softfloat (BSD 3-clause)
//  * fyl2x/fyl2xp1 implementations were derived from fdlibm (Sun license)
//  * fsin/fcos/fsincos/fptan/fpatan implementations were derived from the
//     Cephes math library (MIT license)
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
//=========================================================
//
// softfloat license: (BSD 3-clause)
//
// This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
// Package, Release 3e, by John R. Hauser.
//
// Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice,
//     this list of conditions, and the following disclaimer.
//
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions, and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//  3. Neither the name of the University nor the names of its contributors may
//     be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS", AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, ARE
// DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//==========================================================
//
// fdlibm license:
//
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
//
//==========================================================
//
// Cephes license: (MIT)
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//==========================================================

#include "x87fp64.h"
#include "x87fpext.h"

#include <cstdint>
#include <cmath>
#include <array>

namespace x87
{

#define IC(x) ((int32_t) x)
#define UC(x) ((uint32_t) x)


//===========================================================================
//
// indef / qnan
//
// Shared logic for returning common values.
//
//===========================================================================

//
// return the indefinite value; implicit sets the INVALID exception
//
static inline uint16_t indef(fp64_t &dst, uint16_t flags = 0)
{
    dst = fp64_t::const_indef();
    return flags | X87SW_INVALID_EX;
}

static inline uint16_t indef(fp64_t &dst1, fp64_t &dst2, uint16_t flags = 0)
{
    dst1 = dst2 = fp64_t::const_indef();
    return flags | X87SW_INVALID_EX;
}



//
// convert a NaN to a quiet NaN, setting the INVALID exception if the incoming
// NaN was a signaling NaN
//
static inline uint16_t qnan(fp64_t &dst, uint16_t flags, fp64_t const &src)
{
    x87_assert(src.isnan());

    // if it was a signaling NaNs, set the flag
    if (src.issnan())
        flags |= X87SW_INVALID_EX;

    // make a qNaN from src1
    dst = fp64_t::make_qnan(src);
    return flags;
}

static inline uint16_t qnan(fp64_t &dst1, fp64_t &dst2, uint16_t flags, fp64_t const &src)
{
    uint16_t result = qnan(dst1, flags, src);
    dst2 = dst1;
    return result;
}



//
// convert a NaN to a quiet NaN as above, but for functions with two source inputs;
// in this case, the INVALID exception is signaled if *either* source was a signaling
// NaN, and if both are indeed NaNs, then the result is the greater mantissa
//
static inline uint16_t qnan(fp64_t &dst, uint16_t flags, fp64_t const &src1, fp64_t const &src2)
{
    x87_assert(src1.isnan());

    // if either are signaling NaNs, set the flag
    if (src1.issnan() || src2.issnan())
        flags |= X87SW_INVALID_EX;

    // make a qNaN from src1
    dst = fp64_t::make_qnan(src1);

    // if both are NaNs, the "larger" one is reported
    if (src2.isnan())
    {
        auto man1 = src1.mantissa() & 0x7ffffffffffffull;
        auto man2 = src2.mantissa() & 0x7ffffffffffffull;
        if (man2 > man1 || (man2 == man1 && src1.sign()))
            dst = fp64_t::make_qnan(src2);
    }
    return flags;
}



//
// return a signed infinity
//
static inline uint16_t infinity(fp64_t &dst, uint16_t flags, uint8_t sign)
{
    dst = sign ? fp64_t::const_ninf() : fp64_t::const_pinf();
    return flags;
}

//
// return a signed zero
//
static inline uint16_t zero(fp64_t &dst, uint16_t flags, uint8_t sign)
{
    dst = sign ? fp64_t::const_nzero() : fp64_t::const_zero();
    return flags;
}



//
// Debug helpers for printing out intermediate values
//
#if 0
#include "fmt.h"
void print_val(char const *name, fp64_t const &val)
{
    fmt::print("{} = {:c}{:013X}e{:+05} ({:+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.exponent(), val.as_double());
}
void print_val(char const *name, fp80_t const &val)
{
    fmt::print("{} = {:c}{:016X}e{:+05} ({:+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.exponent(), val.as_double());
}
void print_val(char const *name, fpext64_t const &val)
{
    fmt::print("{} = {:c}{:016X}e{:+05} ({:+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.exponent(), val.as_fp64().as_double());
}
void print_val(char const *name, fpext96_t const &val)
{
    fmt::print("{} = {:c}{:016X}`{:08X}e{:+05} ({:+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.extend(), val.exponent(), val.as_fp64().as_double());
}
#else
#define print_val(x, y) do { } while (0)
#endif



//===========================================================================
//
// poly_eval / poly1_eval
//
// Evaluate a polynomial by iterating over an array of terms. Derived from
// the Cephes math library, found here: https://netlib.org/cephes/
//
//===========================================================================

//
// polynomial evaluator of the form
//    P[0] x^n  +  P[1] x^(n-1)  +  ...  +  P[n]
//
template<typename FpType, size_t Count>
FpType poly_eval(FpType const &x, std::array<FpType, Count> const &terms)
{
    int index = 0;
    FpType dst = terms[0];
    for (int index = 1; index < Count; index++)
        dst = dst * x + terms[index];
    return dst;
}

//
// polynomial evalutaor of the form
//    x^n  +  P[0] x^(n-1)  +  P[1] x^(n-2)  +  ...  +  P[n]
//
template<typename FpType, size_t Count>
FpType poly1_eval(FpType const &x, std::array<FpType, Count> const &terms)
{
    int index = 0;
    FpType dst = x + terms[0];
    for (int index = 1; index < Count; index++)
        dst = dst * x + terms[index];
    return dst;
}



//===========================================================================
//
// x87_fxtract
//
// Extract mantissa and exponent.
//
// fxtract(64) results:
//      128684 matches [200.00%]
//
//===========================================================================

uint16_t fp64_t::x87_fxtract(fp64_t const &src, fp64_t &dst1, fp64_t &dst2)
{
    // denorms always set the flag
    uint16_t flags = 0;
    if (src.isdenorm())
        flags |= X87SW_DENORM_EX;

    // handle special cases
    if (src.ismaxexp())
        goto special;

    // zeros are special
    if (src.iszero())
        goto zero;

    {
        // convert to 80-bit to normalize denorms
        fpext64_t esrc(src);

        // create a new number with 0x3ff exponent, same sign
        dst1 = fp64_t::from_fpbits64((uint64_t(esrc.sign()) << FP64_SIGN_SHIFT) |
                                     (0x3ffull << FP64_EXPONENT_SHIFT) |
                                     ((esrc.mantissa() >> 11) & FP64_MANTISSA_MASK));

        // mantissa in second output
        dst2 = fp64_t(esrc.exponent());
        return flags;
    }

special:
    // NaNs in, NaNs Out
    if (src.isnan())
        return qnan(dst1, dst2, flags, src);
    // infinities return themselves as significand, +infinity as exponent
    dst1 = src;
    dst2 = fp64_t::const_pinf();
    return flags;

zero:
    // zeroes return themselves as significant, -infinity as exponent
    dst1 = src;
    dst2 = fp64_t::const_ninf();
    return X87SW_DIVZERO_EX;
}



//===========================================================================
//
// x87_fscale
//
// Scale value by exponent.
//
// fscale(64) results:
//    165611161 matches [100.00%]
//
//===========================================================================

uint16_t fp64_t::x87_fscale(fp64_t const &src1, fp64_t const &src2, fp64_t &dst)
{
    // denorms always set the flag
    uint16_t flags = 0;
    if (src1.isdenorm() || src2.isdenorm())
        flags |= X87SW_DENORM_EX;

    // handle special cases
    if (src1.ismaxexp())
        goto special1;
    if (src2.ismaxexp())
        goto special2;

    // handle zero significand
    if (src1.iszero())
        goto zero;

    {
        fpext64_t esrc1(src1);

        // extract the exponent
        fp64_t exp;
        if (src2.sign() == 0)
            exp = fp64_t::floor(src2);
        else
            exp = fp64_t::ceil(src2);

        // handle overly large values before converting to int
        if (exp.as_double() >= 32768.0)
            goto overflow;
        if (exp.as_double() <= -32768.0)
            goto underflow;

        // convert to int and check for 0
        int iexp = exp.as_int32();
        if (iexp == 0)
            goto zero;

        // check for underflow
        int newexp = esrc1.exponent() + iexp;
        if (newexp <= -16394)
        {
            // mantissa is considered when setting bits
            auto mantissa = src1.mantissa();
            int thresh = -16394 - ((mantissa == 0) ? 52 : count_trailing_zeros64(mantissa));
            if (newexp <= thresh)
                goto underflow;
        }
        if (newexp <= int(0 - FP64_EXPONENT_BIAS - FP64_MANTISSA_BITS))
            goto underflow2;

        // check for overflow
        if (newexp >= 16384)
            goto overflow;
        if (newexp >= int(FP64_EXPONENT_MAX_BIASED - FP64_EXPONENT_BIAS))
            goto overflow2;

        // all other cases
        dst = fpext64_t::ldexp(esrc1, exp.as_int32()).as_fp64();
        return flags;
    }

underflow:
    // some underflows set the flags, some don't
    flags |= X87CW_MASK_UNDERFLOW_EX | X87CW_MASK_PRECISION_EX;
underflow2:
    return zero(dst, flags, src1.sign());

overflow:
    // some overflows set the flags, some don't
    flags |= X87CW_MASK_OVERFLOW_EX | X87CW_MASK_PRECISION_EX;
overflow2:
    dst = fp64_t::from_fpbits64(src1.sign() ? 0xffefffffffffffffull : 0x7fefffffffffffffull);
    return flags;

special1:
    // NaNs in, NaNs out
    if (src1.isnan())
        return qnan(dst, flags, src1, src2);
    if (src2.isnan())
        return qnan(dst, flags, src2);
    // otherwise, return infinity of our sign
    return infinity(dst, flags, src1.sign());

special2:
    // if NaN, pass it through
    if (src2.isnan())
        return qnan(dst, flags, src2);
    // if we're infinity and src1 is 0, produce a NaN
    if (src1.iszero())
        return indef(dst, flags);
    // otherwise we return an infinity based on the sign of src1
    return infinity(dst, flags, src1.sign());

zero:
    dst = src1;
    return flags;
}



//===========================================================================
//
// x87_fprem / x87_fprem1
//
// Compute fprem(x, y). This code is based off the softfloat 64-bit
// remainder function, with added logic to get an accurate low 3 bits of
// the quotient and handle the fprem1 case.
//
// fprem(64) results:
//    165611161 matches [100.00%]
//
// fprem1(64) results:
//    165611161 matches [100.00%]
//
//===========================================================================

template<bool Rem1>
static uint16_t x87_fprem_core(fp64_t const &src1, fp64_t const &src2, fp64_t &dst)
{
    // denorms always set the flag
    uint16_t flags = 0;
    if (src1.isdenorm() || src2.isdenorm())
        flags |= X87SW_DENORM_EX;

    // handle NaNs and infinities
    if (src1.ismaxexp())
        goto special1;
    if (src2.ismaxexp())
        goto special2;

    // handle divide by zero
    if (src2.iszero())
        return indef(dst, flags);

    {
        fpext64_t esrc1(src1);
        fpext64_t esrc2(src2);
        int32_t dexp = esrc1.exponent() - esrc2.exponent();
        uint64_t rem = esrc1.mantissa() >> 2;
        uint64_t sigb = esrc2.mantissa() >> 2;

        auto factor = (dexp > 63) ? ((dexp - 32) / 32) * 32 : 0;
        dexp -= factor;

        uint64_t altrem;
        uint64_t q;
        if (dexp < 1)
        {
            if (dexp < -1)
            {
                x87_assert(factor == 0);
                dst = src1;
                return flags;
            }
            q = 0;
            if (dexp != 0)
                rem >>= 1;
            else if (sigb <= rem)
                rem -= sigb, q = 1;
        }
        else
        {
            uint32_t recip32 = uint32_t(0x7FFFFFFFFFFFFFFFull / uint32_t(sigb >> 30));
            dexp -= 30;
            uint64_t q64;
            uint64_t qt = 0;
            while (1)
            {
                q64 = uint64_t(uint32_t(rem >> 32)) * recip32;
                if (dexp < 0)
                    break;
                q = (q64 + 0x80000000) >> 32;
                rem <<= 29;
                rem -= q * sigb;
                if (int64_t(rem) < 0)
                {
                    rem += sigb;
                    q--;
                }
                qt = (qt << 29) + q;
                dexp -= 29;
            }
            q = uint32_t(q64 >> 32) >> (~dexp & 31);
            rem = (rem << (dexp + 30)) - q * sigb;
            q = (qt << (dexp + 30)) + q;
            if (int64_t(rem) < 0)
            {
                altrem = rem + sigb;
                goto select;
            }
        }
        do
        {
            altrem = rem;
            q++;
            rem -= sigb;
        } while (int64_t(rem) >= 0);
    select:
        rem = altrem << 2;
        int shift = count_leading_zeros64(rem);
        rem <<= shift;
        fpext64_t res = fpext64_t(rem, 0, (rem == 0) ? fpext64_t::EXPONENT_MIN : esrc2.exponent() - shift + factor, src1.sign());

        // fprem1 returns results from -src2/2..src2/2 instead of 0..src2
        if (Rem1 && factor == 0)
        {
            // if the result is > src/2 or == src/2 and odd, take back one src2
            if (res.exponent() == esrc2.exponent() ||
                (res.exponent() == esrc2.exponent() - 1 && (rem > esrc2.mantissa() || (rem == esrc2.mantissa() && (q & 1) == 0))))
            {
                esrc2.abs();
                if (res.sign() == 0)
                    res -= esrc2;
                else
                    res += esrc2;
                q++;
            }
        }

        dst = res.as_fp64();
        if (factor != 0)
            return flags | X87SW_C2;
        q--;
        return flags | ((q & 1) << X87SW_C1_BIT) | ((q & 2) << (X87SW_C3_BIT - 1)) | ((q & 4) << (X87SW_C0_BIT - 2));
    }

special1:
    // NaNs in, NaNs out
    if (src1.isnan())
        return qnan(dst, flags, src1, src2);
    if (src2.isnan())
        return qnan(dst, flags, src2);
    // infinity is invalid
    return indef(dst, flags);

special2:
    // NaNs in, NaNs out
    if (src2.isnan())
        return qnan(dst, flags, src2);
    // infinite src2 just returns src1
    dst = src1;
    return flags;
}

uint16_t fp64_t::x87_fprem(fp64_t const &src1, fp64_t const &src2, fp64_t &dst)
{
    return x87_fprem_core<false>(src1, src2, dst);
}

uint16_t fp64_t::x87_fprem1(fp64_t const &src1, fp64_t const &src2, fp64_t &dst)
{
    return x87_fprem_core<true>(src1, src2, dst);
}



//===========================================================================
//
// x87_f2xm1
//
// Compute 2^(x - 1). This code is written by Aaron Giles based on the
// algorithm in the paper here:
//
// https://www.researchgate.net/publication/3612479_The_K5_transcendental_functions
//
// f2xm1(64) results:
//        64332 matches [99.98%]
//           10 off by 1 bits [0.02%]
//
//===========================================================================

template<bool Debug>
static uint16_t x87_f2xm1_core(fp64_t const &src, fp64_t &dst)
{
    // special case values outside of defined range
    auto exponent = src.exponent();
    if (exponent >= 0)
        goto special;

    // anything small ends up with a G of 0 and an H == x*ln2
    // this also allows us to avoid dealing with denormals
    if (exponent <= -1000)
        goto tiny;

    static const int LOG_R = 4;
    static const int R = 1 << LOG_R;
    static const int TABLE_SIZE = 2 * R + 1;
    static const int TAYLOR_TERMS = 8;

    static std::array<fpext64_t, TABLE_SIZE> const s_table_g =
    {
        fpext64_t(0x8000000000000000ull, 0x00000000, -1, 1),    // 2^(-16/16) = -0.5l,
        fpext64_t(0xf4aa7930676f09d6ull, 0x746d48e8, -2, 1),    // 2^(-15/16) = -0.47786310878629307983901676063004l,
        fpext64_t(0xe8d47c382ae85232ull, 0x08373af1, -2, 1),    // 2^(-14/16) = -0.45474613366737117039649467211965l,
        fpext64_t(0xdc785918a9dc7993ull, 0xe0524e3f, -2, 1),    // 2^(-13/16) = -0.43060568262165417314808485807924l,
        fpext64_t(0xcf901f5ce48ead21ull, 0x72a5b9d0, -2, 1),    // 2^(-12/16) = -0.40539644249863946664125001471976l,
        fpext64_t(0xc2159b3edcbddca4ull, 0xbeddc1ec, -2, 1),    // 2^(-11/16) = -0.3790710939632579757031612656367l,
        fpext64_t(0xb40252ac9d5d8e2bull, 0xc685013c, -2, 1),    // 2^(-10/16) = -0.35158022267449516703312294110377l,
        fpext64_t(0xa54f822b7abd6a73ull, 0x6cfeae6e, -2, 1),    // 2^( -9/16) = -0.32287222653155363585099262992965l,
        fpext64_t(0x95f619980c4336f7ull, 0x4d04ec99, -2, 1),    // 2^( -8/16) = -0.29289321881345247559915563789515l,
        fpext64_t(0x85eeb8c14fe79282ull, 0xaefdc093, -2, 1),    // 2^( -7/16) = -0.26158692703025034430654625981298l,
        fpext64_t(0xea6357baabe4948bull, 0x0754bcda, -3, 1),    // 2^( -6/16) = -0.22889458729602958819385406895463l,
        fpext64_t(0xc76dcfab81edfc70ull, 0x7729f1c2, -3, 1),    // 2^( -5/16) = -0.1947548340253728459102396663213l,
        fpext64_t(0xa2ec0cd4a58a542full, 0x1965d11a, -3, 1),    // 2^( -4/16) = -0.15910358474628545696887452376679l,
        fpext64_t(0xf999089eab58f777ull, 0xcd3b57dc, -4, 1),    // 2^( -3/16) = -0.12187391981335025844391969031234l,
        fpext64_t(0xa9f9c8c116de3689ull, 0x7e945264, -4, 1),    // 2^( -2/16) = -0.08299595679532876825645840520586l,
        fpext64_t(0xada82eadb7933d38ull, 0x462f3851, -5, 1),    // 2^( -1/16) = -0.04239671930142635306369436485208l,
        fpext64_t(0x0000000000000000ull, 0x00000000, fpext96_t::EXPONENT_MIN, 0), // 0
        fpext64_t(0xb5586cf9890f6298ull, 0xb92b7184, -5, 0),    // 2^( +1/16) = 0.04427378242741384032196647873993l,
        fpext64_t(0xb95c1e3ea8bd6e6full, 0xbe462876, -4, 0),    // 2^( +2/16) = 0.09050773266525765920701065576071l,
        fpext64_t(0x8e1e9b9d588e19b0ull, 0x7eb6c705, -3, 0),    // 2^( +3/16) = 0.13878863475669165370383028384151l,
        fpext64_t(0xc1bf828c6dc54b7aull, 0x356918c1, -3, 0),    // 2^( +4/16) = 0.18920711500272106671749997056048l,
        fpext64_t(0xf7a993048d088d6dull, 0x0488f84f, -3, 0),    // 2^( +5/16) = 0.2418578120734840485936774687266l,
        fpext64_t(0x97fb5aa6c544e3a8ull, 0x72f5fd88, -2, 0),    // 2^( +6/16) = 0.29683955465100966593375411779245l,
        fpext64_t(0xb560fba90a852b19ull, 0x2602a324, -2, 0),    // 2^( +7/16) = 0.3542555469368927282980147401407l,
        fpext64_t(0xd413cccfe7799211ull, 0x65f626ce, -2, 0),    // 2^( +8/16) = 0.4142135623730950488016887242097l,
        fpext64_t(0xf4228e7d6030dafaull, 0xa2047eda, -2, 0),    // 2^( +9/16) = 0.47682614593949931138690748037405l,
        fpext64_t(0x8ace5422aa0db5baull, 0x7c55a193, -1, 0),    // 2^(+10/16) = 0.54221082540794082361229186209073l,
        fpext64_t(0x9c49182a3f0901c7ull, 0xc46b071f, -1, 0),    // 2^(+11/16) = 0.6104903319492543081795206673574l,
        fpext64_t(0xae89f995ad3ad5e8ull, 0x734d1773, -1, 0),    // 2^(+12/16) = 0.68179283050742908606225095246643l,
        fpext64_t(0xc199bdd85529c222ull, 0x0cb12a09, -1, 0),    // 2^(+13/16) = 0.75625216037329948311216061937531l,
        fpext64_t(0xd5818dcfba48725dull, 0xa05aeb67, -1, 0),    // 2^(+14/16) = 0.83400808640934246348708318958829l,
        fpext64_t(0xea4afa2a490d9858ull, 0xf73a18f6, -1, 0),    // 2^(+15/16) = 0.91520656139714729387261127029583l,
        fpext64_t(0x8000000000000000ull, 0x00000000,  0, 0)     // 2^(+16/16) = 1.0
    };
    static std::array<fp64_t, TABLE_SIZE> const s_table_u =
    {
        -16.0/16.0,
        -15.0/16.0,
        -14.0/16.0,
        -13.0/16.0,
        -12.0/16.0,
        -11.0/16.0,
        -10.0/16.0,
         -9.0/16.0,
         -8.0/16.0,
         -7.0/16.0,
         -6.0/16.0,
         -5.0/16.0,
         -4.0/16.0,
         -3.0/16.0,
         -2.0/16.0,
         -1.0/16.0,
          0.0/16.0,
          1.0/16.0,
          2.0/16.0,
          3.0/16.0,
          4.0/16.0,
          5.0/16.0,
          6.0/16.0,
          7.0/16.0,
          8.0/16.0,
          9.0/16.0,
         10.0/16.0,
         11.0/16.0,
         12.0/16.0,
         13.0/16.0,
         14.0/16.0,
         15.0/16.0,
         16.0/16.0
    };
    static std::array<fp64_t, TAYLOR_TERMS - 1> const s_taylor_coeff =
    {
        8.0,
        8.0*7,
        8.0*7*6,
        8.0*7*6*5,
        8.0*7*6*5*4,
        8.0*7*6*5*4*3,
        8.0*7*6*5*4*3*2
    };
    static fp64_t const s_taylor_factorial_inv =
        1.0 / fp64_t(8*7*6*5*4*3*2);  // 1.0/8!

    {
        // round x to the nearest multiple of 1/R by looking at the high bits of the mantissa
        int32_t g_index = 0;

        // anything smaller than -LOG_R - 1 will round to 0, so only do this if above
        if (exponent >= -LOG_R - 1)
        {
            // shift mantissa down (after adding explicit 1) so we just have LOG_R + 1 bits
            auto mantissa = src.mantissa() | (FP64_MANTISSA_MASK + 1);
            g_index = int32_t(mantissa >> (FP64_EXPONENT_SHIFT - LOG_R - exponent - 1));

            // round by adding LSB and shifting to get LOG_R bits
            g_index = (g_index >> 1) + (g_index & 1);

            // if negative, use a negative index
            if (src.sign() != 0)
                g_index = -g_index;
        }

        // compute v = delta from table entry
        fp64_t v = src - s_table_u[g_index + R];

        // multiply v by ln(2) so we can use the e^x Taylor series; do this in
        // extended precision
        fpext64_t w = fpext64_t(v) * fpext64_t::ln2;
        if (Debug) print_val("w", w);

        // Taylor series: this can be done in lower precision; start with h = w + coeff[0]
        fp64_t w64 = w.as_fp64();
        fp64_t h64 = w64 + s_taylor_coeff[0];
        if (Debug) print_val("h1", h64);

        // now compute h = h * w + coeff[term] for terms up through 7
        for (int term = 1; term < TAYLOR_TERMS - 2; term++)
        {
            h64 = h64 * w64 + s_taylor_coeff[term];
            if (Debug) print_val("hn", h64);
        }

        // final term is just times w^2
        h64 *= w64 * w64;
        if (Debug) print_val("h2", h64);

        // then divide by 9!
        h64 = h64 * s_taylor_factorial_inv;

        // back to extended precision for final result; add w for final h value
        fpext64_t h(h64);
        h += w;
        if (Debug) print_val("h3", h);

        // retrieve g from the table
        fpext64_t g = s_table_g[g_index + R];
        if (Debug) print_val("g", g);

        // return g * h + g + h
        if (Debug)
        {
            fpext64_t res = g * h;
            if (Debug) print_val("res", res);
            res += g;
            if (Debug) print_val("res", res);
            res += h;
            if (Debug) print_val("res", res);
        }
        dst = (g * h + g + h).as_fp64();
        return X87SW_PRECISION_EX;
    }

special:
    // return -0.5 for -1
    if (src.as_fpbits64() == 0xbff0000000000000ull)
    {
        dst = fp64_t::from_fpbits64(0xbfe0000000000000ull);
        return X87SW_PRECISION_EX;
    }
    // max exponent cases
    if (src.ismaxexp())
    {
        // NaNs in, NaNs out
        if (src.isnan())
            return qnan(dst, 0, src);
        // return -1 for -inf and set no flags
        if (src.isninf())
        {
            dst = fp64_t::from_fpbits64(0xbff0000000000000ull);
            return 0;
        }
        // return +inf for +inf and set no flags
        if (src.isinf())
            return infinity(dst, 0, 0);
    }
    // for 0 or out-of-range values, just return x
    dst = src;
    return src.iszero() ? 0 : X87SW_PRECISION_EX;

tiny:
    // special case zero
    if (src.iszero())
    {
        dst = src;
        return 0;
    }
    // denorms and other tiny values reduce to a simple multiply
    dst = (fpext64_t(src) * fpext64_t::ln2).as_fp64();
    if (src.isdenorm())
        return X87SW_PRECISION_EX | X87SW_DENORM_EX;
    return X87SW_PRECISION_EX;
}

uint16_t fp64_t::x87_f2xm1(fp64_t const &src, fp64_t &dst)
{
    return x87_f2xm1_core<false>(src, dst);
}



//===========================================================================
//
// x87_fyl2x
//
// Compute y * log2(x). This code is based off the __ieee754_log2
// implementation found in libm.
//
// fyl2x(64) results:
//    162919171 matches [98.37%]
//       186112 off by 1 bits [0.11%]
//      2505878 pseudo infinities [1.51%]
//
//===========================================================================

uint16_t fp64_t::x87_fyl2x(fp64_t const &src1, fp64_t const &src2, fp64_t &dst)
{
    // denorm flag is set regardless
    uint16_t flags = 0;
    if (src1.isdenorm() || src2.isdenorm())
        flags |= X87SW_DENORM_EX;

    // special case max exponent values
    if (src1.ismaxexp())
        goto specialx;
    if (src2.ismaxexp())
        goto specialy;

    // log of negative values produces indefinite
    if (src1.sign() != 0)
        goto invalid;

    // log of 0 is infinity, unless multiplying by 0
    if (src1.iszero())
        goto log0;

    // shortcut if multiplying by 0
    if (src2.iszero())
        goto times0;

    {
        static const fp64_t two54 = fp64_t::from_fpbits64(0x4350000000000000ull); //1.80143985094819840000e+16;
        static const fp64_t Lg1   = fp64_t::from_fpbits64(0x3FE5555555555593ull); //6.666666666666735130e-01;
        static const fp64_t Lg2   = fp64_t::from_fpbits64(0x3FD999999997FA04ull); //3.999999999940941908e-01;
        static const fp64_t Lg3   = fp64_t::from_fpbits64(0x3FD2492494229359ull); //2.857142874366239149e-01;
        static const fp64_t Lg4   = fp64_t::from_fpbits64(0x3FCC71C51D8E78AFull); //2.222219843214978396e-01;
        static const fp64_t Lg5   = fp64_t::from_fpbits64(0x3FC7466496CB03DEull); //1.818357216161805012e-01;
        static const fp64_t Lg6   = fp64_t::from_fpbits64(0x3FC39A09D078C69Full); //1.531383769920937332e-01;
        static const fp64_t Lg7   = fp64_t::from_fpbits64(0x3FC2F112DF3E5244ull); //1.479819860511658591e-01;

        // accuracy/speed results:
        //   fpext52_t: 125778632(0) / 37316356(1) / 2898(2) / 16(3) / 10(4) / 6(5) / 10295(exp), 0.11 ticks
        //   fpext64_t: 162919171(0) / 186112(1), 0.23 ticks
        //   fpext96_t: 162934641(0) / 170642(1), 0.32 ticks
        using fpext_t = fpext64_t;

        static fpext_t const invln2(0xb8aa3b295c17f0bbull, 0xbe87fed0,  0, 0);
        fpext_t src280(src2);
        fpext_t src2invln2 = src280 * invln2;

        if (src1 != fp64_t::const_one())
            flags |= X87SW_PRECISION_EX;

        uint64_t rawsrc = src1.as_fpbits64();
        int32_t hx = int32_t(rawsrc >> 32);

        // handle denorms
        int32_t k = 0;
        fp64_t x = src1;
        if (x.isdenorm())
        {
            k -= 54;
            x *= two54;
            rawsrc = x.as_fpbits64();
            hx = int32_t(rawsrc >> 32);
        }
        k += x.exponent();

        hx &= IC(0x000fffff);
        int32_t i = (hx + IC(0x95f64)) & IC(0x100000);
        x = fp64_t::from_fpbits64((rawsrc & FP64_MANTISSA_MASK) | (uint64_t(i ^ IC(0x3ff00000)) << 32));
        k += (i >> 20);
        fpext_t dk80 = fpext_t(k) * src280;
        fp64_t f = x - 1.0;
        if ((IC(0x000fffff) & (2 + hx)) < 3)
        {                                       // |f| < 2**-20
            if (f == fp64_t::const_zero())
            {
                dst = dk80.as_fp64();
                return flags;
            }
            fp64_t R = f * f * (0.5 - 0.33333333333333333 * f);
            dst = (dk80 - fpext_t(R - f) * src2invln2).as_fp64();
            return flags;
        }
        fp64_t s = f / (2.0 + f);
        fp64_t z = s * s;
        i = hx - IC(0x6147a);
        fp64_t w = z * z;
        int32_t j = IC(0x6b851) - hx;
        fp64_t t1 = w * (Lg2 + w * (Lg4 + w * Lg6));
        fp64_t t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
        i |= j;
        fp64_t R = t2 + t1;
        if (i > 0)
        {
            fp64_t hfsq = 0.5 * f * f;
            dst = (dk80 - fpext_t((hfsq - (s * (hfsq + R))) - f) * src2invln2).as_fp64();
        }
        else
            dst = (dk80 - fpext_t((s * (f - R)) - f) * src2invln2).as_fp64();
        return flags;
    }

specialx:
    // NaNs in, NaNs out
    if (src1.isnan())
        return qnan(dst, flags, src1, src2);
    if (src2.isnan())
        return qnan(dst, flags, src2);
    // log negative infinity or infinity times zero is indefinite
    if (src1.sign() || src2.iszero())
        return indef(dst, flags);
    // log of infinity is infinity of appropriate sign
    return infinity(dst, flags, src2.sign());

specialy:
    // handle special cases for multiplier
    if (src2.isnan())
        return qnan(dst, flags, src2);
    // log of negative values or 1.0 produces indefinite
    if (src1.sign() != 0 || src1 == fp64_t::const_one())
        return indef(dst, flags);
    // multiplying infinity produces infinity of the appropriate sign
    return infinity(dst, flags, (src1.exponent() < 0) ^ src2.sign());

log0:
    // log of 0 is -infinity, unless multiplying by 0
    if (src2.iszero())
        return indef(dst, flags);
    return infinity(dst, flags | X87SW_DIVZERO_EX, src2.sign() ^ 1);

times0:
    return zero(dst, flags, src2.sign() ^ (src1.exponent() < 0));

invalid:
    return indef(dst, flags);
}



//===========================================================================
//
// x87_fyl2xp1
//
// Compute y * log2(x + 1). This code is based off the __ieee754_log1p
// implementation found in libm.
//
// fyl2xp1(64) results:
//    142794485 matches [86.22%]
//     21488338 off by 1 bits [12.98%]
//          126 off by 2 bits [0.00%]
//      1328212 pseudo infinities [0.80%]
//
//===========================================================================

uint16_t fp64_t::x87_fyl2xp1(fp64_t const &src1, fp64_t const &src2, fp64_t &dst)
{
    // denorm flag is set regardless
    uint16_t flags = 0;
    if (src1.isdenorm() || src2.isdenorm())
        flags |= X87SW_DENORM_EX;

    // special case max exponent values
    if (src1.ismaxexp())
        goto specialx;
    if (src2.ismaxexp())
        goto specialy;

    // log of 0 is special
    if (src1 == fp64_t::const_none())
        goto logn1;

    // out of bound values return src1
    if (src1 < fp64_t::const_none())
        goto oob;

    // shortcut if multiplying by 0
    if (src2.iszero())
        goto times0;

    static const fp64_t ln2_hi = fp64_t::from_fpbits64(0x3fe62e42fee00000);     // 6.93147180369123816490e-01;
    static const fp64_t ln2_lo = fp64_t::from_fpbits64(0x3dea39ef35793c76);     // 1.90821492927058770002e-10;
    static const fp64_t two54 =  fp64_t::from_fpbits64(0x4350000000000000);     // 1.80143985094819840000e+16;
    static const std::array<fp64_t, 8> Lp =
    {
        fp64_t::from_fpbits64(0x0000000000000000),    // 0.0,
        fp64_t::from_fpbits64(0x3FE5555555555593),    // 6.666666666666735130e-01,
        fp64_t::from_fpbits64(0x3FD999999997FA04),    // 3.999999999940941908e-01,
        fp64_t::from_fpbits64(0x3FD2492494229359),    // 2.857142874366239149e-01,
        fp64_t::from_fpbits64(0x3FCC71C51D8E78AF),    // 2.222219843214978396e-01,
        fp64_t::from_fpbits64(0x3FC7466496CB03DE),    // 1.818357216161805012e-01,
        fp64_t::from_fpbits64(0x3FC39A09D078C69F),    // 1.531383769920937332e-01,
        fp64_t::from_fpbits64(0x3FC2F112DF3E5244),    // 1.479819860511658591e-01
    };

    {
        // accuracy/speed results:
        //   fpext52_t: fails
        //   fpext64_t: 142794485(0) / 21488338(1) / 126(2), 0.18 ticks
        //   fpext96_t: 142802124(0) / 21480699(1) / 126(2), 0.30 ticks
        using fpext_t = fpext64_t;

        static fpext_t const invln2(0xb8aa3b295c17f0bbull, 0xbe87fed0,  0, 0);
        fpext_t src2invln2 = fpext_t(src2) * invln2;

        if (!src1.iszero())
            flags |= X87SW_PRECISION_EX;

        int32_t hx = int32_t(src1.as_fpbits64() >> 32);
        int32_t ax = hx & IC(0x7fffffff);

        int32_t k = 1;
        fp64_t f, c;
        int32_t hu;
        if (hx < IC(0x3FDA827A))
        {									/* x < 0.41422  */
            if (ax < IC(0x3e200000))
            {								/* |x| < 2**-29 */
                if (ax < IC(0x3c900000))	/* |x| < 2**-54 */
                    dst = (fpext_t(src1) * src2invln2).as_fp64();
                else
                    dst = (fpext_t(src1 - src1 * src1 * 0.5) * src2invln2).as_fp64();
                return flags;
            }
            if (hx > 0 || hx <= IC(0xbfd2bec3))
            {
                k = 0;
                f = src1;
                hu = 1;
                c = 0;
            }								/* -0.2929<x<0.41422 */
        }
        if (k != 0)
        {
            fp64_t u;
            if (hx < IC(0x43400000))
            {
                u = 1.0 + src1;
                hu = int32_t(u.as_fpbits64() >> 32);
                k = u.exponent();
                c = (k > 0) ? 1.0 - (u - src1) : src1 - (u - 1.0);	/* correction term */
                c /= u;
            }
            else
            {
                u = src1;
                hu = int32_t(u.as_fpbits64() >> 32);
                k = u.exponent();
                c = 0;
            }
            hu &= IC(0x000fffff);
            if (hu < IC(0x6a09e))
            {
                u = fp64_t::from_fpbits64(u.mantissa() | 0x3ff0000000000000ull); /* normalize u */
            }
            else
            {
                k += 1;
                u = fp64_t::from_fpbits64(u.mantissa() | 0x3fe0000000000000ull); /* normalize u/2 */
                hu = (IC(0x00100000) - hu) >> 2;
            }
            f = u - 1.0;
        }
        fp64_t hfsq = 0.5 * f * f;
        if (hu == 0)
        {									/* |f| < 2**-20 */
            if (f == fp64_t::const_zero())
            {
                if (k == 0)
                    dst = fp64_t::const_zero();
                else
                {
                    c += k * ln2_lo;
                    dst = (fpext_t(k * ln2_hi + c) * src2invln2).as_fp64();
                }
                return flags;
            }
            fp64_t R = hfsq * (1.0 - 0.66666666666666666 * f);
            if (k == 0)
                dst = (fpext_t(f - R) * src2invln2).as_fp64();
            else
                dst = (fpext_t(k * ln2_hi - ((R - (k * ln2_lo + c)) - f)) * src2invln2).as_fp64();
            return flags;
        }
        fp64_t s = f / (2.0 + f);
        fp64_t z = s * s;
        fp64_t R1 = z * Lp[1];
        fp64_t z2 = z * z;
        fp64_t R2 = Lp[2] + z * Lp[3];
        fp64_t z4 = z2 * z2;
        fp64_t R3 = Lp[4] + z * Lp[5];
        fp64_t z6 = z4 * z2;
        fp64_t R4 = Lp[6] + z * Lp[7];
        fp64_t R = R1 + z2 * R2 + z4 * R3 + z6 * R4;
        if (k == 0)
            dst = (fpext_t(f - (hfsq - s * (hfsq + R))) * src2invln2).as_fp64();
        else
            dst = (fpext_t(k * ln2_hi - ((hfsq - (s * (hfsq + R) + (k * ln2_lo + c))) - f)) * src2invln2).as_fp64();
        return flags;
    }

specialx:
    // NaNs in, NaNs out
    if (src1.isnan())
        return qnan(dst, flags, src1, src2);
    if (src2.isnan())
        return qnan(dst, flags, src2);
    // log negative infinity or infinity times zero is indefinite
    if ((src1.sign() != 0 && src1.exponent() >= 0) || src2.iszero())
        return indef(dst, flags);
    // log of infinity is infinity of appropriate sign
    return infinity(dst, flags, src2.sign());

specialy:
    // handle special cases for multiplier
    if (src2.isnan())
        return qnan(dst, flags, src2);
    // log of values < -1 produces indefinite
    if (src1.iszero() || src1 == fp64_t::const_none())
        return indef(dst, flags);
    // multiplying infinity produces infinity of the appropriate sign
    return infinity(dst, flags, src1.sign() ^ src2.sign());

logn1:
    // log of 0 is infinity, unless multiplying by 0
    if (src2.iszero())
        return indef(dst, flags);
    return infinity(dst, flags, src2.sign());

times0:
    return zero(dst, flags, src2.sign() ^ src1.sign());

oob:
    if (src2.iszero())
        return zero(dst, flags, src2.sign() ^ 1);
    dst = src1;
    return flags | X87SW_PRECISION_EX;
}



//===========================================================================
//
// reduce_trig
//
// Reduce a trigonometric parameter down to a pi/2 quadrant, returning the
// delta in 'delta' and the quadrant index as the return value.
//
//===========================================================================

//
// Intel uses this 66-bit approximation of PI; these constants are computed using
// that as our value so we don't drift from the real x87 behavior
// See: https://www.intel.com/content/dam/develop/external/us/en/documents/x87trigonometricinstructionsvsmathfunctions.pdf
//
template<typename FpType>
static uint32_t reduce_trig(fp64_t src, FpType &delta)
{
    // convert src to FpType in delta
    src = fp64_t::abs(src);
    delta = FpType(src);

    // if < pi/4, return as-is
    if (src < 0.7853981633974483096)
        return 0;

    // srcext is 1.63*2^srcexp
    uint64_t srcman = delta.mantissa();
    int32_t srcexp = delta.exponent();
    x87_assert(int64_t(srcman) < 0);
    x87_assert(srcexp >= -1);
    x87_assert(srcexp < 63);

    // multiply by invpio4, which is a 1.127 value; final result is
    // 2.190; since invpio4 has a 0 exponent, we don't need to adjust
    // for that
    static const uint64_t INV_PIO4_HI = 0xa2f9836e4e44152aull;
    static const uint64_t INV_PIO4_LO = 0x00062bc40da28000ull;
    auto [divmid, divhi] = multiply_64x64(srcman, INV_PIO4_HI);
    auto [divlo, hitemp] = multiply_64x64(srcman, INV_PIO4_LO);
    divmid += hitemp;
    if (divmid < hitemp)
        divhi++;

    // now find the floor; srcexp can only be 62 bits at most, so
    // we only need to worry about masking divhi; the lower parts
    // can be assumed to be 0
    uint64_t result = divhi >> (62 - srcexp);

    // always return an even value
    int evenodd = result & 1;
    result += evenodd;

    // now compute the result times pio4 to high precision; this multiplies
    // the 1.127 pio4 value by a scalar giving a 65.127 result with exponent -1,
    // or effectively a 64.128 value
    static const uint64_t PIO4_HI = 0xc90fdaa22168c234ull;
    static const uint64_t PIO4_LO = 0xc000000000000000ull;
    auto [mulmid, mulhi] = multiply_64x64(result, PIO4_HI);
    auto [mullo, hitemp2] = multiply_64x64(result, PIO4_LO);
    mulmid += hitemp2;
    if (mulmid < hitemp2)
        mulhi++;

    // finally subtract this value from delta, which is a 1.63 value with exponent
    // srcexp; compute the right shift of the mul result above needed to align the
    // mulmid value with src (this can be up to 63 total)
    int shift = 1 + srcexp;
    if (shift != 0)
    {
        mullo = (mullo >> shift) | (mulmid << (64 - shift));
        mulmid = (mulmid >> shift) | (mulhi << (64 - shift));
        x87_assert((evenodd == 0 && (mulhi >> shift) == 0) || (evenodd != 0 && (mulhi >> shift) <= 1));
    }

    // now do the subtraction, tracking the sign; if we divided by pi/4 evenly, then
    // the diff will be src - mult; if we divided oddly, we added one, and the result
    // will be negative, so compute mult - src
    int sign;
    if (evenodd == 0)
    {
        // we subtract an extra 1 because mullo borrows from srcman's implicit 0
        srcman = srcman - mulmid - 1;
        mullo = -int64_t(mullo);
        sign = 0;
    }
    else
    {
        srcman = mulmid - srcman;
        sign = 1;
    }

    // normalize the result
    if (srcman == 0)
    {
        srcman = mullo;
        mullo = 0;
        srcexp -= 64;
    }
    int lz = count_leading_zeros64(srcman);
    if (lz != 0)
    {
        srcman = (srcman << lz) | (mullo >> (64 - lz));
        if (delta.extended())
            mullo <<= lz;
        srcexp -= lz;
    }

    // assemble the result
    delta = FpType(srcman, delta.extended() ? uint32_t(mullo >> 32) : 0, srcexp, sign);
    return uint32_t(result);
}

template<typename FpType>
static uint32_t reduce_trig_alt(fp64_t src, FpType &delta)
{
    using fpext_t = fpext96_t;
    static const fpext_t pi     (0xc90fdaa22168c234ull, 0xc0000000,  1, 0);  // (3.141592653589e+00)
    static const fpext_t pio2   (0xc90fdaa22168c234ull, 0xc0000000,  0, 0);  // (1.570796326794e+00)
    static const fpext_t pio4   (0xc90fdaa22168c234ull, 0xc0000000, -1, 0);  // (7.853981633974e-01)
    static const fpext_t pio4_hi(0xc90fdaa200000000ull, 0x00000000, -1, 0);  // (7.853981633671e-01)
    static const fpext_t pio4_lo(0x85a308d300000000ull, 0x00000000, -35, 0); // (3.038550253152e-11)
    static const fpext_t invpio4(0xa2f9836e4e44152aull, 0x00062bc4,  0, 0);  // (1.273239544735e+00)

    // convert src to FpType in delta
    src = fp64_t::abs(src);
    delta = FpType(src);

    // if < pi/4, return as-is
    if (src < 0.7853981633974483096)
        return 0;

    // compute x mod PIO4
    fpext_t srcext(src);
    srcext.abs();
    uint64_t j;
    fpext_t y = fpext_t::floor_abs_loint(srcext * invpio4, j);

    // map zeros to origin
    if ((j & 1) != 0)
    {
        j += 1;
        y += fpext_t::one;
    }

    // compute delta in 2 stages to preserve max precision
    fpext_t temp1 = y * pio4_hi;
    srcext -= temp1;
    fpext_t temp2 = y * pio4_lo;
    srcext -= temp2;
    delta = FpType(srcext);
    return uint32_t(j);
}



//===========================================================================
//
// x87_fptan
//
// Compute fptan(x). This code is based off the tanl implementation found in
// the 80-bit Cephes library, found here: https://netlib.org/cephes/
//
//===========================================================================

uint16_t fp64_t::x87_fptan(fp64_t const &src, fp64_t &dst1, fp64_t &dst2)
{
    // only works for exponents < 63
    if (src.exponent() >= 63)
        goto oob;

    // accuracy/speed results:
    //   fpext52_t: 107536(0) / 20706(1) /  442(2), 0.11 ticks
    //   fpext64_t: 107020(0) / 19776(1) / 1888(2), 0.59 ticks
    //   fpext96_t: 107020(0) / 19776(1) / 1888(2), 1.11 ticks
    using fpext_t = fpext52_t;

    // constants
    static std::array<fpext_t, 3> const P =
    {
        fpext_t(0xcc96c69279f9bc1cull, 0x3df84886, 13, 1),  // (-1.309369391814e+04)
        fpext_t(0x8ccf652fe4eee5b1ull, 0x4f58e5c3, 20, 0),  // (1.153516648386e+06)
        fpext_t(0x88ff56994c8baf99ull, 0x8b70bfaf, 24, 1),  // (-1.795652519765e+07)
    };
    static std::array<fpext_t, 4> const Q =
    {
        fpext_t(0xd5c52f759b2b8ed3ull, 0xe2c5b9a6, 13, 0),  // (1.368129634707e+04)
        fpext_t(0xa13de2c155e4adcdull, 0x58dfd25f, 20, 1),  // (-1.320892344402e+06)
        fpext_t(0xbecc7e1756c77adfull, 0x21bc5195, 24, 0),  // (2.500838018234e+07)
        fpext_t(0xcd7f01e5f2d186f6ull, 0x1dc3e1c7, 25, 1),  // (-5.386957559295e+07)
    };

    {
        fpext_t z;
        uint32_t j = reduce_trig(src, z);

        auto sign = src.sign();
        uint16_t flags = src.iszero() ? 0 : src.isdenorm() ? (X87SW_PRECISION_EX | X87SW_DENORM_EX) : X87SW_PRECISION_EX;

        fpext_t zz = z * z;
        if (zz.exponent() > -67)
            dst2 = z.as_fp64() + (z * zz * poly_eval(zz, P)).as_fp64() / poly1_eval(zz, Q).as_fp64();
        else
            dst2 = z.as_fp64();

        if ((j & 2) != 0)
            dst2 = -1.0 / dst2;

        if (sign)
            dst2 = fp64_t::chs(dst2);

        dst1 = fp64_t::const_one();
        return flags;
    }

oob:
    // NaNs in, NaNs out
    if (src.isnan())
        return qnan(dst1, dst2, 0, src);
    // infinities return invalid indefinite
    if (src.isinf())
        return indef(dst1, dst2);
    // normal OOBs just return src and set C2
    // for consistency return 0 as the 2nd result here, but this shouldn't be pushed
    dst2 = fp64_t::const_zero();
    dst1 = src;
    return X87SW_C2;
}



//===========================================================================
//
// x87_sin
// x87_cos
// x87_sincos
//
// Compute fsin/fcos/fsincos(x). This code is based off the sinl/cosl
// implementations found in the 80-bit Cephes library, found here:
// https://netlib.org/cephes/
//
//===========================================================================

// accuracy/speed results for sincos:
//   fpext52_t:  76302(0) / 52242(1) / 140(2), 0.34 ticks
//   fpext64_t: 110016(0) / 18668(1) /   0(2), 1.36 ticks
//   fpext96_t: 124936(0) /  3748(1) /   0(2), 2.93 ticks
using fpextsincos_t = fpext52_t;

static std::array<fpextsincos_t, 7> const s_sincoeffs =
{
    fpextsincos_t(0xd5512389e1d64e26ull, 0x9f89cf50, -41, 1),  // (-7.578540409484e-13)
    fpextsincos_t(0xb0904623e70664d7ull, 0x67a8f274, -33, 0),  // (1.605836316732e-10)
    fpextsincos_t(0xd7322946bf3401b0ull, 0xbe53b744, -26, 1),  // (-2.505210488187e-08)
    fpextsincos_t(0xb8ef1d299845c8f6ull, 0xd25b9a66, -19, 0),  // (2.755731921406e-06)
    fpextsincos_t(0xd00d00d00c536514ull, 0x3dde3d85, -13, 1),  // (-1.984126984125e-04)
    fpextsincos_t(0x8888888888885699ull, 0xb8fd9374,  -7, 0),  // (8.333333333333e-03)
    fpextsincos_t(0xaaaaaaaaaaaaaa97ull, 0x2da4d5f5,  -3, 1),  // (-1.666666666667e-01)
};
static std::array<fpextsincos_t, 7> const s_coscoeffs =
{
    fpextsincos_t(0xd55e8c3a6f997436ull, 0x5436d2ee, -45, 0),  // (4.737750796425e-14)
    fpextsincos_t(0xc9c9920f58f42f36ull, 0xfafa14fe, -37, 1),  // (-1.147028484343e-11)
    fpextsincos_t(0x8f76c648659e534full, 0xab5f5d64, -29, 0),  // (2.087675428708e-09)
    fpextsincos_t(0x93f27dbaf5c64d2bull, 0x0e941cac, -22, 1),  // (-2.755731921500e-07)
    fpextsincos_t(0xd00d00d00c6653edull, 0x149dcc8a, -16, 0),  // (2.480158730157e-05)
    fpextsincos_t(0xb60b60b60b607b66ull, 0xd4ce5b04, -10, 1),  // (-1.388888888889e-03)
    fpextsincos_t(0xaaaaaaaaaaaaaa99ull, 0xa9939f52,  -5, 0),  // (4.166666666667e-02)
};

uint16_t fp64_t::x87_fsin(fp64_t const &src, fp64_t &dst)
{
    // only works for exponents < 63
    if (src.exponent() >= 63)
        goto oob;

    {
        using fpext_t = fpextsincos_t;

        auto sign = src.sign();
        uint16_t flags = src.iszero() ? 0 : src.isdenorm() ? (X87SW_PRECISION_EX | X87SW_DENORM_EX) : X87SW_PRECISION_EX;

        fpext_t z;
        uint32_t j = reduce_trig(src, z);

        fpext_t zz = z * z;
        if (((j + 1) & 2) != 0)
            dst = (fpext_t::one - fpext_t::ldexp(zz, -1) + zz * zz * poly_eval(zz, s_coscoeffs)).as_fp64();
        else
            dst = (z + z * zz * poly_eval(zz, s_sincoeffs)).as_fp64();

        if (((sign ^ (j >> 2)) & 1) != 0)
            dst = fp64_t::chs(dst);

        return flags;
    }

oob:
    // NaNs in, NaNs out
    if (src.isnan())
        return qnan(dst, 0, src);
    // infinities return invalid indefinite
    if (src.isinf())
        return indef(dst);
    // normal OOBs just return src and set C2
    dst = src;
    return X87SW_C2;
}

uint16_t fp64_t::x87_fcos(fp64_t const &src, fp64_t &dst)
{
    // only works for exponents < 63
    if (src.exponent() >= 63)
        goto oob;

    {
        using fpext_t = fpextsincos_t;

        fpext_t z;
        uint32_t j = reduce_trig(src, z);

        auto sign = src.sign();
        uint16_t flags = src.iszero() ? 0 : src.isdenorm() ? (X87SW_PRECISION_EX | X87SW_DENORM_EX) : X87SW_PRECISION_EX;

        fpext_t zz = z * z;
        if (((j + 1) & 2) != 0)
            dst = (z + z * zz * poly_eval(zz, s_sincoeffs)).as_fp64();
        else
            dst = (fpext_t::one - fpext_t::ldexp(zz, -1) + zz * zz * poly_eval(zz, s_coscoeffs)).as_fp64();

        if ((((j >> 1) ^ j) & 2) != 0)
            dst = fp64_t::chs(dst);

        return flags;
    }

oob:
    // NaNs in, NaNs out
    if (src.isnan())
        return qnan(dst, 0, src);
    // infinities return invalid indefinite
    if (src.isinf())
        return indef(dst);
    // normal OOBs just return src and set C2
    dst = src;
    return X87SW_C2;
}

uint16_t fp64_t::x87_fsincos(fp64_t const &src, fp64_t &dst1, fp64_t &dst2)
{
    // only works for exponents < 63
    if (src.exponent() >= 63)
        goto oob;

    {
        using fpext_t = fpextsincos_t;

        fpext_t z;
        uint32_t j = reduce_trig(src, z);

        auto sign = src.sign();
        uint16_t flags = src.iszero() ? 0 : src.isdenorm() ? (X87SW_PRECISION_EX | X87SW_DENORM_EX) : X87SW_PRECISION_EX;

        fpext_t zz = z * z;
        fp64_t res1 = (z + z * zz * poly_eval(zz, s_sincoeffs)).as_fp64();
        fp64_t res2 = (fpext_t::one - fpext_t::ldexp(zz, -1) + zz * zz * poly_eval(zz, s_coscoeffs)).as_fp64();
        if (((j + 1) & 2) != 0)
        {
            dst1 = res1;
            dst2 = res2;
        }
        else
        {
            dst1 = res2;
            dst2 = res1;
        }

        if ((((j >> 1) ^ j) & 2) != 0)
            dst1 = fp64_t::chs(dst1);

        if (((sign ^ (j >> 2)) & 1) != 0)
            dst2 = fp64_t::chs(dst2);

        return flags;
    }

oob:
    // NaNs in, NaNs out
    if (src.isnan())
        return qnan(dst1, dst2, 0, src);
    // infinities return invalid indefinite
    if (src.isinf())
        return indef(dst1, dst2);
    // normal OOBs just return src and set C2
    // for consistency return 0 as the 2nd result here, but this shouldn't be pushed
    dst2 = fp64_t::const_zero();
    dst1 = src;
    return X87SW_C2;
}



//===========================================================================
//
// x87_fpatan
//
// Compute atan2(y, x). This code is based off the atanl/atan2l
// implementations found in the 80-bit Cephes library, found here:
// https://netlib.org/cephes/
//
//===========================================================================

uint16_t fp64_t::x87_fpatan(fp64_t const &src1, fp64_t const &src2, fp64_t &dst)
{
    // denorm flag is set regardless
    uint16_t flags = 0;
    if (src1.isdenorm() || src2.isdenorm())
        flags |= X87SW_DENORM_EX;

    // handle special cases
    if (src1.ismaxexp())
        goto specialx;
    if (src2.ismaxexp())
        goto specialy;
    if (src1.iszero())
        goto zerox;
    if (src2.iszero())
        goto zeroy;

    // accuracy/speed results:
    //   fpext52_t: fails all over the place
    //   fpext64_t: 140342188(0) / 25263290(1) / 5683(2), 1.26 ticks
    //   fpext96_t: 141032668(0) / 24572779(1) / 5714(2), 2.60 ticks
    using fpext_t = fpext64_t;

    static std::array<fpext_t, 5> const P =
    {
        fpext_t(0xde5f1266ce538eceull, 0x45933bae, -1, 1),  // (-8.686381817809e-01)
        fpext_t(0xeaefa6bfa06107e6ull, 0x6f351563,  3, 1),  // (-1.468350863318e+01)
        fpext_t(0xffe8557ff29153eeull, 0x47487583,  5, 1),  // (-6.397688865583e+01)
        fpext_t(0xc7fa3f3eeda6f9d5ull, 0xa7a03a0c,  6, 1),  // (-9.998876377727e+01)
        fpext_t(0xcb9393616abcb6c3ull, 0x53e3ffa9,  5, 1),  // (-5.089411689962e+01)
    };
    static std::array<fpext_t, 5> const Q =
    {
        fpext_t(0xb7dae76e894e54d3ull, 0xee74072e,  4, 0),  // (2.298188673359e+01)
        fpext_t(0x8ffdafa27a4676b8ull, 0xd644a00e,  7, 0),  // (1.439909612225e+02)
        fpext_t(0xb4b86beee9c0e3a9ull, 0x5df2ff95,  8, 0),  // (3.614407938615e+02)
        fpext_t(0xc3c9b09850a7abc0ull, 0xb934a367,  8, 0),  // (3.915757017511e+02)
        fpext_t(0x98aeae89100d891bull, 0xd3dd1204,  7, 0),  // (1.526823506989e+02)
    };
    static double const T3P8 = 2.41421356237309504880169;
    static double const TP8 = 4.1421356237309504880169e-1;

    static double const pi64 = 3.1415926535897932384626433832795;
    static double const npi64 = -3.1415926535897932384626433832795;
    static double const pio264 = 1.5707963267948966192313216916398;
    static double const npio264 = -1.5707963267948966192313216916398;
    static double const pio464 = 0.78539816339744830961566084581988;
    static double const npio464 = -0.78539816339744830961566084581988;
    static double const pi3o464 = 2.3561944901923449288469825374596;
    static double const npi3o464 = -2.3561944901923449288469825374596;

    static const fpext_t pio2(0xc90fdaa22168c234ull, 0xc0000000,  0, 0);  // (1.570796326794e+00)
    static const fpext_t pio4(0xc90fdaa22168c234ull, 0xc0000000, -1, 0);  // (7.853981633974e-01)

    {
        fp64_t x = src2 / src1;

        // make argument positive and save the sign
        int sign = 0;
        if (x < 0)
        {
            sign = 1;
            x = fp64_t::chs(x);
        }

        // range reduction
        fpext_t yext;
        fpext_t xext;
        if (x > T3P8)
        {
            yext = pio2;
            xext = fpext_t(-1.0 / x);
        }
        else if (x > TP8)
        {
            yext = pio4;
            xext = fpext_t((x - 1.0) / (x + 1.0));
        }
        else
        {
            yext = fpext_t::zero;
            xext = fpext_t(x);
        }

        fpext_t z = xext * xext;
        yext = yext + poly_eval(z, P).div64(poly1_eval(z, Q)) * z * xext + xext;

        if (sign)
            yext.chs();

        int code = (src1.sign() << 1) | src2.sign();
        dst = yext.as_fp64();

        static fp64_t const s_offsets[4] = { 0.0, 0.0, pi64, npi64 };
        dst += s_offsets[code];

        if (dst == 0 && src2.sign())
            dst = fp64_t::chs(dst);

        return flags | X87SW_PRECISION_EX;
    }

specialx:
    // NaNs in, NaNs out
    if (src1.isnan())
        return qnan(dst, flags, src1, src2);
    if (src2.isnan())
        return qnan(dst, flags, src2);
    // handle infinities
    if (src2.isinf())
    {
        if (src1.sign() == 0)
            dst = (src2.sign() == 0) ? pio464 : npio464;
        else
            dst = (src2.sign() == 0) ? pi3o464 : npi3o464;
    }
    // handle finite numbers
    else
    {
        if (src1.sign() == 0)
            return zero(dst, flags, src2.sign());
        else
            dst = (src2.sign() == 0) ? pi64 : npi64;
    }
    return X87SW_PRECISION_EX | flags;

specialy:
    // NaNs in, NaNs out
    if (src2.isnan())
        return qnan(dst, flags, src2);
    // infinite results depend on sign
    dst = (src2.sign() == 0) ? pio264 : npio264;
    return X87SW_PRECISION_EX | flags;

zerox:
    // handle X and Y both 0
    if (src2.iszero())
    {
        if (src1.sign() == 0)
            return zero(dst, flags, src2.sign());
        else
            dst = (src2.sign() == 0) ? pi64 : npi64;
    }
    // handle other cases
    else
    {
        dst = (src2.sign() == 0) ? pio264 : npio264;
        flags |= X87SW_PRECISION_EX;
    }
    return flags;

zeroy:
    flags |= (src1.sign() ? X87SW_PRECISION_EX : 0);
    if (src1.sign() == 0)
        return zero(dst, flags, src2.sign());
    else
        dst = (src2.sign() == 0) ? pi64 : npi64;
    return flags;
}

}
