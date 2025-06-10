//=========================================================
//  x87common.h
//
//  Common floating-point definitions
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

#ifndef X87COMMON_H
#define X87COMMON_H

//
// enable this to use the standard library's routines for setting/getting the
// current FP rounding mode; when disabled, the code will intrinsics or inline
// assembly to do the same, which may cause some compatibility concerns
//
#define X87_USE_CFENV (0)

//
// all asserts within this code should by x87_assert; define x87_assert at build
// time to map to your own assert macro if you prefer, otherwise we fall back to
// the standard library's assert
//
#ifndef x87_assert
#include <cassert>
#define x87_assert assert
#endif

#include <cstdint>
#include <cmath>
#include <algorithm>

#if X87_USE_CFENV
#include <cfenv>
#elif defined(_WIN32)
#include <float.h>
#endif


//===========================================================================
//
// double_int64_t
// int64_double_t
// float_int32_t
// int32_float_t
//
// Helper unions for converting between floating-point types and raw bits.
// Since initializers apply to the first entry, we have versions with both
// the raw first and the float first.
//
//===========================================================================

namespace x87
{

union double_int64_t { double d; uint64_t i; };
union int64_double_t { uint64_t i; double d; };
union float_int32_t { float d; uint32_t i; };
union int32_float_t { uint32_t i; float d; };

}



//===========================================================================
//
// MANTISSA_MASK/MANTISSA_BITS
// EXPONENT_MASK/EXPONENT_SHIFT/EXPONENT_BIAS/EXPONENT_MIN/EXPONENT_MAX_BIASED
// SIGN_MASK-SIGN_SHIFT
//
// Constants for each of our known floating-point types.
//
//===========================================================================

namespace x87
{

//
// 32-bit FP constants
//
static constexpr int FP32_MANTISSA_BITS = 23;
static constexpr int FP32_EXPONENT_BITS = 8;
static constexpr int FP32_EXPONENT_SHIFT = FP32_MANTISSA_BITS;
static constexpr int FP32_SIGN_SHIFT = FP32_EXPONENT_SHIFT + FP32_EXPONENT_BITS;

static constexpr uint32_t FP32_MANTISSA_MASK = (1 << FP32_MANTISSA_BITS) - 1;
static constexpr uint32_t FP32_EXPONENT_MASK = ((1 << FP32_EXPONENT_BITS) - 1) << FP32_EXPONENT_SHIFT;
static constexpr uint32_t FP32_SIGN_MASK = 1 << FP32_SIGN_SHIFT;
static constexpr uint32_t FP32_ABS_MASK = FP32_MANTISSA_MASK | FP32_EXPONENT_MASK;

static constexpr int32_t FP32_EXPONENT_BIAS = 0x7f;
static constexpr int32_t FP32_EXPONENT_MIN_BIASED = 0;
static constexpr int32_t FP32_EXPONENT_MAX_BIASED = (1 << FP32_EXPONENT_BITS) - 1;

//
// 64-bit FP constants
//
static constexpr int FP64_MANTISSA_BITS = 52;
static constexpr int FP64_EXPONENT_BITS = 11;
static constexpr int FP64_EXPONENT_SHIFT = FP64_MANTISSA_BITS;
static constexpr int FP64_SIGN_SHIFT = FP64_EXPONENT_SHIFT + FP64_EXPONENT_BITS;

static constexpr uint64_t FP64_MANTISSA_MASK = (1ull << FP64_MANTISSA_BITS) - 1;
static constexpr uint64_t FP64_EXPONENT_MASK = ((1ull << FP64_EXPONENT_BITS) - 1) << FP64_EXPONENT_SHIFT;
static constexpr uint64_t FP64_SIGN_MASK = 1ull << FP64_SIGN_SHIFT;
static constexpr uint64_t FP64_ABS_MASK = FP64_MANTISSA_MASK | FP64_EXPONENT_MASK;

static constexpr int32_t FP64_EXPONENT_BIAS = 0x3ff;
static constexpr int32_t FP64_EXPONENT_MIN_BIASED = 0;
static constexpr int32_t FP64_EXPONENT_MAX_BIASED = (1 << FP64_EXPONENT_BITS) - 1;

//
// 80-bit FP constants
//
static constexpr int FP80_MANTISSA_BITS = 63;
static constexpr int FP80_EXPONENT_BITS = 15;
static constexpr int FP80_EXPONENT_SHIFT = 0;
static constexpr int FP80_SIGN_SHIFT = FP80_EXPONENT_SHIFT + FP80_EXPONENT_BITS;

static constexpr uint64_t FP80_MANTISSA_MASK = (1ull << FP80_MANTISSA_BITS) - 1;
static constexpr uint16_t FP80_EXPONENT_MASK = ((1 << FP80_EXPONENT_BITS) - 1) << FP80_EXPONENT_SHIFT;
static constexpr uint16_t FP80_SIGN_MASK = 1 << FP80_SIGN_SHIFT;
static constexpr uint16_t FP80_ABS_MASK = FP80_EXPONENT_MASK;

static constexpr int32_t FP80_EXPONENT_BIAS = 0x3fff;
static constexpr int32_t FP80_EXPONENT_MIN_BIASED = 0;
static constexpr int32_t FP80_EXPONENT_MAX_BIASED = (1 << FP80_EXPONENT_BITS) - 1;

static constexpr uint64_t FP80_EXPLICIT_ONE = 0x8000000000000000ull;

}



//===========================================================================
//
// X87CW_*
// X87SW_*
//
// Constants for x87 control word and status word.
//
//===========================================================================

namespace x87
{

//
// x87 control word values
//
using x87cw_t = uint16_t;
static constexpr x87cw_t X87CW_MASK_INVALID_EX    = 0x0001;
static constexpr x87cw_t X87CW_MASK_DENORM_EX     = 0x0002;
static constexpr x87cw_t X87CW_MASK_DIVZERO_EX    = 0x0004;
static constexpr x87cw_t X87CW_MASK_OVERFLOW_EX   = 0x0008;
static constexpr x87cw_t X87CW_MASK_UNDERFLOW_EX  = 0x0010;
static constexpr x87cw_t X87CW_MASK_PRECISION_EX  = 0x0020;
static constexpr x87cw_t X87CW_MASK_ALL_EX        = X87CW_MASK_INVALID_EX | X87CW_MASK_DENORM_EX | X87CW_MASK_DIVZERO_EX | X87CW_MASK_OVERFLOW_EX | X87CW_MASK_UNDERFLOW_EX | X87CW_MASK_PRECISION_EX;
static constexpr int X87CW_PRECISION_SHIFT        = 8;
static constexpr x87cw_t X87CW_PRECISION_MASK     = 3 << X87CW_PRECISION_SHIFT;
static constexpr x87cw_t X87CW_PRECISION_SINGLE   = 0 << X87CW_PRECISION_SHIFT;
static constexpr x87cw_t X87CW_PRECISION_DOUBLE   = 2 << X87CW_PRECISION_SHIFT;
static constexpr x87cw_t X87CW_PRECISION_EXTENDED = 3 << X87CW_PRECISION_SHIFT;
static constexpr int X87CW_ROUNDING_SHIFT         = 10;
static constexpr x87cw_t X87CW_ROUNDING_MASK      = 3 << X87CW_ROUNDING_SHIFT;
static constexpr x87cw_t X87CW_ROUNDING_NEAREST   = 0 << X87CW_ROUNDING_SHIFT;
static constexpr x87cw_t X87CW_ROUNDING_DOWN      = 1 << X87CW_ROUNDING_SHIFT;
static constexpr x87cw_t X87CW_ROUNDING_UP        = 2 << X87CW_ROUNDING_SHIFT;
static constexpr x87cw_t X87CW_ROUNDING_ZERO      = 3 << X87CW_ROUNDING_SHIFT;
static constexpr x87cw_t X87CW_INFINITY           = 0x1000;

//
// x87 status word values
//
using x87sw_t = uint16_t;
static constexpr x87sw_t X87SW_INVALID_EX    = 0x0001;
static constexpr x87sw_t X87SW_DENORM_EX     = 0x0002;
static constexpr x87sw_t X87SW_DIVZERO_EX    = 0x0004;
static constexpr x87sw_t X87SW_OVERFLOW_EX   = 0x0008;
static constexpr x87sw_t X87SW_UNDERFLOW_EX  = 0x0010;
static constexpr x87sw_t X87SW_PRECISION_EX  = 0x0020;
static constexpr x87sw_t X87SW_STACK_FAULT   = 0x0040;
static constexpr x87sw_t X87SW_ERROR_SUMMARY = 0x0080;
static constexpr int X87SW_C0_BIT            = 8;
static constexpr x87sw_t X87SW_C0            = 1 << X87SW_C0_BIT;
static constexpr int X87SW_C1_BIT            = 9;
static constexpr x87sw_t X87SW_C1            = 1 << X87SW_C1_BIT;
static constexpr int X87SW_C2_BIT            = 10;
static constexpr x87sw_t X87SW_C2            = 1 << X87SW_C2_BIT;
static constexpr int X87SW_TOP_SHIFT         = 11;
static constexpr x87sw_t X87SW_TOP_MASK      = 7 << X87SW_TOP_SHIFT;
static constexpr int X87SW_C3_BIT            = 14;
static constexpr x87sw_t X87SW_C3            = 1 << X87SW_C3_BIT;
static constexpr x87sw_t X87SW_BUSY          = 0x8000;

}



//===========================================================================
//
// fpround_t
//
// Stack-based class to configure the FPU rounding state and restore it
// on destruction.
//
//===========================================================================

namespace x87
{

class fpround_t
{
public:
    //
    // constructor
    //
    explicit fpround_t(x87cw_t round) :
        m_oldmode(fpround_t::get())
    {
        fpround_t::set(round);
    }

    //
    // destructor
    //
    ~fpround_t()
    {
        fpround_t::set(m_oldmode);
    }

    //
    // set the rounding mode on the CPU to match the x87 rounding mode provided
    //
    static void set(x87cw_t round)
    {
#if X87_USE_CFENV
        // see if there's a direct 1:1 mapping between values
        if constexpr(FE_TONEAREST == 0 && FE_UPWARD == 2 * FE_DOWNWARD && FE_TOWARDZERO == 3 * FE_DOWNWARD)
            std::fesetround(FE_DOWNWARD * ((round & X87CW_ROUNDING_MASK) >> X87CW_ROUNDING_SHIFT));
        else
        {
            static int const s_mode[4] = { FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO };
            std::fesetround(s_mode[(round & X87CW_ROUNDING_MASK) >> X87CW_ROUNDING_SHIFT]);
        }
#elif defined(_WIN32)
        static_assert((X87CW_ROUNDING_NEAREST >> 2) == _RC_NEAR);
        static_assert((X87CW_ROUNDING_DOWN >> 2) == _RC_DOWN);
        static_assert((X87CW_ROUNDING_UP >> 2) == _RC_UP);
        static_assert((X87CW_ROUNDING_ZERO >> 2) == _RC_CHOP);
        _controlfp(round >> 2, _MCW_RC);
#elif defined(_M_AMD64) || defined(__amd64__)
        uint32_t temp = round;
        __asm__ volatile(
            "sub $16, %%rsp\n\t"
            "stmxcsr 0(%%rsp)\n\t"
            "shl $3, %0\n\t"
            "andl $0x6000, %0\n\t"
            "andl $0x9fff, 0(%%rsp)\n\t"
            "orl %0, 0(%%rsp)\n\t"
            "ldmxcsr 0(%%rsp)\n\t"
            "add $16, %%rsp\n\t"
            :
            : "r"(temp)
        );
#elif defined(_M_ARM64) || defined(__aarch64__)
        uint32_t temp = round;
        __asm__ volatile(
            "mrs x17, fpcr\n\t"
            "rbit w16, %w0\n\t"
            "ubfx w16, w16, #20, #2\n\t"
            "bfi x17, x16, #22, #2\n\t"
            "msr fpcr, x17\n\t"
            :
            : "r"(temp)
        );
#else
#error Unsupported environment
#endif
    }

    //
    // return the current rounding mode on the CPU
    //
    static x87cw_t get()
    {
#if X87_USE_CFENV
        auto round = std::fegetround();
        if constexpr(FE_TONEAREST == 0 && FE_UPWARD == 2 * FE_DOWNWARD && FE_TOWARDZERO == 3 * FE_DOWNWARD)
            return (round / FE_DOWNWARD) << X87CW_ROUNDING_SHIFT;
        else if (round == FE_TONEAREST)
            return X87CW_ROUNDING_NEAREST;
        else if (round == FE_DOWNWARD)
            return X87CW_ROUNDING_DOWN;
        else if (round == FE_UPWARD)
            return X87CW_ROUNDING_UP;
        else
            return X87CW_ROUNDING_ZERO;
#elif defined(_WIN32)
        return (_controlfp(0, 0) << 2) & X87CW_ROUNDING_MASK;
#elif defined(_M_AMD64) || defined(__amd64__)
        uint32_t mxcsr;
        __asm__ volatile(
            "sub $16, %%rsp\n\t"
            "stmxcsr 0(%%rsp)\n\t"
            "movl 0(%%rsp), %%eax\n\t"
            "add $16, %%rsp\n\t"
            : "=a"(mxcsr)
        );
        return (mxcsr >> 3) & X87CW_ROUNDING_MASK;
#elif defined(_M_ARM64) || defined(__aarch64__)
        uint32_t result;
        __asm__ volatile(
            "mrs x17, fpcr\n\t"
            "rbit w17, w17\n\t"
            "ubfx w17, w17, #8, #2\n\t"
            "lsl %w0, w17, #10\n\t"
            : "=r"(result)
        );
        return result;
#else
#error Unsupported environment
#endif
    }

private:
    //
    // internal state
    //
    x87cw_t m_oldmode;
};

}



//===========================================================================
//
// multiply_64x64
//
// Perform a 64x64-bit multiplication and return the full 128-bit result.
//
//===========================================================================

namespace x87
{

struct result128_t { uint64_t lo, hi; };

#ifdef _MSC_VER

#if defined(_M_X64)

inline result128_t multiply_64x64(uint64_t a, uint64_t b)
{
    result128_t result;
    result.lo = _umul128(a, b, &result.hi);
    return result;
}

#elif defined(_M_ARM64)

extern "C" unsigned __int64 __umulh(unsigned __int64 a, unsigned __int64 b);
#pragma intrinsic(__umulh)
inline result128_t multiply_64x64(uint64_t a, uint64_t b)
{
    return { a * b, __umulh(a, b) };
}

#endif

#else

inline result128_t multiply_64x64(uint64_t a, uint64_t b)
{
    union
    {
        result128_t result;
        unsigned __int128 i128;
    } u;
    u.i128 = (unsigned __int128)a * b;
    return u.result;
}

#endif

}



//===========================================================================
//
// count_leading_zeros64 / count_trailing_zeros64
//
// Count the number of leading or trailing zeros using architecture-specific
// helpers.
//
//===========================================================================

namespace x87
{

//
// count leading zeros in a 64-bit value
//
inline int count_leading_zeros64(uint64_t value)
{
#ifdef _MSC_VER
    unsigned long index;
    return _BitScanReverse64(&index, value) ? (63 - index) : 64;
#else
    return __builtin_clzll(value);
#endif
}


//
// count trailing zeros in a 64-bit value
//
inline int count_trailing_zeros64(uint64_t value)
{
#ifdef _MSC_VER
    unsigned long index;
    return _BitScanForward64(&index, value) ? index : 64;
#else
    return __builtin_ctzll(value);
#endif
}

}

#endif
