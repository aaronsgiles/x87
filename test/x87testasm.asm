; nasm -f win64 x87testasm.asm
    section .data

saved_cw dw 0

    section .text

    bits 64

    global x87consts80
x87consts80:
    fld1
    fstp    tword [rcx + 0*16]
    fldl2t
    fstp    tword [rcx + 1*16]
    fldl2e
    fstp    tword [rcx + 2*16]
    fldpi
    fstp    tword [rcx + 3*16]
    fldlg2
    fstp    tword [rcx + 4*16]
    fldln2
    fstp    tword [rcx + 5*16]
    fldz
    fstp    tword [rcx + 6*16]
    ret

    global x87consts64
x87consts64:
    fld1
    fstp    qword [rcx + 0*8]
    fldl2t
    fstp    qword [rcx + 1*8]
    fldl2e
    fstp    qword [rcx + 2*8]
    fldpi
    fstp    qword [rcx + 3*8]
    fldlg2
    fstp    qword [rcx + 4*8]
    fldln2
    fstp    qword [rcx + 5*8]
    fldz
    fstp    qword [rcx + 6*8]
    ret

    global x87setcw
x87setcw:
    fldcw   [rcx]
    fstcw   [rel saved_cw]
    ret

    global x87getsw
x87getsw:
    fstsw   ax
    movzx   eax, ax
    ret

    global x87test1
x87test1:
    fld     qword [rcx]
    fstp    tword [rdx]
    ret

    global x87test2
x87test2:
    fld     tword [rcx]
    fstp    qword [rdx]
    ret

    global fld8080
fld8080:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    fstsw   ax
    fstp    tword [rdx]
    ret

    global fld6480
fld6480:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rcx]
    fstsw   ax
    fstp    tword [rdx]
    ret

    global fld3280
fld3280:
    finit
    fldcw   [rel saved_cw]
    fld     dword [rcx]
    fstsw   ax
    fstp    tword [rdx]
    ret

    global fild6480
fild6480:
    finit
    fldcw   [rel saved_cw]
    fild    qword [rcx]
    fstsw   ax
    fstp    tword [rdx]
    ret

    global fild3280
fild3280:
    finit
    fldcw   [rel saved_cw]
    fild    dword [rcx]
    fstsw   ax
    fstp    tword [rdx]
    ret

    global fild1680
fild1680:
    finit
    fldcw   [rel saved_cw]
    fild    word [rcx]
    fstsw   ax
    fstp    tword [rdx]
    ret

    global fst8080
fst8080:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    fstp    tword [rdx]
    fstsw   ax
    ret

    global fst8064
fst8064:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    fstp    qword [rdx]
    fstsw   ax
    ret

    global fst8032
fst8032:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    fstp    dword [rdx]
    fstsw   ax
    ret

    global fist8064
fist8064:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    fistp   qword [rdx]
    fstsw   ax
    ret

    global fist8032
fist8032:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    fistp   dword [rdx]
    fstsw   ax
    ret

    global fist8016
fist8016:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    fistp   word [rdx]
    fstsw   ax
    ret

    global fadd80
fadd80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rdx]
    fld     tword [rcx]
    faddp
    fstsw   ax
    fstp    tword [r8]
    ret

    global fsub80
fsub80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rdx]
    fld     tword [rcx]
    fsubp
    fstsw   ax
    fstp    tword [r8]
    ret

    global fmul80
fmul80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rdx]
    fld     tword [rcx]
    fmulp
    fstsw   ax
    fstp    tword [r8]
    ret

    global fdiv80
fdiv80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rdx]
    fld     tword [rcx]
    fdivp
    fstsw   ax
    fstp    tword [r8]
    ret

    global fsqrt80
fsqrt80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    fsqrt
    fstsw   ax
    fstp    tword [rdx]
    ret

    global f2xm180
f2xm180:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    f2xm1
    fstsw   ax
    fstp    tword [rdx]
    ret

    global fyl2x80
fyl2x80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rdx]
    fld     tword [rcx]
    fyl2x
    fstsw   ax
    fstp    tword [r8]
    ret

    global fptan80
fptan80:
    finit
    fldcw   [rel saved_cw]
    fldz
    fld     tword [rcx]
    fptan
    fstsw   ax
    fstp    tword [rdx]
    fstp    tword [r8]
    ret

    global fsincos80
fsincos80:
    finit
    fldcw   [rel saved_cw]
    fldz
    fld     tword [rcx]
    fsincos
    fstsw   ax
    fstp    tword [rdx]
    fstp    tword [r8]
    ret

    global fpatan80
fpatan80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rdx]
    fld     tword [rcx]
    fpatan
    fstsw   ax
    fstp    tword [r8]
    ret

    global fxtract80
fxtract80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    fxtract
    fstsw   ax
    fstp    tword [rdx]
    fstp    tword [r8]
    ret

    global fprem180
fprem180:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rdx]
    fld     tword [rcx]
    fprem1
    fstsw   ax
    fstp    tword [r8]
    ret

    global fprem80
fprem80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rdx]
    fld     tword [rcx]
    fprem
    fstsw   ax
    fstp    tword [r8]
    ret

    global fyl2xp180
fyl2xp180:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rdx]
    fld     tword [rcx]
    fyl2xp1
    fstsw   ax
    fstp    tword [r8]
    ret

    global frndint80
frndint80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    frndint
    fstsw   ax
    fstp    tword [rdx]
    ret

    global fscale80
fscale80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rdx]
    fld     tword [rcx]
    fscale
    fstsw   ax
    fstp    tword [r8]
    ret

    global fsin80
fsin80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    fsin
    fstsw   ax
    fstp    tword [rdx]
    ret

    global fcos80
fcos80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [rcx]
    fcos
    fstsw   ax
    fstp    tword [rdx]
    ret

    global fadd64
fadd64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rdx]
    fld     qword [rcx]
    faddp
    fstsw   ax
    fstp    qword [r8]
    ret

    global fsub64
fsub64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rdx]
    fld     qword [rcx]
    fsubp
    fstsw   ax
    fstp    qword [r8]
    ret

    global fmul64
fmul64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rdx]
    fld     qword [rcx]
    fmulp
    fstsw   ax
    fstp    qword [r8]
    ret

    global fdiv64
fdiv64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rdx]
    fld     qword [rcx]
    fdivp
    fstsw   ax
    fstp    qword [r8]
    ret

    global fsqrt64
fsqrt64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rcx]
    fsqrt
    fstsw   ax
    fstp    qword [rdx]
    ret

    global f2xm164
f2xm164:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rcx]
    f2xm1
    fstsw   ax
    fstp    qword [rdx]
    ret

    global fyl2x64
fyl2x64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rdx]
    fld     qword [rcx]
    fyl2x
    fstsw   ax
    fstp    qword [r8]
    ret

    global fptan64
fptan64:
    finit
    fldcw   [rel saved_cw]
    fldz
    fld     qword [rcx]
    fptan
    fstsw   ax
    fstp    qword [rdx]
    fstp    qword [r8]
    ret

    global fsincos64
fsincos64:
    finit
    fldcw   [rel saved_cw]
    fldz
    fld     qword [rcx]
    fsincos
    fstsw   ax
    fstp    qword [rdx]
    fstp    qword [r8]
    ret

    global fpatan64
fpatan64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rdx]
    fld     qword [rcx]
    fpatan
    fstsw   ax
    fstp    qword [r8]
    ret

    global fxtract64
fxtract64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rcx]
    fxtract
    fstsw   ax
    fstp    qword [rdx]
    fstp    qword [r8]
    ret

    global fprem164
fprem164:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rdx]
    fld     qword [rcx]
    fprem1
    fstsw   ax
    fstp    qword [r8]
    ret

    global fprem64
fprem64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rdx]
    fld     qword [rcx]
    fprem
    fstsw   ax
    fstp    qword [r8]
    ret

    global fyl2xp164
fyl2xp164:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rdx]
    fld     qword [rcx]
    fyl2xp1
    fstsw   ax
    fstp    qword [r8]
    ret

    global frndint64
frndint64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rcx]
    frndint
    fstsw   ax
    fstp    qword [rdx]
    ret

    global fscale64
fscale64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rdx]
    fld     qword [rcx]
    fscale
    fstsw   ax
    fstp    qword [r8]
    ret

    global fsin64
fsin64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rcx]
    fsin
    fstsw   ax
    fstp    qword [rdx]
    ret

    global fcos64
fcos64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [rcx]
    fcos
    fstsw   ax
    fstp    qword [rdx]
    ret
