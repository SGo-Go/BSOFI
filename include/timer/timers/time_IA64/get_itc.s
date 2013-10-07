.text
        .align 16
// C version
//    long clock_tic()
//
        .global clock_tic#
        .proc clock_tic#
clock_tic:
        mov r8 = ar.itc
        br.ret.sptk.many b0
        .endp clock_tic#

        .align 16
// Fortran version
//      integer*8 clock_tic
//
        .global clock_tic_#
        .proc clock_tic_#
clock_tic_:
        mov r8 = ar.itc
        br.ret.sptk.many b0
        .endp clock_tic_#
