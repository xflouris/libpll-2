#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([libpll], [0.3.2], [Tomas.Flouri@h-its.org])
AM_INIT_AUTOMAKE([subdir-objects])
AC_LANG([C])
AC_PROG_CC([gcc icc])
AC_CONFIG_SRCDIR([src/pll.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

# set default gcc options
CFLAGS='-g'

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_LEX
AC_PROG_YACC
AC_PROG_INSTALL

if test "x$LEX" != xflex; then
  AC_MSG_ERROR(could not find required installation of FLEX)
fi

if test "x$YACC" != x"bison -y"; then
  AC_MSG_ERROR(could not find required installation of BISON)
fi

LT_INIT

# Checks for libraries.
AC_CHECK_LIB([m],[exp])

# Checks for header files.
AC_CHECK_HEADERS([assert.h math.h stdio.h stdlib.h string.h ctype.h x86intrin.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_ALLOCA
AC_FUNC_REALLOC
AC_CHECK_FUNCS([asprintf memcpy memset posix_memalign])

have_avx512f=no
have_avx2=no
have_avx=no
have_sse3=no
have_ps2pdf=no

# Detection for ICC/GCC intrinsic flags
my_save_cflags="$CFLAGS"
CFLAGS=-mavx2
AC_MSG_CHECKING([whether CC supports -mavx2])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [AC_SUBST(AM_AVX2_FLAGS, [-mavx2])],
    [AC_MSG_RESULT([no])]
)
CFLAGS=-march=core-avx2
AC_MSG_CHECKING([whether CC supports -march=core-avx2])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [AC_SUBST(AM_AVX2_FLAGS, [-march=core-avx2])],
    [AC_MSG_RESULT([no])]
)
CFLAGS=-mavx512f
AC_MSG_CHECKING([whether CC supports -mavx512f])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [AC_SUBST(AM_AVX512F_FLAGS, [-mavx512f])],
    [AC_MSG_RESULT([no])]
)
CFLAGS=-xCOMMON-AVX512
AC_MSG_CHECKING([whether CC supports -xCOMMON-AVX512])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [AC_SUBST(AM_AVX512F_FLAGS, [-xCOMMON-AVX512])],
    [AC_MSG_RESULT([no])]
)
CFLAGS="$my_save_cflags"

# Compile-time detection of processor features - now disabled
#TODO: No detection for AVX512F yet
#AX_EXT
#
#if test "x${ax_cv_have_fma3_ext}" = "xyes"; then
#  have_avx2=yes
#fi
#
#if test "x${ax_cv_have_avx_ext}" = "xyes"; then
#  have_avx=yes
#fi
#
#if test "x${ax_cv_have_sse3_ext}" = "xyes"; then
#  have_sse3=yes
#fi
#
#AC_ARG_ENABLE(avx2, AS_HELP_STRING([--disable-avx2], [Build without AVX2/FMA support]))
#AS_IF([test "x${ax_cv_have_fma3_ext}" = "xyes"], [
#  have_avx2=yes
#])
#
#AC_ARG_ENABLE(avx, AS_HELP_STRING([--disable-avx], [Build without AVX support]))
#AS_IF([test "x${ax_cv_have_avx_ext}" = "xyes"], [
#  have_avx=yes
#])
#
#AC_ARG_ENABLE(sse, AS_HELP_STRING([--disable-sse],[Build without SSE support]))
#AS_IF([test "x${ax_cv_have_sse3_ext}" = "xyes"], [
#  have_sse3=yes
#])
#

AC_ARG_ENABLE(pdfman, AS_HELP_STRING([--disable-pdfman], [Disable PDF manual creation]))
AS_IF([test "x$enable_pdfman" != "xno"], [
  have_ps2pdf=yes
  AC_CHECK_PROG(HAVE_PS2PDF, ps2pdf, yes, no)
  if test "x$HAVE_PS2PDF" = "xno"; then
    AC_MSG_WARN([*** ps2pdf is required to build a PDF version of the manual])
    have_ps2pdf=no
  fi
])

AC_ARG_ENABLE(sse, AS_HELP_STRING([--disable-sse], [Build without SSE support]))
AS_IF([test "x$enable_sse" != "xno"], [
  have_sse3=yes
  AC_DEFINE([HAVE_SSE3], [1], [Define to 1 to support Streaming SIMD Extensions 3])
])

AC_ARG_ENABLE(avx, AS_HELP_STRING([--disable-avx], [Build without AVX support]))
AS_IF([test "x$enable_avx" != "xno"], [
  have_avx=yes
  AC_DEFINE([HAVE_AVX], [1], [Define to 1 to support Advanced Vector Extensions])
])

AC_ARG_ENABLE(avx2, AS_HELP_STRING([--disable-avx2], [Build without AVX2/FMA support]))
AS_IF([test "x$enable_avx2" != "xno"], [
  have_avx2=yes
  AC_DEFINE([HAVE_AVX2], [1], [Define to 1 to support Advanced Vector Extensions 2])
])

AC_ARG_ENABLE(avx512f, AS_HELP_STRING([--disable-avx512f], [Build without AVX512F support]))
AS_IF([test "x$enable_avx512f" != "xno"], [
  have_avx512f=yes
  AC_DEFINE([HAVE_AVX512F], [1], [Define to 1 to support Advanced Vector Extensions 512 Foundation])
])

AM_CONDITIONAL(HAVE_AVX512F, test "x${have_avx512f}" = "xyes")
AM_CONDITIONAL(HAVE_AVX2, test "x${have_avx2}" = "xyes")
AM_CONDITIONAL(HAVE_AVX, test "x${have_avx}" = "xyes")
AM_CONDITIONAL(HAVE_SSE3, test "x${have_sse3}" = "xyes")
AM_CONDITIONAL(HAVE_PS2PDF, test "x${have_ps2pdf}" = "xyes")
AM_PROG_CC_C_O

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 man/Makefile
                 examples/Makefile])
AC_OUTPUT

AC_MSG_RESULT([
  $PACKAGE $VERSION

    Target:        $host_os $host_cpu
    Compiler:      ${CC}
    CFLAGS:        ${CFLAGS} ${CPPFLAGS}
    LIBS:          ${LIBS} ${LDFLAGS}
    AVX2 Flags:    ${AM_AVX2_FLAGS}
    AVX512F Flags: ${AM_AVX512F_FLAGS}

  Continue with 'make' command
])
