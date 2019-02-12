AC_DEFUN([SWIFT_ADD_LIBRARIES],
[
# Check for the libraries we will need.
AC_CHECK_LIB(m,sqrt,,AC_MSG_ERROR(something is wrong with the math library!))

# Check for METIS.
have_metis="no"
AC_ARG_WITH([metis],
    [AS_HELP_STRING([--with-metis=PATH],
       [root directory where METIS is installed @<:@yes/no@:>@]
    )],
    [with_metis="$withval"],
    [with_metis="no"]
)

METIS_LIBS=""
if test "x$with_metis" != "xno"; then

# Check if we have METIS.
   if test "x$with_metis" != "xyes" -a "x$with_metis" != "x"; then
      METIS_LIBS="-L$with_metis/lib -lmetis"
      METIS_INCS="-I$with_metis/include"
   else
      METIS_LIBS="-lmetis"
      METIS_INCS=""
   fi
   AC_CHECK_LIB([metis],[METIS_PartGraphKway], [have_metis="yes"],
                [have_metis="no"], $METIS_LIBS)
   if test "$have_metis" == "yes"; then
      AC_DEFINE([HAVE_METIS],1,[The METIS library is present.])
   else
      AC_MSG_ERROR("Failed to find a METIS library")
   fi
fi

AC_SUBST([METIS_LIBS])
AC_SUBST([METIS_INCS])
AM_CONDITIONAL([HAVEMETIS],[test -n "$METIS_LIBS"])

# Check for ParMETIS note we can have both as ParMETIS uses METIS.
have_parmetis="no"
AC_ARG_WITH([parmetis],
    [AS_HELP_STRING([--with-parmetis=PATH],
       [root directory where ParMETIS is installed @<:@yes/no@:>@]
    )],
    [with_parmetis="$withval"],
    [with_parmetis="no"]
)

if test "x$with_parmetis" != "xno"; then

# Check if we have ParMETIS.
   if test "x$with_parmetis" != "xyes" -a "x$with_parmetis" != "x"; then
      PARMETIS_LIBS="-L$with_parmetis/lib -lparmetis"
      PARMETIS_INCS="-I$with_parmetis/include"
   else
      PARMETIS_LIBS="-lparmetis"
      PARMETIS_INCS=""
   fi
   AC_CHECK_LIB([parmetis],[ParMETIS_V3_RefineKway], [have_parmetis="yes"],
                [have_parmetis="no"], $PARMETIS_LIBS)
   if test "$have_parmetis" == "no"; then

# A build may use an external METIS library, check for that.

      if test "x$with_parmetis" != "xyes" -a "x$with_parmetis" != "x"; then
         PARMETIS_LIBS="-L$with_parmetis/lib -lparmetis -lmetis"
         PARMETIS_INCS="-I$with_parmetis/include"
      else
         PARMETIS_LIBS="-lparmetis -lmetis"
         PARMETIS_INCS=""
      fi
      AC_CHECK_LIB([parmetis],[ParMETIS_V3_RefineKway], [have_parmetis="yes"],
                   [have_parmetis="no"], [$METIS_LIBS $PARMETIS_LIBS])

   fi
   if test "$have_parmetis" == "yes"; then
      AC_DEFINE([HAVE_PARMETIS],1,[The ParMETIS library is present.])
   else
      AC_MSG_ERROR("Failed to find a ParMETIS library")
   fi
fi

AC_SUBST([PARMETIS_LIBS])
AC_SUBST([PARMETIS_INCS])
AM_CONDITIONAL([HAVEPARMETIS],[test -n "$PARMETIS_LIBS"])

# METIS fixed width integer printing can require this, so define. Only needed
# for some non C99 compilers, i.e. C++ pre C++11.
AH_VERBATIM([__STDC_FORMAT_MACROS],
            [/* Needed to get PRIxxx macros from stdint.h when not using C99 */
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS 1
#endif])

# Check for FFTW. We test for this in the standard directories by default,
# and only disable if using --with-fftw=no or --without-fftw. When a value
# is given FFTW must be found.
# If FFTW is found, we check whether this is the threaded version.
have_fftw="no"
AC_ARG_WITH([fftw],
    [AS_HELP_STRING([--with-fftw=PATH],
       [root directory where fftw is installed @<:@yes/no@:>@]
    )],
    [with_fftw="$withval"],
    [with_fftw="test"]
)
if test "x$with_fftw" != "xno"; then

   # Was FFTW's location specifically given?
   if test "x$with_fftw" != "xyes" -a "x$with_fftw" != "xtest" -a "x$with_fftw" != "x"; then
      FFTW_LIBS="-L$with_fftw/lib -lfftw3"
      FFTW_INCS="-I$with_fftw/include"
   else
      FFTW_LIBS="-lfftw3"
      FFTW_INCS=""
   fi

   #  FFTW is not specified, so just check if we have it.
   if test "x$with_fftw" = "xtest"; then
      AC_CHECK_LIB([fftw3],[fftw_malloc],[have_fftw="yes"],[have_fftw="no"],$FFTW_LIBS)
      if test "x$have_fftw" != "xno"; then
      	 AC_DEFINE([HAVE_FFTW],1,[The FFTW library appears to be present.])
      fi
   # FFTW was specified, check that it was a valid location.
   else
      AC_CHECK_LIB([fftw3],[fftw_malloc],
         AC_DEFINE([HAVE_FFTW],1,[The FFTW library appears to be present.]),
         AC_MSG_ERROR(something is wrong with the FFTW library!), $FFTW_LIBS)
      have_fftw="yes"
   fi

   # FFTW was requested not to be used.
   if test "$have_fftw" = "no"; then
      FFTW_LIBS=""
      FFTW_INCS=""
   fi

   # Now, check whether we have the threaded version of FFTW
   if test "x$have_fftw" = "xyes"; then

      # Was FFTW's location specifically given?
      if test "x$with_fftw" != "xyes" -a "x$with_fftw" != "xtest" -a "x$with_fftw" != "x"; then
        FFTW_THREADED_LIBS="-L$with_fftw/lib -lfftw3_threads -lfftw3"
        FFTW_THREADED_INCS="-I$with_fftw/include"
      else
        FFTW_THREADED_LIBS="-lfftw3_threads -lfftw3"
        FFTW_THREADED_INCS=""
      fi

      # Verify that the library is threaded
      AC_CHECK_LIB([fftw3],[fftw_init_threads],[have_threaded_fftw="yes"],
		   [have_threaded_fftw="no"], $FFTW_THREADED_LIBS)

      # If found, update things
      if test "x$have_threaded_fftw" = "xyes"; then
         AC_DEFINE([HAVE_THREADED_FFTW],1,[The threaded FFTW library appears to be present.])
         FFTW_LIBS=$FFTW_THREADED_LIBS
         FFTW_INCS=$FFTW_THREADED_INCS
	 have_fftw="yes - threaded"
      fi
   fi
fi
AC_ARG_WITH([arm-fftw],
    [AS_HELP_STRING([--with-arm-fftw=PATH],
      [root directory where arm fft library is installed @<:@yes/no@:>@]
    )],
    [with_arm_fftw="$withval"],
    [with_arm_fftw=no]
)
if test "x$with_arm_fftw" != "xno"; then

   # Was FFTW's location specifically given?
   if test "x$with_arm_fftw" != "xyes" -a "x$with_arm_fftw" != "xtest" -a "x$with_arm_fftw" != "x"; then
      FFTW_LIBS="-L$with_arm_fftw/lib -larmpl_lp64"
      FFTW_INCS="-I$with_arm_fftw/include"
   else
      FFTW_LIBS="-larmpl_lp64"
      FFTW_INCS=""
   fi

   #  FFTW is not specified, so just check if we have it.
   if test "x$with_arm_fftw" = "xtest"; then
      AC_CHECK_LIB([armpl_lp64],[fftw_malloc],[have_fftw="yes"],[have_fftw="no"],$FFTW_LIBS)
      if test "x$have_arm_fftw" != "xno"; then
      	 AC_DEFINE([HAVE_FFTW],1,[The FFTW library appears to be present.])
	 have_fftw="yes - ARM"
      fi
   # FFTW was specified, check that it was a valid location.
   else
      AC_CHECK_LIB([armpl_lp64],[fftw_malloc],
         AC_DEFINE([HAVE_FFTW],1,[The FFTW library appears to be present.]),
         AC_MSG_ERROR(something is wrong with the FFTW library!), $FFTW_LIBS)
      have_fftw="yes - ARM"
   fi

   # FFTW was requested not to be used.
   if test "$have_arm_fftw" = "no"; then
      FFTW_LIBS=""
      FFTW_INCS=""
   fi

   # Now, check whether we have the threaded version of FFTW
   if test "x$have_arm_fftw" = "xyes"; then

      # Was FFTW's location specifically given?
      if test "x$with_arm_fftw" != "xyes" -a "x$with_arm_fftw" != "xtest" -a "x$with_arm_fftw" != "x"; then
        FFTW_THREADED_LIBS="-L$with_arm_fftw/lib -larmpl_lp64_threads -larmpl_lp64"
        FFTW_THREADED_INCS="-I$with_arm_fftw/include"
      else
        FFTW_THREADED_LIBS="-larmpl_lp64_threads -larmpl_lp64"
        FFTW_THREADED_INCS=""
      fi

      # Verify that the library is threaded
      AC_CHECK_LIB([armpl_lp64],[fftw_init_threads],[have_threaded_fftw="yes"],
                  [have_threaded_fftw="no"], $FFTW_THREADED_LIBS)

      # If found, update things
      if test "x$have_threaded_fftw" = "xyes"; then
         AC_DEFINE([HAVE_THREADED_FFTW],1,[The threaded FFTW library appears to be present.])
         FFTW_LIBS=$FFTW_THREADED_LIBS
         FFTW_INCS=$FFTW_THREADED_INCS
         have_fftw="yes - ARM - threaded"
      fi
   fi
fi
AC_SUBST([FFTW_LIBS])
AC_SUBST([FFTW_INCS])
AM_CONDITIONAL([HAVEFFTW],[test -n "$FFTW_LIBS"])

# Check for GSL. We test for this in the standard directories by default,
# and only disable if using --with-gsl=no or --without-gsl. When a value
# is given GSL must be found.
have_gsl="no"
AC_ARG_WITH([gsl],
    [AS_HELP_STRING([--with-gsl=PATH],
       [root directory where GSL is installed @<:@yes/no@:>@]
    )],
    [with_gsl="$withval"],
    [with_gsl="test"]
)
if test "x$with_gsl" != "xno"; then
   if test "x$with_gsl" != "xyes" -a "x$with_gsl" != "xtest" -a "x$with_gsl" != "x"; then
      GSL_LIBS="-L$with_gsl/lib -lgsl -lgslcblas"
      GSL_INCS="-I$with_gsl/include"
   else
      GSL_LIBS="-lgsl -lgslcblas"
      GSL_INCS=""
   fi
   #  GSL is not specified, so just check if we have it.
   if test "x$with_gsl" = "xtest"; then
      AC_CHECK_LIB([gslcblas],[cblas_dgemm],[have_gsl="yes"],[have_gsl="no"],$GSL_LIBS)
      if test "x$have_gsl" != "xno"; then
         AC_DEFINE([HAVE_LIBGSLCBLAS],1,[The GSL CBLAS library appears to be present.])
         AC_CHECK_LIB([gsl],[gsl_integration_qag],
            AC_DEFINE([HAVE_LIBGSL],1,[The GSL library appears to be present.]),
            [have_gsl="no"],$GSL_LIBS)
      fi
   else
      AC_CHECK_LIB([gslcblas],[cblas_dgemm],
         AC_DEFINE([HAVE_LIBGSLCBLAS],1,[The GSL CBLAS library appears to be present.]),
         AC_MSG_ERROR(something is wrong with the GSL CBLAS library!), $GSL_LIBS)
      AC_CHECK_LIB([gsl],[gsl_integration_qag],
         AC_DEFINE([HAVE_LIBGSL],1,[The GSL library appears to be present.]),
         AC_MSG_ERROR(something is wrong with the GSL library!), $GSL_LIBS)
      have_gsl="yes"
   fi
   if test "$have_gsl" = "no"; then
      GSL_LIBS=""
      GSL_INCS=""
   fi
fi
AC_SUBST([GSL_LIBS])
AC_SUBST([GSL_INCS])
AM_CONDITIONAL([HAVEGSL],[test -n "$GSL_LIBS"])

# If available check for NUMA as well. There is a problem with the headers of
# this library, mainly that they do not pass the strict prototypes check when
# installed outside of the system directories. So we actually do this check
# in two phases. The basic ones first (before strict-prototypes is added to CFLAGS).
have_numa="no"
AC_ARG_WITH([numa],
    [AS_HELP_STRING([--with-numa=PATH],
       [Directory where the NUMA library exists @<:@yes/no@:>@]
    )],
    [with_numa="$withval"],
    [with_numa="yes"]
)
if test "$ac_cv_func_pthread_setaffinity_np" = "yes" -a "x$with_numa" != "xno"; then

    if test "x$with_numa" != "xyes" -a "x$with_numa" != "x"; then
        NUMA_LIBS="-L$with_numa/lib -lnuma"
        NUMA_INCS="-I$with_numa/include"
    else
        NUMA_LIBS="-lnuma"
        NUMA_INCS=""
    fi

    #  Test for header file.
    old_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $NUMA_INCS"
    AC_CHECK_HEADER([numa.h])
    CPPFLAGS="$old_CPPFLAGS"
    if test "$ac_cv_header_numa_h" = "yes"; then

        #  If NUMA location is specified check if we have it.
        if test "x$with_numa" != "xyes" -a "x$with_numa" != "x"; then
            AC_CHECK_LIB([numa],[numa_available],
                AC_DEFINE([HAVE_LIBNUMA],1,[The NUMA library appears to be present.]),
                AC_MSG_ERROR(something is wrong with the NUMA library!), $NUMA_LIBS)
            have_numa="yes"
        else
            AC_CHECK_LIB([numa],[numa_available],[have_numa="yes"],[have_numa="no"],$NUMA_LIBS)
            if test "x$have_numa" != "xno"; then
                AC_DEFINE([HAVE_LIBNUMA],1,[The NUMA library appears to be present.])
            fi
        fi
    fi

    #  We can live without this.
    if test "$have_numa" = "no"; then
       NUMA_LIBS=""
    fi
fi
AC_SUBST([NUMA_LIBS])

# Check for Intel and PowerPC intrinsics header optionally used by vector.h.
AC_CHECK_HEADERS([immintrin.h], [], [],
[#ifdef HAVE_IMMINTRIN_H
# include <immintrin.h>
#endif
])
AC_CHECK_HEADERS([altivec.h], [], [],
[#ifdef HAVE_ALTIVEC_H
# include <altivec.h>
#endif
])

# Check for timing functions needed by cycle.h.
AC_HEADER_TIME
AC_CHECK_HEADERS([sys/time.h c_asm.h intrinsics.h mach/mach_time.h])
AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t
is defined in <sys/time.h>])],,
[#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif])
AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time])
AC_MSG_CHECKING([for _rtc intrinsic])
rtc_ok=yes
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#ifdef HAVE_INTRINSICS_H
#include <intrinsics.h>
#endif]],
[[_rtc()]])],
[AC_DEFINE(HAVE__RTC,1,[Define if you have the UNICOS _rtc() intrinsic.])],[rtc_ok=no])
AC_MSG_RESULT($rtc_ok)

# Check whether we have any of the ARM v8.1 tick timers
AX_ASM_ARM_PMCCNTR
AX_ASM_ARM_CNTVCT

# Special timers for the ARM v7 platforms (taken from FFTW-3 to match their cycle.h)
AC_ARG_ENABLE(armv7a-cntvct, [AC_HELP_STRING([--enable-armv7a-cntvct],[enable the cycle counter on Armv7a via the CNTVCT register])], have_armv7acntvct=$enableval)
if test "$have_armv7acntvct"x = "yes"x; then
	AC_DEFINE(HAVE_ARMV7A_CNTVCT,1,[Define if you have enabled the CNTVCT cycle counter on ARMv7a])
fi

AC_ARG_ENABLE(armv7a-pmccntr, [AC_HELP_STRING([--enable-armv7a-pmccntr],[enable the cycle counter on Armv7a via the PMCCNTR register])], have_armv7apmccntr=$enableval)
if test "$have_armv7apmccntr"x = "yes"x; then
	AC_DEFINE(HAVE_ARMV7A_PMCCNTR,1,[Define if you have enabled the PMCCNTR cycle counter on ARMv7a])
fi


# Second part of the NUMA library checks. We now decide if we need to use
# -isystem to get around the strict-prototypes problem. Assumes isystem
# is available when strict-prototypes is.
if test "$have_numa" != "no"; then
    if test "x$with_numa" != "xyes" -a "x$with_numa" != "x"; then
        case "$CFLAGS" in
            *strict-prototypes*)
                NUMA_INCS="-isystem$with_numa/include"
                # This may still fail if CPATH is used, so we check if the
                # headers are usable.
                AS_UNSET(ac_cv_header_numa_h)
                old_CPPFLAGS="$CPPFLAGS"
                CPPFLAGS="$CPPFLAGS $NUMA_INCS"
                numa_failed="no"
                AC_CHECK_HEADER([numa.h],[numa_failed="no"],
                                [numa_failed="yes"])
                if test "$numa_failed" = "yes"; then
                    AC_MSG_ERROR([Failed to compile the numa.h header file: you may need to set --enable-compiler-warnings to yes or no])
                fi
                CPPFLAGS="$old_CPPFLAGS"
            ;;
            *)
                NUMA_INCS="-I$with_numa/include"
            ;;
        esac
   fi
fi
AC_SUBST([NUMA_INCS])
])
