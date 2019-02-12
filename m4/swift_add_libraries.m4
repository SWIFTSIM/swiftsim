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
AC_CHECK_HEADERS([immintrin.h])
AC_CHECK_HEADERS([altivec.h])

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

# Special timers for the ARM v7 and ARM v8 platforms (taken from FFTW-3 to match their cycle.h)
AC_ARG_ENABLE(armv8-pmccntr-el0, [AC_HELP_STRING([--enable-armv8-pmccntr-el0],[enable the cycle counter on ARMv8 via the PMCCNTR_EL0 register])], have_armv8pmccntrel0=$enableval)
if test "$have_armv8pmccntrel0"x = "yes"x; then
	AC_DEFINE(HAVE_ARMV8_PMCCNTR_EL0,1,[Define if you have enabled the PMCCNTR_EL0 cycle counter on ARMv8])
fi

AC_ARG_ENABLE(armv8-cntvct-el0, [AC_HELP_STRING([--enable-armv8-cntvct-el0],[enable the cycle counter on ARMv8 via the CNTVCT_EL0 register])], have_armv8cntvctel0=$enableval)
if test "$have_armv8cntvctel0"x = "yes"x; then
	AC_DEFINE(HAVE_ARMV8_CNTVCT_EL0,1,[Define if you have enabled the CNTVCT_EL0 cycle counter on ARMv8])
fi

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


# Check for special allocators
have_special_allocator="no"

#  Check for tcmalloc a fast malloc that is part of the gperftools.
have_tcmalloc="no"
AC_ARG_WITH([tcmalloc],
   [AS_HELP_STRING([--with-tcmalloc=PATH],
      [use tcmalloc library or specify the directory with lib @<:@yes/no@:>@]
   )],
   [with_tcmalloc="$withval"],
   [with_tcmalloc="no"]
)
if test "x$with_tcmalloc" != "xno" -a "x$have_special_allocator" != "xno"; then
   AC_MSG_ERROR("Cannot activate more than one alternative malloc library")
fi

if test "x$with_tcmalloc" != "xno"; then
   if test "x$with_tcmalloc" != "xyes" -a "x$with_tcmalloc" != "x"; then
      tclibs="-L$with_tcmalloc -ltcmalloc"
   else
      tclibs="-ltcmalloc"
   fi
   AC_CHECK_LIB([tcmalloc],[tc_cfree],[have_tcmalloc="yes"],[have_tcmalloc="no"],
                $tclibs)

   #  Could just have the minimal version.
   if test "$have_tcmalloc" = "no"; then
      if test "x$with_tcmalloc" != "xyes" -a "x$with_tcmalloc" != "x"; then
         tclibs="-L$with_tcmalloc -ltcmalloc_minimal"
      else
         tclibs="-ltcmalloc_minimal"
      fi
      AC_CHECK_LIB([tcmalloc],[tc_cfree],[have_tcmalloc="yes"],[have_tcmalloc="no"],
                   $tclibs)
   fi

   if test "$have_tcmalloc" = "yes"; then
      TCMALLOC_LIBS="$tclibs"

      AC_DEFINE([HAVE_TCMALLOC],1,[The tcmalloc library appears to be present.])

      have_special_allocator="tcmalloc"

      # Prevent compilers that replace the calls with built-ins (GNU 99) from doing so.
      case "$ax_cv_c_compiler_vendor" in
        intel | gnu | clang)
             CFLAGS="$CFLAGS -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free"
          ;;
      esac

   else
      TCMALLOC_LIBS=""
   fi
fi
AC_SUBST([TCMALLOC_LIBS])
AM_CONDITIONAL([HAVETCMALLOC],[test -n "$TCMALLOC_LIBS"])

#  Check for jemalloc another fast malloc that is good with contention.
have_jemalloc="no"
AC_ARG_WITH([jemalloc],
   [AS_HELP_STRING([--with-jemalloc=PATH],
      [use jemalloc library or specify the directory with lib @<:@yes/no@:>@]
   )],
   [with_jemalloc="$withval"],
   [with_jemalloc="no"]
)
if test "x$with_jemalloc" != "xno" -a "x$have_special_allocator" != "xno"; then
   AC_MSG_ERROR("Cannot activate more than one alternative malloc library")
fi

if test "x$with_jemalloc" != "xno"; then
   if test "x$with_jemalloc" != "xyes" -a "x$with_jemalloc" != "x"; then
      jelibs="-L$with_jemalloc -ljemalloc"
   else
      jelibs="-ljemalloc"
   fi
   AC_CHECK_LIB([jemalloc],[malloc_usable_size],[have_jemalloc="yes"],[have_jemalloc="no"],
                $jelibs)

   if test "$have_jemalloc" = "yes"; then
      JEMALLOC_LIBS="$jelibs"

      AC_DEFINE([HAVE_JEMALLOC],1,[The jemalloc library appears to be present.])

      have_special_allocator="jemalloc"

      # Prevent compilers that replace the regular calls with built-ins (GNU 99) from doing so.
      case "$ax_cv_c_compiler_vendor" in
        intel | gnu | clang)
             CFLAGS="$CFLAGS -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free"
          ;;
      esac

   else
      JEMALLOC_LIBS=""
   fi
fi
AC_SUBST([JEMALLOC_LIBS])
AM_CONDITIONAL([HAVEJEMALLOC],[test -n "$JEMALLOC_LIBS"])

#  Check for tbbmalloc, Intel's fast and parallel allocator
have_tbbmalloc="no"
AC_ARG_WITH([tbbmalloc],
   [AS_HELP_STRING([--with-tbbmalloc=PATH],
      [use tbbmalloc library or specify the directory with lib @<:@yes/no@:>@]
   )],
   [with_tbbmalloc="$withval"],
   [with_tbbmalloc="no"]
)
if test "x$with_tbbmalloc" != "xno" -a "x$have_special_allocator" != "xno"; then
   AC_MSG_ERROR("Cannot activate more than one alternative malloc library")
fi

if test "x$with_tbbmalloc" != "xno"; then
   if test "x$with_tbbmalloc" != "xyes" -a "x$with_tbbmalloc" != "x"; then
      tbblibs="-L$with_tbbmalloc -ltbbmalloc_proxy -ltbbmalloc"
   else
      tbblibs="-ltbbmalloc_proxy -ltbbmalloc"
   fi
   AC_CHECK_LIB([tbbmalloc],[scalable_malloc],[have_tbbmalloc="yes"],[have_tbbmalloc="no"],
                $tbblibs)

   if test "$have_tbbmalloc" = "yes"; then
      TBBMALLOC_LIBS="$tbblibs"

      AC_DEFINE([HAVE_TBBMALLOC],1,[The TBBmalloc library appears to be present.])

      have_special_allocator="TBBmalloc"

      # Prevent compilers that replace the calls with built-ins (GNU 99) from doing so.
      case "$ax_cv_c_compiler_vendor" in
        intel | gnu | clang)
             CFLAGS="$CFLAGS -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free"
          ;;
      esac

   else
      TBBMALLOC_LIBS=""
   fi
fi
AC_SUBST([TBBMALLOC_LIBS])
AM_CONDITIONAL([HAVETBBMALLOC],[test -n "$TBBMALLOC_LIBS"])

# Check for HDF5. This is required.
AX_LIB_HDF5
if test "$with_hdf5" != "yes"; then
    AC_MSG_ERROR([Could not find a working HDF5 library])
fi

# We want to know if this HDF5 supports MPI and whether we should use it.
# The default is to use MPI support if it is available, i.e. this is
# a parallel HDF5.
have_parallel_hdf5="no"
if test "$with_hdf5" = "yes"; then
    AC_ARG_ENABLE([parallel-hdf5],
       [AS_HELP_STRING([--enable-parallel-hdf5],
         [Enable parallel HDF5 library MPI functions if available. @<:@yes/no@:>@]
       )],
       [enable_parallel_hdf5="$enableval"],
       [enable_parallel_hdf5="yes"]
    )

    if test "$enable_parallel_hdf5" = "yes"; then
        AC_MSG_CHECKING([for HDF5 parallel support])

	# Check if the library is capable, the header should define H5_HAVE_PARALLEL.
        old_CPPFLAGS="$CPPFLAGS"
        CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
        AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
        #include "hdf5.h"
        #ifndef H5_HAVE_PARALLEL
        # error macro not defined
        #endif
        ]])], [parallel="yes"], [parallel="no"])
        if test "$parallel" = "yes"; then
            have_parallel_hdf5="yes"
            AC_DEFINE([HAVE_PARALLEL_HDF5],1,[HDF5 library supports parallel access])
        fi
        AC_MSG_RESULT($parallel)
        CPPFLAGS="$old_CPPFLAGS"
    fi
fi
AM_CONDITIONAL([HAVEPARALLELHDF5],[test "$have_parallel_hdf5" = "yes"])

# Check for grackle.
have_grackle="no"
AC_ARG_WITH([grackle],
    [AS_HELP_STRING([--with-grackle=PATH],
       [root directory where grackle is installed @<:@yes/no@:>@]
    )],
    [with_grackle="$withval"],
    [with_grackle="no"]
)
if test "x$with_grackle" != "xno"; then

   if test "x$with_grackle" != "xyes" -a "x$with_grackle" != "x"; then
      GRACKLE_LIBS="-L$with_grackle/lib -lgrackle"
      GRACKLE_INCS="-I$with_grackle/include"
   else
      GRACKLE_LIBS="-lgrackle"
      GRACKLE_INCS=""
   fi

   have_grackle="yes"

   echo $GRACKLE_LIBS

   AC_CHECK_LIB(
      [grackle],
      [initialize_chemistry_data],
      [AC_DEFINE([HAVE_GRACKLE],1,[The GRACKLE library appears to be present.])
        AC_DEFINE([CONFIG_BFLOAT_8],1,[Use doubles in grackle])
      ],
      [AC_MSG_ERROR(Cannot find grackle library!)],
      [$GRACKLE_LIBS])
fi
AC_SUBST([GRACKLE_LIBS])
AC_SUBST([GRACKLE_INCS])
AM_CONDITIONAL([HAVEGRACKLE],[test -n "$GRACKLE_LIBS"])

# Check for VELOCIraptor.
have_velociraptor="no"
AC_ARG_WITH([velociraptor],
    [AS_HELP_STRING([--with-velociraptor=PATH],
       [Directory where velociraptor library exists @<:@yes/no@:>@]
    )],
    [with_velociraptor="$withval"],
    [with_velociraptor="no"]
)
if test "x$with_velociraptor" != "xno"; then
   if test "x$with_velociraptor" != "xyes" -a "x$with_velociraptor" != "x"; then
      VELOCIRAPTOR_LIBS="-L$with_velociraptor -lvelociraptor -lmpi -lstdc++ -lhdf5_cpp"
      CFLAGS="$CFLAGS -fopenmp"
   else
      VELOCIRAPTOR_LIBS=""
   fi

   have_velociraptor="yes"

   AC_CHECK_LIB(
      [velociraptor],
      [InitVelociraptor],
      [AC_DEFINE([HAVE_VELOCIRAPTOR],1,[The VELOCIraptor library appears to be present.])],
      [AC_MSG_ERROR(Cannot find VELOCIraptor library at $with_velociraptor)],
      [$VELOCIRAPTOR_LIBS $HDF5_LDFLAGS $HDF5_LIBS $GSL_LIBS]
   )
fi
AC_SUBST([VELOCIRAPTOR_LIBS])
AM_CONDITIONAL([HAVEVELOCIRAPTOR],[test -n "$VELOCIRAPTOR_LIBS"])

# Check for dummy VELOCIraptor.
AC_ARG_ENABLE([dummy-velociraptor],
    [AS_HELP_STRING([--enable-dummy-velociraptor],
       [Enable dummy velociraptor compilation @<:@yes/no@:>@]
    )],
    [enable_dummy_velociraptor="$enableval"],
    [enable_dummy_velociraptor="no"]
)

if test "$enable_dummy_velociraptor" = "yes"; then
  have_velociraptor="yes"

  AC_DEFINE(HAVE_VELOCIRAPTOR,1,[The VELOCIraptor library appears to be present.])
  AC_DEFINE(HAVE_DUMMY_VELOCIRAPTOR,1,[The dummy VELOCIraptor library is present.])
fi


])
