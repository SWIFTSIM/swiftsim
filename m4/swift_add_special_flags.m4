# SYNOPSIS
#
#   SWIFT_ADD_SPECIAL_FLAGS
#
# DESCRIPTION
#
# This macro defines a serie of special flags used in SWIFT.
# It includes the logger, MPI, debugging flags and vectorization.
#
# LICENSE
#
# This file is part of SWIFT.
# Copyright (C) 2012 pedro.gonnet@durham.ac.uk.
#               2016 p.w.draper@durham.ac.uk.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


AC_DEFUN([SWIFT_ADD_SPECIAL_FLAGS],
[
# If debug is selected then we also define SWIFT_DEVELOP_MODE to control
# any developer code options.
if test "x$ax_enable_debug" != "xno"; then
   AC_DEFINE([SWIFT_DEVELOP_MODE],1,[Enable developer code options])
fi

# logger
AC_ARG_ENABLE([logger],
	[AS_HELP_STRING([--enable-logger],
		[enable the particle logger]
	)],
	[with_logger="${enableval}"],
	[with_logger="no"]
)

if test "$with_logger" = "yes"; then
   AC_DEFINE([WITH_LOGGER], 1, [logger enabled])
fi

# Interprocedural optimization support. Needs special handling for linking and
# archiving as well as compilation with Intels, needs to be done before
# libtool is configured (to use correct LD).
AC_ARG_ENABLE([ipo],
   [AS_HELP_STRING([--enable-ipo],
     [Enable interprocedural optimization @<:@no/yes@:>@]
   )],
   [enable_ipo="$enableval"],
   [enable_ipo="no"]
)

if test "$enable_ipo" = "yes"; then
   if test "$ax_cv_c_compiler_vendor" = "intel"; then
      CFLAGS="$CFLAGS -ip -ipo"
      LDFLAGS="$LDFLAGS -ipo"
      : ${AR="xiar"}
      : ${LD="xild"}
      AC_MSG_RESULT([added Intel interprocedural optimization support])
   elif test "$ax_cv_c_compiler_vendor" = "gnu"; then
      CFLAGS="$CFLAGS -flto"
      LDFLAGS="$LDFLAGS -flto"
      AX_COMPARE_VERSION($ax_cv_c_compiler_version, [ge], [5.0.0],
                          [
      : ${AR="gcc-ar"}
      : ${RANLIB="gcc-ranlib"}
                          ], [:] )
      AC_MSG_RESULT([added GCC interprocedural optimization support])
   elif test "$ax_cv_c_compiler_vendor" = "clang"; then
      CFLAGS="$CFLAGS -flto"
      : ${AR="llvm-ar"}
      : ${LD="llvm-ld"}
      : ${RANLIB="llvm-ranlib"}
      AC_MSG_RESULT([added LLVM interprocedural optimization support])
   else
      AC_MSG_WARN([Compiler does not support interprocedural optimization])
   fi
fi

# Check for MPI. Need to do this before characterising the compiler (C99 mode),
# as this changes the compiler.
# We should consider using AX_PROG_CC_MPI to replace AC_PROG_CC when compiling
# whole applications. There are issues with mixing compilers when using this
# macro. See
# http://lists.gnu.org/archive/html/autoconf-archive-maintainers/2011-05/msg00004.html.
AC_ARG_ENABLE([mpi],
    [AS_HELP_STRING([--enable-mpi],
      [Compile with functionality for distributed-memory parallelism using MPI @<:@yes/no@:>@]
    )],
    [enable_mpi="$enableval"],
    [enable_mpi="yes"]
)
good_mpi="yes"
if test "$enable_mpi" = "yes"; then
    AX_MPI([CC="$MPICC" AC_DEFINE(HAVE_MPI, 1, [Define if you have the MPI library.]) ])
    MPI_LIBRARY="Unknown MPI"

    # Various MPI implementations require additional libraries when also using
    # threads. Use mpirun (on PATH) as that seems to be only command with
    # version flag, allow MPIRUN to override for systems that insist on
    # a non-standard name (PRACE).
    : ${MPIRUN='mpirun'}
    if test "$MPIRUN" = "mpirun"; then
       AC_PATH_PROG([MPIRUN],[mpirun],[notfound])
    fi
    if test "$MPIRUN" = "notfound"; then
       # This may not be fatal (some systems do not allow mpirun on
       # development nodes)., so push on.
       AC_MSG_WARN([Cannot find mpirun command on PATH, thread support may not be correct])
    else
       # Special options we know about.
       # Intel: -mt_mpi
       # PLATFORM: -lmtmpi
       # OpenMPI: nothing, but library should be built correctly.
       # Set MPI_THREAD_LIBS and add to linker commands as necessary.
       AC_MSG_CHECKING([MPI threads options])
       version=`$MPIRUN -version 2>&1`
       case "$version" in
         *Intel*MPI*)
            MPI_THREAD_LIBS="-mt_mpi"
            MPI_LIBRARY="Intel MPI"
            AC_MSG_RESULT([Intel MPI])
         ;;
         *Platform*)
            MPI_THREAD_LIBS="-lmtmpi"
            MPI_LIBRARY="PLATFORM MPI"
            AC_MSG_RESULT([PLATFORM MPI])
         ;;
         *"Open MPI"*)
            MPI_THREAD_LIBS=""
            MPI_LIBRARY="Open MPI"
            AC_MSG_RESULT([Open MPI])
            #  OpenMPI should be 1.8.6 or later, if not complain.
            #  Version is last word on first line of -version output.
            revision=`mpirun -version 2>&1 | grep "Open MPI" | awk '{print $NF}'`
            AX_COMPARE_VERSION( $revision, [ge], [1.8.6],,[good_mpi="no"] )
            if test "$good_mpi" = "no"; then
                AC_MSG_WARN([
    Open MPI version should be at least 1.8.6 (is $revision)])
                enable_mpi="yes (but with warning)"
            fi
         ;;
         *)
            MPI_THREAD_LIBS=""
            AC_MSG_RESULT([unknown])
         ;;
       esac
       AC_SUBST([MPI_THREAD_LIBS])
    fi
    AC_DEFINE_UNQUOTED([SWIFT_MPI_LIBRARY], ["$MPI_LIBRARY"], [The MPI library name, if known.])
fi
AM_CONDITIONAL([HAVEMPI],[test $enable_mpi = "yes"])

# Indicate that MPIRUN can be modified by an environment variable
AC_ARG_VAR(MPIRUN, Path to the mpirun command if non-standard)
# If debugging try to show inlined functions.
if test "x$enable_debug" = "xyes"; then
   #  Show inlined functions.
   if test "$ax_cv_c_compiler_vendor" = "gnu"; then
      # Would like to use -gdwarf and let the compiler pick a good version
      # but that doesn't always work.
      AX_CHECK_COMPILE_FLAG([-gdwarf -fvar-tracking-assignments],
        [inline_EXTRA_FLAGS="-gdwarf -fvar-tracking-assignments"],
        [inline_EXTRA_FLAGS="-gdwarf-2 -fvar-tracking-assignments"])
      CFLAGS="$CFLAGS $inline_EXTRA_FLAGS"
   elif test "$ax_cv_c_compiler_vendor" = "intel"; then
      CFLAGS="$CFLAGS -debug inline-debug-info"
   fi
fi

# Check if task debugging is on.
AC_ARG_ENABLE([task-debugging],
   [AS_HELP_STRING([--enable-task-debugging],
     [Store extra information for generating task dump files @<:@yes/no@:>@]
   )],
   [enable_task_debugging="$enableval"],
   [enable_task_debugging="no"]
)
if test "$enable_task_debugging" = "yes"; then
   AC_DEFINE([SWIFT_DEBUG_TASKS],1,[Enable task debugging])
fi

# Check if threadpool debugging is on.
AC_ARG_ENABLE([threadpool-debugging],
   [AS_HELP_STRING([--enable-threadpool-debugging],
     [Store threadpool mapper timing information and generate threadpool dump files @<:@yes/no@:>@]
   )],
   [enable_threadpool_debugging="$enableval"],
   [enable_threadpool_debugging="no"]
)
if test "$enable_threadpool_debugging" = "yes"; then
   AC_DEFINE([SWIFT_DEBUG_THREADPOOL],1,[Enable threadpool debugging])
   LDFLAGS="$LDFLAGS -rdynamic"
fi

# Check if the general timers are switched on.
AC_ARG_ENABLE([timers],
   [AS_HELP_STRING([--enable-timers],
     [Activate the basic timers @<:@yes/no@:>@]
   )],
   [enable_timers="$enableval"],
   [enable_timers="no"]
)
if test "$enable_timers" = "yes"; then
   AC_DEFINE([SWIFT_USE_TIMERS],1,[Enable individual timers])
fi

# Check if expensive debugging is on.
AC_ARG_ENABLE([debugging-checks],
   [AS_HELP_STRING([--enable-debugging-checks],
     [Activate expensive consistency checks @<:@yes/no@:>@]
   )],
   [enable_debugging_checks="$enableval"],
   [enable_debugging_checks="no"]
)
if test "$enable_debugging_checks" = "yes"; then
   AC_DEFINE([SWIFT_DEBUG_CHECKS],1,[Enable expensive debugging])
fi

# Check if using our custom icbrtf is enalbled.
AC_ARG_ENABLE([custom-icbrtf],
   [AS_HELP_STRING([--enable-custom-icbrtf],
     [Use SWIFT's custom icbrtf function instead of the system cbrtf @<:@yes/no@:>@]
   )],
   [enable_custom_icbrtf="$enableval"],
   [enable_custom_icbrtf="no"]
)
if test "$enable_custom_icbrtf" = "yes"; then
   AC_DEFINE([WITH_ICBRTF],1,[Enable custom icbrtf])
fi

# Check whether we want to default to naive cell interactions
AC_ARG_ENABLE([naive-interactions],
   [AS_HELP_STRING([--enable-naive-interactions],
     [Activate use of naive cell interaction functions @<:@yes/no@:>@]
   )],
   [enable_naive_interactions="$enableval"],
   [enable_naive_interactions="no"]
)
if test "$enable_naive_interactions" = "yes"; then
   AC_DEFINE([SWIFT_USE_NAIVE_INTERACTIONS],1,[Enable use of naive cell interaction functions])
fi

# Check if gravity force checks are on for some particles.
AC_ARG_ENABLE([gravity-force-checks],
   [AS_HELP_STRING([--enable-gravity-force-checks],
     [Activate expensive brute-force gravity checks for a fraction 1/N of all particles @<:@N@:>@]
   )],
   [gravity_force_checks="$enableval"],
   [gravity_force_checks="no"]
)
if test "$gravity_force_checks" == "yes"; then
   AC_MSG_ERROR(Need to specify the fraction of particles to check when using --enable-gravity-force-checks!)
elif test "$gravity_force_checks" != "no"; then
   AC_DEFINE_UNQUOTED([SWIFT_GRAVITY_FORCE_CHECKS], [$enableval] ,[Enable gravity brute-force checks])
fi

# Check whether we want to switch on glass making
AC_ARG_ENABLE([glass-making],
   [AS_HELP_STRING([--enable-glass-making],
     [Activate the glass-making procedure by reversing the sign of gravity @<:@yes/no@:>@]
   )],
   [gravity_glass_making="$enableval"],
   [gravity_glass_making="no"]
)
if test "$gravity_glass_making" == "yes"; then
   AC_DEFINE([SWIFT_MAKE_GRAVITY_GLASS], 1, [Make the code run in a way to produce a glass file for gravity/cosmology])
fi

# Only optimize if allowed, otherwise assume user will set CFLAGS as
# appropriate.
AC_ARG_ENABLE([optimization],
   [AS_HELP_STRING([--enable-optimization],
     [Enable compile time optimization flags for host @<:@yes/no@:>@]
   )],
   [enable_opt="$enableval"],
   [enable_opt="yes"]
)

#  Disable vectorisation for known compilers. This switches off optimizations
#  that could be enabled above, so in general should be appended. Slightly odd
#  implementation as want to describe as --disable-vec, but macro is enable
#  (there is no enable action).
AC_ARG_ENABLE([vec],
   [AS_HELP_STRING([--disable-vec],
     [Disable vectorization]
   )],
   [enable_vec="$enableval"],
   [enable_vec="yes"]
)

#  Disable hand written vectorisation. Slightly odd implementation as want
# to describe as --disable-hand-vec, but macro is enable (there is no enable action).
AC_ARG_ENABLE([hand-vec],
   [AS_HELP_STRING([--disable-hand-vec],
     [Disable intrinsic vectorization]
   )],
   [enable_hand_vec="$enableval"],
   [enable_hand_vec="yes"]
)

HAVEVECTORIZATION=0

if test "$enable_opt" = "yes" ; then

   # Add code optimisation flags and tuning to host. This is a funny macro
   # that does not like CFLAGS being already set. Work around that as we have
   # at least set it to "", so it is set.
   ac_test_CFLAGS="no"
   old_CFLAGS="$CFLAGS"
   AX_CC_MAXOPT
   ac_test_CFLAGS="yes"
   CFLAGS="$old_CFLAGS $CFLAGS"

   # Check SSE & AVX support (some overlap with AX_CC_MAXOPT).
   # Don't use the SIMD_FLAGS result with Intel compilers. The -x<code>
   # value from AX_CC_MAXOPT should be sufficient.
   AX_EXT
   if test "$SIMD_FLAGS" != ""; then
       if test "$ax_cv_c_compiler_vendor" != "intel"; then
           CFLAGS="$CFLAGS $SIMD_FLAGS"
       fi
   fi

   if test "$enable_vec" = "no"; then
      if test "$ax_cv_c_compiler_vendor" = "intel"; then
      	 CFLAGS="$CFLAGS -no-vec -no-simd"
      	 AC_MSG_RESULT([disabled Intel vectorization])
      elif test "$ax_cv_c_compiler_vendor" = "gnu"; then
      	 CFLAGS="$CFLAGS -fno-tree-vectorize"
      	 AC_MSG_RESULT([disabled GCC vectorization])
      elif test "$ax_cv_c_compiler_vendor" = "clang"; then
         CFLAGS="$CFLAGS -fno-vectorize -fno-slp-vectorize"
         AC_MSG_RESULT([disabled clang vectorization])
      else
         AC_MSG_WARN([Do not know how to disable vectorization for this compiler])
      fi
   elif test "$enable_hand_vec" = "yes"; then
      AC_DEFINE([WITH_VECTORIZATION],1,[Enable hand-written vectorization])
      HAVEVECTORIZATION=1
   fi
fi
AM_CONDITIONAL([HAVEVECTORIZATION],[test -n "$HAVEVECTORIZATION"])


# Add address sanitizer options to flags, if requested. Only useful for GCC
# version 4.8 and later and clang.
AC_ARG_ENABLE([sanitizer],
   [AS_HELP_STRING([--enable-sanitizer],
     [Enable memory error detection using address sanitizer @<:@no/yes@:>@]
   )],
   [enable_san="$enableval"],
   [enable_san="no"]
)

if test "$enable_san" = "yes"; then
   if test "$ax_cv_c_compiler_vendor" = "gnu"; then
      AX_COMPARE_VERSION( $ax_cv_c_compiler_version, [ge], [4.8.0],
                          [enable_san="yes"], [enable_san="no"] )
   elif test "$ax_cv_c_compiler_vendor" = "clang"; then
      AX_COMPARE_VERSION( $ax_cv_c_compiler_version, [ge], [3.2.0],
                          [enable_san="yes"], [enable_san="no"] )
   fi
   if test "$enable_san" = "yes"; then
      CFLAGS="$CFLAGS -fsanitize=address -fno-omit-frame-pointer"
      AC_MSG_RESULT([added address sanitizer support])
   else
      AC_MSG_WARN([Compiler does not support address sanitizer option])
   fi
fi

# Add the undefined sanitizer option to flags. Only useful for GCC
# version 4.9 and later and clang to detected undefined code behaviour
# such as integer overflow and memory alignment issues.
AC_ARG_ENABLE([undefined-sanitizer],
   [AS_HELP_STRING([--enable-undefined-sanitizer],
     [Enable detection of code that causes undefined behaviour @<:@no/yes@:>@]
   )],
   [enable_ubsan="$enableval"],
   [enable_ubsan="no"]
)

if test "$enable_ubsan" = "yes"; then
   if test "$ax_cv_c_compiler_vendor" = "gnu"; then
      AX_COMPARE_VERSION( $ax_cv_c_compiler_version, [ge], [4.9.0],
                          [enable_ubsan="yes"], [enable_ubsan="no"] )
   elif test "$ax_cv_c_compiler_vendor" = "clang"; then
      AX_COMPARE_VERSION( $ax_cv_c_compiler_version, [ge], [3.7.0],
                          [enable_ubsan="yes"], [enable_ubsan="no"] )
   fi
   if test "$enable_ubsan" = "yes"; then
      CFLAGS="$CFLAGS -fsanitize=undefined"
      AC_MSG_RESULT([added undefined sanitizer support])
   else
      AC_MSG_WARN([Compiler does not support undefined sanitizer option])
   fi
fi

# Check for pthreads.
AX_PTHREAD([LIBS="$PTHREAD_LIBS $LIBS" CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
    CC="$PTHREAD_CC" LDFLAGS="$LDFLAGS $PTHREAD_LIBS $LIBS"],
    AC_MSG_ERROR([Could not find a working version of
    the pthread library. Make sure you have the library and header files installed
    or use CPPFLAGS and LDFLAGS if the library is installed in a
    non-standard location.]))

# Check whether POSIX thread barriers are implemented (e.g. OSX does not have them)
have_pthread_barrier="no"
AC_CHECK_LIB(pthread, pthread_barrier_init,
	     have_pthread_barrier="yes",
	     AC_MSG_WARN(POSIX implementation does not have barriers. SWIFT will use home-made ones.))
if test "x$have_pthread_barrier" == "xyes"; then
  AC_DEFINE([HAVE_PTHREAD_BARRIERS], [1], [The posix library implements barriers])
fi

# Check whether POSIX file allocation functions exist (e.g. OSX does not have them)
AC_CHECK_LIB(pthread, posix_fallocate,
	     AC_DEFINE([HAVE_POSIX_FALLOCATE], [1], [The posix library implements file allocation functions.]),
	     AC_MSG_WARN(POSIX implementation does not have file allocation functions.))

#  Check for -lprofiler usually part of the gperftools along with tcmalloc.
have_profiler="no"
AC_ARG_WITH([profiler],
   [AS_HELP_STRING([--with-profiler=PATH],
      [use cpu profiler library or specify the directory with lib @<:@yes/no@:>@]
   )],
   [with_profiler="$withval"],
   [with_profiler="no"]
)
if test "x$with_profiler" != "xno"; then
   if test "x$with_profiler" != "xyes" -a "x$with_profiler" != "x"; then
      proflibs="-L$with_profiler -lprofiler"
   else
      proflibs="-lprofiler"
   fi
   AC_CHECK_LIB([profiler],[ProfilerFlush],
    [have_profiler="yes"
      AC_DEFINE([WITH_PROFILER],1,[Link against the gperftools profiling library.])],
    [have_profiler="no"], $proflibs)

   if test "$have_profiler" = "yes"; then
      PROFILER_LIBS="$proflibs"
   else
      PROFILER_LIBS=""
   fi
fi
AC_SUBST([PROFILER_LIBS])
AM_CONDITIONAL([HAVEPROFILER],[test -n "$PROFILER_LIBS"])

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

# Check for floating-point execeptions
AC_CHECK_FUNC(feenableexcept, AC_DEFINE([HAVE_FE_ENABLE_EXCEPT],[1],
    [Defined if the floating-point exception can be enabled using non-standard GNU functions.]))

# Check for setaffinity.
AC_CHECK_FUNC(pthread_setaffinity_np, AC_DEFINE([HAVE_SETAFFINITY],[1],
    [Defined if pthread_setaffinity_np exists.]) )
AM_CONDITIONAL(HAVESETAFFINITY,
    [test "$ac_cv_func_pthread_setaffinity_np" = "yes"])

# Add warning flags by default, if these can be used. Option =error adds
# -Werror to GCC, clang and Intel.  Note do this last as compiler tests may
# become errors, if that's an issue don't use CFLAGS for these, use an AC_SUBST().
AC_ARG_ENABLE([compiler-warnings],
   [AS_HELP_STRING([--enable-compiler-warnings],
     [Enable compile time warning flags, if compiler is known @<:@error/no/yes)@:>@]
   )],
   [enable_warn="$enableval"],
   [enable_warn="error"]
)
if test "$enable_warn" != "no"; then

    # AX_CFLAGS_WARN_ALL does not give good warning flags for the Intel compiler
    # We will do this by hand instead and only default to the macro for unknown compilers
    case "$ax_cv_c_compiler_vendor" in
          gnu | clang)
             CFLAGS="$CFLAGS -Wall -Wextra -Wno-unused-parameter -Wshadow"
          ;;
	  intel)
             CFLAGS="$CFLAGS -w2 -Wunused-variable -Wshadow"
          ;;
	  *)
	     AX_CFLAGS_WARN_ALL
	  ;;
    esac

    # Add a "choke on warning" flag if it exists
    if test "$enable_warn" = "error"; then
       case "$ax_cv_c_compiler_vendor" in
          intel | gnu | clang)
             CFLAGS="$CFLAGS -Werror"
          ;;
       esac
    fi

    # We want strict-prototypes, but this must still work even if warnings
    # are an error.
    AX_CHECK_COMPILE_FLAG([-Wstrict-prototypes],[CFLAGS="$CFLAGS -Wstrict-prototypes"],
                          [CFLAGS="$CFLAGS"],[$CFLAGS],[AC_LANG_SOURCE([int main(void){return 0;}])])
fi

])
