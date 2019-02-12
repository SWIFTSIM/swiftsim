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

# Init the project.
AC_INIT([SWIFT],[0.8.0],[https://gitlab.cosma.dur.ac.uk/swift/swiftsim])
swift_config_flags="$*"

#  We want to stop when given unrecognised options. No subdirs so this is safe.
enable_option_checking=${enable_option_checking:-fatal}
if test -n "$ac_unrecognized_opts"; then
    case $enable_option_checking in
        no)
        ;;
        fatal)
            { $as_echo "$as_me: error: unrecognized options: $ac_unrecognized_opts" >&2
              { (exit 1); exit 1; }; }
        ;;
        *)
            $as_echo "$as_me: WARNING: unrecognized options: $ac_unrecognized_opts" >&2
        ;;
    esac
fi

AC_COPYRIGHT
AC_CONFIG_SRCDIR([src/space.c])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([subdir-objects])

# Add local macro collection.
AC_CONFIG_MACRO_DIR([m4])

# Stop default CFLAGS from anyone except the environment.
: ${CFLAGS=""}

# Generate header file.
AM_CONFIG_HEADER(config.h)

# Find and test the compiler.
AX_CHECK_ENABLE_DEBUG
AC_PROG_CC
AM_PROG_CC_C_O

# Enable POSIX and platform extension preprocessor macros.
AC_USE_SYSTEM_EXTENSIONS

# Check for compiler version and vendor.
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

#  Restrict support.
AC_C_RESTRICT


# Add libtool support (now that CC is defined).
LT_INIT

# Need C99 and inline support.
AC_PROG_CC_C99
AC_C_INLINE

# Need fortran for grackle and velociraptor
AC_PROG_FC
AC_FC_LIBRARY_LDFLAGS


# Define HAVE_POSIX_MEMALIGN if it works.
AX_FUNC_POSIX_MEMALIGN

# Autoconf stuff.
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_HEADER_STDC

SWIFT_ADD_LIBRARIES

SWIFT_ADD_SPECIAL_FLAGS

# this one should be at the beginning
# it overwrites all the required parameters
SWIFT_ADD_SUBGRID_MODELS



SWIFT_ADD_GRAVITY

SWIFT_ADD_HYDRO

SWIFT_ADD_STARS

SWIFT_ADD_KERNEL

SWIFT_ADD_DIMENSION

SWIFT_ADD_EQUATION_OF_STATE

SWIFT_ADD_CHEMISTRY

SWIFT_ADD_COOLING

SWIFT_ADD_RIEMANN

SWIFT_ADD_TRACER




# Check for git, needed for revision stamps.
AC_PATH_PROG([GIT_CMD], [git])
AC_SUBST([GIT_CMD])

# Make the documentation. Add conditional to handle disable option.
DX_INIT_DOXYGEN(libswift,doc/Doxyfile,doc/)
AM_CONDITIONAL([HAVE_DOXYGEN], [test "$ac_cv_path_ac_pt_DX_DOXYGEN" != ""])

# Check if using EAGLE cooling
AM_CONDITIONAL([HAVEEAGLECOOLING], [test $with_cooling = "EAGLE"])

# Handle .in files.
AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile examples/CoolingRates/Makefile doc/Makefile doc/Doxyfile tests/Makefile])
AC_CONFIG_FILES([argparse/Makefile tools/Makefile])
AC_CONFIG_FILES([tests/testReading.sh], [chmod +x tests/testReading.sh])
AC_CONFIG_FILES([tests/testActivePair.sh], [chmod +x tests/testActivePair.sh])
AC_CONFIG_FILES([tests/test27cells.sh], [chmod +x tests/test27cells.sh])
AC_CONFIG_FILES([tests/test27cellsPerturbed.sh], [chmod +x tests/test27cellsPerturbed.sh])
AC_CONFIG_FILES([tests/test27cellsStars.sh], [chmod +x tests/test27cellsStars.sh])
AC_CONFIG_FILES([tests/test27cellsStarsPerturbed.sh], [chmod +x tests/test27cellsStarsPerturbed.sh])
AC_CONFIG_FILES([tests/test125cells.sh], [chmod +x tests/test125cells.sh])
AC_CONFIG_FILES([tests/test125cellsPerturbed.sh], [chmod +x tests/test125cellsPerturbed.sh])
AC_CONFIG_FILES([tests/testPeriodicBC.sh], [chmod +x tests/testPeriodicBC.sh])
AC_CONFIG_FILES([tests/testPeriodicBCPerturbed.sh], [chmod +x tests/testPeriodicBCPerturbed.sh])
AC_CONFIG_FILES([tests/testInteractions.sh], [chmod +x tests/testInteractions.sh])
AC_CONFIG_FILES([tests/testParser.sh], [chmod +x tests/testParser.sh])
AC_CONFIG_FILES([tests/testSelectOutput.sh], [chmod +x tests/testSelectOutput.sh])
AC_CONFIG_FILES([tests/testFormat.sh], [chmod +x tests/testFormat.sh])

# Save the compilation options
AC_DEFINE_UNQUOTED([SWIFT_CONFIG_FLAGS],["$swift_config_flags"],[Flags passed to configure])

# Make sure the latest git revision string gets included, when we are
# working in a checked out repository.
test -d ${srcdir}/.git && touch ${srcdir}/src/version.c

#  Need to define this, instead of using fifth argument of AC_INIT, until
#  2.64. Defer until now as this redefines PACKAGE_URL, which can emit a
#  compilation error when testing with -Werror.
AC_DEFINE([PACKAGE_URL],["www.swiftsim.com"], [Package web pages])

# Generate output.
AC_OUTPUT

# Report general configuration.
AC_MSG_RESULT([
 ------- Summary --------

   $PACKAGE_NAME v.$PACKAGE_VERSION

   Compiler             : $CC
    - vendor            : $ax_cv_c_compiler_vendor
    - version           : $ax_cv_c_compiler_version
    - flags             : $CFLAGS
   MPI enabled          : $enable_mpi
   HDF5 enabled         : $with_hdf5
    - parallel          : $have_parallel_hdf5
   METIS/ParMETIS       : $have_metis / $have_parmetis
   FFTW3 enabled        : $have_fftw
   GSL enabled          : $have_gsl
   libNUMA enabled      : $have_numa
   GRACKLE enabled      : $have_grackle
   Special allocators   : $have_special_allocator
   CPU profiler         : $have_profiler
   Pthread barriers     : $have_pthread_barrier
   VELOCIraptor enabled : $have_velociraptor
   Particle Logger      : $with_logger

   Hydro scheme       : $with_hydro
   Dimensionality     : $with_dimension
   Kernel function    : $with_kernel
   Equation of state  : $with_eos
   Adiabatic index    : $with_gamma
   Riemann solver     : $with_riemann

   Gravity scheme      : $with_gravity
   Multipole order     : $with_multipole_order
   No gravity below ID : $no_gravity_below_id
   Make gravity glass  : $gravity_glass_making
   External potential  : $with_potential

   Entropy floor        : $with_entropy_floor
   Cooling function     : $with_cooling
   Chemistry            : $with_chemistry
   Tracers              : $with_tracers
   Stellar model        : $with_stars
   Star formation model : $with_star_formation
   Feedback model       : $with_feedback

   Individual timers           : $enable_timers
   Task debugging              : $enable_task_debugging
   Threadpool debugging        : $enable_threadpool_debugging
   Debugging checks            : $enable_debugging_checks
   Interaction debugging       : $enable_debug_interactions
   Stars interaction debugging : $enable_debug_interactions_stars
   Naive interactions          : $enable_naive_interactions
   Gravity checks              : $gravity_force_checks
   Custom icbrtf               : $enable_custom_icbrtf

 ------------------------])
