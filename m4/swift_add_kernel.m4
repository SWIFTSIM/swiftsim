# SYNOPSIS
#
#   SWIFT_ADD_KERNEL
#
# DESCRIPTION
#
# This macro defines the kernel.
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

AC_DEFUN([SWIFT_ADD_KERNEL],
[
# SPH Kernel function
AC_ARG_WITH([kernel],
   [AS_HELP_STRING([--with-kernel=<kernel>],
      [Kernel function to use @<:@cubic-spline, quartic-spline, quintic-spline, wendland-C2, wendland-C4, wendland-C6 default: cubic-spline@:>@]
   )],
   [with_kernel="$withval"],
   [with_kernel="cubic-spline"]
)
case "$with_kernel" in
   cubic-spline)
      AC_DEFINE([CUBIC_SPLINE_KERNEL], [1], [Cubic spline kernel])
   ;;
   quartic-spline)
      AC_DEFINE([QUARTIC_SPLINE_KERNEL], [1], [Quartic spline kernel])
   ;;
   quintic-spline)
      AC_DEFINE([QUINTIC_SPLINE_KERNEL], [1], [Quintic spline kernel])
   ;;
   wendland-C2)
      AC_DEFINE([WENDLAND_C2_KERNEL], [1], [Wendland-C2 kernel])
   ;;
   wendland-C4)
      AC_DEFINE([WENDLAND_C4_KERNEL], [1], [Wendland-C4 kernel])
   ;;
   wendland-C6)
      AC_DEFINE([WENDLAND_C6_KERNEL], [1], [Wendland-C6 kernel])
   ;;
   *)
      AC_MSG_ERROR([Unknown kernel function: $with_kernel])
   ;;
esac

])
