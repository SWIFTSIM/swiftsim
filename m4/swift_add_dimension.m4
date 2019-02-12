# SYNOPSIS
#
#   SWIFT_ADD_DIMENSION
#
# DESCRIPTION
#
# This macro defines the number of dimension to use.
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


AC_DEFUN([SWIFT_ADD_DIMENSION],
[
#  Dimensionality of the hydro scheme.
AC_ARG_WITH([hydro-dimension],
   [AS_HELP_STRING([--with-hydro-dimension=<dim>],
      [dimensionality of problem @<:@3/2/1 default: 3@:>@]
   )],
   [with_dimension="$withval"],
   [with_dimension="3"]
)
case "$with_dimension" in
   1)
      AC_DEFINE([HYDRO_DIMENSION_1D], [1], [1D solver])
   ;;
   2)
      AC_DEFINE([HYDRO_DIMENSION_2D], [2], [2D solver])
   ;;
   3)
      AC_DEFINE([HYDRO_DIMENSION_3D], [3], [3D solver])
   ;;
   *)
      AC_MSG_ERROR([Dimensionality must be 1, 2 or 3])
   ;;
esac

])
