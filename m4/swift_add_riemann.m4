# SYNOPSIS
#
#   SWIFT_ADD_RIEMANN
#
# DESCRIPTION
#
# This macro defines the riemann solver.
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

AC_DEFUN([SWIFT_ADD_RIEMANN],
[
#  Riemann solver
AC_ARG_WITH([riemann-solver],
   [AS_HELP_STRING([--with-riemann-solver=<solver>],
      [riemann solver (gizmo-sph only) @<:@none, exact, trrs, hllc, default: none@:>@]
   )],
   [with_riemann="$withval"],
   [with_riemann="none"]
)
case "$with_riemann" in
   none)
      AC_DEFINE([RIEMANN_SOLVER_NONE], [1], [No Riemann solver])
   ;;
   exact)
      AC_DEFINE([RIEMANN_SOLVER_EXACT], [1], [Exact Riemann solver])
   ;;
   trrs)
      AC_DEFINE([RIEMANN_SOLVER_TRRS], [1], [Two Rarefaction Riemann Solver])
   ;;
   hllc)
      AC_DEFINE([RIEMANN_SOLVER_HLLC], [1], [Harten-Lax-van Leer-Contact Riemann solver])
   ;;
   *)
      AC_MSG_ERROR([Unknown Riemann solver: $with_riemann])
   ;;
esac

])
