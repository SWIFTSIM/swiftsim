# SYNOPSIS
#
#   SWIFT_ADD_COOLING
#
# DESCRIPTION
#
# This macro defines the cooling model to use.
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

AC_DEFUN([SWIFT_ADD_COOLING],
[
#  Cooling function
AC_ARG_WITH([cooling],
   [AS_HELP_STRING([--with-cooling=<function>],
      [cooling function @<:@none, const-du, const-lambda, EAGLE, grackle, grackle1, grackle2, grackle3 default: none@:>@]
   )],
   [with_cooling="$withval"],
   [with_cooling="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_cooling" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-cooling together])
   else
      with_cooling="$with_subgrid_cooling"
   fi
fi

case "$with_cooling" in
   none)
      AC_DEFINE([COOLING_NONE], [1], [No cooling function])
   ;;
   const-du)
      AC_DEFINE([COOLING_CONST_DU], [1], [Const du/dt cooling function])
   ;;
   const-lambda)
      AC_DEFINE([COOLING_CONST_LAMBDA], [1], [Const Lambda cooling function])
   ;;
   compton)
      AC_DEFINE([COOLING_COMPTON], [1], [Compton cooling off the CMB])
   ;;
   grackle)
      AC_DEFINE([COOLING_GRACKLE], [1], [Cooling via the grackle library])
      AC_DEFINE([COOLING_GRACKLE_MODE], [0], [Grackle chemistry network, mode 0])
   ;;
   grackle1)
      AC_DEFINE([COOLING_GRACKLE], [1], [Cooling via the grackle library])
      AC_DEFINE([COOLING_GRACKLE_MODE], [1], [Grackle chemistry network, mode 1])
   ;;
   grackle2)
      AC_DEFINE([COOLING_GRACKLE], [1], [Cooling via the grackle library])
      AC_DEFINE([COOLING_GRACKLE_MODE], [2], [Grackle chemistry network, mode 2])
   ;;
   grackle3)
      AC_DEFINE([COOLING_GRACKLE], [1], [Cooling via the grackle library])
      AC_DEFINE([COOLING_GRACKLE_MODE], [3], [Grackle chemistry network, mode 3])
   ;;
   EAGLE)
      AC_DEFINE([COOLING_EAGLE], [1], [Cooling following the EAGLE model])
   ;;
   *)
      AC_MSG_ERROR([Unknown cooling function: $with_cooling])
   ;;
esac

])
