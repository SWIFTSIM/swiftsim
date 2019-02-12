# SYNOPSIS
#
#   SWIFT_ADD_EQUATION_OF_STATE
#
# DESCRIPTION
#
# This macro defines the equation of state to use.
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

AC_DEFUN([SWIFT_ADD_EQUATION_OF_STATE],
[
#  Equation of state
AC_ARG_WITH([equation-of-state],
   [AS_HELP_STRING([--with-equation-of-state=<EoS>],
      [equation of state @<:@ideal-gas, isothermal-gas, planetary default: ideal-gas@:>@]
   )],
   [with_eos="$withval"],
   [with_eos="ideal-gas"]
)
case "$with_eos" in
   ideal-gas)
      AC_DEFINE([EOS_IDEAL_GAS], [1], [Ideal gas equation of state])
   ;;
   isothermal-gas)
      AC_DEFINE([EOS_ISOTHERMAL_GAS], [1], [Isothermal gas equation of state])
   ;;
   planetary)
      AC_DEFINE([EOS_PLANETARY], [1], [All planetary equations of state])
   ;;
   *)
      AC_MSG_ERROR([Unknown equation of state: $with_eos])
   ;;
esac

#  Adiabatic index
AC_ARG_WITH([adiabatic-index],
   [AS_HELP_STRING([--with-adiabatic-index=<gamma>],
      [adiabatic index @<:@5/3, 7/5, 4/3, 2 default: 5/3@:>@]
   )],
   [with_gamma="$withval"],
   [with_gamma="5/3"]
)
case "$with_gamma" in
   5/3)
      AC_DEFINE([HYDRO_GAMMA_5_3], [5./3.], [Adiabatic index is 5/3])
   ;;
   7/5)
      AC_DEFINE([HYDRO_GAMMA_7_5], [7./5.], [Adiabatic index is 7/5])
   ;;
   4/3)
      AC_DEFINE([HYDRO_GAMMA_4_3], [4./3.], [Adiabatic index is 4/3])
   ;;
   2)
      AC_DEFINE([HYDRO_GAMMA_2_1], [2.], [Adiabatic index is 2])
   ;;
   *)
      AC_MSG_ERROR([Unknown adiabatic index: $with_gamma])
   ;;
esac

])
