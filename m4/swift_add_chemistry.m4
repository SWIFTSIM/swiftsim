# SYNOPSIS
#
#   SWIFT_ADD_CHEMISTRY
#
# DESCRIPTION
#
# This macro defines the chemistry model to use.
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

AC_DEFUN([SWIFT_ADD_CHEMISTRY],
[
#  chemistry function
AC_ARG_WITH([chemistry],
   [AS_HELP_STRING([--with-chemistry=<function>],
      [chemistry function @<:@none, GEAR, EAGLE default: none@:>@]
   )],
   [with_chemistry="$withval"],
   [with_chemistry="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_chemistry" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-chemistry together])
   else
      with_chemistry="$with_subgrid_chemistry"
   fi
fi

case "$with_chemistry" in
   none)
      AC_DEFINE([CHEMISTRY_NONE], [1], [No chemistry function])
   ;;
   GEAR)
      AC_DEFINE([CHEMISTRY_GEAR], [1], [Chemistry taken from the GEAR model])
   ;;
   EAGLE)
      AC_DEFINE([CHEMISTRY_EAGLE], [1], [Chemistry taken from the EAGLE model])
   ;;
   *)
      AC_MSG_ERROR([Unknown chemistry function: $with_chemistry])
   ;;
esac
])
