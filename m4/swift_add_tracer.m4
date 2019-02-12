# SYNOPSIS
#
#   SWIFT_ADD_TRACER
#
# DESCRIPTION
#
# This macro defines the tracer options.
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

AC_DEFUN([SWIFT_ADD_TRACER],
[
#  Particle tracers
AC_ARG_WITH([tracers],
   [AS_HELP_STRING([--with-tracers=<function>],
      [chemistry function @<:@none, EAGLE default: none@:>@]
   )],
   [with_tracers="$withval"],
   [with_tracers="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_tracers" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-tracers together])
   else
      with_tracers="$with_subgrid_tracers"
   fi
fi

case "$with_tracers" in
   none)
      AC_DEFINE([TRACERS_NONE], [1], [No tracers function])
   ;;
   EAGLE)
      AC_DEFINE([TRACERS_EAGLE], [1], [Tracers taken from the EAGLE model])
   ;;
   *)
      AC_MSG_ERROR([Unknown tracers choice: $with_tracers])
   ;;
esac

])
