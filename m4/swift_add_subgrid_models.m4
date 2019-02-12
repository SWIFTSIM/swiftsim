# SYNOPSIS
#
#   SWIFT_ADD_SUBGRID_MODELS
#
# DESCRIPTION
#
# This macro defines the subgrid models to use.
# It overwrites the default value of different models and
# therefore should be called before any other specific
# subgrid models.
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

AC_DEFUN([SWIFT_ADD_SUBGRID_MODELS],
[
# Various package configuration options.

# Master subgrid options
# If you add a restriction (e.g. no cooling, chemistry or hydro)
# you will need to check for overwrite after reading the additional options.
# As an example for this, see the call to AC_ARG_WITH for cooling.
AC_ARG_WITH([subgrid],
	[AS_HELP_STRING([--with-subgrid=<subgrid>],
		[Master switch for subgrid methods. Inexperienced user should start from here @<:@none, GEAR, EAGLE default: none@:>@]
	)],
	[with_subgrid="$withval"],
	[with_subgrid=none]
)

# Default values
with_subgrid_cooling=none
with_subgrid_chemistry=none
with_subgrid_tracers=none
with_subgrid_entropy_floor=none
with_subgrid_stars=none
with_subgrid_star_formation=none
with_subgrid_feedback=none

case "$with_subgrid" in
   yes)
      AC_MSG_ERROR([Invalid option. A subgrid model must be chosen.])
   ;;
   none)
   ;;
   GEAR)
	with_subgrid_cooling=grackle
	with_subgrid_chemistry=GEAR
	with_subgrid_tracers=none
	with_subgrid_entropy_floor=none
	with_subgrid_stars=GEAR
	with_subgrid_star_formation=none
	with_subgrid_feedback=thermal
   ;;
   EAGLE)
	with_subgrid_cooling=EAGLE
	with_subgrid_chemistry=EAGLE
	with_subgrid_tracers=EAGLE
	with_subgrid_entropy_floor=EAGLE
	with_subgrid_stars=EAGLE
	with_subgrid_star_formation=EAGLE
	with_subgrid_feedback=none
   ;;
   *)
      AC_MSG_ERROR([Unknown subgrid choice: $with_subgrid])
   ;;
esac
])
