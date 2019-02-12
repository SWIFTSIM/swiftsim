# SYNOPSIS
#
#   SWIFT_ADD_STARS
#
# DESCRIPTION
#
# This macro defines the stars options.
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

AC_DEFUN([SWIFT_ADD_STARS],
[
# Check if debugging interactions stars is switched on.
AC_ARG_ENABLE([debug-interactions-stars],
   [AS_HELP_STRING([--enable-debug-interactions-stars],
     [Activate interaction debugging for stars, logging a maximum of @<:@N@:>@ neighbours. Defaults to 256 if no value set.]
   )],
   [enable_debug_interactions_stars="$enableval"],
   [enable_debug_interactions_stars="no"]
)
if test "$enable_debug_interactions_stars" != "no"; then
    AC_DEFINE([DEBUG_INTERACTIONS_STARS],1,[Enable interaction debugging for stars])
    if test "$enable_debug_interactions_stars" == "yes"; then
      AC_DEFINE([MAX_NUM_OF_NEIGHBOURS_STARS],256,[The maximum number of particle neighbours to be logged for stars])
      [enable_debug_interactions_stars="yes (Logging up to 256 neighbours)"]
    else
      AC_DEFINE_UNQUOTED([MAX_NUM_OF_NEIGHBOURS_STARS], [$enableval] ,[The maximum number of particle neighbours to be logged for stars])
      [enable_debug_interactions_stars="yes (Logging up to $enableval neighbours)"]
    fi
fi

# Stellar model.
AC_ARG_WITH([stars],
   [AS_HELP_STRING([--with-stars=<model>],
      [Stellar model to use @<:@none, EAGLE, GEAR, debug default: none@:>@]
   )],
   [with_stars="$withval"],
   [with_stars="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_stars" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-stars together])
   else
      with_stars="$with_subgrid_stars"
   fi
fi

case "$with_stars" in
   EAGLE)
      AC_DEFINE([STARS_EAGLE], [1], [EAGLE stellar model])
   ;;
   GEAR)
      AC_DEFINE([STARS_GEAR], [1], [GEAR stellar model])
   ;;
   none)
      AC_DEFINE([STARS_NONE], [1], [None stellar model])
   ;;

   *)
      AC_MSG_ERROR([Unknown stellar model: $with_stars])
   ;;
esac

# Feedback model
AC_ARG_WITH([feedback],
   [AS_HELP_STRING([--with-feedback=<model>],
      [Feedback model to use @<:@none, thermal, debug default: none@:>@]
   )],
   [with_feedback="$withval"],
   [with_feedback="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_feedback" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-feedback together])
   else
      with_feedback="$with_subgrid_feedback"
   fi
fi

case "$with_feedback" in
   thermal)
      AC_DEFINE([FEEDBACK_THERMAL], [1], [Thermal Blastwave])
   ;;
   none)
   ;;

   *)
      AC_MSG_ERROR([Unknown feedback model: $with_feedback])
   ;;
esac


#  Star formation
AC_ARG_WITH([star-formation], 
    [AS_HELP_STRING([--with-star-formation=<sfm>],
       [star formation @<:@none, EAGLE, default: none@:>@] 
    )],
    [with_star_formation="$withval"],
    [with_star_formation="none"]
)
if test "$with_subgrid" != "none"; then
   if test "$with_star_formation" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-star-formation together])
   else
      with_star_formation="$with_subgrid_star_formation"
   fi
fi

case "$with_star_formation" in
   none)
      AC_DEFINE([STAR_FORMATION_NONE], [1], [No star formation])
   ;;
   EAGLE)
      AC_DEFINE([STAR_FORMATION_EAGLE], [1], [EAGLE star formation model (Schaye and Dalla Vecchia (2008))])
   ;;
   *)
      AC_MSG_ERROR([Unknown star formation model])
   ;;
esac 

])
