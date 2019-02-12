# SYNOPSIS
#
#   SWIFT_ADD_HYDRO
#
# DESCRIPTION
#
# This macro defines the hydro options.
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

AC_DEFUN([SWIFT_ADD_HYDRO],
[
# Hydro scheme.
AC_ARG_WITH([hydro],
   [AS_HELP_STRING([--with-hydro=<scheme>],
      [Hydro dynamics to use @<:@gadget2, minimal, pressure-entropy, pressure-energy, pressure-energy-monaghan, default, gizmo-mfv, gizmo-mfm, shadowfax, planetary, anarchy-pu debug default: gadget2@:>@]
   )],
   [with_hydro="$withval"],
   [with_hydro="gadget2"]
)

case "$with_hydro" in
   gadget2)
      AC_DEFINE([GADGET2_SPH], [1], [Gadget-2 SPH])
   ;;
   minimal)
      AC_DEFINE([MINIMAL_SPH], [1], [Minimal SPH])
   ;;
   pressure-entropy)
      AC_DEFINE([HOPKINS_PE_SPH], [1], [Pressure-Entropy SPH])
   ;;
   pressure-energy)
      AC_DEFINE([HOPKINS_PU_SPH], [1], [Pressure-Energy SPH])
   ;;
   pressure-energy-monaghan)
      AC_DEFINE([HOPKINS_PU_SPH_MONAGHAN], [1], [Pressure-Energy SPH with M&M Variable A.V.])
   ;;
   default)
      AC_DEFINE([DEFAULT_SPH], [1], [Default SPH])
   ;;
   gizmo-mfv)
      AC_DEFINE([GIZMO_MFV_SPH], [1], [GIZMO MFV SPH])
   ;;
   gizmo-mfm)
      AC_DEFINE([GIZMO_MFM_SPH], [1], [GIZMO MFM SPH])
   ;;
   shadowfax)
      AC_DEFINE([SHADOWFAX_SPH], [1], [Shadowfax SPH])
   ;;
   planetary)
      AC_DEFINE([PLANETARY_SPH], [1], [Planetary SPH])
   ;;
   anarchy-pu)
      AC_DEFINE([ANARCHY_PU_SPH], [1], [ANARCHY (PU) SPH])
   ;;


   *)
      AC_MSG_ERROR([Unknown hydrodynamics scheme: $with_hydro])
   ;;
esac

# Check if debugging interactions is switched on.
AC_ARG_ENABLE([debug-interactions],
   [AS_HELP_STRING([--enable-debug-interactions],
     [Activate interaction debugging, logging a maximum of @<:@N@:>@ neighbours. Defaults to 256 if no value set.]
   )],
   [enable_debug_interactions="$enableval"],
   [enable_debug_interactions="no"]
)
if test "$enable_debug_interactions" != "no"; then
  if test "$with_hydro" = "gadget2"; then
      AC_DEFINE([DEBUG_INTERACTIONS_SPH],1,[Enable interaction debugging])
    if test "$enable_debug_interactions" == "yes"; then
      AC_DEFINE([MAX_NUM_OF_NEIGHBOURS],256,[The maximum number of particle neighbours to be logged])
      [enable_debug_interactions="yes (Logging up to 256 neighbours)"]
    else
      AC_DEFINE_UNQUOTED([MAX_NUM_OF_NEIGHBOURS], [$enableval] ,[The maximum number of particle neighbours to be logged])
      [enable_debug_interactions="yes (Logging up to $enableval neighbours)"]
    fi
  else
    [enable_debug_interactions="no (only available for gadget2 hydro scheme)"]
  fi
fi


#  Entropy floor
AC_ARG_WITH([entropy-floor], 
    [AS_HELP_STRING([--with-entropy-floor=<floor>],
       [entropy floor @<:@none, EAGLE, default: none@:>@] 
    )],
    [with_entropy_floor="$withval"],
    [with_entropy_floor="none"]
)
if test "$with_subgrid" != "none"; then
   if test "$with_entropy_floor" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-entropy-floor together])
   else
      with_entropy_floor="$with_subgrid_entropy_floor"
   fi
fi

case "$with_entropy_floor" in
   none)
      AC_DEFINE([ENTROPY_FLOOR_NONE], [1], [No entropy floor])
   ;;
   EAGLE)
      AC_DEFINE([ENTROPY_FLOOR_EAGLE], [1], [EAGLE entropy floor])
   ;;
   *)
      AC_MSG_ERROR([Unknown entropy floor model])
   ;;
esac 

])
