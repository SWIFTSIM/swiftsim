# SYNOPSIS
#
#   SWIFT_ADD_GRAVITY
#
# DESCRIPTION
#
# This macro defines the gravity options.
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

AC_DEFUN([SWIFT_ADD_GRAVITY],
[

# Gravity scheme.
AC_ARG_WITH([gravity],
   [AS_HELP_STRING([--with-gravity=<scheme>],
      [Gravity scheme to use @<:@default, with-potential, default: default@:>@]
   )],
   [with_gravity="$withval"],
   [with_gravity="default"]
)

case "$with_gravity" in
   with-potential)
      AC_DEFINE([POTENTIAL_GRAVITY], [1], [Gravity scheme with potential calculation])
   ;;
   default)
      AC_DEFINE([DEFAULT_GRAVITY], [1], [Default gravity scheme])
   ;;
   *)
      AC_MSG_ERROR([Unknown gravity scheme: $with_gravity])
   ;;
esac

#  External potential
AC_ARG_WITH([ext-potential],
   [AS_HELP_STRING([--with-ext-potential=<pot>],
      [external potential @<:@none, point-mass, point-mass-ring, point-mass-softened, isothermal, softened-isothermal, nfw, hernquist, disc-patch, sine-wave, default: none@:>@]
   )],
   [with_potential="$withval"],
   [with_potential="none"]
)

case "$with_potential" in
   none)
      AC_DEFINE([EXTERNAL_POTENTIAL_NONE], [1], [No external potential])
   ;;
   point-mass)
      AC_DEFINE([EXTERNAL_POTENTIAL_POINTMASS], [1], [Point-mass external potential])
   ;;
   isothermal)
      AC_DEFINE([EXTERNAL_POTENTIAL_ISOTHERMAL], [1], [Isothermal external potential])
   ;;
   hernquist)
      AC_DEFINE([EXTERNAL_POTENTIAL_HERNQUIST], [1], [Hernquist external potential])
   ;;
   nfw)
      AC_DEFINE([EXTERNAL_POTENTIAL_NFW], [1], [Navarro-Frenk-White external potential])
   ;;
   disc-patch)
      AC_DEFINE([EXTERNAL_POTENTIAL_DISC_PATCH], [1], [Disc-patch external potential])
   ;;
   sine-wave)
      AC_DEFINE([EXTERNAL_POTENTIAL_SINE_WAVE], [1], [Sine wave external potential in 1D])
   ;;
   point-mass-ring)
      AC_DEFINE([EXTERNAL_POTENTIAL_POINTMASS_RING], [1], [Point mass potential for Keplerian Ring (Hopkins 2015).])
   ;;
   point-mass-softened)
      AC_DEFINE([EXTERNAL_POTENTIAL_POINTMASS_SOFT], [1], [Softened point-mass potential with form 1/(r^2 + softening^2).])
   ;;
   *)
      AC_MSG_ERROR([Unknown external potential: $with_potential])
   ;;
esac

#  Gravity multipole order
AC_ARG_WITH([multipole-order],
   [AS_HELP_STRING([--with-multipole-order=<order>],
      [order of the multipole and gravitational field expansion @<:@ default: 4@:>@]
   )],
   [with_multipole_order="$withval"],
   [with_multipole_order="4"]
)
AC_DEFINE_UNQUOTED([SELF_GRAVITY_MULTIPOLE_ORDER], [$with_multipole_order], [Multipole order])


# Check if we want to zero the gravity forces for all particles below some ID.
AC_ARG_ENABLE([no-gravity-below-id],
   [AS_HELP_STRING([--enable-no-gravity-below-id],
     [Zeros the gravitational acceleration of all particles with an ID smaller than @<:@N@:>@]
   )],
   [no_gravity_below_id="$enableval"],
   [no_gravity_below_id="no"]
)
if test "$no_gravity_below_id" == "yes"; then
   AC_MSG_ERROR(Need to specify the ID below which particles get zero forces when using --enable-no-gravity-below-id!)
elif test "$no_gravity_below_id" != "no"; then
   AC_DEFINE_UNQUOTED([SWIFT_NO_GRAVITY_BELOW_ID], [$enableval] ,[Particles with smaller ID than this will have zero gravity forces])
fi


])
