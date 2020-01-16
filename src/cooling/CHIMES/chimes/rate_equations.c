/****************************************************************************
 * This file is part of CHIMES.
 * Copyright (c) 2020 Alexander Richings (alexander.j.richings@durham.ac.uk)
 *
 * CHIMES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h> 
#include <sundials/sundials_types.h>
#include <sundials/sundials_dense.h> 
#include "chimes_vars.h"
#include "chimes_proto.h"

/** 
 * @brief Checks the constraint equations. 
 * 
 * Checks the constraint equations that ensure that the 
 * sum on all species containing a particular element 
 * matches the total element abundance, and that the 
 * net charge of all ions is zero. If any of these 
 * constraints are not met to within 1 per cent, the 
 * abundances of all species involved in that constraint 
 * are re-scaled accordingly. This routine also enforces 
 * that all abundances are non-negative. 
 * 
 * @param myGasVars The #gasVariables struct. 
 * @param myGlobalVars The #globalVariables struct. 
 */ 
void check_constraint_equations(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars)
{
  ChimesFloat x;
  int i;

  for (i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    myGasVars->abundances[i] = chimes_max(myGasVars->abundances[i], 0.0); 

  /* Helium */
  if (myGasVars->element_abundances[0] > METALS_MINIMUM_THRESHOLD)
    {
      x = myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeI]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeII]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeIII]];

      if (x <= METALS_MINIMUM_THRESHOLD) 
	myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeI]] = myGasVars->element_abundances[0]; 
      else if (fabs((x - myGasVars->element_abundances[0]) / myGasVars->element_abundances[0]) > 0.01)
	{
	  for (i = 0; i < 3; i++)
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeI + i]] *= myGasVars->element_abundances[0] / x;
	}
    }
  else 
    {
      for (i = 0; i < 3; i++)
	myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeI + i]] = 0.0;
    }

  /* Nitrogen */
  if (myGlobalVars->element_included[1] == 1)
    {
      if (myGasVars->element_abundances[2] > METALS_MINIMUM_THRESHOLD)
	{
	  x = 0.0;
	  for (i = 0; i < 8; i++)
	    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_NI + i]];

	  if (x <= METALS_MINIMUM_THRESHOLD) 
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_NI]] = myGasVars->element_abundances[2]; 
	  else if (fabs((x - myGasVars->element_abundances[2]) / myGasVars->element_abundances[2]) > 0.01)
	    {
	      for (i = 0; i < 8; i++)
		myGasVars->abundances[myGlobalVars->speciesIndices[sp_NI + i]] *= myGasVars->element_abundances[2] / x;
	    }
	}
      else 
	{
	  for (i = 0; i < 8; i++)
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_NI + i]] = 0.0; 
	}
    }

  /* Neon */
  if (myGlobalVars->element_included[3] == 1)
    {
      if (myGasVars->element_abundances[4] > METALS_MINIMUM_THRESHOLD) 
	{
	  x = 0.0;
	  for (i = 0; i < 11; i++)
	    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_NeI + i]];

	  if (x <= METALS_MINIMUM_THRESHOLD) 
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_NeI]] = myGasVars->element_abundances[4]; 
	  else if (fabs((x - myGasVars->element_abundances[4]) / myGasVars->element_abundances[4]) > 0.01)
	    {
	      for (i = 0; i < 11; i++)
		myGasVars->abundances[myGlobalVars->speciesIndices[sp_NeI + i]] *= myGasVars->element_abundances[4] / x;
	    }
	}
      else 
	{
	  for (i = 0; i < 11; i++)
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_NeI + i]] = 0.0;
	}
    }

  /* Magnesium */
  if (myGlobalVars->element_included[4] == 1)
    {
      if (myGasVars->element_abundances[5] > METALS_MINIMUM_THRESHOLD) 
	{
	  x = 0.0;
	  for (i = 0; i < 13; i++)
	    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_MgI + i]];

	  if (x <= METALS_MINIMUM_THRESHOLD) 
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_MgI]] = myGasVars->element_abundances[5]; 
	  else if (fabs((x - myGasVars->element_abundances[5]) / myGasVars->element_abundances[5]) > 0.01)
	    {
	      for (i = 0; i < 13; i++)
		myGasVars->abundances[myGlobalVars->speciesIndices[sp_MgI + i]] *= myGasVars->element_abundances[5] / x;
	    }
	}
      else 
	{
	  for (i = 0; i < 13; i++)
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_MgI + i]] = 0.0;
	}
    }

  /* Silicon */
  if (myGlobalVars->element_included[5] == 1)
    {
      if (myGasVars->element_abundances[6] > METALS_MINIMUM_THRESHOLD) 
	{
	  x = 0.0;
	  for (i = 0; i < 15; i++)
	    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_SiI + i]];

	  if (x <= METALS_MINIMUM_THRESHOLD) 
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_SiI]] = myGasVars->element_abundances[6]; 
	  else if (fabs((x - myGasVars->element_abundances[6]) / myGasVars->element_abundances[6]) > 0.01)
	    {
	      for (i = 0; i < 15; i++)
		myGasVars->abundances[myGlobalVars->speciesIndices[sp_SiI + i]] *= myGasVars->element_abundances[6] / x;
	    }
	}
      else 
	{
	  for (i = 0; i < 15; i++)
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_SiI + i]] = 0.0;
	}
    }
  /* Sulphur */
  if (myGlobalVars->element_included[6] == 1)
    {
      if (myGasVars->element_abundances[7] > METALS_MINIMUM_THRESHOLD) 
	{ 
	  x = 0.0;
	  for (i = 0; i < 17; i++)
	    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_SI + i]];

	  if (x <= METALS_MINIMUM_THRESHOLD) 
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_SI]] = myGasVars->element_abundances[7]; 
	  else if (fabs((x - myGasVars->element_abundances[7]) / myGasVars->element_abundances[7]) > 0.01)
	    {
	      for (i = 0; i < 17; i++)
		myGasVars->abundances[myGlobalVars->speciesIndices[sp_SI + i]] *= myGasVars->element_abundances[7] / x;
	    }
	}
      else 
	{
	  for (i = 0; i < 17; i++)
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_SI + i]] = 0.0;
	}
    }

  /* Calcium */
  if (myGlobalVars->element_included[7] == 1)
    {
      if (myGasVars->element_abundances[8] > METALS_MINIMUM_THRESHOLD) 
	{
	  x = 0.0;
	  for (i = 0; i < 21; i++)
	    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_CaI + i]];

	  if (x <= METALS_MINIMUM_THRESHOLD) 
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_CaI]] = myGasVars->element_abundances[8]; 
	  else if (fabs((x - myGasVars->element_abundances[8]) / myGasVars->element_abundances[8]) > 0.01)
	    {
	      for (i = 0; i < 21; i++)
		myGasVars->abundances[myGlobalVars->speciesIndices[sp_CaI + i]] *= myGasVars->element_abundances[8] / x;
	    }
	}
      else 
	{
	  for (i = 0; i < 21; i++)
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_CaI + i]] = 0;
	}
    }

  /* Iron */
  if (myGlobalVars->element_included[8] == 1)
    {
      if (myGasVars->element_abundances[9] > METALS_MINIMUM_THRESHOLD) 
	{
	  x = 0.0;
	  for (i = 0; i < 27; i++)
	    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_FeI + i]];

	  if (x <= METALS_MINIMUM_THRESHOLD) 
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_FeI]] = myGasVars->element_abundances[9]; 
	  else if (fabs((x - myGasVars->element_abundances[9]) / myGasVars->element_abundances[9]) > 0.01)
	    {
	      for (i = 0; i < 27; i++)
		myGasVars->abundances[myGlobalVars->speciesIndices[sp_FeI + i]] *= myGasVars->element_abundances[9] / x;
	    }
	}
      else 
	{
	  for (i = 0; i < 27; i++)
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_FeI + i]] = 0.0; 
	}
    }

  /* Carbon */
  if (myGlobalVars->element_included[0] == 1)
    {
      if (myGasVars->element_abundances[1] > METALS_MINIMUM_THRESHOLD) 
	{
	  x = 0.0;
	  for (i = 0; i < 8; i++)   /* Includes Cm */
	    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_CI + i]];
	  x += 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_C2]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH2]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH3p]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_CHp]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH2p]];
	  if (myGlobalVars->element_included[2] == 1)
	    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_HCOp]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_COp]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_HOCp]];

	  if (x <= METALS_MINIMUM_THRESHOLD) 
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_CI]] = myGasVars->element_abundances[1]; 
	  else if (fabs((x - myGasVars->element_abundances[1]) / myGasVars->element_abundances[1]) > 0.01)
	    {
	      for (i = 0; i < 8; i++)
		myGasVars->abundances[myGlobalVars->speciesIndices[sp_CI + i]] *= myGasVars->element_abundances[1] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_C2]] *= myGasVars->element_abundances[1] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH]] *= myGasVars->element_abundances[1] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH2]] *= myGasVars->element_abundances[1] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH3p]] *= myGasVars->element_abundances[1] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_CHp]] *= myGasVars->element_abundances[1] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH2p]] *= myGasVars->element_abundances[1] / x;
	      if (myGlobalVars->element_included[2] == 1)
		{
		  myGasVars->abundances[myGlobalVars->speciesIndices[sp_HCOp]] *= myGasVars->element_abundances[1] / x;
		  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] *= myGasVars->element_abundances[1] / x;
		  myGasVars->abundances[myGlobalVars->speciesIndices[sp_COp]] *= myGasVars->element_abundances[1] / x;
		  myGasVars->abundances[myGlobalVars->speciesIndices[sp_HOCp]] *= myGasVars->element_abundances[1] / x;
		}
	    }
	}
      else 
	{
	  for (i = 0; i < 8; i++)
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_CI + i]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_C2]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH2]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH3p]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CHp]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH2p]] = 0.0; 
	  if (myGlobalVars->element_included[2] == 1)
	    {
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_HCOp]] = 0.0;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] = 0.0;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_COp]] = 0.0; 
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_HOCp]] = 0.0; 
	    }
	}
    }

  /* Oxygen */
  if (myGlobalVars->element_included[2] == 1)
    {
      if (myGasVars->element_abundances[3] > METALS_MINIMUM_THRESHOLD) 
	{
	  x = 0.0;
	  for (i = 0; i < 10; i++)   /* Includes Om */
	    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_OI + i]]; 
	  x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_OH]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_O2]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_OHp]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2Op]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_O2p]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_H3Op]];
	  if (myGlobalVars->element_included[0] == 1)
	    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_HCOp]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_COp]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_HOCp]];

	  if (x <= METALS_MINIMUM_THRESHOLD) 
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_OI]] = myGasVars->element_abundances[3]; 
	  else if (fabs((x - myGasVars->element_abundances[3]) / myGasVars->element_abundances[3]) > 0.01)
	    {
	      for (i = 0; i < 10; i++)
		myGasVars->abundances[myGlobalVars->speciesIndices[sp_OI + i]] *= myGasVars->element_abundances[3] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_OH]] *= myGasVars->element_abundances[3] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] *= myGasVars->element_abundances[3] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_O2]] *= myGasVars->element_abundances[3] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_OHp]] *= myGasVars->element_abundances[3] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2Op]] *= myGasVars->element_abundances[3] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_O2p]] *= myGasVars->element_abundances[3] / x;
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_H3Op]] *= myGasVars->element_abundances[3] / x;
	      if (myGlobalVars->element_included[0] == 1)
		{
		  myGasVars->abundances[myGlobalVars->speciesIndices[sp_HCOp]] *= myGasVars->element_abundances[3] / x;
		  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] *= myGasVars->element_abundances[3] / x;
		  myGasVars->abundances[myGlobalVars->speciesIndices[sp_COp]] *= myGasVars->element_abundances[3] / x;
		  myGasVars->abundances[myGlobalVars->speciesIndices[sp_HOCp]] *= myGasVars->element_abundances[3] / x;
		}
	    }
	}
      else 
	{ 
	  for (i = 0; i < 10; i++)
	    myGasVars->abundances[myGlobalVars->speciesIndices[sp_OI + i]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_OH]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_O2]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_OHp]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2Op]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_O2p]] = 0.0; 
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_H3Op]] = 0.0; 
	  if (myGlobalVars->element_included[0] == 1)
	    {
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_HCOp]] = 0.0; 
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] = 0.0; 
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_COp]] = 0.0; 
	      myGasVars->abundances[myGlobalVars->speciesIndices[sp_HOCp]] = 0.0; 
	    }
	}
    }
  
  /* Hydrogen */
  x = myGasVars->abundances[myGlobalVars->speciesIndices[sp_HI]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_HII]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_Hm]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2p]] + 3.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H3p]];

  if (myGlobalVars->element_included[0] == 1)
    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH2]] + 3.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH3p]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_CHp]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH2p]];
  if (myGlobalVars->element_included[2] == 1)
    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_OH]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]]+ myGasVars->abundances[myGlobalVars->speciesIndices[sp_OHp]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2Op]] + 3.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H3Op]];
  if (myGlobalVars->element_included[0] == 1 && myGlobalVars->element_included[2] == 1)
    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_HCOp]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_HOCp]];

  if (x <= METALS_MINIMUM_THRESHOLD) 
    myGasVars->abundances[myGlobalVars->speciesIndices[sp_HI]] = 1.0; 
  else if (fabs(x - 1.0) > 0.01)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_HI]] /= x;
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_HII]] /= x;
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_Hm]] /= x;
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]] /= x;
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2p]] /= x;
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_H3p]] /= x;
      if (myGlobalVars->element_included[0] == 1)
	{
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH]] /= x;
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH2]] /= x;
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH3p]] /= x;
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CHp]] /= x;
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH2p]] /= x;
	}
      if (myGlobalVars->element_included[2] == 1)
	{
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_OH]] /= x;
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] /= x;
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_OHp]] /= x;
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2Op]] /= x;
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_H3Op]] /= x;
	}
      if (myGlobalVars->element_included[0] == 1 && myGlobalVars->element_included[2] == 1)
	{
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_HCOp]] /= x;
	  myGasVars->abundances[myGlobalVars->speciesIndices[sp_HOCp]] /= x;
	}
    }

  /* Electrons */
  x = myGasVars->abundances[myGlobalVars->speciesIndices[sp_HII]];
  x -= myGasVars->abundances[myGlobalVars->speciesIndices[sp_Hm]];
  x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeII]];
  x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeIII]] * 2.0;

  if (myGlobalVars->element_included[0] == 1)
    {
      for (i = 1; i <= 6; i++)
	x += i * myGasVars->abundances[myGlobalVars->speciesIndices[sp_CI + i]];
      x -= myGasVars->abundances[myGlobalVars->speciesIndices[sp_Cm]];
      x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH3p]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_CHp]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_CH2p]];
    }

  if (myGlobalVars->element_included[1] == 1)
    {
      for (i = 1; i <= 7; i++)
	x += i * myGasVars->abundances[myGlobalVars->speciesIndices[sp_NI + i]];
    }

  if (myGlobalVars->element_included[2] == 1)
    {
      for (i = 1; i <= 8; i++)
	x += i * myGasVars->abundances[myGlobalVars->speciesIndices[sp_OI + i]];
      x -= myGasVars->abundances[myGlobalVars->speciesIndices[sp_Om]];
      x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_OHp]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2Op]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_H3Op]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_O2p]];
    }

  if (myGlobalVars->element_included[0] == 1 && myGlobalVars->element_included[2] == 1)
    x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_HCOp]] +  + myGasVars->abundances[myGlobalVars->speciesIndices[sp_COp]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_HOCp]];

  if (myGlobalVars->element_included[3] == 1)
    {
      for (i = 1; i <= 10; i++)
	x += i * myGasVars->abundances[myGlobalVars->speciesIndices[sp_NeI + i]];
    }

  if (myGlobalVars->element_included[4] == 1)
    {
      for (i = 1; i <= 12; i++)
	x += i * myGasVars->abundances[myGlobalVars->speciesIndices[sp_MgI + i]];
    }

  if (myGlobalVars->element_included[5] == 1)
    {
      for (i = 1; i <= 14; i++)
	x += i * myGasVars->abundances[myGlobalVars->speciesIndices[sp_SiI + i]];
    }

  if (myGlobalVars->element_included[6] == 1)
    {
      for (i = 1; i <= 16; i++)
	x += i * myGasVars->abundances[myGlobalVars->speciesIndices[sp_SI + i]];
    }

  if (myGlobalVars->element_included[7] == 1)
    {
      for (i = 1; i <= 20; i++)
	x += i * myGasVars->abundances[myGlobalVars->speciesIndices[sp_CaI + i]];
    }

  if (myGlobalVars->element_included[8] == 1)
    {
      for (i = 1; i <= 26; i++)
	x += i * myGasVars->abundances[myGlobalVars->speciesIndices[sp_FeI + i]];
    }

  x += myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2p]] + myGasVars->abundances[myGlobalVars->speciesIndices[sp_H3p]];

  if (fabs((x - myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]]) / chimes_max(myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]], 1.0e-100)) > 0.01)
    myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] = x;
}

/** 
 * @brief Defines the right-hand side function. 
 * 
 * Defines the system of differential equations that make 
 * up the right-hand side function, which will be integrated 
 * by CVode. 
 * 
 * @param t Current time. 
 * @param y Vector containing the variables to be integrated. 
 * @param ydot Vector containing the time derivatives of the variables. 
 * @param user_data The #UserData struct containing the input data. 
 */
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  int i, j;		
  struct UserData *data;
  
  data = (struct UserData *) user_data;
  int indices[data->network_size];  /* We will use this array to relate the enum types of                                         
				     * each (non-eq) species to their position in y */  

  /* First, loop through the enum types of all
   * non-eq species. If they are included in 
   * the network then their abundance is in
   * the vector y. */

  i = 0;	/* We use this to keep track of where we are in the vector y */
  for (j = 0; j < data->myGlobalVars->totalNumberOfSpecies; j++)
    {
      if (data->species[j].include_species == 1)
	{
	  data->myGasVars->abundances[j] = (ChimesFloat) NV_Ith_S(y, i);
	  indices[i] = j;
	  i++;
	}
    }
	
  /* If Thermal Evolution is switched on, the final element in the
   * vector y is the internal energy (per unit volume). Use this 
   * to update the temperature, and also the rates that depend on T */
  if (data->myGasVars->ThermEvolOn == 1)
    data->myGasVars->temperature = chimes_max(((ChimesFloat) NV_Ith_S(y, data->network_size)) / (1.5 * calculate_total_number_density(data->myGasVars->abundances, data->myGasVars->nH_tot, data->myGlobalVars) * BOLTZMANNCGS), 10.1); /* The rates are not defined below ~10 K */

  // Update rates 
  update_rate_coefficients(data->myGasVars, data->myGlobalVars, *data, data->myGasVars->ThermEvolOn); 
  update_rates(data->myGasVars, data->myGlobalVars, *data); 

  // Zero all species rates 
  for (i = 0; i < data->network_size; i++)
    {
      data->species[indices[i]].creation_rate = 0.0;
      data->species[indices[i]].destruction_rate = 0.0;
    }

  // Compute creation and destruction rates 
  update_rate_vector(data->species, data->myGasVars, data->myGlobalVars, *data); 
  
  
  /* Now set the output ydot vector for the chemical abundances */
  for (i = 0; i < data->network_size; i++)
    NV_Ith_S(ydot, i) = (realtype) (data->species[indices[i]].creation_rate - data->species[indices[i]].destruction_rate);
		
  // Finally, if Thermal Evolution is switched on, calculate the cooling rate 
  if (data->myGasVars->ThermEvolOn == 1)
    {
      if (data->myGasVars->temperature > data->myGasVars->TempFloor)
	NV_Ith_S(ydot, data->network_size) = (realtype) -calculate_total_cooling_rate(data->myGasVars, data->myGlobalVars, *data, 0);		/* Note that network_size is the number of chemcial species, hence No. of eqns = network_size + 1 when ThermEvol is on */
      else
  	NV_Ith_S(ydot, data->network_size) = (realtype) chimes_max(-calculate_total_cooling_rate(data->myGasVars, data->myGlobalVars, *data, 0), 0.0);  /* Once T falls below T_floor, set T_dot >= 0 */ 
    }

  return 0;
}
