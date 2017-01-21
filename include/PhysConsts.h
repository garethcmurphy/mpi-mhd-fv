/***************************************************************************
 *   Copyright (C) 2017 by Gareth Murphy   *
 *   gmurphy@cp.dias.ie *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/** 
\brief Physical Constant Database
 */
#ifndef _PHYSCONST_
#define _PHYSCONST_
namespace PhysConsts
{
/** The adiabatic index or ratio of specific heats */
  const double gamma = 5.0 / 3.0;
/** Speed of light in a vacuum*/
  const double c = 3e10;
/** Boltzmann's Constant */
  const double kb = 1.380658e-16;
/** Mass of a proton in grams */
  const double mp = 1.6726231e-24;
/** Planck's constant */
  const double h = 1.67e-27;
/** Pi */
  const double pi = 3.1415926535897932384626433832795028841971693993751058;
/** Universal Gravitational Constant */
  const double G = 6.67259e-8;
/** Mass of a Star */
  const double Mstar = 4e33;
/** Mass of a Sun */
  const double Msun = 1.98892e33;
/** Astronomical Unit */
  const double AU = 1.49598e13;
}
#endif
