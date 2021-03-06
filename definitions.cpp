/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2014  Alex Luschan <alexander.luschan@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

#include "definitions.h"

/*
 * ΔvapH = A exp(-αTr) (1 − Tr)β
 *  ΔvapH = Enthalpy of vaporization (at saturation pressure) (J/mol)
 *  Tr = reduced temperature (T / Tc)
 */
double getFuelVaporizationEnthaly(double T, double Tmin, double Tmax, double A, double a, double b, double Tc){
	if(T < Tmin) T = Tmin;
	if(T > Tmax) T = Tmax;
	T=T/Tc;
	return A*exp(-a*T)*pow(1.0-T, b);
}

/*
 * Antoine Equation
 * returns: p [Pa];
 * params: T [K], relHum [%];
 */
double getPressureFromRelHum(double T, double relHum){
	return pow(10, 3 + WaterAntoinePars[0] - WaterAntoinePars[1]/(WaterAntoinePars[2]+T)) * relHum; // [Pa] = 10^5 [bar] * [%] / 100
}

/*
 * "project" any number to the operating range of a 4-stroke combustion engine: [0.0 .. 4*pi[
 */
double normalizeCranckAngle(double phi) {
	double result = phi;
	if(phi < 0.0) {
		result += 4*M_PI;
	}else if(phi >= 4*M_PI){
		result -= 4*M_PI;
	}
	return result;
}
