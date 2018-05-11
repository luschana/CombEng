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

#ifndef GASCOMPONENT_H
#define GASCOMPONENT_H

#include "definitions.h"

double getPressureFromRelHum(double T, double relHum);

class GasComponent
{
public:
	GasComponent();
	GasComponent(double V, double T, double p);
	GasComponent(double V, double T, double p, const double nu[defs::Fuel+1]);
	GasComponent(double V, double T, double p, const double nu[defs::Fuel+1], bool isContainer);
	/*
	 * remove mols from higher pressured component and mix it to the other
	 */
	void calcGasExchange(double A_crosssection, GasComponent *pgc);
	void calcStateChange(double cmpFactor, double H_cooling, double n_Fuel);//, const GasComponent &inlet, GasComponent &exhaust);
//setter methods
	void setCombustionStarted(bool combustionStarted);
	//void setMols(double n);
//getter methods
	double getSpecHeatCapacity() const;
	double getEnthalpy() const;
	double getMolareWeight() const;
	double getMols() const;
	const double* getNu() const;
	double getP() const;
	double getT() const;
	double getspecV() const;
	double getV() const;
	bool isCombustionStarted() const;

protected:
	double _n_g, _nu[defs::Fuel+1];
	double _T,_p,_v;
	double _V,_H;
	double _MW,_cp;
	bool _combustionStarted;
	bool _isContainer;

	void transferFrom(double n, GasComponent &gc);
	void normalizeMols();


private:
	double _n_chemR[4];
	double isentropicStateChange(double cmpFactor);
	double injection(double n_Fuel);
	double chemReaction();
	double calcMolareWeight();
	double calcMols();
};


#endif // GASCOMPONENT_H
