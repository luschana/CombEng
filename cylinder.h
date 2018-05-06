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

#ifndef CYLINDER_H
#define CYLINDER_H

#include "definitions.h"
#include "gascomponent.h"
#include "ecu.h"
#include "environment.h"
#include "oil.h"

class Cylinder
{
public:
	Cylinder();
    Cylinder(double phi0, GasComponent *pinlet, GasComponent *pexhaust, Ecu *pecu, Oil *poil);


    void setT_CW(double T_CW);

    void calcCylinder(double dPhi);

	double getH_Cool() const;
	double getM_G() const;
	double getM_P() const;

	double getX_p() const;
	double getdX_p() const;
	double getV_p() const;
	double getdV_p() const;

	double getT_Cyl() const;
	double getPhi() const;
	double getHxGas() const;
	double getPressure() const;
	double getTemperature() const;

private:
    double _phi, _x_p, _dx_p, _v_p, _dv_p; //cs angle, abs. pos of piston, delta pos, abs. velocity, delta in velocity[m/s];
    double _M_p, _M_g, _F_fr;
    double _H_cooling, _T_CW, _T_cyl, _H_hx_gas; // to water [J], [K], [K]
    double _n_Fuel;
    GasComponent _gc;
    GasComponent *_pinlet, *_pexhaust;
    Ecu *_pecu;
    Oil *_poil;

    double getHeatExchangeEnthalpy();
    double getCylArea(double x_pos) const;
    double getCmpFactor() const;
    bool passedAngle(double alpha, double dphi) const;
    double calcFuelInj();
};

#endif // CYLINDER_H
