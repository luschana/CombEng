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

#ifndef ENGINE_H
#define ENGINE_H

#include "definitions.h"
#include "cylinder.h"
#include "ecu.h"
#include "oil.h"
#include "environment.h"

class CombustionEngine
{
public:
	static CombustionEngine *getInst();
	/** Method: run(w, thrPos)
	 * Parameters:
	 * w ... speed of rotation in [rad/s];
	 * thrPos ... throttle position, absolute value ==> [0..1]
	 */
	void run(double w, double thrPos);
	void setPhiSpark(double sparkAngle);
	void setPhiInjection(double injectionAngle);
	void setTempCoolingWater(double T_CW);
	const GasComponent& getExhaust() const;
	const GasComponent& getIntake() const;
	double getThrPos() const;
	double getW() const;
	double getMShaft() const;
	double getPhi() const;
	double getp_Cyl(int cylNum) const;
	double getT_Cyl(int cylNum) const;
	double getCyl_T(int cylNum) const;
	double getFuelConsumption() const;
	double getH_Cooling() const;
	unsigned long getCnt() const;

protected:
private:
	CombustionEngine();
	static CombustionEngine* _pInst;

	Ecu _ecu;
	Oil _oil;
	double _w;
	Cylinder _cyl[Ncyl];
	double _phiCyl[Ncyl]; // angle offset off cylinders
	GasComponent _intake;
	GasComponent _exhaust;
	double _T_CW;
	double _M_Shaft;
	unsigned long _cnt;

};

#endif // ENGINE_H
