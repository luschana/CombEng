/*
 * ecu.h
 *
 *  Created on: 23.11.2014
 *      Author: alex
 */

#ifndef ECU_H_
#define ECU_H_

#include "definitions.h"
#include "valve.h"

class Ecu{

public:
	Ecu();
	double getPhiValveInOpen();
	double getPhiValveInClose();
	double getPhiValveOutOpen();
	double getPhiValveOutClose();

	void setPhiValveInClose(double phiValveInClose);
	void setPhiValveInOpen(double phiValveInOpen);
	void setPhiValveOutClose(double phiValveOutClose);
	void setPhiValveOutOpen(double phiValveOutOpen);

	double getPhiSpark() const;
	void setPhiSpark(double sparkAngle);
	double getPhiInjection() const;
	void setPhiInjection(double phiInjection);
	double getThrottlePosition() const;
	void setThrottlePosition(double throttlePosition, double w);
	double getFuelConsumed() const; // get mols of fuel that were consumed
	double fillInjector(double p, double T); // consumes the fuel

private:
	double _phiSpark;
    double _phiValveInOpen, _phiValveInClose, _phiValveOutOpen, _phiValveOutClose;
    double _phiInjection;
    double _throttlePosition;
	double _w;
    double _n_fuel_consumed;

};



#endif /* ECU_H_ */
