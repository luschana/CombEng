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
	double getValveOut_A(double phi);
	double getValveIn_A(double phi);
//	double getPhiValveInClose() const;
	void setPhiValveInClose(double phiValveInClose);
//	double getPhiValveInOpen() const;
	void setPhiValveInOpen(double phiValveInOpen);
//	double getPhiValveOutClose() const;
	void setPhiValveOutClose(double phiValveOutClose);
//	double getPhiValveOutOpen() const;
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
    Valve _valve_exh, _valve_inl;
    double _n_fuel_consumed;

};



#endif /* ECU_H_ */
