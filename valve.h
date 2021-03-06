/*
 * valve.h
 *
 *  Created on: 31.01.2016
 *      Author: alex
 */

#ifndef VALVE_H_
#define VALVE_H_

#include "definitions.h"
#include "gascomponent.h"

class Valve{

public:
	Valve();
	Valve(double A_Valves, int num_Valves, double dphi_open, double dphi_close);
	Valve(double A_Valves);
	double getPhiM();
	double getDPhi();
	void setPhiM(double phi_M);
	void setDPhi(double d_Phi);
	void calcFlow(double phi, GasComponent *pIn, GasComponent *pOut);
	const GasComponent *getGasComponent();

private:
	double _phi_m;		// angle of middle position of camshaft (in means of crankshaft working cycle)
	double _dphi;		// angle open<->middle and middle<->close position
	double _r;			// radius of a single valve
	double _stroke;		// stroke of valve(s)
	int _num;			// number of valves
	GasComponent *_pgc;	// GasComponent to hold the floating trough gas
	bool _isAperture;	// specify this valve to have a const A_cr
	double _A_Aperture;		// value of cross section area

	double getCrosssection(double phi);
	double getActStroke(double phi);
};



#endif /* VALVE_H_ */
