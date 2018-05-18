/*
 * valve.h
 *
 *  Created on: 31.01.2016
 *      Author: alex
 */

#ifndef VALVE_H_
#define VALVE_H_

#include "definitions.h"

class Valve{

public:
	Valve();
	Valve(double A_Valves, int num_Valves, double dphi_open, double dphi_close);
	double getPhiM();
	double getDPhi();
	void setPhiM(double phi_M);
	void setDPhi(double d_Phi);
	double getCrosssection(double phi);

private:
	double _phi_m;		// angle of middle position of camshaft (in means of crankshaft working cycle)
	double _dphi;		// angle open<->middle and middle<->close position
	double _r;			// radius of a single valve
	double _stroke;		// stroke of valve(s)
	int _num;			// number of valves


	double getActStroke(double phi);
};



#endif /* VALVE_H_ */
