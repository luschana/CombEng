/*
 * valve.cpp
 *
 *  Created on: 31.01.2016
 *      Author: alex
 */

#include "valve.h"

/**
 * constructor for a valve
 */
Valve::Valve(double A_Valves){
	_phi_m = 0.0;
	_dphi = 0.0;
	_num = 1;
	_r = 0.0;
	_stroke = 0.0;
	_pgc = new GasComponent();
	_isAperture = true;
	_A_Aperture = A_Valves;
}

/**
 * constructor for a valve
 */
Valve::Valve(double A_Valves, int num_Valves, double phi_open, double phi_close){
	_phi_m = (phi_open+phi_close)/2;
	_dphi = (phi_close-phi_open)/2;
	_num = num_Valves;
	_r = sqrt(A_Valves/((double)num_Valves * M_PI));
	_stroke = _r*2;
	_pgc = new GasComponent();
	_isAperture = false;
	_A_Aperture = 0.0;
}

double Valve::getPhiM() {
	return _phi_m;
}

double Valve::getDPhi() {
	return _dphi;
}

void Valve::setPhiM(double phi_M) {
	_phi_m = phi_M;
}

void Valve::setDPhi(double d_Phi) {
	_dphi = d_Phi;
}

/**
 * init the calculation of the transfer component -- do this for closed valve too...
 */
void Valve::calcFlow(double phi, GasComponent* pIn, GasComponent* pOut) {
	_pgc->calcFlow(getCrosssection(phi), pIn, pOut);
}

const GasComponent* Valve::getGasComponent() {
	return _pgc;
}

double Valve::getActStroke(double phi){
	double result = 0.0;
	if (!( _phi_m + _dphi > 4*M_PI && phi < M_PI)) { // closing might happen in between 0.. (<M_PI)
		result = (1 - fabs(_phi_m - phi)/_dphi)*_stroke;
	}else{
		result = (1 -fabs(_phi_m - phi - 4*M_PI)/_dphi) * _stroke;
	}
	return result;
}

double Valve::getCrosssection(double phi) {
	double result = _A_Aperture;
	if(!_isAperture){
		result = getActStroke(phi);
		if(result > 0.0){
			if(result > _r/2.0){ 		// fully opened
				result = _r*_r*M_PI*_num;
			}else{						// partially opened
				result *= 2.0*_r*M_PI*_num;
			}
		}
	}
	return result;
}

Valve::Valve(){
	_phi_m = 0.0;
	_dphi = 0.0;
	_num = 1;
	_r = 0.0;
	_stroke = 0.0;
	_pgc = new GasComponent();
	_isAperture = true;
	_A_Aperture = 0.0;
}
