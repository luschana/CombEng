/*
 * ecu.cpp
 *
 *  Created on: 23.11.2014
 *      Author: alex
 */

#include "ecu.h"

Ecu::Ecu() {
	_phiSpark = (1.0 - 10.0/180.0)*M_PI;
	_phiInjection = (1.0 - 10.0/180.0)*M_PI;
	_phiValveOutOpen = (2.0	-	10.0/180.0)*M_PI;
	_phiValveOutClose = (3.0+	10.0/180.0)*M_PI;
	_phiValveInOpen = (3.0	-	10.0/180.0)*M_PI;
	_phiValveInClose = (4.0	+	10.0/180.0)*M_PI;
	_valve_exh = Valve(A_Valve_out, num_Valve, 	_phiValveOutOpen, _phiValveOutClose);
	_valve_inl = Valve(A_Valve_in, num_Valve, 	_phiValveInOpen, _phiValveInClose);
	_throttlePosition = 0.0;
	_w = 0.0;
	_n_fuel_consumed = 0.0;
}

double Ecu::getValveOut_A(double phi){
	return _valve_exh.getCrosssection(phi);
}

double Ecu::getValveIn_A(double phi){
	return _valve_inl.getCrosssection(phi);
}

void Ecu::setPhiValveInClose(double phiValveInClose) {
	_phiValveInClose = phiValveInClose;
}

void Ecu::setPhiValveInOpen(double phiValveInOpen) {
	_phiValveInOpen = phiValveInOpen;
}

void Ecu::setPhiValveOutClose(double phiValveOutClose) {
	_phiValveOutClose = phiValveOutClose;
}

void Ecu::setPhiValveOutOpen(double phiValveOutOpen) {
	_phiValveOutOpen = phiValveOutOpen;
}


/**
 * valve class should do this!!!

double Ecu::getPhiValveInClose() const {
	return _phiValveInClose;
}

double Ecu::getPhiValveInOpen() const {
	return _phiValveInOpen;
}

double Ecu::getPhiValveOutClose() const {
	return _phiValveOutClose;
}

double Ecu::getPhiValveOutOpen() const {
	return _phiValveOutOpen;
}*/


double Ecu::getPhiSpark() const {
	return _phiSpark;
}

void Ecu::setPhiSpark(double sparkAngle) {
	if(sparkAngle >= -40.0 && sparkAngle <= 20.0){
		_phiSpark = (1.0 + sparkAngle / 180.0)*M_PI;
	}
}

double Ecu::getPhiInjection() const {
	return _phiInjection;
}

void Ecu::setPhiInjection(double phiInjection) {
	if(phiInjection >= -80.0 && phiInjection <= 20.0){
		_phiInjection = (1.0 + phiInjection / 180.0)*M_PI;
	}
}

double Ecu::getThrottlePosition() const {
	return _throttlePosition;
}

void Ecu::setThrottlePosition(double throttlePosition, double w) {
	_w = w;
	if(throttlePosition >= 0.0 && throttlePosition <= 1.0) {
		_throttlePosition = throttlePosition;
	}
	if( _w < w_engine_min || _w > w_engine_max){
		_throttlePosition = 0.0;
	}
}

double Ecu::getFuelConsumed() const{
	return _n_fuel_consumed;
}

/**
 * mols of O2 in the "default filled cylinder"...
 * _n_Fuel = n_O2 / (Fuel_O2_req * 1.05) * _pecu->getThrottlePosition();
 */
double Ecu::fillInjector(double p, double T){
	double result = p*Vcyl/(R*T)*nu_Air[defs::O2]/ (Fuel_O2_req * 1.05) * _throttlePosition;
	_n_fuel_consumed += result;
	return result;
}
