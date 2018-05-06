/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2014  Alex Luschan <a.luschan@inode.at>
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

#include "cylinder.h"

Cylinder::Cylinder(){
	_phi = 0.0;
	_x_p = 0.0;
	_dx_p = 0.0;
	_v_p = 0.0;
	_M_p = 0.0;
	_M_g = 0.0;
	_F_fr = 0.0;
	_n_Fuel = 0.0;
	_H_cooling = 0.0;
	_H_hx_gas = 0.0;
	_T_cyl = T_ref;
	_T_CW = T_ref;
	_gc = GasComponent(Vcyl*(1+1/chi), T_ref, p_ref);
	_pinlet = NULL;
	_pexhaust = NULL;
	_pecu = NULL;
}

Cylinder::Cylinder(double phi0, GasComponent *pinlet, GasComponent *pexhaust, Ecu *pecu, Oil *poil){
	_phi = phi0;
	_x_p = 0.0;
	_dx_p = 0.0;
	_v_p = 0.0;
	_M_p = 0.0;
	_M_g = 0.0;
	_F_fr = 0.0;
	_n_Fuel = 0.0;
	_H_cooling = 0.0;
	_H_hx_gas = 0.0;
	_T_cyl = T_ref;
	_T_CW = T_ref;
	_gc = GasComponent(Vcyl*(1+1/chi), T_ref, p_ref);
	_pinlet = pinlet;
	_pexhaust = pexhaust;
	_pecu = pecu;
	_poil = poil;
}

void Cylinder::setT_CW(double T_CW) {
	_T_CW = T_CW;
}

double Cylinder::getH_Cool() const {
	return _H_cooling;
}

double Cylinder::getX_p() const {
	return _x_p;
}

double Cylinder::getdX_p() const {
	return _dx_p;
}

double Cylinder::getV_p() const {
	return _v_p;
}

double Cylinder::getdV_p() const {
	return _dv_p;
}

double Cylinder::getM_G() const {
	return _M_g;
}

double Cylinder::getM_P() const {
	return _M_p;
}

double Cylinder::getT_Cyl() const {
	return _T_cyl;
}

double Cylinder::getPhi() const {
	return _phi;
}

double Cylinder::getPressure() const {
	return _gc.getP();
}

double Cylinder::getTemperature() const {
	return _gc.getT();
}

void Cylinder::calcCylinder(double dPhi){
	_phi += dPhi;
	if(_phi >= 4.0*M_PI){
		_phi -= 4.0*M_PI;
	}
	if(_phi < 0.0){
		_phi += 4.0*M_PI;
	}
	_dx_p = r_cs*( sqrt(pow(l2r,2.0) - pow(sin(_phi),2.0)) - cos(_phi) -(l2r-1.0)) -_x_p; // change in x pos (from 0)
	_x_p += _dx_p;
	_dv_p = _v_p - _dx_p/Ts;
	_v_p = _dx_p/Ts;
	//force of friction: F = eta(T) * A * v/d; M = F*r(phi)
	_F_fr = _poil->getEta(_T_cyl)* 2*r_cs*M_PI*h_Piston * _v_p / d_Piston;
	_M_p = m_Piston*_dv_p / Ts * r_cs * sin(_phi)- fabs(_F_fr * r_cs*sin(_phi)); // speed dep. && friction
	if(passedAngle(_pecu->getPhiInjection(), dPhi)){
		_n_Fuel = _pecu->fillInjector(_pinlet->getP(), _pinlet->getT());
	}
	if(passedAngle(_pecu->getPhiSpark(), dPhi)){
		_gc.setCombustionStarted(true);
	}
	_gc.calcGasExchange(_pecu->getValveOut_A(_phi), _pexhaust);
	_pinlet->calcGasExchange(_pecu->getValveIn_A(_phi), &_gc);
	_gc.calcStateChange(getCmpFactor() , getHeatExchangeEnthalpy(), calcFuelInj());//, *_pinlet, *_pexhaust);
	_M_g = -ACyl*(_gc.getP() - Environment::getInst()->getAmbientAir()->getP())*r_cs*sin(_phi);
}

double Cylinder::getCylArea(double x_pos) const {
	return 2*ACyl*(1.0 + (hCyl-_x_p));
}

/*
 * Heat exchange
 *  gas <-> cyl. wall <-> cooling water
 */
double Cylinder::getHeatExchangeEnthalpy(){
	double hx_factor = _gc.getspecV()/(R*T_ref/p_ref) * getCylArea(hCyl - _x_p); // v_ref/v(p,T) * A [m3/mol / m3/mol * m2]
	_H_cooling = hx_a_CW * (_T_cyl - _T_CW)*Ts;
	_H_hx_gas = hx_a_CG* hx_factor * (_T_cyl - _gc.getT())*Ts;
	_T_cyl += (-_H_hx_gas - _H_cooling + fabs(_F_fr*_dx_p))/hx_C_W; // gas hx && friction
	return _H_hx_gas;
}

double Cylinder::getCmpFactor() const {
	return 1.0/ (1.0 + _dx_p/(hCyl-_x_p));
}

bool Cylinder::passedAngle(double alpha, double dphi) const {
	return (alpha < _phi && _phi <= alpha + dphi);
}

double Cylinder::getHxGas() const {
	return _H_hx_gas;
}
/**
 * calc the mols of fuel to be injected
 */
double Cylinder::calcFuelInj() {
	double result = 0.0;
	if(_n_Fuel > EPSILON && _phi >= _pecu->getPhiInjection()){
		if(_n_Fuel > Fuel_n_Inject){
			result = Fuel_n_Inject;
			_n_Fuel -= Fuel_n_Inject;
		}else{
			result = _n_Fuel;
			_n_Fuel = 0.0;
		}
	}
	return result;
}

