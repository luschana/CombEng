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

#include "gascomponent.h"

/*
 * returns: p [Pa];
 * params: T [K], relHum [%];
 */
double getPressureFromRelHum(double T, double relHum){
	return pow(10, 3 + Water_AntoinePars[0] - Water_AntoinePars[1]/(Water_AntoinePars[2]+T)) * relHum; // [Pa] = 10^5 [bar] * [%] / 100
}

GasComponent::GasComponent(){
	_T=T_ref;
	_p=p_ref;
	_V=1.0;
	_n_g = calcMols();
	_v=_V/_n_g;
	int i = 0;
	for (i = 0; i < defs::Fuel+1; i++) {
		_nu[i]=nu_Air[i];
	}
	_MW = calcMolareWeight();
	_cp = Shomate::getInst()->getHeatCapacity(_T, _nu);
	_H = _n_g * _cp * _T;
	_combustionStarted = false;
	_isContainer = false;
}

GasComponent::GasComponent(double V, double T, double p){
	_T=T;
	_p=p;
	_V=V;
	_n_g = calcMols();
	_v=_V/_n_g;
	int i = 0;
	for (i = 0; i < defs::Fuel+1; i++) {
		_nu[i]=nu_Air[i];
	}
	_MW = calcMolareWeight();
	_cp = Shomate::getInst()->getHeatCapacity(_T, _nu);
	_H = _n_g * _cp * T;
	_combustionStarted = false;
	_isContainer = false;
}

GasComponent::GasComponent(double V, double T, double p, const double nu[defs::Fuel+1]){
	_T=T;
	_p=p;
	_V=V;
	_n_g = calcMols();
	_v=_V/_n_g;
	int i = 0;
	for (i = 0; i < defs::Fuel+1; i++) {
		_nu[i]=nu[i];
	}
	normalizeMols();
	_MW = calcMolareWeight();
	_cp = Shomate::getInst()->getHeatCapacity(_T, _nu);
	_H = _n_g * _cp * T;
	_combustionStarted = false;
	_isContainer = false;
}

GasComponent::GasComponent(double V, double T, double p, const double nu[defs::Fuel+1], bool isContainer){
	_T=T;
	_p=p;
	_V=V;
	_n_g = calcMols();
	_v=_V/_n_g;
	int i = 0;
	for (i = 0; i < defs::Fuel+1; i++) {
		_nu[i]=nu[i];
	}
	normalizeMols();
	_MW = calcMolareWeight();
	_cp = Shomate::getInst()->getHeatCapacity(_T, _nu);
	_H = _n_g * _cp * T;
	_combustionStarted = false;
	_isContainer = isContainer;
}


void GasComponent::calcGasExchange(double A_crosssection, GasComponent *pgc){
	if(A_crosssection > 0 && fabs(_p - pgc->_p)>EPSILON){
		double deltaN = 0.0;
		if(_p > pgc->_p){
			deltaN = A_crosssection * pow( (_p*(_p - pgc->_p))/(2.0*_MW*R*_T) , 0.5)*Ts;
			pgc->transferFrom(deltaN, *this);
		}else{
			deltaN = A_crosssection * pow( (pgc->_p*(pgc->_p - _p))/(2.0*pgc->_MW*R*pgc->_T) , 0.5)*Ts;
			transferFrom(deltaN, *pgc);
		}
	}
}

void GasComponent::calcStateChange(double cmpFactor, double H_cooling, double n_Fuel){
		//,const GasComponent &inlet, GasComponent &exhaust){
	double deltaH = H_cooling;
	deltaH += isentropicStateChange(cmpFactor);
	deltaH += injection(n_Fuel);
	deltaH += chemReaction();
	double dT_est = _T*deltaH/_H; // == dH/(n*cp)
	_cp = Shomate::getInst()->getHeatCapacity(_T + dT_est, _nu);
	_MW = getMolareWeight();
	_H += deltaH;
	_T = _H/(_n_g * _cp);
	if(_T > Fuel_T_Autoignition) {
		_combustionStarted = true;
	}
	_p = R*_T/_v;

}

void GasComponent::setCombustionStarted(bool combustionStarted) {
	_combustionStarted = combustionStarted;
}

bool GasComponent::isCombustionStarted() const {
	return _combustionStarted;
}

double GasComponent::getSpecHeatCapacity() const {
	return _cp;
}

double GasComponent::getEnthalpy() const {
	return _H;
}

double GasComponent::getMolareWeight() const {
	return _MW;
}

double GasComponent::getMols() const {
	return _n_g;
}

const double* GasComponent::getNu() const {
	return _nu;
}

double GasComponent::getP() const {
	return _p;
}

double GasComponent::getT() const {
	return _T;
}

double GasComponent::getspecV() const {
	return _v;
}

double GasComponent::getV() const {
	return _V;
}

// --- private methods

/*
 * cmpFactor: V_i/V_(i-1)
 */
double GasComponent::isentropicStateChange(double cmpFactor) {
	double deltaH = 0.0;
	if(fabs(cmpFactor-1.0)>EPSILON){
		_V*=cmpFactor;
		_v*=cmpFactor;
		deltaH = _H*(pow(cmpFactor, 1.0/(1.0 - _cp/R)) - 1.0);
		if(cmpFactor > 1.0){// expansion
			deltaH *= eta_is;
		}else{ // compression
			deltaH /= eta_is;
		}
	}
	return deltaH;
}

/*double GasComponent::isochoricStateChange(double deltaH) {
	return deltaH;
}*/

void GasComponent::transferFrom(double dn, GasComponent &gc) {
	if(dn > EPSILON && gc._n_g > dn){ //do not take it from an "near empty" component
		// remove gas from gc
		double dH = dn* gc._cp * gc._T;
		double cmpFactor = 1 - dn/gc._n_g;
		gc._H -= dH;
		gc._n_g -= dn;

		if(!gc._isContainer){ //isentropic expansion (V const)
			gc._H = gc._H * pow(cmpFactor, 1.0/(gc._cp/R - 1.0));
			gc._T = gc._H/(gc._n_g * gc._cp);
			gc._v = gc._V/gc._n_g;
			gc._p = R*gc._T / gc._v;
		}else{ //isobaric expansion (p,T,v const)
			gc._V *= cmpFactor;
		}

		// add to 'this'
		_H += dH;
		cmpFactor = dn/_n_g;
		int i = 0;
		for (i = 0; i < defs::Fuel+1; i++) {
			_nu[i] = (_nu[i]+cmpFactor*gc._nu[i])/(1+cmpFactor);
		}
		_n_g += dn;
		_MW = calcMolareWeight();

		if(!_isContainer){
			_H = _H * pow( (1+cmpFactor) , 1.0/(_cp/R - 1.0));
			_T = _H/(_n_g * _cp);
			_v = _V/_n_g;
			_p = R*_T / _v;
		}else{
			_V *= (1+cmpFactor);
		}
		_combustionStarted &= gc._combustionStarted;

	}
}

double GasComponent::injection(double n_Fuel){
	double deltaH = 0.0;
	double k = 1.0;
	if(n_Fuel > EPSILON){
		deltaH = n_Fuel * (-Fuel_deltaH_vap + Shomate::getInst()->getFuelHeatCapacity(T_ref)*T_ref - Shomate::getInst()->getFuelHeatCapacity(_T)*_T);
		k = n_Fuel/_n_g;
		_nu[defs::Fuel] += k;
		_p *= (1+k);
		normalizeMols();
		_v = _V / _n_g;
	}
	return deltaH;
}

double GasComponent::chemReaction(){
	double deltaH = 0.0;
	if(_combustionStarted){
		double k_Tpt =  _p/p_ref * exp( _T/T_ref - 1.0)*Ts;
		if(k_Tpt < 0) k_Tpt = 0.0;
		double k_O2 = sqrt(_nu[defs::O2]);
		// Fuel ==> 8*C + 9*H2
		_n_chemR[0] = ChemRectionRate[0]*k_Tpt*_nu[defs::Fuel];
		// H2 + 1/2*O2 ==> H20
		_n_chemR[1] = ChemRectionRate[1]*k_Tpt*_nu[defs::H2]*k_O2;
		// C + 1/2*O2 ==> C0
		_n_chemR[2] = ChemRectionRate[2]*k_Tpt*_nu[defs::C]*k_O2;
		// CO + 1/2*O2 ==> C02
		_n_chemR[3] = ChemRectionRate[3]*k_Tpt*_nu[defs::CO]*k_O2;
		deltaH =  _n_chemR[0] * Shomate::getInst()->getEnthalpyOfFormation(defs::Fuel);
		deltaH -= _n_chemR[1] * Shomate::getInst()->getEnthalpyOfFormation(defs::H2O);
		deltaH -= _n_chemR[2] * Shomate::getInst()->getEnthalpyOfFormation(defs::CO);
		deltaH -= _n_chemR[3] * (Shomate::getInst()->getEnthalpyOfFormation(defs::CO2) - Shomate::getInst()->getEnthalpyOfFormation(defs::CO));
		deltaH *= _n_g;
		_nu[defs::Fuel] -= _n_chemR[0];
		_nu[defs::H2]   += _n_chemR[0]*(Fuel_n_C+1.0) - _n_chemR[1];
		_nu[defs::H2O]  += _n_chemR[1];
		_nu[defs::C]    += _n_chemR[0]*Fuel_n_C - _n_chemR[2];
		_nu[defs::CO]   += _n_chemR[1] - _n_chemR[3];
		_nu[defs::CO2]  += _n_chemR[3];
		_nu[defs::O2]   -= 0.5*(_n_chemR[1] + _n_chemR[2] + _n_chemR[3]);
		normalizeMols();
		_v = _V / _n_g;
		if(_nu[defs::O2] < EPSILON ||
				( _nu[defs::H2] < EPSILON && _nu[defs::C] < EPSILON && _nu[defs::CO] < EPSILON && _nu[defs::Fuel] < EPSILON)){
			_combustionStarted = false;
		}
	}
	return deltaH;
}

// helper methods

double GasComponent::calcMolareWeight(){
	double MW=0.0;
	int i=0;
	for (i = 0; i < defs::Fuel+1; i++) {
		MW += _nu[i] * MolWeights[i];
	}
	return MW;
}

double GasComponent::calcMols(){
	return _p*_V/(_T*R);
}

/**
 * recalcs the num of mols
 * sum(nu_i) == 1.0 && nu_i >= 0.0;
 */
void GasComponent::normalizeMols(){
	double k = 0.0;
	int i = 0;
	for (i = 0; i < defs::Fuel+1; i++) {
		if(_nu[i] > EPSILON){
			k += _nu[i];
		}else{
			_nu[i] = 0.0;
		}
	}
	if( fabs(k-1.0) > EPSILON ){
		_n_g *= k;
		for (i = 0; i < defs::Fuel+1; i++) {
			_nu[i] /= k;
		}
	}
	_MW = calcMolareWeight();
}


