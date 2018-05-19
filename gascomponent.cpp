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
	return pow(10, 3 + WaterAntoinePars[0] - WaterAntoinePars[1]/(WaterAntoinePars[2]+T)) * relHum; // [Pa] = 10^5 [bar] * [%] / 100
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

GasComponent* GasComponent::getFuelComponent() {
	GasComponent *result = new GasComponent();
	for (int i = 0; i < defs::Fuel+1; i++) {
		result->_nu[i]=0.0;
	}
	result->_nu[defs::Fuel] = 1.0;
	result->_n_g = 0.0;
	result->_H = 0.0;
	result->_cp = Shomate::getInst()->getFuelHeatCapacity(result->_T);
	result->_MW = MolWeights[defs::Fuel];
	return result;
}

/*
 * set moles to be n_Fuel and H to the enthalpy of its vaporization
 */
void GasComponent::setFuelComponent(double n_Fuel) {
	_n_g = n_Fuel;
	_H = - _n_g * Fuel_deltaH_vap;
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

/*
 * calculate the "exchanged gas component"; Bernoulli in its basic form
 */
void GasComponent::calcFlow(double A_crosssection, GasComponent* pIn, GasComponent* pOut){
	_n_g=0.0;
	if(A_crosssection > EPSILON){
		GasComponent *pSrc = pIn;
		if(pIn->_p < pOut->_p){
			pSrc = pOut;
		}
		_n_g = A_crosssection * pow( (2.0*(pIn->_p - pOut->_p)*pSrc->_MW/pSrc->_v) , 0.5)*Ts; // neg. flow ok
		_p = pSrc->_p;// orig version: _p = pDest->_p; the kinetic energy is recuperated, so...
		_T = pSrc->_T;
		setNu(pSrc);
		_cp = pSrc->_cp;
		_MW = pSrc->_MW;
		_v = pSrc->_v; //R*T/p...
		_H = _n_g*_cp*_T;
		_V = _n_g*_v;
	}
}
void GasComponent::calcStateChange(bool add, const GasComponent *pgc){
	double dH = 0.0;
	if(fabs(pgc->_n_g)>EPSILON){
		if(add){
			dH=addGC(pgc);
		}else{
			dH=removeGC(pgc);
		}
	}
}

/*
 * keep track of total enthalpy and all the mole based values; state vals are calculated at the end...
 */
void GasComponent::calcStateChange(double cmpFactor, double H_cooling, const GasComponent *pFuel, const GasComponent *pIntake, const GasComponent *pExhaust){
	if(_T > Fuel_T_Autoignition) {
		_combustionStarted = true;
	}
	double deltaH = H_cooling;
	deltaH += isentropicStateChange(cmpFactor); // changes v&V
	if(pExhaust->_n_g > 0){
		deltaH += removeGC(pExhaust);
	}else{
		deltaH += addGC(pExhaust);
	}
	if(pIntake->_n_g > 0){
		deltaH += addGC(pIntake);
	}else{
		deltaH += removeGC(pIntake);
	}
	deltaH += injection(pFuel);
	deltaH += chemReaction();


	double dT_est = deltaH/(_n_g*_cp); // == dH/(n*cp)
	if(_T + dT_est < 200) {
		dT_est = 200 - _T; //debugging / Fangnetz -- should never occur!!
		deltaH = dT_est * _cp * _n_g; // changed T&p to be nan instead of -0;
	}
	_cp = Shomate::getInst()->getHeatCapacity(_T + dT_est, _nu);
	_MW = getMolareWeight();
	_H += deltaH;
	_T = _H/(_n_g * _cp);
	_p = R*_T/_v;
}

// --- private methods

/*
 * there is an isentropic/adiabatic expansion to the free space left by the removed element
 */
double GasComponent::removeGC(const GasComponent *pgc){
	double dH = 0.0;
	double dn = fabs(pgc->_n_g);
	if(dn>0){
		if(_n_g > dn){
			_n_g -= dn;
			dH = -dn*_cp*_T;
		}else{ // it should be more removed than there's available?!?
			dH = -(_n_g-EPSILON)*_cp*_T;
			_n_g = EPSILON;
		}
	}
	return dH;
}

/*
 * calc the missing enthalpy to bring the gc to gas temp and mix the moles
 */
double GasComponent::addGC(const GasComponent *pgc){
	double dH = 0.0;
	double dn = fabs(pgc->_n_g);
	if(dn>0){
		// nonsense: adds just the enthalpy of gc at its current temp // bring gc to _T (H(T) = n cp(T) T = H(T_old) + dH -- dH was missing ==> neg. sign
		dH = pgc->_H;// - dn*Shomate::getInst()->getHeatCapacity(_T, pgc->_nu)*_T;
		for (int i = 0; i<defs::Fuel+1; i++) {
			_nu[i] = (_nu[i]*_n_g + pgc->_nu[i]*dn)/(_n_g + dn);
		}
	}
	return dH;
}

/*
 * calc the enthalpy needed for liquid->gas state change;
 */
double GasComponent::injection(const GasComponent * pFuel){
	double deltaH = 0.0;
	double k = 1.0;
	if(pFuel->_n_g > EPSILON){
		deltaH = -pFuel->_n_g * Fuel_deltaH_vap;//  + Shomate::getInst()->getFuelHeatCapacity(T_ref)*T_ref - Shomate::getInst()->getFuelHeatCapacity(_T)*_T);
		addGC(pFuel);
	}
	return deltaH;
}


/* calculates the enthalpy difference but does not change the state values
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

/*
 * evaluate state change from change of enthalpy;
 * V = const; p/T = const;
 * _cp needs to be "upToDate"!
 */
void GasComponent::isochoricStateChange(double deltaH) {
	if(fabs(deltaH)>0){
		_p *= (1.0 + deltaH/_H); // results from: p1/p2 = T1/T2 = cp1*T1/cp2*T2;
		_H +=deltaH;
		// T = f(deltaH)...
		double dT_est = deltaH/(_n_g*_cp);
		_cp = Shomate::getInst()->getHeatCapacity(_T+dT_est, _nu); // cp close to the new temperature
		_T = _H/(_n_g*_cp);
	}
}

/*
 * isentropic state change and enthalpy taken from itself
 */
void GasComponent::adiabaticStateChange(double cmpFactor) {
	isochoricStateChange(isentropicStateChange(cmpFactor));
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

/*
 * no checks, just write it...
 */
void GasComponent::setNu(const GasComponent *pSrc){
	for (int i = 0; i < defs::Fuel+1; i++) {
		_nu[i] = pSrc->_nu[i];
	}
}

