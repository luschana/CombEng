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
	_dirtyFlag = false;
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
	_dirtyFlag = false;
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
	_dirtyFlag = false;
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
	_dirtyFlag = false;
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
	result->_dirtyFlag = false;
	return result;
}

/*
 * set moles to be n_Fuel and H to the enthalpy of its vaporization
 */
void GasComponent::setFuelComponent(double n_Fuel) {
	_n_g = n_Fuel;
	_H = - _n_g * Fuel_deltaH_vap;
	_dirtyFlag = false;
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
	if(A_crosssection > EPSILON && fabs(pIn->_p-pOut->_p) > EPSILON*p_ref){
		GasComponent *pSrc = pIn;
		GasComponent *pDest = pOut;
		bool reverse = (pIn->_p < pOut->_p);
		if(reverse){
			pSrc = pOut;
			pDest = pIn;
		}
		_n_g = A_crosssection * sqrt( 2.0*fabs(pIn->_p - pOut->_p)*pSrc->_MW/pSrc->_v )*Ts; // sqrt of delta p!!!
		if((pSrc->_n_g - _n_g) < EPSILON) {
			_n_g = pSrc->_n_g - EPSILON;
		}
		_p = pDest->_p;
		_T = pSrc->_T;
		setNu(pSrc);
		_cp = pSrc->_cp;
		_MW = pSrc->_MW;
		_v = pSrc->_v; //R*T/p...
		_H = _n_g*_cp*_T;
		_V = _n_g*_v;
		if(reverse) _n_g *= -1; // neg. flow ok
		_dirtyFlag = pSrc->_dirtyFlag;
	}
}

/*
 * get all Valve GCs for intake/exhaust at once -- all cylinders + "environment valve"
 */
void GasComponent::calcStateChange(bool *add, const GasComponent **pgc){
	for (int i = 0; i < Ncyl+1; i++) {
		if(!add[i] && fabs(pgc[i]->_n_g)>EPSILON*Ts){
			removeGC(pgc[i]);
		}
	}
	for (int i = 0; i < Ncyl+1; i++) {
		if(add[i] && fabs(pgc[i]->_n_g)>EPSILON*Ts){
			addGC(pgc[i]);
		}
	}

	double T_est =_H/(_n_g*_cp); // (_H_act - _H_old)/(_n_g*_cp)
	cleanUpStateChange(T_est);
}

/*
 * keep track of total enthalpy and all the mole based values; state vals are calculated at the end...
 */
void GasComponent::calcStateChange(double cmpFactor, double H_cooling, const GasComponent *pFuel, const GasComponent *pIntake, const GasComponent *pExhaust){
	double deltaH = _H; //used as storage in the first place
	_H += isentropicStateChange(cmpFactor); // changes v&V
	// remove components first, than add
	if(pIntake->_n_g < 0){
		removeGC(pIntake);
	}
	if(pExhaust->_n_g > 0){
		removeGC(pExhaust);
	}else{
		addGC(pExhaust);
	}
	if(pIntake->_n_g > 0){
		addGC(pIntake);
	}
	addGC(pFuel);
	_H += chemReaction();
	_H += H_cooling;

	deltaH = _H - deltaH; //now it's the delta

	double T_est = _T + deltaH/(_n_g*_cp);
	cleanUpStateChange(T_est);
}

// --- private methods
/*
 * "harmonize" what got unsorted by calc* methodes
 */
void GasComponent::cleanUpStateChange(double T_estimate){
	if(!_isContainer){
		if(_dirtyFlag){
			_MW = calcMolareWeight();
			_v = _V/_n_g;
		}
		if(_dirtyFlag || fabs(T_estimate - _T)>EPSILON){
			_cp = Shomate::getInst()->getHeatCapacity(T_estimate, _nu);
			_T = _H/(_n_g * _cp);
			_p = R*_T/_v;
			_dirtyFlag = false;
		}
	}else{
		if(_dirtyFlag || fabs(T_estimate - _T)>EPSILON){
			_cp = Shomate::getInst()->getHeatCapacity(T_estimate, _nu);
			_T = _H/(_n_g * _cp);
		}
		if(_dirtyFlag){
			_MW = calcMolareWeight();
			_v = R*_T/_p;
			_V = _v * _n_g;
			_dirtyFlag = false;
		}

	}
}

/*
 * take away the moles and enthalpy
 */
void GasComponent::removeGC(const GasComponent *pgc){
	double dn = fabs(pgc->_n_g);
	if(dn>EPSILON*Ts){
		if(_n_g > dn){
			_n_g -= dn;
			_H -= dn*_cp*_T;
		}else{ // more should be removed than there's available?!?
			_H -= (_n_g-EPSILON)*_cp*_T;
			_n_g = EPSILON;
		}
		_H += isentropicEnthalpyChange((_n_g + dn)/_n_g); // v1/v0 = n0/n1
		_dirtyFlag = true;
	}
}

/*
 * molare mixture and add mole/enthalpy;
 */
void GasComponent::addGC(const GasComponent *pgc){
	double dn = fabs(pgc->_n_g);
	if(dn>EPSILON*Ts){
		_H += pgc->_H;
		_H += isentropicEnthalpyChange(_n_g/(_n_g + dn)); // V0/V1 ~ n1/n0
		for (int i = 0; i<defs::Fuel+1; i++) {
			_nu[i] = (_nu[i]*_n_g + pgc->_nu[i]*dn)/(_n_g + dn);
		}
		_n_g += dn;
		_dirtyFlag = true;
	}
}

/* calculates the enthalpy difference
 * cmpFactor: v1/v0
 * H1/H0~T1/T0~(v0/v1)^(k-1) = (v1/v0)^(1-k)
 */
double GasComponent::isentropicEnthalpyChange(double cmpFactor) {
	double deltaH = 0.0;
	if(fabs(cmpFactor-1.0)>EPSILON){
		deltaH = _H*(pow(cmpFactor, 1.0/(1.0 - _cp/R)) - 1.0);
	}
	return deltaH;
}

/* calculates the enthalpy difference but does not change the state values
 * cmpFactor: V_i/V_(i-1)
 */
double GasComponent::isentropicStateChange(double cmpFactor) {
	double deltaH = isentropicEnthalpyChange(cmpFactor);
	if(fabs(cmpFactor-1.0)>EPSILON){
		_V*=cmpFactor;
		_v*=cmpFactor;
		if(cmpFactor > 1.0){// expansion
			deltaH *= eta_is;
		}else{ // compression
			deltaH /= eta_is;
		}
		_dirtyFlag = true;
	}
	return deltaH;
}


double GasComponent::chemReaction(){
	double deltaH = 0.0;
	if(_T > Fuel_T_Autoignition) {
		_combustionStarted = true;
	}
	// check if combustion is to be turned off
	if(_combustionStarted &&
			((_nu[defs::O2] < EPSILON && _nu[defs::Fuel] < EPSILON) ||
			( _nu[defs::H2] < EPSILON && _nu[defs::C] < EPSILON && _nu[defs::CO] < EPSILON && _nu[defs::Fuel] < EPSILON))){
		_combustionStarted = false;
	}

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
		_dirtyFlag = true;
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
		_MW = calcMolareWeight();
		_dirtyFlag = true;
	}
}

/*
 * no checks, just write it...
 */
void GasComponent::setNu(const GasComponent *pSrc){
	for (int i = 0; i < defs::Fuel+1; i++) {
		_nu[i] = pSrc->_nu[i];
	}
	_dirtyFlag = true;
}
