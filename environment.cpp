/*
 * Enviroment.cpp
 *
 *  Created on: 09.11.2014
 *      Author: alex
 */

#include "environment.h"

Environment* Environment::_pinst = NULL;

Environment::Environment() {
	_pAmbientAir = new GasComponent(pow(10.0,9.0), T_ref, p_ref, nu_Air, true);
	_pExhaustGas = new GasComponent(EPSILON, T_ref, p_ref, nu_Air, true);
}

GasComponent* Environment::getAmbientAir(){
	return _pAmbientAir;
}
GasComponent* Environment::getExhaustGas(){
	return _pExhaustGas;
}


Environment* Environment::getInst() {
	if(!_pinst){
		_pinst = new Environment;
	}
	return _pinst;
}

