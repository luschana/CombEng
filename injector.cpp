/*
 * Injector.cpp
 *
 *  Created on: May 19, 2018
 *      Author: Alex Luschan <alexander.luschan@gmail.com>
 */

#include "injector.h"


Injector::Injector() {
	_pFuel = GasComponent::getFuelComponent();
	_pEcu = NULL;
	_n_Fuel = 0.0;
}

Injector::Injector(Ecu* pEcu) {
	_pFuel = GasComponent::getFuelComponent();
	_pEcu = pEcu;
	_n_Fuel = 0.0;
}

const GasComponent* Injector::getFuel(double phi) {
    _pFuel->setFuelComponent(0.0);
	if(_n_Fuel > EPSILON*Fuel_n_Inject && phi >= _pEcu->getPhiInjection()){
		if(_n_Fuel > Fuel_n_Inject){
			_pFuel->setFuelComponent(Fuel_n_Inject);
			_n_Fuel -= Fuel_n_Inject;
		}else{
			_pFuel->setFuelComponent(_n_Fuel);
			_n_Fuel = 0.0;
		}
	}
	return _pFuel;
}

void Injector::fill(double n_Fuel) {
	_n_Fuel = n_Fuel;
}
