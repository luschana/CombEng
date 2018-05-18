/*
 * Time slice based calculation of a combustion engine;
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

#include "combustionengine.h"

CombustionEngine* CombustionEngine::_pInst = NULL;

CombustionEngine* CombustionEngine::getInst(){
	if (_pInst == NULL){
		_pInst = new CombustionEngine();
	}
	return _pInst;
}

CombustionEngine::CombustionEngine() {
	int i = 0;
	_oil = Oil();
	_ecu = Ecu();
	_w = 0.0;
	_cnt = 0;
	_intake = GasComponent(V_intake, Environment::getInst()->getAmbientAir()->getT(),
			Environment::getInst()->getAmbientAir()->getP(), Environment::getInst()->getAmbientAir()->getNu());
	_exhaust = GasComponent(V_exhaust, Environment::getInst()->getAmbientAir()->getT(),
			Environment::getInst()->getAmbientAir()->getP(), Environment::getInst()->getAmbientAir()->getNu());
	_pValveIntake = new Valve(A_intake, num_Valve, 0.0, 4*M_PI);
	_pValveExhaust = new Valve(A_exhaust, num_Valve, 0.0, 4*M_PI);
//	_thrPos = 0.0;
	_T_CW = 340.0;
	_M_Shaft = 0.0;
	for (i = 0; i < Ncyl; i++) {
		_phiCyl[i] = -4*M_PI/Ncyl*(double)i;
		_cyl[i] = Cylinder(_phiCyl[i], &_intake, &_exhaust, &_ecu, &_oil);
	}

}
/**
 * Parameters:
 * w ... speed of rotation in [rad/s];
 * thrPos ... throttle position, absolute value ==> [0..1]
 */
void CombustionEngine::run(double w, double thrPos) {
	int i = 0;
	_ecu.setThrottlePosition(thrPos, w);
	_M_Shaft = 0.0;
	double dphi = (_w + w)/2.0 * Ts;
	_w = w;
	_cnt++;

	_pValveIntake->calcFlow(1.0, Environment::getInst()->getAmbientAir(), &_intake);
	_pValveExhaust->calcFlow(1.0, &_exhaust, Environment::getInst()->getExhaustGas());

	for (i = 0; i < Ncyl; i++) {
		_cyl[i].run(dphi);
		_M_Shaft += _cyl[i].getM_G() + _cyl[i].getM_P();
	}

	_intake.calcStateChange(1.0 , 0.0, 0.0, _pValveIntake->getGasComponent(), _pValveExhaust->getGasComponent());

	//Environment::getInst()->getAmbientAir()->calcGasExchange(A_intake, &_intake);
	//_exhaust.calcGasExchange(A_exhaust, Environment::getInst()->getExhaustGas());
	/* use amb for tests
	Environment::getInst()->getExhaustGas()->calcGasExchange(A_exhaust, &_exhaust);
	*/
//	Environment::getInst()->getAmbientAir()->calcGasExchange(A_exhaust, &_exhaust);
}

void CombustionEngine::setPhiSpark(double sparkAngle) {
	_ecu.setPhiSpark(sparkAngle);
}

void CombustionEngine::setPhiInjection(double injectionAngle){
	_ecu.setPhiInjection(injectionAngle);
}


const GasComponent& CombustionEngine::getExhaust() const {
	return _exhaust;
}

const GasComponent& CombustionEngine::getIntake() const {
	return _intake;
}

double CombustionEngine::getThrPos() const {
	return _ecu.getThrottlePosition();
}

double CombustionEngine::getW() const {
	return _w;
}

double CombustionEngine::getMShaft() const {
	return _M_Shaft;
}

double CombustionEngine::getPhi() const {
	return _cyl[0].getPhi();
}

double CombustionEngine::getp_Cyl(int cylNum) const {
	return _cyl[cylNum].getPressure();
}

double CombustionEngine::getT_Cyl(int cylNum) const {
	return _cyl[cylNum].getTemperature();
}

double CombustionEngine::getCyl_T(int cylNum) const {
	return _cyl[cylNum].getT_Cyl();
}

double CombustionEngine::getFuelConsumption() const {
	return _ecu.getFuelConsumed();
}

double CombustionEngine::getH_Cooling() const {
	double result = 0.0;
	for (int i = 0; i < Ncyl; i++) {
		result += _cyl[i].getH_Cool();
	}
	return result;
}

unsigned long CombustionEngine::getCnt() const {
		return _cnt;
	}

void CombustionEngine::setTempCoolingWater(double T_CW){
	_T_CW = T_CW;
	for (int i = 0; i < Ncyl; i++) {
		_cyl[i].setT_CW(_T_CW);
	}
}


