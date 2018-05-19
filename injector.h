/*
 * Injector.h
 *
 *  Created on: May 19, 2018
 *      Author: Alex Luschan <alexander.luschan@gmail.com>
 */

#ifndef INJECTOR_H_
#define INJECTOR_H_

#include "definitions.h"
#include "gascomponent.h"
#include "ecu.h"
/*
 *
 */
class Injector {
public:
	Injector();
	Injector(Ecu *pEcu);
	const GasComponent *getFuel(double phi);
	void fill(double n_Fuel);

private:
	double _n_Fuel;
	Ecu *_pEcu;
	GasComponent *_pFuel;
};

#endif /* INJECTOR_H_ */
