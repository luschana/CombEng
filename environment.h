/*
 * Enviroment.h
 *
 *  Created on: 09.11.2014
 *      Author: alex
 */

#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

#include "definitions.h"
#include "gascomponent.h"

class Environment {
public:
	static Environment* getInst();
	GasComponent* getAmbientAir();
	GasComponent* getExhaustGas();
/*	const GasComponent*& getAmbientAir() const;
	const GasComponent*& getExhaustGas() const;
*/
private:
	Environment();
	static Environment* _pinst;

	GasComponent* _pAmbientAir;
	GasComponent* _pExhaustGas;
};



#endif /* ENVIRONMENT_H_ */
