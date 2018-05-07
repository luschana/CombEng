/*
 * Oil.h
 *
 *  Created on: Apr 15, 2018
 *      Author: alex
 */

#ifndef OIL_H_
#define OIL_H_

#include "definitions.h"

class Oil {
public:
	Oil();
	~Oil();
	double getEta(double T) const;

private:
	double _eta_base, _eta_exp, _k;
};

#endif /* OIL_H_ */
