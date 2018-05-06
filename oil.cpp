/*
 * Oil.cpp
 *
 *  Created on: Apr 15, 2018
 *      Author: alex
 */

#include "oil.h"

Oil::Oil() {
	_eta_base = Oil_p[3]/Oil_p[1];
	_eta_exp = Oil_p[1];
	_k = Oil_p[2]/(Oil_p[0]-Oil_p[2]);
}

Oil::~Oil() {
}

/**
 * implementation of the Arrhenius eq.
 * eta = eta_0 exp(k/T);
 * with the 2-point definition the bases changes to eta1/eta2...
 *
 * eta = eta1 * (eta2/eta1) ^ ( (1/T - 1/T1) / (1/T2 - 1/T1) )
 *     = eta1 * _eta_base^(-_k) * _eta_base^(_k*T1/T)
 */
double Oil::getEta(double T) const {
	double eta = Oil_p[5];
	if (T < Oil_p[4]){
		eta = Oil_p[1]*pow(_eta_base, -_k)*pow(_eta_base, _k*Oil_p[0]/T);
	}else{
		if (T < 250){
			eta = 1;
		}
	}
	return eta;
}
