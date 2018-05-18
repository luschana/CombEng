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
#include "shomate.h"
#include "definitions.h"

TRange::TRange(){
  T_min_ = 0.0;
  T_max_ = 0.0;
}

TRange::TRange(double Tmin, double Tmax){
  TRange();
  set(Tmin, Tmax);
}

void TRange::set(double Tmin, double Tmax){
  T_min_ = Tmin;
  T_max_ = Tmax;
}



ShParDef::ShParDef(){
}

ShParDef::ShParDef(double A, double B, double C, double D, double E, double F, double G, double Hf){
  set(A, B, C, D, E, F, G, Hf);
}

void ShParDef::set(double A, double B, double C, double D, double E, double F, double G, double Hf){
  data_[shomate::A] = A;
  data_[shomate::B] = B;
  data_[shomate::C] = C;
  data_[shomate::D] = D;
  data_[shomate::E] = E;
  data_[shomate::F] = F;
  data_[shomate::G] = G;
  data_[shomate::Hf]= Hf;
}

void ShParDef::add(const double factor, const ShParDef * pars) {
	for (int i = 0; i <= shomate::Hf; i++) {
		if(pars != 0){ // && (pars->data_[i]) > 0.0){
			data_[i] += (pars->data_[i])* factor;
		}
	}
}

void ShParDef::reset(){
  set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}

ShDataEntry::ShDataEntry(double Tmin, double Tmax, double A, double B, double C, double D, double E, double F, double G, double Hf){
	  T_ = TRange(Tmin, Tmax);
	  pars_= ShParDef(A, B, C, D, E, F, G, Hf);
}

ShData::ShData(){
	defs_ = 0;
}

ShData::ShData(int dataSets){
  if(dataSets > 0 && dataSets <= 5){
	  defs_ = dataSets;
  }
}

ShData::ShData(ShDataEntry shEntry){
	defs_ = 1;
	T_[0] = shEntry.T_;
	pars_[0] = shEntry.pars_;
}

ShData::ShData(int dataSets, const ShDataEntry shEntry[]){
	defs_ = 0;
	if(dataSets > 0 && dataSets <= 5){
		  defs_ = dataSets;
		  int i = 0;
		  for (i = 0; i < defs_; i++) {
			  T_[i] = shEntry[i].T_;
			  pars_[i] = shEntry[i].pars_;
		}
	  }
}

void ShData::set(int dataSet, double Tmin, double Tmax, double A, double B, double C, double D, double E, double F, double G, double Hf){
  if(dataSet >=0 && dataSet < defs_){
	  T_[dataSet].set(Tmin, Tmax);
	  pars_[dataSet].set(A, B, C, D, E, F, G, Hf);
  }
}

const ShParDef * ShData::getParams(double T) const  {
  const ShParDef *result = NULL;
  if(defs_ > 0 ){
	  if(defs_ == 1 || T <= T_[0].T_max_){
		result = &pars_[0];
	  } else if(defs_ == 2 ){
		result = &pars_[1];
	  } else{
		for(int i = 1; i<defs_;i++){
		  if(T <= T_[i].T_max_){
			result = &pars_[i];
			break;
		  }
		}
	  }
  }
  return result;
}

ShData::~ShData(){
}


Shomate* Shomate::_pInst = NULL;

Shomate::Shomate(){
	_pActShParams = new ShParDef();
	_pActShParams->reset();
}

double Shomate::getHeatCapacity(double T, const double *pn_def){
  double cp = 0.0;
  if (T < Shomate_T_min) T=Shomate_T_min;
  if (T > Shomate_T_max) T=Shomate_T_max;
  double t = T/1000.0;
  int i = 0;
  _pActShParams->reset();
  for(i=0;i<=defs::Fuel;i++){
	  _pActShParams->add( *(pn_def+i), ShDataDB[i].getParams(T));
  }
  cp= _pActShParams->data_[shomate::A] + _pActShParams->data_[shomate::B]*t + _pActShParams->data_[shomate::C]*pow(t,2.0)
		  + _pActShParams->data_[shomate::D]*pow(t,3.0) + _pActShParams->data_[shomate::E]/pow(t,2.0);
  return cp;
}

double Shomate::getFuelHeatCapacity(double T){
	//return ShDataDB[defs::Fuel].pars_[0].data_[shomate::A] + ShDataDB[defs::Fuel].pars_[0].data_[shomate::B] * T/1000.0
	//		+ ShDataDB[defs::Fuel].pars_[0].data_[shomate::C] + ShDataDB[defs::Fuel].pars_[0].data_[shomate::D] * T/1000.0;
	  double cp = 0.0;
	  if (T < Shomate_T_min) T=Shomate_T_min;
	  if (T > Shomate_T_max) T=Shomate_T_max;
	  double t = T/1000.0;
	  const ShParDef *pfuel = ShDataDB[defs::Fuel].getParams(T);
	  cp= pfuel->data_[shomate::A] + pfuel->data_[shomate::B]*t + pfuel->data_[shomate::C]*pow(t,2.0) + pfuel->data_[shomate::D]*pow(t,3.0) + pfuel->data_[shomate::E]/pow(t,2.0);
	  return cp;
}

double Shomate::getEnthalpyOfFormation(enum defs::eChemList substance){
	return ShDataDB[substance].pars_[0].data_[shomate::Hf];
}

Shomate* Shomate::getInst() {
	if (!_pInst){
			_pInst = new Shomate();
		}
		return _pInst;
}





