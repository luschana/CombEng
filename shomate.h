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

#ifndef SHOMATE_H
#define SHOMATE_H
#include <iostream>
#include <stdio.h>
#include <stdlib.h>


// chemical substances
namespace defs{
  enum eChemList {N2=0, O2=1, H2O=2, CO2=3, CO=4, H2=5, C=6, Fuel=7};
}

namespace shomate{
	enum shName{A ,B ,C ,D ,E ,F ,G ,Hf};
}

class TRange
{
public:
  double T_min_, T_max_;
  TRange();
  TRange(double Tmin, double Tmax);
  void set(double Tmin, double Tmax);
};

/* unit of Hf: [J/mol]*/
class ShParDef{
public:
  ShParDef();
  ShParDef(double A,double B,double C,double D,double E,double F,double G,double Hf);
  void set(double A,double B,double C,double D,double E,double F,double G,double Hf);
  //void add(double *factor, ShParDef * pars);
  void add(const double factor, const ShParDef *pars);
  void reset();
  double data_[shomate::Hf+1];
};

class ShDataEntry{
public:
	TRange T_;
	ShParDef pars_;
	ShDataEntry(double Tmin, double Tmax, double A,double B,double C,double D,double E,double F,double G,double Hf);
};

class ShData
{
public:
	ShData();
	ShData(int dataSets);
	ShData(ShDataEntry shEntry);
	ShData(int dataSets, const ShDataEntry shEntry[]);
	~ShData();
	void set(int dataSet, double Tmin, double Tmax, double A,double B,double C,double D,double E,double F,double G,double Hf);
	TRange T_[5];
	ShParDef pars_[5];
	const ShParDef* getParams(double T)  const ;
private:
	int defs_;
};

class Shomate{
public:
	static Shomate* getInst();

	double getHeatCapacity(double T, const double *pn_def);
	double getFuelHeatCapacity(double T);
	double getEnthalpyOfFormation(enum defs::eChemList substance);

private:
	Shomate();
	static Shomate* _pInst;

	ShParDef* _pActShParams;
};

#endif // SHOMATE_H
