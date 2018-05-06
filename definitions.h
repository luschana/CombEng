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

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

//#include <stdio.h>
#include <math.h>
#include "shomate.h"

// calc consts
const double Ts = 10.0*pow(10.0,-6.0); //[s] cycle time
const double EPSILON = pow(10.0,-6.0); // max deviation of double numbers

// engine consts
// mechanical
const int Ncyl = 1; // number of cylinders
const double Vcyl = 0.5*pow(10.0,-3.0); //[m3] cylinder volume
const double chi = 18.0; // compression factor
const double w_engine_min = 150.0*M_PI/30.0; // minimum engine speed -- injection starting point
const double w_engine_max = 12000.0*M_PI/30.0; // maximum engine speed -- injection end point
const double l2r = 2.5; //[m/m] pleuel to crank shaft radius
const double r_cs = pow(Vcyl/(2.0*M_PI),1.0/3.0);// crank shaft radius =(Vcyl/2pi)^1/3: V_h(r^2 pi 2r);
const double m_Piston = 0.35; // [kg] piston mass
const double h_Piston = 2.0*r_cs; // height of piston (https://www.thm.de/me/images/user/herzog-91/Kolbenmaschinen/Kolbenmaschinen_6_Konstruktionselemente.pdf)
const double d_Piston = pow(10.0,-4.0); // distance between piston and cylinder [m];
const double ACyl = pow(r_cs,2.0) * M_PI;
const double hCyl = 2.0*r_cs*(1.0+1.0/chi); //[m] height from lower dead center to cylinder head
// fluid dynamics
const double k_Aval = 0.2; //ratio A_valves/A_cyl
const double A_Valve_out = k_Aval * ACyl; // cross section of outlet valves
const double A_Valve_in = k_Aval * ACyl; // cross section of inlet valves
const int num_Valve = 2; // number of inlet and outlet valve per cylinder;

const double A_intake = A_Valve_in * (1.0+Ncyl/2.0); // cross section of intake
const double A_exhaust = A_intake; // (min.) cross section of the exhaust gas pipe
const double V_intake = Vcyl*(10.0+Ncyl*2.0); //[m3] volume of intake manifold
const double V_exhaust = V_intake*2.0; //[m3] volume of exhaust pipe
// heat transfer
const double hx_a_CW = 400.0; // [J/(K s)] heat transfer coefficient cooling water <-> wall
const double hx_a_CG = 2000.0; // [J/(K m2 s)] heat transfer coefficient wall <-> gas (delta T ~300K; 10kW thermal flow rate)
const double hx_C_W = 15000.0; // [J/K] heat capacity of wall
//const double hx_gamma = hx_C_W / hx_a_CW; // ratio hx_C_W / hx_a_CW
//const double hx_alpha = hx_a_CG / hx_a_CW; // ratio hx_a_WG / hx_a_CW

// standard definitions
const double R = 8.314462175; // [J/mol K]
const double p_ref = pow(10.0,5.0); // [Pa]
const double T_ref = 295.15;// [K]
const double nu_Air[defs::Fuel+1]={0.7825, 0.2099, 0.0076, 0.0, 0.0, 0.0, 0.0, 0.0}; // [mol/mol] N2, O2, H2O, CO2, CO
const double eta_is = 0.95; // isentropic efficiency
const double MolWeights[defs::Fuel+1] = {0.028, 0.032, 0.018, 0.044, 0.028, 0.002, 0.012, 0.114};//molare weight [kg/mol]

// oil definitions
// https://de.wikipedia.org/wiki/Viskosit%C3%A4t#Typische_Werte; last 2 values are "burnt oil" -- blocking eng?!
const double Oil_p[6] = {273.15 + 25.0, 0.1, 273.15 + 150.0, 0.003, 500.0, 10}; // T1, eta1, T2, eta2, end of life...

// fuel definitions
//static const double Hf_ref = 50.0*1000.0*1000.0; // [J/kg]
  //defined in ShData_Fuel: static const double H_fuel_form = -208700; // [J/mol] enthalpy of formation (http://webbook.nist.gov/cgi/cbook.cgi?ID=C111659&Units=SI&Mask=1)
const double Fuel_n_C = 10.0;
const double Fuel_O2_req = Fuel_n_C * 1.5 + 0.5; //[mol_O2/mol_Fuel] mols of O2 for stoichiometric reaction of an alkane
// injection
const double Fuel_n_Inject = nu_Air[defs::O2]*(p_ref*Vcyl/(R*T_ref))/Fuel_O2_req*(Ts/(2*pow(10.0,-3.0))); //[mol/s] molare amount injected per sample (def.: duration = 2 ms for "std filled" cyl)
const double Fuel_T_Autoignition = 273.15 + 255.0; //https://de.wikipedia.org/wiki/Z%C3%BCndtemperatur

// phys/chem/math consts [Fuel, H2O, CO, CO2]
//const double ChemRectionRate[4] = {12.0, 40.0, 10.0, 5.0}; // was too slow at ~5000rpm
const double ChemRectionRate[4] = {50.0, 40.0, 10.0, 3.0}; // reaction rate at T&p_ref

//const double ChemRectionRate[4] = {0.40, 0.026, 0.026, 0.026}; // reaction rate at T&p_ref
//const double ChemRectionRate[4] = {30.0, 30.0, 30.00, 15.0}; // reaction rate at T&p_ref

const double AntoinePars[3] = {5.40221, 1838.675, -31.737}; // Antoine pars for 273 <= T <= 303 K

const ShDataEntry ShData_N2[1] = {
		ShDataEntry(298.0, 6000.0,26.092, 8.218801, -1.976141, 0.159274, 0.044434, -7.98923, 221.02, 0.0)};
const ShDataEntry ShData_O2[1] = {
		ShDataEntry(298.0, 6000.0, 29.659, 6.137261, -1.186521, 0.09578, -0.219663, -9.861391, 237.948, 0.0)};
const ShDataEntry ShData_H2O[2] = {
		ShDataEntry(500.0, 1700.0, 30.092, 6.832514, 6.793435, -2.53448, 0.082139, -250.881, 223.3967, -241826.4),
		ShDataEntry(1700.0, 6000.0, 41.96426, 8.622053, -1.49978, 0.098119, -11.15764, -272.1797, 219.7809, -241826.4)};
const ShDataEntry ShData_CO2[2] = {
		ShDataEntry(298.0, 1200.0, 24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431, -393522.4),
		ShDataEntry(1200.0, 6000.0, 58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, 263.6125, -393522.4)};
const ShDataEntry ShData_CO[2] = {
		ShDataEntry(298.0, 1300.0, 25.56759, 6.09613, 4.054656, -2.671301, 0.131021, -118.0089, 227.3665, -110527.1),
		ShDataEntry(1300.0, 6000.0, 35.1507, 1.300095, -0.205921, 0.01355, -3.28278, -127.8375, 231.712, -110527.1)};
const ShDataEntry ShData_Fuel[1] = {
		ShDataEntry(298.0, 6000.0, 351.455, 279.288, 0.0, 0.0, 0.0, 0.0, 0.0, -249700.0)}; // Decane; only Hf for C10H22, others for Octane
		//ShDataEntry(298.0, 6000.0, 351.455, 279.288, 0.0, 0.0, 0.0, 0.0, 0.0, -208700.0)}; //Octane


const ShData ShDataDB[defs::Fuel+1] = {ShData(1, ShData_N2), ShData(1, ShData_O2),
		ShData(2, ShData_H2O), ShData(2, ShData_CO2), ShData(2, ShData_CO),
		ShData(), ShData(), ShData(1, ShData_Fuel)
};

#endif // DEFINITIONS_H
