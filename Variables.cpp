#include "Variables.h"

consVariables calcConserve(const natVariables& natvar, const Gas& gas)
{
	consVariables out;
	out.rho = natvar.rho;
	out.rhoU = natvar.rho * natvar.u;
	double Cv = gas.Cp / gas.gamma;
	double E = Cv * calcTemperature(natvar, gas) + natvar.u * natvar.u / 2.0;
	out.rhoE = natvar.rho * E;
	return out;
}

natVariables calcNatural(const consVariables& conserv, const Gas& gas)
{
	natVariables out;
	out.rho = conserv.rho;
	out.u = conserv.rhoU / conserv.rho;
	double E = conserv.rhoE / conserv.rho,
		Cv = gas.Cp / gas.gamma,
		Rm = gas.R / gas.molMass;
	double T = (E - out.u * out.u / 2.0) / Cv;
	out.P = out.rho * Rm * T;
	return out;
}

Fluxes calcFluxes(const natVariables& natvar, const Gas& gas)
{
	Fluxes out;
	out.rhoU = natvar.rho * natvar.u;
	out.rhoUU_P = out.rhoU * natvar.u + natvar.P;
	double H = gas.Cp * calcTemperature(natvar, gas) + natvar.u * natvar.u / 2.0;
	out.rhoUH = out.rhoU * H;
	return out;
}

double calcTemperature(const natVariables& natvar, const Gas& gas)
{
	double Rm = gas.R / gas.molMass;
	return natvar.P / Rm / natvar.rho;
}

double calcTemperature(const double p, const double rho, const Gas& gas)
{
	double Rm = gas.R / gas.molMass;
	return p / Rm / rho;
}

double calcEnergy(const natVariables& natvar, const Gas& gas)
{
	return natvar.P / (gas.gamma - 1.0) / natvar.rho;
}