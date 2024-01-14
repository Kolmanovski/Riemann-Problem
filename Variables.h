#pragma once
#ifndef VARIABLES_H
#define VARIABLES_H
#include "InputData.h"

struct consVariables
{
	double rho = 0;
	double rhoU = 0;
	double rhoE = 0;

	consVariables(double rho, double rhoU, double rhoE) :rho(rho), rhoU(rhoU), rhoE(rhoE) {}

	consVariables() {}

	consVariables operator*(const double v) const
	{
		return consVariables(this->rho * v, this->rhoU * v, this->rhoE * v);
	}

	consVariables operator-(const consVariables& other) const
	{
		return consVariables(this->rho - other.rho, this->rhoU - other.rhoU, this->rhoE - other.rhoE);
	}

};

struct natVariables
{
	double rho = 0;
	double u = 0;
	double P = 0;

	natVariables operator-(const natVariables& other) const
	{
		return natVariables{ this->rho - other.rho, this->u - other.u, this->P - other.P };
	}
	natVariables operator/(const natVariables& other) const
	{
		return natVariables{ this->rho / other.rho, this->u / other.u, this->P / other.P };
	}
	natVariables operator+(const natVariables& other) const
	{
		return natVariables{ this->rho + other.rho, this->u + other.u, this->P + other.P };
	}
	natVariables operator*(const natVariables& other) const
	{
		return natVariables{ this->rho * other.rho, this->u * other.u, this->P * other.P };
	}
	natVariables operator*(const double other) const
	{
		return natVariables{ this->rho * other, this->u * other, this->P * other };
	}

	static natVariables getMin(const natVariables& v1, const natVariables& v2)
	{
		return natVariables{ v1.rho < v2.rho ? v1.rho : v2.rho, v1.u < v2.u ? v1.u : v2.u, v1.P < v2.P ? v1.P : v2.P };
	}
};

struct Fluxes
{
	double rhoU = 0;
	double rhoUU_P = 0;
	double rhoUH = 0;

	Fluxes(double rhoU, double rhoUU_P, double rhoUH) :rhoU(rhoU), rhoUU_P(rhoUU_P), rhoUH(rhoUH) {}

	Fluxes() {}

	
	Fluxes operator-(const Fluxes& other)
	{
		return Fluxes(this->rhoU - other.rhoU, this->rhoUU_P - other.rhoUU_P, this->rhoUH - other.rhoUH);
	}

	Fluxes operator*(const double v)
	{
		return Fluxes(this->rhoU*v, this->rhoUU_P*v, this->rhoUH*v);
	}

	Fluxes operator+(const consVariables& other)
	{
		return Fluxes(this->rhoU + other.rho, this->rhoUU_P + other.rhoU, this->rhoUH + other.rhoE);
	}

	Fluxes operator/(const double v)
	{
		return Fluxes(this->rhoU / v, this->rhoUU_P / v, this->rhoUH / v);
	}
};

consVariables calcConserve(const natVariables& natvar, const Gas& gas);

natVariables calcNatural(const consVariables& conserv, const Gas& gas);

Fluxes calcFluxes(const natVariables& natvar, const Gas& gas);

double calcTemperature(const natVariables& natvar, const Gas& gas);

double calcTemperature(const double p, const double rho, const Gas& gas);

double calcEnergy(const natVariables& natvar, const Gas& gas);

#endif // !VARIABLES_H