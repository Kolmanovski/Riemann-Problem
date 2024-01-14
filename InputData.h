#pragma once
#ifndef INPUTDATA_H
#define INPUTDATA_H

struct Gas
{
	double Cp = 0;
	double molMass = 0;
	double gamma = 0;
	const double R = 8.3144262;
};

struct inputData
{
	int NX = 0;
	double L = 0;
	double CFL = 0;			// } Свойства(или поля) структуры inputData
	double calcTime = 0;
	double rhoL, uL, PL = 0;
	double rhoR, uR, PR = 0;
	Gas gasData;
	int scheme = 0;
	int TVD = 0;
};

#endif // !INPUTDATA_H

