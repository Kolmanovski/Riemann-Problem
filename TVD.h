#pragma once
#include "Variables.h"

class TVD
{
public:
	enum Side
	{
		left = 0,
		right = 1
	};
	natVariables calcVariables(const natVariables& natvarLeft, const natVariables& natvarCenter, const natVariables& natvarRight, const Side& side)
	{
		natVariables arg;
		natVariables temp;
		natVariables znam;
		natVariables limit;
		
		if (side == left)
		{
			temp = (natvarRight - natvarCenter) / (natvarCenter - natvarLeft);
			znam = natvarCenter - natvarLeft;
			arg = natVariables{ znam.rho == 0 ? 0 : temp.rho, znam.u == 0 ? 0 : temp.u, znam.P == 0 ? 0 : temp.P };
			arg = natVariables{ arg.rho <= 0 ? 0 : temp.rho, arg.u <= 0 ? 0 : temp.u, arg.P <= 0 ? 0 : temp.P };
			limit = (arg * arg + arg) / (arg * arg + natVariables{ 1.0,1.0,1.0 }); // psi(r)
			return natvarCenter + limit * (natvarCenter - natvarLeft) * 0.5;
		}
		else
		{
			temp = (natvarCenter - natvarLeft) / (natvarRight - natvarCenter);
			znam = natvarRight - natvarCenter;
			arg = natVariables{ znam.rho == 0 ? 0 : temp.rho, znam.u == 0 ? 0 : temp.u, znam.P == 0 ? 0 : temp.P };
			arg = natVariables{ arg.rho <= 0 ? 0 : temp.rho, arg.u <= 0 ? 0 : temp.u, arg.P <= 0 ? 0 : temp.P };
			limit = (arg * arg + arg) / (arg * arg + natVariables{ 1.0,1.0,1.0 }); // psi(r)
			return natvarCenter - limit * (natvarRight - natvarCenter) * 0.5;
		}
	} 
};
//limit = arg * 2.0 / (arg * arg + natVariables{ 1.0,1.0,1.0 }); // fi(r)
//return natvarCenter - limit * (natvarRight - natvarLeft) * 0.25;
// minmod -> natVariables::getMin(natVariables{ 1.0,1.0,1.0 }, arg);