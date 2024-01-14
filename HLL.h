#pragma once
#include <iostream>
#include <math.h>
#include "Variables.h"
#include "Inputdata.h"

class HLL
{
private:

    double SL{ 0 }, SR{ 0 };

    Gas gas;

public:

    HLL(Gas gas) : gas(gas) {}

    Fluxes CalculateFlux(const natVariables& varLeft, const natVariables& varRight)
    {
        calcSpeed(varLeft, varRight);

        if (SL > 0)
        {
            return POTOK(varLeft);
        }
        else if (SR < 0)
        { 
            return POTOK(varRight);
        }
        else
        {
            return Fstar(varLeft, varRight);
        }
    }
private:

    Fluxes Fstar(const natVariables& left, const natVariables& right)
    {
        Fluxes fluxStar;
        
        fluxStar = (POTOK(left) * SR - POTOK(right) * SL + (calcConserve(right, gas) - calcConserve(left, gas)) * SL * SR) / (SR - SL);
        return fluxStar;
    }

    Fluxes POTOK(const natVariables& natvar)
    {
     
        double rhoU = natvar.rho * natvar.u;
        double rhoUU_P = rhoU * natvar.u + natvar.P;
        double Rm = gas.R / gas.molMass;
        double H = gas.Cp * natvar.P / Rm / natvar.rho + (natvar.u * natvar.u) / 2.0;
        double rhoUH = rhoU * H;
        return Fluxes(rhoU, rhoUU_P, rhoUH);
    }
    void calcSpeed(const natVariables& left, const natVariables& right)
    {
        double AK1, AK2;
        double CL, CR;
        AK1 = gas.gamma;
        AK2 = gas.gamma;

        CL = sqrt(AK1 * left.P / left.rho);
        CR = sqrt(AK1 * right.P / right.rho);

        SL = std::min({ left.u - CL, right.u - CR });
        SR = std::min({ left.u + CL, right.u + CR });
    }
};
