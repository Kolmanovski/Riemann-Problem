#pragma once
#include "Variables.h"
#include "InputData.h"
class Schemes
{
public:
    virtual Fluxes CalculateFlux(const natVariables& varLeft, const natVariables& varRight) = 0;
    Gas gas;
    Schemes(Gas gas) : gas(gas) {}
};

class HLL :public Schemes
{
private:

    double SL{ 0 }, SR{ 0 };

public:

    HLL(Gas gas) : Schemes(gas) {}

    Fluxes CalculateFlux(const natVariables& varLeft, const natVariables& varRight) override
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
            //if (SL < 0 && SR > 0)
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
        CR = sqrt(AK2 * right.P / right.rho);

        SL = std::min({ left.u - CL, right.u - CR });
        SR = std::max({ left.u + CL, right.u + CR }); // Изменил min на max
    }
};

class Godunov :public Schemes
{
private:
    double PBIG{ 0 }, UBIG{ 0 }, UDOT{ 0 };
    double RBIG1{ 0 }, RBIG2{ 0 };
    double DL1{ 0 }, DL2{ 0 };
    double DP1{ 0 }, DP2{ 0 };

    natVariables varFace;

public:

    Godunov(Gas gas) : Schemes(gas) {}

    Fluxes CalculateFlux(const natVariables& varLeft, const natVariables& varRight) override
    {
        RASPAD(varLeft, varRight);
        POTOK(varLeft, varRight);

        Fluxes flux;
        flux.rhoU = varFace.rho * varFace.u;
        flux.rhoUU_P = flux.rhoU * varFace.u + varFace.P;
        double Rm = gas.R / gas.molMass;
        double H = gas.Cp * varFace.P / Rm / varFace.rho + (varFace.u * varFace.u) / 2.0;
        flux.rhoUH = flux.rhoU * H;
        return flux;
        std::printf, flux;
    }
private:

    void POTOK(const natVariables& left, const natVariables& right)
    {
        double AK1, AK2;
        double CZJM, CZJP, CZT;
        AK1 = gas.gamma;
        AK2 = gas.gamma;
        varFace.u = UBIG;
        varFace.P = PBIG;

        if (UBIG - UDOT > 0.)
        {
            if ((DL1 - DL2) == 0.0)
            {
                if (DL1 - UDOT >= 0.)
                {
                    varFace.rho = left.rho;
                    varFace.u = left.u;
                    varFace.P = left.P;
                }
                else
                {
                    varFace.rho = RBIG1;
                }
            }
            else
            {
                if (DL2 - UDOT <= 0.)
                    varFace.rho = RBIG1;
                else
                {
                    if (DL1 - UDOT >= 0.)
                    {
                        varFace.rho = left.rho;
                        varFace.u = left.u;
                        varFace.P = left.P;
                    }
                    else
                    {
                        CZJM = sqrt(AK1 * left.P / left.rho);
                        CZT = ((AK1 - 1.) * (left.u - UDOT) + 2. * CZJM) / (AK1 + 1.);
                        varFace.u = CZT + UDOT;
                        varFace.P = left.P * pow(CZT / CZJM, 2. * AK1 / (AK1 - 1.));
                        varFace.rho = AK1 * varFace.P / (CZT * CZT);
                    }
                }
            }
        }
        else
        {
            if (DP1 - DP2 == 0.0)
            {
                if (DP1 - UDOT <= 0.)
                {
                    varFace.rho = right.rho;
                    varFace.u = right.u;
                    varFace.P = right.P;
                }
                else
                    varFace.rho = RBIG2;
            }
            else
            {
                if (DP1 - UDOT >= 0.)
                    varFace.rho = RBIG2;
                else
                {
                    if (DP2 - UDOT <= 0.)
                    {
                        varFace.rho = right.rho;
                        varFace.u = right.u;
                        varFace.P = right.P;
                    }
                    else
                    {
                        CZJP = sqrt(AK2 * right.P / right.rho);
                        CZT = (-(AK2 - 1.) * (right.u - UDOT) + 2. * CZJP) / (AK2 + 1.);
                        varFace.u = -CZT + UDOT;
                        varFace.P = right.P * pow(CZT / CZJP, 2. * AK2 / (AK2 - 1.));
                        varFace.rho = AK2 * varFace.P / (CZT * CZT);
                    }
                }
            }
        }
        return;
    }


    double FGOD(const double P, const double R, const double AK)
    {
        double PI = PBIG / P;
        double ST = 0.5 * (AK - 1.0) / AK;

        if (PI >= 1.0)
            return sqrt(P / R) * (PI - 1.0) / sqrt(0.5 * (AK + 1.0) * PI + 0.5 * (AK - 1.0));
        else
            return 2.0 / (AK - 1.0) * sqrt(AK * P / R) * (pow(PI, ST) - 1.0);

    }

    double DIVFGOD(const double P, const double R, const double AK)
    {
        double PI = PBIG / P;
        double ST = 0.5 * (AK - 1.0) / AK;

        if (PI >= 1.0)

            return ((AK + 1.0) * PI + (3.0 * AK - 1.0)) / (4.0 * AK * R * sqrt(AK * P / R) *
                sqrt(pow(0.5 * (AK + 1.0) / AK * PI + 0.5 * (AK - 1.0) / AK, 3)));
        else

            return 1.0 / AK / PBIG * sqrt(AK * P / R) * pow(PI, ST);
    }

    void RASPAD(const natVariables& left, const natVariables& right)
    {
        bool FLAG;
        double AK1, AK2, C1, C2, ST1, ST2, A1, A2, C2Z, C1Z, PI, EPSP, EPS;
        double PBIGN;

        EPSP = 1.0E-8;
        EPS = 1.0E-8;

        AK1 = gas.gamma;
        AK2 = gas.gamma;
        C1 = sqrt(AK1 * left.P / left.rho);
        C2 = sqrt(AK2 * right.P / right.rho);

        if ((left.u - right.u + 2.0 * (C1 / (AK1 - 1.0) + C2 / (AK2 - 1.0))) < 0.0)
        {
            PBIG = 1.0e-05;
            RBIG1 = 0.0;
            RBIG2 = 0.0;
            return;
        }

        FLAG = false;
        PBIG = (left.P * right.rho * C2 + right.P * left.rho * C1 + (left.u - right.u) * left.rho * C1 * right.rho * C2) /
            (left.rho * C1 + right.rho * C2);

        if (PBIG < 0.0)
            PBIG = 1.0E-05;

    CalculatePBIG:
        PBIGN = PBIG - (FGOD(left.P, left.rho, AK1) + FGOD(right.P, right.rho, AK2) - (left.u - right.u)) /
            (DIVFGOD(left.P, left.rho, AK1) + DIVFGOD(right.P, right.rho, AK2));

        if (PBIGN < 0.0)
        {
            std::cout << PBIG << "\t" << FLAG << "\n";
            PBIG = 1.0E-05;
            RBIG1 = 0.0;
            RBIG2 = 0.0;
            if (FLAG) return;
            FLAG = true;
            goto CalculatePBIG;
        }
        if (abs(PBIGN / PBIG - 1.0) > EPS)
        {
            PBIG = PBIGN;
            goto CalculatePBIG;
        }

        ST1 = 0.5 * (AK1 - 1.0) / AK1;
        ST2 = 0.5 * (AK2 - 1.0) / AK2;
        if (PBIG >= (left.P - EPS))
            A1 = sqrt(left.rho * (0.5 * (AK1 + 1.0) * PBIG + 0.5 * (AK1 - 1.0) * left.P));
        else
        {
            PI = PBIG / left.P;
            A1 = 0.5 * (AK1 - 1.0) / AK1 * left.rho * C1 * (1.0 - PI) / (1.0 - pow(PI, ST1));
        }


        if (PBIG >= (right.P - EPS))
            A2 = sqrt(right.rho * (0.5 * (AK2 + 1.0) * PBIG + 0.5 * (AK2 - 1.0) * right.P));
        else
        {
            PI = PBIG / right.P;
            A2 = 0.5 * (AK2 - 1.0) / AK2 * right.rho * C2 * (1.0 - PI) / (1.0 - pow(PI, ST2));
        }

        UBIG = (A1 * left.u + A2 * right.u + left.P - right.P) / (A1 + A2);

        if ((left.P - EPS) < PBIG)
        {
            RBIG1 = left.rho * A1 / (A1 - left.rho * (left.u - UBIG));
            DL1 = left.u - A1 / left.rho;
            DL2 = DL1;
        }
        else
        {
            C1Z = C1 + 0.5 * (AK1 - 1.0) * (left.u - UBIG);
            RBIG1 = AK1 * PBIG / (C1Z * C1Z);
            DL1 = left.u - C1;
            DL2 = UBIG - C1Z;
        }


        if ((right.P - EPS) < PBIG)
        {
            RBIG2 = right.rho * A2 / (A2 + right.rho * (right.u - UBIG));
            DP1 = right.u + A2 / right.rho;
            DP2 = DP1;
        }
        else
        {
            C2Z = C2 - 0.5 * (AK2 - 1.0) * (right.u - UBIG);
            RBIG2 = AK2 * PBIG / (C2Z * C2Z);
            DP2 = right.u + C2;
            DP1 = UBIG + C2Z;
        }
        return;
    }
};