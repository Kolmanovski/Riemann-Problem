#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include "InputData.h"
#include "Variables.h"
#include <vector>
#include <string>
#include <fstream>

class Solver
{
public:
    explicit Solver(const inputData& inData);
    void createMesh(); // функци€ создани€ сетки
    void initialize();
    void solve();
    void calcBC();
    std::vector<double> getMesh() const; // метод, который что-то получает от класса Solver
    std::vector<natVariables> getSolution() const;
    void monPoint(std::ofstream &filename, double time, const natVariables &v) const;

private:
    inputData m_inData; // inputData - структура
    std::vector<double> m_mesh;
    std::vector<double> m_face; // дл€ расчета площади грани
    std::vector<double> m_pos_face; // дл€ определени€ позиции грани
    std::vector<double> m_volume; // дл€ расчета объема (площади) €чейки
    std::vector<natVariables> m_var;
   
    double m_dh = 0;
    int index = 0;
};

#endif // !SOLVER_H

//std::vector<Variables> getInitialize() const;