#include <iostream>
#include "Solver.h"
#include <algorithm>
#include "Schemes.h"
#include "TVD.h"
//#include "Godunov.h"
//#include "HLL.h"

Solver::Solver(const inputData& inData) : m_inData(inData) {} // Конструктор класса Solver

void Solver::createMesh()
{
	m_mesh.resize(m_inData.NX + 2);
	m_face.resize(m_inData.NX + 1);
	m_pos_face.resize(m_inData.NX + 1);
	m_volume.resize(m_inData.NX);
	m_var.resize(m_inData.NX + 2);

	//Лямбда-функция:
	auto Square = [](const double x) {return 1.0; }; // rectangle
	double
		x0 = 0.0,
		x1 = 0.3,
		x2 = 0.7,
		xL = m_inData.L;
	double k = 0.2;
	//auto Square = [&](const double x) 
	//{
	//	if (x >= x0 && x <= x1)
	//		return k * x + 0.2;
	//	else if (x >= x1 && x <= x2)
	//		return 0.26;
	//	else if (x >= x2 && x <= xL)
	//		return -k * (x - x2) + 0.26;
	//}; // trapeze

	m_dh = m_inData.L / m_inData.NX; // Шаг по пространству

	m_mesh[0] = -m_dh / 2.0; // Заграничное значение
	m_pos_face[0] = 0; // Позиция грани в нуле
	m_face[0] = Square(m_pos_face[0]);

	for (int i = 1; i < m_mesh.size() - 1; ++i)
	{
		m_mesh[i] = m_mesh[i - 1] + m_dh;
		m_pos_face[i] = m_pos_face[i - 1] + m_dh;
		m_face[i] = Square(m_pos_face[i]);
	}

	m_mesh[m_inData.NX + 1] = m_mesh[m_inData.NX] + m_dh;

	for (int i = 0; i < m_volume.size(); ++i)
	{
		m_volume[i] = 0.5 * (m_face[i] + m_face[i + 1]) * (m_pos_face[i + 1] - m_pos_face[i]);
	}
}

void Solver::initialize()
{
	const double pos = m_inData.L / 2.0;
	for (int i = 0; i < m_mesh.size(); ++i)
	{
		if (m_mesh.at(i) < pos)
		{
			m_var[i].rho = m_inData.rhoL;
			m_var[i].u = m_inData.uL;
			m_var[i].P = m_inData.PL;
			index = i;
		}
		else
		{
			m_var[i].rho = m_inData.rhoR;
			m_var[i].u = m_inData.uR;
			m_var[i].P = m_inData.PR;
		}
	}
	calcBC();
}

/* РЕШЕНИЕ ОСНОВНЫХ УРАВНЕНИЙ */

void Solver::solve()
{
	std::ofstream fout("C:\\STUDY\\kolesnik-3sem\\program_vs\\test\\praktika\\1D_CFL\\result_point_1D_Godunov_1_od_CFL=0.85.plt");
	fout << "Variables = \"Time\", \"rho\", \"U\", \"P\", \"T\", \"E\"\n";
	//Godunov g(m_inData.gasData);
	//HLL g(m_inData.gasData);

	Schemes* scheme; // Создаю абстрактный объект класса
	TVD tvd;
	Fluxes fL, fR;

	if (m_inData.scheme) // Если в input.txt "0" - схема Godunov, иначе "1" - схема HLL
	{
		scheme = new HLL(m_inData.gasData);
	}
	else
	{
		scheme = new Godunov(m_inData.gasData);
	}

	double T = calcTemperature(std::max({ m_inData.PL, m_inData.PR }), std::max({ m_inData.rhoL, m_inData.rhoR }), m_inData.gasData);
	double a = std::sqrt(m_inData.gasData.gamma * m_inData.gasData.R / m_inData.gasData.molMass * T);
	double dt = m_inData.CFL * m_dh / a;
	int NT = m_inData.calcTime / dt;

	std::cout << "Total time: " << m_inData.calcTime << " s" << "\t"
		<< "Time step = " << dt << " s" << "\t"
		<< "Time iterations: " << NT << "\n";

	for (int n = 0; n < NT; ++n) // Цикл по времени
	{
		std::vector<consVariables> W(m_var.size()), Wnew; // Задаю вектор U и переменную Unew, которой буду переприсваивать значения
		// На каждом шаге по времени преобразую m_var. в W
		std::transform(m_var.begin(), m_var.end(), W.begin(),
			[&](const natVariables& natvar) { return calcConserve(natvar, m_inData.gasData); });

		Wnew = W;
		for (int i = 1; i < m_var.size() - 1; ++i)
		{
			//Fluxes fR = scheme->CalculateFlux(m_var[i], m_var[i + 1]);
			//Fluxes fL = scheme->CalculateFlux(m_var[i - 1], m_var[i]);
			if (m_inData.TVD != 0 && i > 1 && i < m_var.size() - 2) // Если не "0", то TVD
			{
				natVariables right = tvd.calcVariables(m_var[i - 1], m_var[i], m_var[i + 1], TVD::right);
				natVariables left = tvd.calcVariables(m_var[i - 2], m_var[i - 1], m_var[i], TVD::left);
				fL = scheme->CalculateFlux(left, right);

				left = tvd.calcVariables(m_var[i - 1], m_var[i], m_var[i + 1], TVD::left);
				right = tvd.calcVariables(m_var[i], m_var[i + 1], m_var[i + 2], TVD::right);
				fR = scheme->CalculateFlux(left, right);
			}
			else
			{
				fR = scheme->CalculateFlux(m_var[i], m_var[i + 1]);
				fL = scheme->CalculateFlux(m_var[i - 1], m_var[i]);
			}
			double Pc = calcNatural(W[i], m_inData.gasData).P;
			Wnew[i].rho = W[i].rho - dt * (fR.rhoU * m_face[i] - fL.rhoU * m_face[i - 1]) / m_volume[i - 1];
			Wnew[i].rhoU = W[i].rhoU - dt * ((fR.rhoUU_P * m_face[i] - fL.rhoUU_P * m_face[i - 1]) - Pc * (m_face[i] - m_face[i - 1])) / m_volume[i - 1];
			Wnew[i].rhoE = W[i].rhoE - dt * (fR.rhoUH * m_face[i] - fL.rhoUH * m_face[i - 1]) / m_volume[i - 1];
		}
		// На каждом шаге по времени преобразую W в m_var.
		std::transform(Wnew.begin(), Wnew.end(), m_var.begin(),
			[&](const consVariables& conserv) { return calcNatural(conserv, m_inData.gasData); });

		calcBC(); // Граничные условия

		std::cout << dt * (n + 1) << std::endl;
		monPoint(fout, dt * (n + 1), m_var.at(index)); // Точка мониторинга
	}
}

void Solver::monPoint(std::ofstream& filename, double time, const natVariables& v) const
{
	filename << time << "\t" << v.rho << "\t" << v.u << "\t" << v.P << "\t" << calcTemperature(v, m_inData.gasData) << "\t" << calcEnergy(v, m_inData.gasData) << std::endl;
}

std::vector<double> Solver::getMesh() const
{
	return m_mesh;
}

std::vector<natVariables> Solver::getSolution() const
{
	return m_var;
}

void Solver::calcBC()
{
	std::size_t bIndx = 0, lastIndx = bIndx + 1;

	m_var[bIndx].rho = m_var[lastIndx].rho;
	m_var[bIndx].u = -m_var[lastIndx].u; // wall: m_var[bIndx].u = -m_var[lastIndx].u;
	m_var[bIndx].P = m_var[lastIndx].P;

	bIndx = m_var.size() - 1;
	lastIndx = bIndx - 1;

	m_var[bIndx].rho = m_var[lastIndx].rho;
	m_var[bIndx].u = -m_var[lastIndx].u; // wall: m_var[bIndx].u = -m_var[lastIndx].u;
	m_var[bIndx].P = m_var[lastIndx].P;
}

//std::vector<Variables> Solver::getInitialize() const // Работает!
//{
//	return m_var;
//}