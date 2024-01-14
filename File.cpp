#include <iostream>
#include "File.h"
#include <map>

inputData File::readInputData(const std::string& fileName) // –аботает!
{
	std::ifstream file(fileName);
	inputData inData; // —оздаем объект (переменную) структуры inputData

	if (!file.is_open())
	{
		std::cout << "FILE IS NOT OPENING\n";
	}
	else
	{
		std::cout << "FILE IS OPENED!\n";
	}

	std::map <std::string, std::string> inMap;
	std::string line;

	while (!file.eof())
	{
		std::getline(file, line);
		if (line.empty()) // если обнаружил пустую строку, все равно продолжай считывать данные
			continue;
		int sepIndx = line.find("\t");
		std::string key = line.substr(0, sepIndx);
		key.erase(std::remove_if(key.begin(), key.end(), ::isspace), key.end()); // key - это NX, например
		std::string val = line.substr(sepIndx);
		val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end()); // val - это значение NX, т.е., 10
		if (key.empty() || val.empty())
			continue;
		inMap[key] = val;
	}
    file.close();

    inData.NX = std::stoi(inMap["NX"]);
    inData.L = std::stod(inMap["L"]);
	inData.CFL = std::stod(inMap["CFL"]);
	inData.calcTime = std::stod(inMap["Time"]);
	inData.rhoL = std::stod(inMap["rhoL"]);
    inData.uL = std::stod(inMap["uL"]);
	inData.PL = std::stod(inMap["PL"]);
	inData.rhoR = std::stod(inMap["rhoR"]);
    inData.uR = std::stod(inMap["uR"]);
	inData.PR = std::stod(inMap["PR"]);
	inData.gasData.Cp = std::stod(inMap["Cp"]);
	inData.gasData.molMass = std::stod(inMap["molMass"]);
	inData.gasData.gamma = std::stod(inMap["gamma"]);
	inData.scheme = std::stoi(inMap["scheme"]);
	inData.TVD = std::stoi(inMap["TVD"]);
    return inData;
}

void File::outputData(const std::vector<double>& x, const std::vector<natVariables>& natvar, const Gas& gas, const std::string& fileName)
{
    std::ofstream file(fileName);
	file << "Variables = \"X\", \"rho\", \"U\", \"P\", \"T\", \"E\"\n";
	for (int i = 1; i < x.size() - 1; ++i)
	{
		const natVariables& v = natvar.at(i);
		file << x.at(i) << "\t" << v.rho << "\t" << v.u << "\t" << v.P << "\t" << calcTemperature(v, gas) << "\t" << calcEnergy(v, gas) << std::endl;
	}
    file.close();
}