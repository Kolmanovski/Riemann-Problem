#pragma once
#ifndef FILE_H
#define FILE_H

#include <fstream>
#include <string>
#include <vector>
#include "Inputdata.h"
#include "Variables.h"

class File
{
public:
    static inputData readInputData(const std::string& file);
    static void outputData(const std::vector<double>& x, const std::vector<natVariables>& y, const Gas& gas, const std::string& fileName);
};

#endif // FILE_H
