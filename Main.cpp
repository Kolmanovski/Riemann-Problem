#include <iostream>
#include "Solver.h"
#include "InputData.h"
#include "File.h"

int main(int argc, char* argv[])
{
	inputData inData = File::readInputData("C:\\STUDY\\kolesnik-3sem\\program_vs\\test\\praktika\\input.txt");
	Solver solver(inData);
	solver.createMesh();
	solver.initialize();
	solver.solve();
	File::outputData(solver.getMesh(), solver.getSolution(), inData.gasData,
		"C:\\STUDY\\kolesnik-3sem\\program_vs\\test\\praktika\\1D_CFL\\1D_Godunov_1_od_CFL=0.85.plt");
	//File::monitorPoint(,"C:\\STUDY\\kolesnik-3sem\\program_vs\\test\\praktika\\result_point.plt");
	return 0;
}