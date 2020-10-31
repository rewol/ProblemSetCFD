/****************************************************/
/* Refik Alper Tuncer, AÐUSTOS2020, CFD Projesi-1	*/
/* Anderson, Chapter-10 problem						*/
/****************************************************/

#include <iostream>
#include "Field2D.h"
#include "Solver.h"
#include<algorithm>
#include<chrono>

int main()
{
	// Set timer
	auto start = std::chrono::steady_clock::now();	
	
	Field2D field(121, 121, 1e-5, 4, 1, 0);
	Field2D* f_ptr = &field;
	Solver s(f_ptr, 1e-3, 1e-3);
	s.MAC();
	field.printPrimitive();

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "Elapsed Time: " << elapsed_seconds.count() << " sn" << "\n";

	return 0;
}
