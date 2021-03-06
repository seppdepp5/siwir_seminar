#include <iostream>
#include <sys/time.h>
#include "MGSolver.hh"
#include "Smoother.hh"
#include "Timer.hh"
#include "Stencil.hh"


int main(int argc, char **args)
{
	if(argc != 2){
		std::cout << "Usage: ./mgsolve <number_of_levels> <number_of_V-cycles>" << std::endl;
		return 0;
	}

	int l;
	int n;
	n = 10;

	std::istringstream iss(args[1]);
	if(!(iss >> l)){
		std::cerr << "Could not parse number of level argument: " << args[1] << std::endl;
		return 1;
	}
	iss.str("");
	iss.clear();


	Smoother smoother;
	MGSolver solver(l, smoother);

	solver.initialize_assignment_01();
//	std::cout << "writing init.dat" << std::endl;
	solver.saveToFile("init.dat");

//	gettimeofday(&t0, NULL);

	std::cout<<"Your Alias: "<< "David&Thomas" <<std::endl;
	struct timeval t0, t;
	gettimeofday(&t0, NULL);

	solver.v_cycle(2, 1, n);

	gettimeofday(&t, NULL);
	std::cout << "Wall clock time of MG execution: " << ((int64_t) (t.tv_sec - t0.tv_sec) * (int64_t)1000000 + (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3 << " ms" << std::endl;

//	std::cout << "writing solution.dat" << std::endl;
	solver.saveToFile("solution.dat");

	return 0;
}

