
#include "Smoother.hh"
#include "Types.hh"
#include "Array.hh"
#include "MGSolver.hh"
#include <omp.h>
#include <assert.h>

void Smoother::smooth_red_black_gauss_seidel_2d ( Array & u,    // modify this array
                                                  Array & f,    // rhs
                                                  int times,     // number of sweeps
                                                  real h,        // spacing         
                                                  bool finest_grid
						)
{

//if(finest_grid) omp_set_num_threads(32);
	
	int width  = u.getSize(DIM_1D);
	int height = u.getSize(DIM_2D);


	real h_2     = h * h;
	real h_2_inv = 1.0 / h_2;
	real factor  = h_2 * 0.25;
	real w = 1.141421356237;

	for (int iter = 0; iter < times; iter++)
	{
		// red points
#pragma omp parallel for
		for (int j = 1; j < height-1; j++)
		{
	//		std::cout << omp_get_num_threads() << std::endl;
			for (int i = 1; i < width-1; i+=2)
			{
				u(i,j) = (1-w)*u(i,j) + w*factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
			}

			j++;
			if(j == height-1) continue;
			for (int i = 2; i < width-2; i+=2)
			{
				if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
				u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
			}

		}

		// black points
#pragma omp parallel for
		for (int j = 1; j < height-1; j++)
		{
			for (int i = 2; i < width-2; i+=2)
			{
				u(i,j) = (1-w)*u(i,j) + w*factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
			}
			
			j++;
			if(j == height-1) continue;
			for (int i = 1; i < width-1; i+=2)
			{
				if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
				u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
			}
		}
	} 
//if(finest_grid) omp_set_num_threads(8);
}
