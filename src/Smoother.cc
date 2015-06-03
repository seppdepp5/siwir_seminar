
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

//	if(finest_grid) omp_set_num_threads(32);

	int width  = u.getSize(DIM_1D);
	int height = u.getSize(DIM_2D);

//TODO
//da ist doch die if oder??

	real h_2     = h * h;
	real h_2_inv = 1.0 / h_2;
	real factor  = h_2 * 0.25;

	for (int iter = 0; iter < times; iter++)
	{
		// red points
#pragma omp parallel for
		for (int j = 1; j < height-1; j++)
		{
			//std::cout << omp_get_num_threads() << std::endl;
			for (int i = 1; i < width-1; i++)
			{
				// inner domain
				// i+j gerade
				if( ((i + j) % 2) == 0){
					if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
				}
			}
/*
			j++;
			for (int i = 2; i < width-2; i+=2)
			{
				// inner domain
				// i+j gerade
				//if( ((i + j) % 2) == 0){
					if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
				//}
			}
*/
		}

		// black points
#pragma omp parallel for
		for (int j = 1; j < height-1; j++)
		{
			for (int i = 1; i < width-1; i++)
			{
				// inner domain
				// i+j ungerade
				if( ((i + j) % 2) == 1){
					if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
				}
			}
/*			
			j++;
			for (int i = 1; i < width-1; i+=2)
			{
				// inner domain
				// i+j ungerade
				//if( ((i + j) % 2) == 1){
				if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
				//}
			}
*/

		}
	} 
}
