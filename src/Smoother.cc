#include "Smoother.hh"
#include "Types.hh"
#include "Array.hh"
#include "MGSolver.hh"
#include <omp.h>

void Smoother::smooth_red_black_gauss_seidel_2d ( Array & u,    // modify this array
		Array & f,    // rhs
		int times,     // number of sweeps
		real h,       // spacing         
		bool finest_grid
		)
{

	(void) finest_grid;

	int width  = u.getSize(DIM_1D);
	int height = u.getSize(DIM_2D);

	real h_2     = h * h;
	real h_2_inv = 1.0 / h_2;
	real factor  = h_2 * 0.25;

	for (int iter = 0; iter < times; iter++)
	{
		// red points
#pragma omp parallel for
		for (int j = 1; j < height-1; j++)
		{
			for (int i = 1; i < width-1; i+=2)
			{
			//	if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
				// inner domain
				// i+j gerade
				u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
			}
			j++;
			for (int i = 2; i < width-1; i+=2)
			{
				if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
				// inner domain
				// i+j gerade
				u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
			}

		}

#pragma omp parallel for
		// black points
		for (int j = 1; j < height-1; j++)
		{
			for (int i = 2; i < width-1; i+=2)
			{
			//	if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
				// inner domain
				// i+j ungerade
				u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
			}
			j++;
			for (int i = 1; i < width-1; i+=2)
			{
				if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
				// inner domain
				// i+j ungerade
				u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
			}

		}
	}
}
