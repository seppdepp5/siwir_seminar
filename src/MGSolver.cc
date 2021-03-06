#include <cmath>
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <omp.h>

#include "MGSolver.hh"
#include "Array.hh"
#include "Smoother.hh"
#include "Stencil.hh"

#ifndef PI
#define PI (M_PI)
#endif

#define PRINT_ERROR (1)
#define PRINT_RESIDUAL (1)
#define LERROR 9.18e-5

	MGSolver::MGSolver ( int levels, Smoother & smoother )
	: levels_ (levels)
	, v_grids_(levels, NULL)
	, r_grids_(levels, NULL)
	, tmp_grids_(levels, NULL)
	, h_intervals_(levels, 0)
	  , smoother_(smoother)
{

	//v_grids_(levels, 0);
	//r_grids_(levels, 0);

	// initialize solution
	solution_ = new Array(std::pow(2, levels) + 1, std::pow(2, levels) + 1);

	// allocate memory for v and r
	for (int i = 1; i <= levels; i++)
	{

		v_grids_[i-1] = new Array ( std::pow(2, i) + 1, std::pow(2, i) + 1);
		r_grids_[i-1] = new Array ( std::pow(2, i) + 1, std::pow(2, i) + 1);
		tmp_grids_[i-1] = new Array ( std::pow(2, i) + 1, std::pow(2, i) + 1);
		h_intervals_[i-1] = 1.0 / std::pow(2, i);


#if 0
		std::cout << h_intervals_[i-1] << std::endl;
		v_grids_[i-1]->print();
#endif

	}	

}

MGSolver::~MGSolver()
{
	delete solution_;
}

void MGSolver::initialize_assignment_01 ()
{

	// initialize bc on the unknown array boundary
	// bc is sin(x*pi) sinh(y*pi) which is only != 0 at the top boundary (y = 1)

	Array * finest_grid = v_grids_.back();
	real h		= h_intervals_.back()*2;
	int xleft	= (finest_grid->getSize(DIM_1D)) * (-0.5);
	int xright 	= (finest_grid->getSize(DIM_1D)) * 0.5;
	int ydown 	= (finest_grid->getSize(DIM_2D)) * (-0.5);
	int yup 	= (finest_grid->getSize(DIM_2D)) * 0.5;
	int xsize 	= (finest_grid->getSize(DIM_1D)) * 0.5;
	int ysize 	= (finest_grid->getSize(DIM_2D)) * 0.5;

	//bottom and upper
#pragma omp parallel for
	for (int col = -1;col >= xleft;col--)
	{
		finest_grid->operator()(col + xsize, 0) = 				-sqrt(sqrt(1 + col*h*col*h)) * sin(0.5*atan2(-1,(col*h)));
		finest_grid->operator()(col + xsize, finest_grid->getSize(DIM_2D)-1) =  sqrt(sqrt(1 + col*h*col*h)) * sin(0.5*atan2(1,(col*h)));

	}
#pragma omp parallel for
	for (int col = 1;col <= xright;col++)
	{
		finest_grid->operator()(col + xsize, 0) = 				-sqrt(sqrt(1 + col*h*col*h)) * sin(0.5*atan2(-1,(col*h)));
		finest_grid->operator()(col + xsize, finest_grid->getSize(DIM_2D)-1) =  sqrt(sqrt(1 + col*h*col*h)) * sin(0.5*atan2(1,(col*h)));
	}

	//left and right
#pragma omp parallel for
	for (int row = 0;row >= ydown;row--)
	{
		finest_grid->operator()(0,row + ysize) = 				-sqrt(sqrt(row*row*h*h + 1)) * sin(0.5*atan2(row*h,-1));
		finest_grid->operator()(finest_grid->getSize(DIM_1D)-1,row + ysize) =   -sqrt(sqrt(row*row*h*h + 1)) * sin(0.5*atan2(row*h,1));
	}

#pragma omp parallel for
	for (int row = 0;row <= yup;row++)
	{
		finest_grid->operator()(0,row + ysize) = 				sqrt(sqrt(row*row*h*h + 1)) * sin(0.5*atan2(row*h,-1));
		finest_grid->operator()(finest_grid->getSize(DIM_1D)-1,row + ysize) =   sqrt(sqrt(row*row*h*h + 1)) * sin(0.5*atan2(row*h,1));
	}

#pragma omp parallel for
	for(int col = 1; col <= xright; col++)
	{
		finest_grid->operator()(col + xsize, ysize) = 0.0;
	}
	finest_grid->operator()(xsize, ysize) = 0.0;
	finest_grid->operator()(xsize, finest_grid->getSize(DIM_2D)-1) = 1/sqrt(2);
	finest_grid->operator()(xsize, 0) = 1/sqrt(2);
/*
	// initialize solution
	// initialize solution
	int solleft	= (solution_->getSize(DIM_1D)-1) * (-0.5);
	int solright	= (solution_->getSize(DIM_1D)-1) * 0.5;
	int soldown	= (solution_->getSize(DIM_2D)-1) * (-0.5);
	int solup 	= (solution_->getSize(DIM_2D)-1) * 0.5;
	int xsol 	= (solution_->getSize(DIM_1D)-1) * 0.5;
	int ysol 	= (solution_->getSize(DIM_2D)-1) * 0.5;

	//solution  im array
	for (int row = soldown; row <= solup; row++)
	{
		for (int col = solleft; col <= solright; col++)
		{
			if(col == 0 && row == 0){ solution_->operator()(col+xsize,row+ysize) = 0.0; continue;}
			solution_->operator()(col + xsol, row + ysol) = sqrt(sqrt(row*h2*row*h2 + col*h2*col*h2)) * sin(0.5*atan2(row*h2,col*h2));//(1 - (col*h2/sqrt(col*h2*col*h2 + row*h2*row*h2))));

			//sin(PI * (real) col * h) * sinh(PI * (real) row * h);	
		}
	}
*/
}

void MGSolver::v_cycle( int pre_smooth, int post_smooth, int times)
{

	//	real error = error_L2 ( * v_grids_.back(), * solution_,	h_intervals_.back());
	for(int i = 1; i <= 9; i++)
	{
		v_cycle_pvt (pre_smooth, post_smooth, levels_);
	}

#if PRINT_RESIDUAL
//		real residual = residual_2d ( * v_grids_.back(), * r_grids_.back(), h_intervals_.back());
//		std::cout << "Residual (cylcle no " << i + 1 << "):  " << residual << std::endl;
#endif
//		real error = error_L2 ( * v_grids_.back(), * solution_, h_intervals_.back());
		//		i++;
#if PRINT_ERROR
//		std::cout << "Error: "  << error << std::endl;
#endif

}

// solve when level == 1
void MGSolver::v_cycle_pvt ( int pre_smooth, int post_smooth, int level)
{

	(void) pre_smooth;
	(void) post_smooth;
	(void) level;

	if (level == 1)
	{
		// SOLVE
		v_grids_[level-1]->operator()(0, 0) = r_grids_[level-1]->operator()(0, 0) * h_intervals_[level-1] * h_intervals_[level-1] * 0.25;
		return;
	}

	v_grids_[level-2]->fill(0.0);
	tmp_grids_[level-1]->fill(0.0);

	// 1. perform pre_smooth gauss seidel iterations on v
	smoother_.smooth_red_black_gauss_seidel_2d( *v_grids_[level-1], *r_grids_[level-1], pre_smooth, h_intervals_[level-1], level == (int)v_grids_.size());	// true if finest grid

#if 1

	// 2. calculate coarser right hand side
	compose_right_hand_side ( *v_grids_[level-1], *r_grids_[level-1], *r_grids_[level-2], level, h_intervals_[level-1]);

	// 3.recursive call
	v_cycle_pvt( pre_smooth, post_smooth, level-1 );

	// 4. correction
	error_correction( *v_grids_[level-1], *v_grids_[level-2], level);

	// 5. post smoothing
	smoother_.smooth_red_black_gauss_seidel_2d( *v_grids_[level-1], *r_grids_[level-1], post_smooth, h_intervals_[level-1], level == (int)v_grids_.size());	// true if finest grid
#endif

}

void MGSolver::error_correction(Array &u, Array &e_2h, int current_level)
{
	Array &e_h = *tmp_grids_[current_level-1];
	e_h.fill(0.0);

	int width  = e_2h.getSize(DIM_1D);
	int height = e_2h.getSize(DIM_2D);

	Stencil rest(0.25, 0.5, 0.25, 0.5, 1.0, 0.5, 0.25, 0.5, 0.25);
	
	// calculate I * e_2h (error to finer grid)
#pragma omp parallel for schedule(static)
	for (int j = 1; j < height-1; j++)
	{
		for (int i = 1; i < width-1; i++)
		{   
			if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
			int mid_i = 2 * i;	
			int mid_j = 2 * j;	

			e_h(mid_i - 1, mid_j + 1) += rest.getw1() * e_2h(i, j);
			e_h(mid_i    , mid_j + 1) += rest.getw2() * e_2h(i, j);
			e_h(mid_i + 1, mid_j + 1) += rest.getw3() * e_2h(i, j);
			e_h(mid_i - 1, mid_j    ) += rest.getw4() * e_2h(i, j);
			e_h(mid_i    , mid_j    ) += rest.getw5() * e_2h(i, j);
			e_h(mid_i + 1, mid_j    ) += rest.getw6() * e_2h(i, j);
			e_h(mid_i - 1, mid_j - 1) += rest.getw7() * e_2h(i, j);
			e_h(mid_i    , mid_j - 1) += rest.getw8() * e_2h(i, j);
			e_h(mid_i + 1, mid_j - 1) += rest.getw9() * e_2h(i, j);

		}
	}

	// add to the solution

	width  = u.getSize(DIM_1D);
	height = u.getSize(DIM_2D);

#pragma omp parallel for schedule(static)
	for (int j = 1; j < height-1; j++) {
		for (int i = 1; i < width-1; i++)
		{   
			if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
			u(i, j) += e_h(i, j);
		}
	}

}

void MGSolver::compose_right_hand_side( Array &u, Array &f, Array &r_2h, int current_level, real h)
{

	// store f - Au in temporary grid storage
	Array & res = * tmp_grids_[current_level-1];

	real h_2_inv = 1.0 / (h*h);

	int width  = u.getSize(DIM_1D);
	int height = u.getSize(DIM_2D);

	// calculate f - Au
#pragma omp parallel for schedule(static)
	for (int j = 1; j < height-1; j++)
	{
		for (int i = 1; i < width-1; i++)
		{   
			if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
			res(i, j) = f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j) );
		}   
	}

	Stencil rest(0.0625, 0.125, 0.0625, 0.125, 0.25, 0.125, 0.0625, 0.125, 0.0625);

	restrict_2d(rest, res, r_2h);

}

void MGSolver::restrict_2d (Stencil &rest, Array &u, Array &u_2h)
{

	//	rest.print();

	int width  = u_2h.getSize(DIM_1D);
	int height = u_2h.getSize(DIM_2D);

	// restrict to coarser domain
#pragma omp parallel for schedule(static)
	for(int j = 1; j < height-1; j++)
	{
		for(int i = 1; i < width-1; i++)
		{   
			if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
			int mid_i = 2*i;
			int mid_j = 2*j;

			u_2h(i, j) = rest.getw1() * u(mid_i - 1, mid_j + 1) +
				rest.getw1() * u(mid_i    , mid_j + 1) +
				rest.getw3() * u(mid_i + 1, mid_j + 1) +
				rest.getw4() * u(mid_i - 1, mid_j    ) +
				rest.getw5() * u(mid_i    , mid_j    ) +
				rest.getw6() * u(mid_i + 1, mid_j    ) +
				rest.getw7() * u(mid_i - 1, mid_j - 1) +
				rest.getw8() * u(mid_i    , mid_j - 1) +
				rest.getw9() * u(mid_i + 1, mid_j - 1);
		}
	}
}



real MGSolver::residual_2d( Array &u, Array &f, real h)
{

	real sum = 0.0;	
	real h_2_inv = 1.0 / (h*h);

	int width  = u.getSize(DIM_1D);
	int height = u.getSize(DIM_2D);

	// add up squares of the entries of the residual
#pragma omp parallel for schedule(static) reduction(+:sum)
	for (int j = 1; j < height-1; j++) {
		for (int i = 1; i < width-1; i++)
		{   
			if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
			sum +=  (f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j))) *
				(f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j)));
		}   
	}


	return sqrt(sum / (real) ((width-2) * (height-2)));
}

real MGSolver::error_L2( Array &approximation, Array &solution, real h)
{
	real sum = 0.0;	

	int width  = approximation.getSize(DIM_1D);
	int height = approximation.getSize(DIM_2D);

	// add up squares of the entries of the residual
#pragma omp parallel for schedule(static) reduction(+:sum)
	for (int j = 1; j < height-1; j++) {
		for (int i = 1; i < width-1; i++)
		{   
			if(j == (height-1)*0.5 && i >= (width-1)*0.5) continue;
			sum +=  (approximation(i,j) - solution(i,j)) *
				(approximation(i,j) - solution(i,j));
		}   
	}

	return sqrt(sum / (real) ((width-2) * (height-2)));


}

int MGSolver::saveToFile(std::string filename) const
{
	Array *u = v_grids_.back();
	//	Array *u = solution_;
	//	std::cout << "width: " << u->getSize(DIM_1D) << std::endl;

	std::ofstream gnuFile(filename);
	if (gnuFile.is_open())
	{
	//	gnuFile << "# x y u(x,y)" << "\n";
		for (int j = -0.5*(u->getSize(DIM_2D) - 1); j < 0.5*(u->getSize(DIM_2D)); j++)
		{
			for (int i = -0.5*(u->getSize(DIM_1D) - 1); i < 0.5*(u->getSize(DIM_1D)) ; i++)
			{
				gnuFile << (double) 2*i/(u->getSize(DIM_1D)-1) << " " << (double) 2*j/(u->getSize(DIM_2D)-1) << " " << u->operator()(i + 0.5*(u->getSize(DIM_1D) - 1) , j + 0.5*(u->getSize(DIM_2D) - 1) ) << "\n";
			}
//			gnuFile << "\n";
		}
		gnuFile.close();
		return 0;
	}
	else
	{
		return 1;
	}
}
