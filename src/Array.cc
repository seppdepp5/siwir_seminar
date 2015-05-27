#include <cstdlib>
#include <iostream>
#include <cmath>
#include "Debug.hh"
#include "Array.hh"


	// CONSTRUCTORS

Array::Array(int xSize, int ySize)
{

	xSize_ = xSize;
	ySize_ = ySize;
	size_ = xSize*ySize;
	if(xSize < 0 || ySize < 0)
	{
		std::cerr << "xSize and ySize must be positive" << std::endl;
	}
	ar = new real[xSize*ySize];

}

Array::Array(const Array &s)
{

	xSize_ = s.getSize(DIM_1D);
	ySize_ = s.getSize(DIM_2D);
	size_ = xSize_*ySize_;
	
	ar = new real[xSize_*ySize_];
	real *p_s = s.getArray();

	for(int i = 0; i < size_; i++)
	{
		ar[i] = p_s[i];
	}

}

Array& Array::operator=(const Array& s)
{
        real *local_ar = new real[s.getSize()];
        real *s_ar = s.getArray();

        for (int i = 0; i < s.getSize(); i++)
        {
                local_ar[i] = s_ar[i];
        }

        this->ar = local_ar;
        this->xSize_ = s.getSize(DIM_1D);
        this->ySize_ = s.getSize(DIM_2D);

        return *this;
}

//DESTRUCTOR
Array::~Array()
{

	if(NULL != ar)
	{
		delete[] ar;
	}

}


	//CONVENIENCE FUNCTIONS

//initialize the whole array with a constant value
void Array::fill( real value )
{

        for (int i = 0; i < getSize(); i++)
        {
                ar[i] = value;
        }

}

real * Array::getArray() const
{
        return ar;
}

void Array::print()
{

                for (int i = ySize_ - 1; i >= 0; i--)
                {
                        for (int j = 0; j < xSize_; j++)
                        {
                                std::cout << Array::operator ()(j,i) << " ";
                        }
                        std::cout << std::endl;
                }
}

int Array::getSize( int dimension ) const
{

        int ret = 0;
        switch (dimension)
        {
        case DIM_1D:
                ret = xSize_;
                break;
        case DIM_2D:
                ret = ySize_;
                break;
	}

	return ret;

}

int Array::getSize() const
{

   return size_;

}

real& Array::operator ()(int i,int j)
{

   //static double dummy;
//    CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
//    CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
//    CHECK_MSG(dimension_ == DIM_2D, "Wrong dimension.");
   return ar[xSize_*j + i];

}

