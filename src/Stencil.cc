#include <cstdlib>
#include <iostream>
#include "Stencil.hh"
#include "Types.hh"


	//CONSTRUCTORS

Stencil::Stencil(real w1, real w2, real w3,
		 real w4, real w5, real w6,
		 real w7, real w8, real w9)
{

	w1_ = w1;
	w2_ = w2;
	w3_ = w3;
	w4_ = w4;
	w5_ = w5;
	w6_ = w6;
	w7_ = w7;
	w8_ = w8;
	w9_ = w9;

}

Stencil::Stencil(const Stencil &s)
{
	
	w1_ = s.getw1();
	w2_ = s.getw2();
	w3_ = s.getw3();
	w4_ = s.getw4();
	w5_ = s.getw5();
	w6_ = s.getw6();
	w7_ = s.getw7();
	w8_ = s.getw8();
	w9_ = s.getw9();

}

Stencil::~Stencil()
{
	
	w1_ = 0.0;
	w2_ = 0.0;
	w3_ = 0.0;
	w4_ = 0.0;
	w5_ = 0.0;
	w6_ = 0.0;
	w7_ = 0.0;
	w8_ = 0.0;
	w9_ = 0.0;

}

	//CLASS FUNCTIONS

real Stencil::getw1() const
{
	return w1_;
}
real Stencil::getw2() const
{
	return w2_;
}
real Stencil::getw3() const
{
	return w3_;
}
real Stencil::getw4() const
{
	return w4_;
}
real Stencil::getw5() const
{
	return w5_;
}
real Stencil::getw6() const
{
	return w6_;
}
real Stencil::getw7() const
{
	return w7_;
}
real Stencil::getw8() const
{
	return w8_;
}
real Stencil::getw9() const
{
	return w9_;
}

void Stencil::print()
{

		std::cout << "w1: " << Stencil::getw1() << std::endl;
		std::cout << "w2: " << Stencil::getw2() << std::endl;
		std::cout << "w3: " << Stencil::getw3() << std::endl;
		std::cout << "w4: " << Stencil::getw4() << std::endl;
		std::cout << "w5: " << Stencil::getw5() << std::endl;
		std::cout << "w6: " << Stencil::getw6() << std::endl;
		std::cout << "w7: " << Stencil::getw7() << std::endl;
		std::cout << "w8: " << Stencil::getw8() << std::endl;
		std::cout << "w9: " << Stencil::getw9() << std::endl;

}

/*
struct Stencil
{
	
	// Constructor
	explicit inline Stencil ();

	// Member variables
	real C; // The central stencil entry .
	real N; // The stencil entry for access
	real S; // The stencil entry for access
	real W; // The stencil entry for access
	real E; // The stencil entry for access
	real NW; // The stencil entry for access
	real NE; // The stencil entry for access
	real SW; // The stencil entry for access
	real SE; // The stencil entry for access

};
*/
