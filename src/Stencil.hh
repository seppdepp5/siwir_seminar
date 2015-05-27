#ifndef STENCIL_HH
#define STENCIL_HH

#include <cstdlib>
#include <iostream>
#include "Types.hh"

class Stencil
{
public:
	//CONSTRUCTORS

	Stencil(real w1, real w2, real w3, real w4, real w5, real w6, real w7, real w8, real w9);

	~Stencil();

	Stencil(const Stencil &s);

	//CLASS FUNCTIONS

	void print();
	real getw1() const;
	real getw2() const;
	real getw3() const;
	real getw4() const;
	real getw5() const;
	real getw6() const;
	real getw7() const;
	real getw8() const;
	real getw9() const;

private:

	real w1_;
	real w2_;
	real w3_;
	real w4_;
	real w5_;
	real w6_;
	real w7_;
	real w8_;
	real w9_;

};
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

#endif //STENCIL_HH
