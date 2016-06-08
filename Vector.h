#ifndef __VECTOR_H_INCLUDED__
#define __VECTOR_H_INCLUDED__

#include <iostream>
#include <Interpolation.h>
#include <math.h>

using namespace std;

class Vector{

	protected:

		double x;
		double y;

	public:
		Vector (Point P, Point Q);
		Vector (double _x, double _y);

		double getX();
		double getY();

		//vetor normal unitario orientacao positiva ao vetor da classe
		Vector* normal_unit_pos_vector();	
		
		//modulo (intensidade) do vetor	
		double module();
};


#endif