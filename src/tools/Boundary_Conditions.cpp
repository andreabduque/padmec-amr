/*
 * Boundary_Conditions.cpp
 *
 *  Created on: Jan 29, 2015
 *      Author: rogerio
 */

#include "Boundary_Conditions.h"

double Benchmark3D_case1__ES(double x, double y, double z){
	return 1.0 + sin(pi*x) * sin(pi*(y+.5)) * sin(pi*(z+.33333333333));
}

double Benchmark3D_case1__SST(double x, double y, double z){
	double A = double(sin(pi*x));
	double B = double(sin(pi*(y+.5)));
	double C = double(sin(pi*(z+.3333333333)));
	double D = double(cos(pi*x));
	double E = double(cos(pi*(y+.5)));
	double F = double(cos(pi*(z+.3333333333)));
	double MINUS_DIV_K_GRAD_U = pi*pi*( 3.0*A*B*C - D*E*C - A*E*F );
	return MINUS_DIV_K_GRAD_U;
}

double Benchmark3D_case5__ES(double x, double y, double z){
	double alpha;
	if (y<=.5 && z<=.5){
		alpha = 0.1;
	}
	else if (y>.5 && z<=.5){
		alpha = 10.;
	}
	else if (y>.5 && z>.5){
		alpha = 100.;
	}
	else if (y<=.5 && z>.5){
		alpha = 0.01;
	}
	else{
		throw Exception(__LINE__,__FILE__,"Alpha could not be calculated.");
	}

	double A = (double)sin(2.*pi*x);
	double B = (double)sin(2.*pi*y);
	double C = (double)sin(2.*pi*z);
	double excsol = (double)(alpha*A*B*C);
	if (fabs(excsol)<1e-8){
		excsol = .0;
	}
	return excsol;
}

double Benchmark3D_case5__SST(double x, double y, double z){
	double alpha, ax, ay, az;
	if (y<=.5 && z<=.5){
		alpha = 0.1;
		ax = 1.;
		ay = 10.;
		az = 0.01;
	}
	else if (y>.5 && z<=.5){
		alpha = 10.;
		ax = 1;
		ay = 0.1;
		az = 100.;
	}
	else if (y>.5 && z>.5){
		alpha = 100.;
		ax = 1.;
		ay = 0.01;
		az = 10.;
	}
	else if (y<=.5 && z>.5){
		alpha = 0.01;
		ax = 1.;
		ay = 100.;
		az = 0.1;
	}
	return (ax + ay+ az)*(4.*pi*pi*alpha*sin(2.*pi*x)*sin(2.*pi*y)*sin(2.*pi*z));
}

