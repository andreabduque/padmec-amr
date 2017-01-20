/*
 * Boundary_Conditions.h
 *
 *  Created on: Jan 29, 2015
 *      Author: rogerio
 */

#ifndef BOUNDARY_CONDITIONS_H_
#define BOUNDARY_CONDITIONS_H_

#include "auxiliar.h"

/*
 * When prescribed value (Dirichlet) and prescribed flux (Neumann) were necessary, you MUST define them here.
 */

enum BENCHMARK {CASE_1, CASE_5};

/* Case 1: */
/* Exact Solution  : */ double Benchmark3D_case1__ES(double x, double y, double z);
/* Source/sink term: */ double Benchmark3D_case1__SST(double x, double y, double z);

/* Case 5: */
/* Exact Solution  : */ double Benchmark3D_case5__ES(double x, double y, double z);
/* Source/sink term: */ double Benchmark3D_case5__SST(double x, double y, double z);


#endif /* BOUNDARY_CONDITIONS_H_ */
