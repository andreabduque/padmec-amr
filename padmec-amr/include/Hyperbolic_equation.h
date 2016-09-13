/*
 * Hyperbolic_equation.h
 *
 *  Created on: 09/01/2009
 *      Author: rogerio
 */

#ifndef HYPERBOLIC_EQUATION_H_
#define HYPERBOLIC_EQUATION_H_

/**
 * Base class for hypebolic solvers.
 */

#include "OilProductionManagement.h"
#include "exportVTK.h"

namespace PRS
{
	class Hyperbolic_equation{
	public:

		Hyperbolic_equation() {}
		virtual ~Hyperbolic_equation(){}
		virtual double solver(pMesh,double&)=0;
	};
}
#endif /* HYPERBOLIC_EQUATION_H_ */
