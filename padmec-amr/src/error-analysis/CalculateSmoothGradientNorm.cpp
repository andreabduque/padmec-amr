/*
 * CalculateSmoothGradientNorm.cpp
 *
 *  Created on: 24/04/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis.h"

/*
 * Calculate smooth gradient Norm (SGN):
 * 	SNG = sum((|gradSw_I|+|gradSw_J|+|gradSw_K|)/3.0)^2
 * 	where sum is for each triangle element with nodes I, J, K.
 */
void ErrorAnalysis::calculate_SmoothedGradientNorm(SimulatorParameters *pSimPar, GeomData* pGCData, void(*pFunc_getGrad)(FIELD,int,int,int,const double*&), FIELD field){
	int i;
	int dim = pGCData->getMeshDim();			// mesh dimension
	int pos = dim+1;							// Auxiliary variable
	const int* indices = NULL;					// indices: elements connectivities

	// arrays of pointers for arrays: gradients and delta gradients
	const double* grad[4] = {NULL,NULL,NULL,NULL};
//	for (i=0; i<dim+1; i++){
//		grad[i] = new double[3];
//	}

	// loop over domains
	double sgn_tmp, SGN = 0;
	for (int dom=0; dom<pSimPar->getNumDomains(); dom++){
		int nelements = pGCData->getNumElemPerDomain(dom);

		// loop over domains' faces
		for (int row=0; row<nelements; row++){
			pGCData->getElement(dom,row,indices);

			// get nodal gradients for one element
			sgn_tmp = 0;
			for (i=0; i<pos; i++){
				pFunc_getGrad(field,dom,indices[i],indices[i+pos],grad[i]);
				sgn_tmp += inner_product(grad[i],grad[i],dim);;
			}
			SGN += (double)(sgn_tmp/pos);
		}
	}
	setSmoothedGradNorm(sqrt(SGN));
//	for (i=0; i<dim+1; i++){
//		delete[] grad[i]; grad[i] = 0;
//	}
	//cout << "SGN (all elements) = " << sqrt(SGN) << endl;
}

void ErrorAnalysis::calculate_SmoothedGradientNorm_Singularity(SimulatorParameters *pSimPar, GeomData* pGCData, void(*pFunc_getGrad)(FIELD,int,int,int,const double*&), FIELD field){
	int i;
	int dim = pGCData->getMeshDim();			// mesh dimension
	int pos = dim+1;							// Auxiliary variable
	const int* indices = NULL;					// indices: elements connectivities

	// arrays of pointers for arrays: gradients and delta gradients
	const double* grad[4] = {NULL,NULL,NULL,NULL};;
//	for (i=0; i<dim+1; i++){
//		grad[i] = new double[3];
//	}

	int ith_elem = 0;
	int counter = 0;

	// loop over domains
	double sgn_tmp, SGN = 0;
	for (int dom=0; dom<pSimPar->getNumDomains(); dom++){
		int nelements = pGCData->getNumElemPerDomain(dom);

		// loop over domains' faces
		for (int row=0; row<nelements; row++){

			// exclude elements flagged as singular
			if ( !isSingular(ith_elem) ){
				counter++;
				pGCData->getElement(dom,row,indices);

				// get nodal gradients for one element
				sgn_tmp = 0;
				for (i=0; i<pos; i++){
					pFunc_getGrad(field,dom,indices[i],indices[i+pos],grad[i]);
					sgn_tmp += inner_product(grad[i],grad[i],dim);;
				}
				SGN += (double)(sgn_tmp/pos);
			}
			ith_elem++;
		}
	}
	setSmoothedGradNorm_Singularity( sqrt(SGN) );

//	for (i=0; i<dim+1; i++){
//		delete[] grad[i]; grad[i] = 0;
//	}
//	cout << "SGN (singular) = " << sqrt(SGN) << endl;
//	cout << "counter = " << counter << endl;
}
