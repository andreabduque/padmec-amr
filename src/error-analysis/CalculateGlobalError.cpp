/*
 * CalculateGlobalError.cpp
 *
 *  Created on: 28/02/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis.h"

void ErrorAnalysis::calculate_GlobalError(GeomData* pGCData){
	double error_sum = calculate_ErrorSum(pGCData,false);
	double smooth_gradNorm = getSmoothedGradNorm();
	double global_error = sqrt(fabs(error_sum))/smooth_gradNorm;
	setGlobalError(global_error);

//	cout << "\n\nerror_sum:       " << error_sum << endl;
//	cout << "smooth_gradNorm: " << smooth_gradNorm << endl;
//	cout << "global_error:    " << global_error << "\n\n";
}

void ErrorAnalysis::calculate_GlobalError_Singularity(GeomData* pGCData){
	double error_sum = calculate_ErrorSum(pGCData,true);
	double smooth_gradNorm = getSmoothedGradNorm_Singularity();
	double global_error_singular = sqrt(fabs(error_sum))/smooth_gradNorm;
	setGlobalError_Singularity(global_error_singular);

//	cout << "\n\nerror_sum:       " << error_sum << endl;
//	cout << "smooth_gradNorm: " << smooth_gradNorm << endl;
//	cout << "global_error:    " << global_error_singular << "\n\n";
}

// Define an average error. It's an distributed error over the mesh. * 	//eq. 4.24 e 3.2.1 (p�g 32) do algoritmo - tese de Filipe
void ErrorAnalysis::calculate_AvgError(int nelem, double tol, bool singularity){
	tol = .01;
	double SGN = (singularity)?getSmoothedGradNorm_Singularity():getSmoothedGradNorm();
	countElements(nelem,singularity);
	double numElements = (double)(singularity)?getNumElements_Singularity():nelem;
	double avgError = (numElements)?(tol*( SGN/sqrt(numElements))):.0;
	if (singularity){
		setAverageError_Singularity(avgError);
	}
	else{
		setAverageError(avgError);
	}

//	cout << "\nSGN:       " << SGN << endl;
//	cout << "singularity:    " << singularity << "\n";
//	cout << "numElements: " << numElements << endl;
//	cout << "tol: " << tol << endl;
//	cout << "avgError:    " << avgError << "\n\n";
}
