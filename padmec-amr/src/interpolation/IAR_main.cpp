/*
 * IAR_main.cpp
 *
 *  Created on: 14/02/2013
 *      Author: rogsoares
 *
 *      IAR - Interpolation Adaptative Remeshing
 *      Here are called all interpolation function for remeshing procedure
 */


#include "interpolation.h"

//void Linear(InterpolationDataStruct* pIData){
//#ifdef TRACKING_PROGRAM_STEPS
//	cout << "TRACKING_PROGRAM_STEPS: Interpolation (Linear Method)\tIN\n";
//#endif
//	int dim = pIData->m1->getDim();
//	pIData->theOctree = OctreeCreate2<iterall>(pIData->m2->beginall(dim),pIData->m2->endall(dim),dim);
//	calculate_GeometricCoefficients(pIData,dim);
//	calculate_LinearInterpolation(pIData,dim);
//#ifdef TRACKING_PROGRAM_STEPS
//	cout << "TRACKING_PROGRAM_STEPS: Interpolation (Linear Method)\tOUT\n";
//#endif
//}
//
//void Quadratic(InterpolationDataStruct* pIData){
//#ifdef TRACKING_PROGRAM_STEPS
//	cout << "TRACKING_PROGRAM_STEPS: Interpolation (Quadratic Method)\tIN\n";
//#endif
//	initialize(pIData); cout << __LINE__ << endl;
//	Linear(pIData); cout << __LINE__ << endl;
//	for (int field=0; field<pIData->numFields; field++){
//		calculate_Gradients(pIData,field); cout << __LINE__ << endl;
//		calculate_DerivativesError(pIData); cout << __LINE__ << endl;
//		calculate_QuadraticInterpolation(pIData); cout << __LINE__ << endl;
//	}
//	finalize(pIData); cout << __LINE__ << endl;
//#ifdef TRACKING_PROGRAM_STEPS
//	cout << "TRACKING_PROGRAM_STEPS: Interpolation (Quadratic Method)\tIN\n";
//#endif
//}

//void Adaptative(InterpolationDataStruct* pIData){
//	throw Exception(__LINE__,__FILE__,"Under construction!");
//}
//
//void Conservative(InterpolationDataStruct* pIData){
//	throw Exception(__LINE__,__FILE__,"Under construction!");
//}
//
//void PureInjection(InterpolationDataStruct* pIData){
//	throw Exception(__LINE__,__FILE__,"Under construction!");
//}
//
//void HalfWeighting(InterpolationDataStruct* pIData){
//	throw Exception(__LINE__,__FILE__,"Under construction!");
//}
//
//void FullWighting(InterpolationDataStruct* pIData){
//	throw Exception(__LINE__,__FILE__,"Under construction!");
//}

