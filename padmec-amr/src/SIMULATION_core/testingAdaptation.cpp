///*
// * testingAdaptation.cpp
// *
// *  Created on: 15/05/2012
// *      Author: rogsoares
// */
//
//#include "SIMULATION_core.h"
//
//namespace PRS{
//
//
//
//	void initializeParameters(pMesh theMesh,ErrorAnalysis* pErrorAnalysis, MeshAdaptation* pMeshAdapt){
//		FIter fit = M_faceIter (theMesh);
//		while (pEntity face = FIter_next(fit)){
//			if (theMesh->getRefinementDepth(face)){
//				pErrorAnalysis->setLevelOfRefinement(face,0);
//				pMeshAdapt->setElementAsSpecial(face,false);
//			}
//		}
//		FIter_delete(fit);
//	}
//
//
//	void createRefinementLevel(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, double x_center, double step){
//		const double delta_1 = 0.03*x_center;
//		const double delta_2 = 0.07*x_center;
//		const double delta_3 = 0.11*x_center;
//		double fcenter[3];
//
//		FIter fit = M_faceIter (theMesh);
//		while (pEntity face = FIter_next(fit)){
//			if (!theMesh->getRefinementDepth(face)){
//				getFCenter(face,fcenter);
////				if (fabs(fcenter[0]-x_center) <= delta_1){
////					pErrorAnalysis->setLevelOfRefinement(face,3);
////				}
//
//				if (fabs(fcenter[0]-x_center) <= delta_3){
//					pErrorAnalysis->setLevelOfRefinement(face,2);
//				}
////				else
////				if (fabs(fcenter[0]-x_center) <= delta_2){
////					pErrorAnalysis->setLevelOfRefinement(face,2);
////				}
//			}
//		}
//	}
//
//	void createUnrefinementLevel(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, double x_center, double step){
////			const double x_center = 0.5;
//			const double delta_1 = 0.03*x_center;
//			const double delta_2 = 0.07*x_center;
//			const double delta_3 = 0.11*x_center;
//			double fcenter[3];
//
//			FIter fit = M_faceIter (theMesh);
//			while (pEntity face = FIter_next(fit)){
//				if (!theMesh->getRefinementDepth(face) && theMesh->getRefinementLevel(face)==2){
//					getFCenter(face,fcenter);
//					if (fabs(fcenter[0]-x_center) <= delta_3){
//						pErrorAnalysis->setLevelOfRefinement(face,-2);
//					}
////					else if (fabs(fcenter[0]-x_center) <= delta_2){
////						pErrorAnalysis->setLevelOfRefinement(face,2);
////					}
////					else if (fabs(fcenter[0]-x_center) <= delta_3){
////						pErrorAnalysis->setLevelOfRefinement(face,1);
////					}
//				}
//			}
//		}
//}
