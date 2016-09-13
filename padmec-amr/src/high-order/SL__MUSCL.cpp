///*
// * SL__MUSCL.cpp
// *
// *  Created on: 08/02/2011
// *      Author: rogerio
// */
//
///*
// * SlopeLimiterFunctions.cpp
// *
// *  Created on: 07/07/2009
// *      Author: rogerio
// */
//
//#include "SL__functions.h"
//
//namespace PRS{
//
//	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//	 * 							MUSCL node slope limiter
//	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//
//	Node_MUSCL::Node_MUSCL(){}
//
//	Node_MUSCL::Node_MUSCL(void *p){
//		pStruct = (PointerStruct*)(p);
//	}
//
//	Node_MUSCL::~Node_MUSCL(){}
//
//	void Node_MUSCL::calculateNodeSlopeLimiters(pEntity I, pEntity J, double &SLII, double &SLJJ){
//		SLII = SLJJ =1.0;
//	}
//
//
//	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//	 * 							MUSCL edge slope limiter
//	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//	Edge_MUSCL::Edge_MUSCL(){}
//
//	Edge_MUSCL::Edge_MUSCL(void *p){
//		pStruct = (PointerStruct*)(p);
//	}
//
//	Edge_MUSCL::~Edge_MUSCL(){
//	}
//
//	void Edge_MUSCL::calculateEdgeSlopeLimiters(pEntity edge, const double &Sw_I,const double &Sw_J,double& SLII,double& SLJJ, double& DSwII,double& DSwJJ, int dim){
//		// get nodes I and J
//		pEntity I = (pVertex)edge->get(0,0);
//		pEntity J = (pVertex)edge->get(0,1);
//
//		if (EN_id(edge->get(0,0)) > EN_id(edge->get(0,1))){
//			std::swap(I,J);
//		}
//
//		double coord1[3]; V_coord(I,coord1);
//		double coord2[3]; V_coord(J,coord2);
//		double edIJ[3] = {.0,.0,.0};
//		for (int i=0; i<dim; i++){
//			edIJ[i] = (coord2[i]-coord1[i]);
//		}
//
//		double Sw_grad_I[3] = {.0,.0,.0}, Sw_grad_J[3] = {.0,.0,.0};
//		pStruct->pPPData->get_Sw_Grad(I,0,0,Sw_grad_I);
//		pStruct->pPPData->get_Sw_Grad(J,0,0,Sw_grad_J);
//		double delta_Sw = (Sw_J - Sw_I);
//
//		// Upwind-biased Interpolations
//		DSwII = 2.*inner_product(Sw_grad_I,edIJ,dim) - delta_Sw;
//		DSwJJ = 2.*inner_product(Sw_grad_J,edIJ,dim) - delta_Sw;
//
//		double ratioI = (2.*DSwII*delta_Sw + qsi) / ( DSwII*DSwII + delta_Sw*delta_Sw + qsi);
//		double ratioJ = (2.*DSwJJ*delta_Sw + qsi) / ( DSwJJ*DSwJJ + delta_Sw*delta_Sw + qsi);
//
////		switch (pStruct->pSimPar->getEdgeSlopeLimitFunc()){
////			// MUSCL-MUSCL � o Van Leer
////			case MUSCL:{
//				SLII = (ratioI + fabs(ratioI) + qsi)/(1. + fabs(ratioI) + qsi);
//				SLJJ = (ratioJ + fabs(ratioJ) + qsi)/(1. + fabs(ratioJ) + qsi);
////				break;
////			}
////			case VAN_ALBADA:{
////				SLII = max(.0,ratioI);
////				SLJJ = max(.0,ratioJ);
////				break;
////			}
////			case SUPERBEE:{
////				double AuxLI = max(.0, min(1.,2.*ratioI));
////				double AuxLJ = max(.0, min(1.,2.*ratioJ));
////				SLII = max(AuxLI, min(2.,ratioI));
////				SLJJ = max(AuxLJ, min(2.,ratioJ));
////				break;
////			}
////			case MINMOD:{
////				SLII = max(.0, min(1.,ratioI));
////				SLJJ = max(.0, min(1.,ratioJ));
////				break;
////			}
////			case OSHER:{
////				SLII = max(.0, min(2.,ratioI));
////				SLJJ = max(.0, min(2.,ratioJ));
////				break;
////			}
////			default:
////				throw Exception(__LINE__,__FILE__,"ERROR: unknown edge limiter.\n");
////		}
//	}
//}
