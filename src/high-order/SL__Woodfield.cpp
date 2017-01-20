///*
// * SL__Woodfield.cpp
// *
// *  Created on: 08/02/2011
// *      Author: rogerio
// */
//
///* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * 							WoodField node slope limiter                       *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//
//#include "SL__functions.h"
//
//namespace PRS{
//
//	Node_WoodField::Node_WoodField(){}
//
//	Node_WoodField::Node_WoodField(void *p){
//		pStruct = (PointerStruct*)(p);
//		firstNode_I = 0;
//		firstNode_J = 0;
//	}
//
//	Node_WoodField::~Node_WoodField(){}
//
//	/*
//	 * Initialize max and min nodal saturation
//	 */
//	double Node_WoodField::initializeMaxMin_Sw(){
//		pEntity node;
//		double startt = MPI_Wtime();
//		VIter vit = M_vertexIter(pStruct->theMesh);
//		while ( (node = VIter_next(vit)) ){
//			double SwI = pStruct->pPPData->getSaturation(node);
//			pStruct->pPPData->setSw_max(node,.0);
//			pStruct->pPPData->setSw_min(node,1e+6);
//		}
//		VIter_delete(vit);
//		double endt = MPI_Wtime();
//		return endt-startt;
//	}
//
//	// WoodField node limiter
//	void Node_WoodField::calculateNodeSlopeLimiters(pEntity I, pEntity J, double &SLII, double &SLJJ){
//		SLII = pStruct->pPPData->getS_Limit(I);
//		SLJJ = pStruct->pPPData->getS_Limit(J);
//	}
//
//	/**
//	 * Define slope limiters for all nodes for domain dom.
//	 * It must be called before calculate fluxes across volume
//	 * control surfaces.
//	 */
//	double Node_WoodField::defineSlopeLimiters(){
//		double startt = MPI_Wtime();
//
//		// set parameters to find the min and max local saturation
//		initializeMaxMin_Sw();
//
//		const double DELTA = pStruct->pSimPar->get_WoodfieldDelta();
//		pEntity edge;
//
//		double Sw_min, Sw_max;
//
//		// loop over all edges
//		EIter eit = M_edgeIter(pStruct->theMesh);
//		while ( (edge = EIter_next(eit)) ){
//			pEntity I = (pVertex)edge->get(0,0);
//			pEntity J = (pVertex)edge->get(0,1);
//
//			// Finding Maximum and Minimum Saturations From Closest Neighbours
//			double SwI = pStruct->pPPData->getSaturation(I);
//			double SwJ = pStruct->pPPData->getSaturation(J);
//
//			// Finding Minimum Saturation From Closest Neighbors
//			if (SwJ < pStruct->pPPData->getSw_min(I) ) pStruct->pPPData->setSw_min(I,SwJ);
//			if (SwI < pStruct->pPPData->getSw_min(J) ) pStruct->pPPData->setSw_min(J,SwI);
//
//			// Finding Maximum Saturation From Closest Neighbors
//			if (SwJ > pStruct->pPPData->getSw_max(I) ) pStruct->pPPData->setSw_max(I,SwJ);
//			if (SwI > pStruct->pPPData->getSw_max(J) ) pStruct->pPPData->setSw_max(J,SwI);
//		}
//		EIter_delete(eit);
//
//		double S_Limit;
//		pEntity node;
//
//		VIter vit = M_vertexIter(pStruct->theMesh);
//		while ( (node = VIter_next(vit)) ){
//			double Sw = pStruct->pPPData->getSaturation(node);
//			Sw_max = pStruct->pPPData->getSw_max(node);
//			Sw_min = pStruct->pPPData->getSw_min(node);
//			double delta_Sw_max_min =  Sw_max - Sw_min;
//
//			// To avoid Division By Zero for defining Gama(I)
//			if (delta_Sw_max_min <= 1.0e-20)
//				S_Limit = 1.0;
//			else{
//				double GamaI = (Sw-Sw_min)/delta_Sw_max_min;
//				// Nodal Slope Limiter
//				if (GamaI >= 1.0 || GamaI <= .0)
//					S_Limit = .0;
//				else if ((GamaI >= DELTA) && (GamaI <= (1.0 - DELTA)))
//					S_Limit = 1.0;
//				else if ( (GamaI > .0) && (GamaI < DELTA) )
//					S_Limit = GamaI/DELTA;
//				else if ( (GamaI > 1. - DELTA) && (GamaI < 1.) )
//					S_Limit = (1. - GamaI)/DELTA;
//			}
//			// update S_Limit
//			pStruct->pPPData->setS_Limit(node,S_Limit);
//		}
//		VIter_delete(vit);
//
//		double endt = MPI_Wtime();
//		return endt-startt;
//	}
//
//	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//	 * 							WoodField edge slope limiter
//	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//
//	Edge_WoodField::Edge_WoodField(){}
//	Edge_WoodField::Edge_WoodField(void *p){
//		pStruct = (PointerStruct*)(p);
//		pNode_WF = new Node_WoodField(p);
//		firstEdge = 0;
//	}
//
//	Edge_WoodField::~Edge_WoodField(){
//	}
//
//	void Edge_WoodField::calculateEdgeSlopeLimiters(pEntity edge, const double &Sw_I,
//			                                        const double &Sw_J,double& SLII,
//			                                        double& SLJJ, double& DSwII,
//			                                        double& DSwJJ,int dim){
//
//		if (!firstEdge) firstEdge = edge;
//		if (edge==firstEdge) pNode_WF->defineSlopeLimiters();
//
//		int dom = 3300;
//
//		// get nodes I and J
//		pEntity I = (pVertex)edge->get(0,0);
//		pEntity J = (pVertex)edge->get(0,1);
//
//		/*
//		 * edIJ is a versor that points from I to J. I corresponds to smallest
//		 * edge node global ID
//		 */
//		double coord1[3]; V_coord(I,coord1);
//		double coord2[3]; V_coord(J,coord2);
//		dblarray edIJ(dim);
//		for (int i=0;i<dim;i++) edIJ[i] = coord1[i] - coord2[i];
//
//		/*
//		 * Take saturation gradients from nodes I and J
//		 */
//		dblarray Sw_grad_I(dim,.0), Sw_grad_J(dim,.0);
////		pStruct->pPPData->get_Sw_Grad(I,dom,Sw_grad_I);
////		pStruct->pPPData->get_Sw_Grad(J,dom,Sw_grad_J);
//
//		double delta_Sw = Sw_J - Sw_I;
//
//		double slimit_I = pStruct->pPPData->getS_Limit(I);
//		double slimit_J = pStruct->pPPData->getS_Limit(J);
//
//		// edge limiter for node I
//		double a = 1.0;
//		double inner = Sw_grad_I[0]*edIJ[0]+Sw_grad_I[1]*edIJ[1];
//		double b = 2.*delta_Sw/(slimit_I*inner + qsi);
//		double theta_I = std::min(a,b);
//		SLII = std::max(.0,theta_I);
//
//		// edge limiter for node J
//		a = 1.0;
//		inner = Sw_grad_J[0]*edIJ[0]+Sw_grad_J[1]*edIJ[1];
//		b = -2.*delta_Sw/(slimit_J*inner + qsi);
//		double theta_J = std::min(a,b);
//		SLJJ = std::max(.0,theta_J);
//	}
//}
