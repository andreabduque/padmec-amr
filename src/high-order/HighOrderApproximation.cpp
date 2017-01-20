//#include "HighOrderApproximation.h"
//
//namespace PRS{
//
//	HighOrderApproximation::HighOrderApproximation(){
//	}
//
//	HighOrderApproximation::HighOrderApproximation(void *p){
//		pStruct = (PointerStruct*)(p);
//		pNodeSL = getNodeSlopeLimitFunc();
//		pEdgeSL = getEdgeSlopeLimitFunc();
//		koef = pStruct->pSimPar->get_koef();
//	}
//
//	HighOrderApproximation::~HighOrderApproximation(){
//	}
//
//	/*
//	 * Initialize a pointer to node slope limiter
//	 */
//	NodeSlopeLimiter* HighOrderApproximation::getNodeSlopeLimitFunc(){
//		int val = pStruct->pSimPar->getNodeSlopeLimitFunc();
//		if (val==node_MUSCL || val==node_Superbee || val== node_Minmod || val==node_Osher || val==node_Van_Albada)
//			return new Node_MUSCL(pStruct);
////		else if (val==node_WoodField)
////			return new Node_WoodField(pStruct);
//		else
//			throw Exception(__LINE__,__FILE__,"Unknown high order methods.\n");
//	}
//
//	EdgeSlopeLimiter* HighOrderApproximation::getEdgeSlopeLimitFunc(){
//		int val = pStruct->pSimPar->getNodeSlopeLimitFunc();
//		if (val==MUSCL || val==SUPERBEE || val== MINMOD || val==OSHER || val==VAN_ALBADA)
//			return new Edge_MUSCL(pStruct);
////		else if (val==WOODFIELD)
////			return new Edge_WoodField(pStruct);
//		else
//			throw Exception(__LINE__,__FILE__,"Unknown high order methods.\n");
//	}
//
//	double HighOrderApproximation::getSw_HighOrderApproximation(pEntity edge, int dom, double &Sw_I, double &Sw_J, int dim){
//		// start CPU time clock
//		double startt = MPI_Wtime();
//
//		// get nodes I and J
//		pEntity I = (pVertex)edge->get(0,0);
//		pEntity J = (pVertex)edge->get(0,1);
//
//		double DSwII, DSwJJ;
//		double delta_Sw_JI = Sw_J - Sw_I;
//
//		double SLII,SLJJ;			// edge slope limiter
//		double slimit_I, slimit_J;	// nodal slope limiter
//
//		pNodeSL->calculateNodeSlopeLimiters(I,J,slimit_I,slimit_J);
//		pEdgeSL->calculateEdgeSlopeLimiters(edge,Sw_I,Sw_J,SLII,SLJJ,DSwII,DSwJJ,dim);
//
//		// final limiter function
//		const double SLI = SLII*slimit_I;
//		const double SLJ = SLJJ*slimit_J;
//
//		/*
//		 * Modified Taylor expansion series extrapolating nodal saturation value
//		 * on volume control interfaces. Injection wells are excluded.
//		 */
//		if ( !pStruct->pSimPar->isInjectionWell(GEN_tag(I->getClassification())) )
//			Sw_I = Sw_I + (SLI/4.)*((1.-koef)*DSwII + (1.+koef)*delta_Sw_JI);
//		if ( !pStruct->pSimPar->isInjectionWell(GEN_tag(J->getClassification())) )
//			Sw_J = Sw_J - (SLJ/4.)*((1.-koef)*DSwJJ + (1.+koef)*delta_Sw_JI);
//
//		// finish CPU time clock. Return CPU time consumption (seconds)
//		double endt = MPI_Wtime();
//		return endt-startt;
//	}
//}
