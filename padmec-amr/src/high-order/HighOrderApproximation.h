///*
// * highOrderApproximation.h
// *
// *  Created on: 14/05/2009
// *      Author: rogerio
// */
//
//#ifndef HIGHORDERSATURATION_H_
//#define HIGHORDERSATURATION_H_
//
//#include "SL__functions.h"
//
//namespace PRS{
//
//	class HighOrderApproximation{
//
//	public:
//
//		HighOrderApproximation();
//		HighOrderApproximation(void *);
//		~HighOrderApproximation();
//
//		/**
//		 * To calculate fluxes with high order approximation for saturation field.
//		 * Arguments: edge pointer, domain flag, Sw_I, Sw_J, respectively.
//		 */
//		virtual double getSw_HighOrderApproximation(pEntity, int, double&, double&,int);
//
//		NodeSlopeLimiter* getNodeSL_Ptr(){
//			return pNodeSL;
//		}
//
//	private:
//
//		PointerStruct *pStruct;
//		SimulatorParameters *pSimPar;
//
//		/*
//		 * Slope limiter function pointers
//		 */
//		EdgeSlopeLimiter* pEdgeSL;
//		NodeSlopeLimiter* pNodeSL;
//
//		NodeSlopeLimiter* getNodeSlopeLimitFunc();
//		EdgeSlopeLimiter* getEdgeSlopeLimitFunc();
//
//		double koef;
//	};
//}
//
//#endif /* HIGHORDERAPPROXIMATION_H_ */
