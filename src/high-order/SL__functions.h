///*
// * SlopeLimiterFunctions.h
// *
// *  Created on: 07/07/2009
// *      Author: rogerio
// */
//
//#ifndef SLOPELIMITERFUNCTIONS_H_
//#define SLOPELIMITERFUNCTIONS_H_
//
//#include "PointerStruct.h"
//
//namespace PRS{
//
//	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//	 * The following two classes (NodeSlopeLimiter and EdgeSlopeLimiter) are base
//	 * classes with pure virtual functions and intended to be function interfaces.
//	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//
//	/*
//	 * Nodal slope limiter base class.
//	 */
//
//	class NodeSlopeLimiter{
//
//	public:
//
//		NodeSlopeLimiter(){}
//		~NodeSlopeLimiter(){}
//		virtual void calculateNodeSlopeLimiters(pEntity,pEntity,double&,double&)=0;
//
//		/**
//		 * It performs some calculations that must be done every time step before
//		 * applying a high order approximation for fluxes.
//		 */
//		virtual double defineSlopeLimiters()=0;
//
//	protected:
//		/**
//		 * It points to a set of pointers commonly used by main classes
//		 * It help to reduce the number of passing arguments
//		 */
//		PointerStruct *pStruct;
//	};
//
//	/*
//	 * Edge slope limiter base class.
//	 */
//	class EdgeSlopeLimiter{
//	public:
//		EdgeSlopeLimiter(){}
//		~EdgeSlopeLimiter(){}
//		virtual void calculateEdgeSlopeLimiters(pEntity,const double&,const double&,
//				double&,double&,double&,double&,int)=0;
//
//	protected:
//		PointerStruct *pStruct;
//	};
//
//	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//	 * The classes below are implementations for high order approximations.
//	 * They are derived from the classes above.
//	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//
//	/*
//	 * Nodal slope limiter MUSCL class.
//	 */
//	class Node_MUSCL : public NodeSlopeLimiter{
//	public:
//		Node_MUSCL();
//		Node_MUSCL(void *);
//		~Node_MUSCL();
//		void calculateNodeSlopeLimiters(pEntity,pEntity,double&,double&);
//
//		/**
//		 * MUSCL approximation does not need any previous calculation to approximate fluxes.
//		 */
//		double defineSlopeLimiters(){
//			return 0;
//		}
//	};
//
//	/*
//	 * Edge slope limiter MUSCL class.
//	 */
//	class Edge_MUSCL : public EdgeSlopeLimiter{
//	public:
//		Edge_MUSCL();
//		Edge_MUSCL(void*);
//		~Edge_MUSCL();
//		void calculateEdgeSlopeLimiters(pEntity,const double&,const double&,
//				                        double&,double&, double&, double&,int);
//	};
//
//
//	/*
//	 * Nodal slope limiter WoodField class.
//	 */
//	class Node_WoodField : public NodeSlopeLimiter{
//	public:
//		Node_WoodField();
//		Node_WoodField(void *p);
//		~Node_WoodField();
//		void calculateNodeSlopeLimiters(pEntity,pEntity,double&,double&);
//
//		double defineSlopeLimiters();
//
//	private:
//		double initializeMaxMin_Sw();
//		pEntity firstNode_I;
//		pEntity firstNode_J;
//	};
//
//	/*
//	 * Edge slope limiter WoodField class.
//	 */
//	class Edge_WoodField : public EdgeSlopeLimiter{
//	public:
//		Edge_WoodField();
//		Edge_WoodField(void *p);
//		~Edge_WoodField();
//		void calculateEdgeSlopeLimiters(pEntity,const double&,const double&,
//				double&,double&, double&, double&, int);
//	private:
//
//		Node_WoodField* pNode_WF;
//		pEntity firstEdge;
//	};
//}
//#endif /* SLOPELIMITERFUNCTIONS_H_ */
