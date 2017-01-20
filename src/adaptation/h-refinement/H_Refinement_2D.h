///*
// * MeshAdaptation__h-type.h
// *
// *  Created on: 24/01/2013
// *      Author: rogsoares
// */
//
//#ifndef MESHADAPTATION__H_TYPE_H_
//#define MESHADAPTATION__H_TYPE_H_
//
//#include "H_Refinement.h"
//
//class H_Refinement_2D : public H_Refinement{
//
//public:
//	H_Refinement_2D();
//	H_Refinement_2D(pMesh);
//	virtual ~H_Refinement_2D();
//
//	void removeSpecialElements(pMesh, ErrorAnalysis *);
//	void createSpecialElementsList(pMesh, ErrorAnalysis*);
//	void eliminateHangNodes(pMesh);
//	void unevenElements(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, std::list<pEntity>&, int, int);
//	void refineMesh(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, int maxLevel);
//	int checkMaxDORdifference(ErrorAnalysis* pErrorAnalysis, pMesh theMesh);
//	void meshRegularization(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, int min, std::list<pEntity>& faceList);
//	bool leavesCorrectionMain(ErrorAnalysis*, pMesh, int, std::list<pEntity>&);
//	void firstLeavesCorrection(ErrorAnalysis*, pMesh, pEntity, int, bool&);
//	void secondLeavesCorrection(ErrorAnalysis*, pMesh, pEntity, int depth, bool&);
//	void refreshDataStructure(pMesh);
//	void elementsFlagsCorrection(ErrorAnalysis*, pMesh);
//	void rodar(ErrorAnalysis  *pErrorAnalysis, pMesh & theMesh);
//
//	// Inline functions (do not use it!)
//	int getNumElements(pMesh theMesh){
//		return -1;
//	}
//};
//#endif /* MESHADAPTATION__H_TYPE_H_ */
