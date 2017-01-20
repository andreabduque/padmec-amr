///*
// * H_Refinement_2D.cpp
// *
// *  Created on: 06/06/2012
// *      Author: rogsoares
// */
//
//#include "H_Refinement_2D.h"
//
//H_Refinement_2D::H_Refinement_2D(){
//}
//
//H_Refinement_2D::H_Refinement_2D(pMesh theMesh){
//	pCallback = new Adaptation_SplitCallBack(theMesh);
//}
//
//H_Refinement_2D::~H_Refinement_2D(){
//}
//
///* STRATEGY:
// * step 1: element's flags correction (unity unevennes and leaves correction)
// * step 2: mesh adaptation
// * step 3: reconstruct mesh database
// * step 4: eliminate hangnodes
// * step 5: reconstruct mesh database
// */
//
//void H_Refinement_2D::rodar(ErrorAnalysis *pErrorAnalysis, pMesh & theMesh){
//	if (theMesh->getDim() != 2){
//		throw Exception(__LINE__,__FILE__,"Adaptation procedure has been implemented only for 2D meshes!\n");
//	}
//	removeSpecialElements(theMesh,pErrorAnalysis);
//	refreshDataStructure(theMesh);							// step 5:
//	elementsFlagsCorrection(pErrorAnalysis,theMesh);		// step 1:
//	adaptation(theMesh);									// step 2:
//	refreshDataStructure(theMesh);							// step 3:
//	createSpecialElementsList(theMesh,pErrorAnalysis);
//	eliminateHangNodes(theMesh);							// step 4:
//	refreshDataStructure(theMesh);							// step 5:
//
//}
//
//void H_Refinement_2D::refreshDataStructure(pMesh theMesh){
//	theMesh->modifyState(1,2,1);
//
//}
//
//void H_Refinement_2D::elementsFlagsCorrection(ErrorAnalysis* pErrorAnalysis, pMesh theMesh){
//	//pErrorAnalysis->calculate_MaxMinErrorData(theMesh);
//	int maxDepth = pErrorAnalysis->getMaxDepth();
//	int max      = pErrorAnalysis->getMaxRefinementFlag();
//	int min      = pErrorAnalysis->getMinRefinementFlag();
//
//	list<pEntity> faceList;
//
//	//desnUnit(max, min, theMesh, faceList);
//	cout << __LINE__ << endl;
//	unevenElements(pErrorAnalysis,theMesh,faceList,max,min);
//	cout << __LINE__ << endl;
//
//	//Correcao entre irmaos
//	//correction = correcIrmao (maxDepth, theMesh);
//	bool correction = leavesCorrectionMain(pErrorAnalysis,theMesh,maxDepth,faceList);
//
//	faceList.clear();
//
//	//Desnivel e regularizacao (se houver alteracao em correcao entre irmaos)
//	if (correction){
//		//desnUnit(max, min, theMesh, faceList);
//		cout << __LINE__ << endl;
//		unevenElements(pErrorAnalysis,theMesh,faceList,max,min);
//		cout << __LINE__ << endl;
//		correction = false;
//		if(!faceList.empty())
//			//correction = correctIrmao (theMesh, maxDepth, faceList);
//			correction = leavesCorrectionMain(pErrorAnalysis,theMesh,maxDepth,faceList);
//
//		faceList.clear();
//		while (correction){
//			//desnUnit(max, min, theMesh, faceList);
//			cout << __LINE__ << endl;
//			unevenElements(pErrorAnalysis,theMesh,faceList,max,min);
//			cout << __LINE__ << endl;
//
//			correction = false;
//			if(!faceList.empty())
////				correction = correctIrmao (theMesh, maxDepth, faceList);
//				correction = leavesCorrectionMain(pErrorAnalysis,theMesh,maxDepth,faceList);
//
//			faceList.clear();
//		}
//	}
//	this->setEFC(true);	// OK. Now mesh adaptation can be done.
//}
