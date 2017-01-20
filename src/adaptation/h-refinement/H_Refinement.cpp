///*
// * H_Refinement.cpp
// *
// *  Created on: 04/06/2012
// *      Author: rogsoares
// */
//
//#include "H_Refinement.h"
//
//
///* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// *function getLeaves: get the leaves of a parent
// *function getAllSubTree: get the tree (include the parent itself)
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//void H_Refinement::getChildren (pMesh theMesh, pEntity parent, std::vector<pFace> &fchildren){
////	list<mEntity*>::iterator lit;
////	list<mEntity*> leaves;
////	//parent was refined only once
////	if (theMesh->getRefinementDepth(parent) == 1){
////		parent->getLeaves(leaves);
////		for (lit = leaves.begin(); lit!=leaves.end(); lit++)
////			fchildren.push_back(*lit);
////	}
////	//parent was refined more than once
////	else if (theMesh->getRefinementDepth(parent) > 1){
////		parent->getAllSubTree(leaves);
////		for (lit = leaves.begin(); lit != leaves.end(); lit++){
////			//get only the parent children (difference between level is 1)
////			int levelDiff = theMesh->getRefinementLevel(*lit) - theMesh->getRefinementLevel(parent);
////			if ( levelDiff == 1 ){
////				fchildren.push_back(*lit);
////			}
////		}
////	}
////	leaves.clear();
//}
//
//void H_Refinement::adaptation(pMesh theMesh){
//////	if (!this->EFC_isDone()){
//////		throw Exception(__LINE__,__FILE__,"You MUST call elementsFlagsCorrection() to guarantee the correct element flags\n");
//////	}
////
////	/*
////	 * First unrefine (eliminate unnecessaries elements)
////	 * Then, refine what is required.
////	 */
////	pCallback->setUnrefine();
////	theMesh->refunref(*pCallback);
////	pCallback->setRefine();
////	theMesh->refunref(*pCallback);
////	this->setEFC(false);
//}
