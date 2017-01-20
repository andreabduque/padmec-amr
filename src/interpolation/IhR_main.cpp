///*
// *  Universidade Federal de Pernambuco
// *  IhR.cpp
// *  IhR - Interpolation h-Refinement
// *  Created by Rogerio Soares on 15/02/2013.
// */
//
//#include "interpolation.h"
//
//void interpolate(pMesh theMesh, std::list<mEntity*> &leaves, int parentDepth, GetDblFunction* pArray_GetFunctions, SetDblFunction* pArray_SetFunctions, int numFields);
//
//void hRefinement(InterpolationDataStruct* pIntpData){
//	std::cout<<"Data interpolation begins...";
//	list<mEntity*> leaves;
//
//	FIter fit = M_faceIter(pIntpData->m1);
//	while (pEntity face = FIter_next(fit)){
//		// Do not look for special elements (with hang-nodes) and leaves elements (elements without children)
//		int levelToRefine = pIntpData->getLevelOfRefinement(face);				// how many time element must be refined
//		int depth = pIntpData->m1->getRefinementDepth(face);					// how many time element was refined
//		int ref = pIntpData->m1->getRefinementLevel(face);						// how many time element was refined
//		bool isSpecial = pIntpData->isElementSpecial(face);
//
//		if ( !isSpecial ){
//			if ( depth!=0 && !ref){
//				if (levelToRefine==depth || levelToRefine == depth-1){
//					face->getAllSubTree(leaves);
//					interpolate(pIntpData->m1,leaves,depth,pIntpData->pGetDblFunctions,pIntpData->pSetDblFunctions,pIntpData->numFields);
//				}
//			}
//		}
//	}
//
//	FIter_delete(fit);
//	//unflagEdgesAsInterpolated(pIntpData->m1);
//	std::cout<<"  done!\n";
//}
//
//void interpolate(pMesh theMesh, std::list<mEntity*> &leaves, int parentDepth, GetPFuncScalar* pArray_GetFunctions, SetPFuncScalar* pArray_SetFunctions, int numFields){
//	pMeshDataId EAI = MD_newMeshDataId("EAI"); 				// EAI = Edge Already Interpolated
//	for (list<mEntity*>::iterator lit = leaves.begin(); lit != leaves.end(); lit++){
//		pEntity face = *lit;
//		int childDepth = theMesh->getRefinementDepth(face);
//		if ( childDepth == parentDepth - 1){
//			// get child face's edges
//			// for each edge, take its parent (if it exist)
//			for (int i=0; i<3; i++){
//				pEdge edge = face->get(1,i);
//				pEdge edgeParent = edge->parent();
//				if ( edgeParent ){
//					int isEAI = 0;
//					EN_getDataInt(edgeParent,MD_lookupMeshDataId("EAI"),&isEAI);
//					if (!isEAI){
//						EN_attachDataInt(edgeParent,MD_lookupMeshDataId("EAI"),1);
//						// find node to receive an interpolated value
//						int ID0 = EN_id(edgeParent->get(0,0));				// edge parent IDs
//						int ID1 = EN_id(edgeParent->get(0,1));				// edge parent IDs
//						int id0 = EN_id(edge->get(0,0));					// edge child IDs
//						int id1 = EN_id(edge->get(0,1));					// edge child IDs
//						int ID = (id0 == ID0 || id0 == ID1)?id1:id0;		// check what is the new node ID
//						pVertex newNode = theMesh->getVertex(ID);			// get new node
//
//						// take node solution from parent nodes I and J
//						// interpolation done for all fields
//						for (int j=0; j<numFields; j++){
//							double sol_I = pArray_GetFunctions[j]((pVertex)edgeParent->get(0,0));
//							double sol_J = pArray_GetFunctions[j]((pVertex)edgeParent->get(0,1));
//							double v = (double)0.5*(sol_I + sol_J);				// interpolate
//							pArray_SetFunctions[j](newNode,v);								// set interpolated value to new node
//						}
//					}
//				}
//			}
//		}
//	}
//	// go to next depth: depth = depth - 1
//	if (--parentDepth){
//		interpolate(theMesh,leaves,parentDepth,pArray_GetFunctions,pArray_SetFunctions,numFields);
//	}
//}
//
////void unflagEdgesAsInterpolated(pMesh theMesh){
////	//pMeshDataId EAI = MD_newMeshDataId("EAI"); 				// EAI = Edge Already Interpolated
////	EIter eit = M_edgeIter(theMesh);
////	while (pEntity edge = EIter_next(eit)){
////		EN_attachDataInt(edge,MD_lookupMeshDataId("EAI"),0);
////	}
////	EIter_delete(eit);
////}
//
//
//
