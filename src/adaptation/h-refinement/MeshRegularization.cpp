///*
// * MeshRegularization.cpp
// *
// *  Created on: 03/06/2012
// *      Author: rogsoares
// */
//
//#include "H_Refinement_2D.h"
//
//void H_Refinement_2D::meshRegularization(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, int min, std::list<pEntity>& faceList){
//	pFace faces[3];																				// faces around face
//	int numFaces;																				// Number of face's neighbors
//	list<mEntity*> subTree;																		// stores all face leaves
//	list<pEntity>::iterator sit;
//	bool flagModification = false;
//
//	int c = 0;
//	FIter fit = M_faceIter(theMesh);
//	while (pEntity face = FIter_next(fit)){
//		//get elements without children: Depth = 0
//		if (!theMesh->getRefinementDepth(face)){
//			int count = 0;																		// what is 'count' for????
//			int faceLevel = pErrorAnalysis->getLevelOfRefinement(face);							// level to refine of face
//			getFacesAroundFace(face,faces,numFaces);											// get faces around face
//			for(int i=0; i<numFaces; i++){
//				int adjElemLevel = pErrorAnalysis->getLevelOfRefinement(faces[i]);				// neighbor's face ref level
//				if (  faceLevel == adjElemLevel - 1 ) {
//					count++;
//				}
//			}
//			if (count>1){
//				pErrorAnalysis->setLevelOfRefinement(face,faceLevel+1);
//				flagModification = true;
//				c++;
//			}
//		}
//	}
//	FIter_delete(fit);
//
//	//cout << "Num c: " << c << endl;
//	//throw 1;
//	/*
//	 * If any element flag for refinement has been modified, then the uneven procedure must be performed again.
//	 */
//	if (flagModification){
//		unevenElements(pErrorAnalysis,theMesh,faceList,4,min);
//	}
//}
//
//bool H_Refinement_2D::leavesCorrectionMain(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, int maxDepth, list<pEntity> &faceList){
//	list<pEntity>::iterator iter;
//	int depth = 0;
//	bool correction = false;
//
//	//mark from children to parent
//	while (depth < maxDepth){
//		for (iter = faceList.begin(); iter != faceList.end(); iter++)
//			firstLeavesCorrection(pErrorAnalysis,theMesh,*iter,depth,correction);
//		depth++;
//	}
//	depth = maxDepth;
//
//	//mark from parent to children
//	while (depth > 0){
//		for (iter = faceList.begin(); iter != faceList.end(); iter++)
//			secondLeavesCorrection(pErrorAnalysis,theMesh,*iter,depth,correction);
//		depth--;
//	}
//	return correction;
//}
//
//void H_Refinement_2D::firstLeavesCorrection(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, pEntity parent, int depth, bool& correction){
//	//EN_getDataInt ((pFace)parent , MD_lookupMeshDataId ("mRef"), &parentLevel);
//	int parentLevel = pErrorAnalysis->getLevelOfRefinement(parent);
//	if ( !theMesh->getRefinementLevel(parent) && !theMesh->getRefinementDepth(parent) && parentLevel<0 ){
//		//EN_attachDataInt ( (pFace) parent, MD_lookupMeshDataId ("mRef"),0);
//		pErrorAnalysis->setLevelOfRefinement(parent,0);
//		correction = true;
//	}
//	else if (theMesh->getRefinementDepth(parent) == (depth+1)){
//		std::vector<pFace> fchildren;
//		//get children (2 children -> special element from hangnode treatment
//		//or 4 children->FMDB refinement)
//		//---------------------------------------------------------
//		getChildren(theMesh,parent,fchildren);
////		int *array = new int[fchildren.size()];
////		for (int i=0; i<fchildren.size(); i++){
////			int nivel;
////			EN_getDataInt ((pFace)fchildren[i] , MD_lookupMeshDataId ("mRef"), &nivel);
////			array[i] = nivel;
////		}
////		int maxLevel = *max_element(array,array+fchildren.size());
//
//		int fc_size = (int)fchildren.size();
//		int maxLevel = -1000;
//		int minLevel = 1000;
//		int fchildrenLevels[fc_size];
//		for (int i=0; i<fc_size; i++){
//			fchildrenLevels[i] = pErrorAnalysis->getLevelOfRefinement(fchildren[i]);
//			maxLevel = std::max(fchildrenLevels[i],maxLevel);
//			minLevel = std::min(fchildrenLevels[i],minLevel);
//		}
//
//		//case 1: all brothers with negative mark
//		//----------------------------------------------
//		if (maxLevel < 0){
//			int childLevel = std::min(abs(maxLevel),theMesh->getRefinementLevel(fchildren[0]));
//			//Loop nos filhos e marcacao de filhos vai receber nivel min(nivel filho,level do parent)
//			//----------------------------------------------
//			for (int i=0; i<fchildren.size(); i++){
//				//EN_attachDataInt ( (pFace)fchildren[i], MD_lookupMeshDataId ("mRef"), -childLevel );
//				pErrorAnalysis->setLevelOfRefinement(fchildren[i],-childLevel);
////				int nivelP = (1 - childLevel);
////				EN_attachDataInt ( (pFace) parent, MD_lookupMeshDataId ("mRef"), nivelP );
//				pErrorAnalysis->setLevelOfRefinement(parent,1 - childLevel);
//				if ( fabs(fchildrenLevels[i]) != childLevel )
//					correction = true;
//			}
//		}
//		else{
////			// int nivel_min = *min_element (fchildrenLevels,fchildrenLevels+fchildren.size());
//			//case 2: brothers with negative and positive mark
//			if (minLevel < 0){
//				for (int i=0; i<fc_size; i++){
//					int childLevel = std::max(fchildrenLevels[i],0);
//					//EN_attachDataInt ((pFace) fchildren[i], MD_lookupMeshDataId ("mRef"), nivelSon );
//					pErrorAnalysis->setLevelOfRefinement(fchildren[i],childLevel);
//					//EN_attachDataInt ((pFace) parent, MD_lookupMeshDataId ("mRef"), 1);
//					pErrorAnalysis->setLevelOfRefinement(parent,1);
//					if ( fchildrenLevels[i] != childLevel )
//						correction = true;
//				}
//			}
//			//case 3: all brothers with positive mark
//			else{
//				//EN_attachDataInt ((pFace) parent, MD_lookupMeshDataId ("mRef"),depth+1);
//				pErrorAnalysis->setLevelOfRefinement(parent,depth+1);
//			}
//		}
//		fchildren.clear();
//	}
//}
//
///*
// *	Pass parent's flag to its leaves
// */
//void H_Refinement_2D::secondLeavesCorrection(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, pEntity parent, int depth, bool& correction){
//	if (theMesh->getRefinementDepth(parent) == depth){
//		std::vector<pFace> fchildren;
//		getChildren (theMesh, parent,fchildren);
//		int parentLevel;
//		parentLevel = pErrorAnalysis->getLevelOfRefinement(parent);
//		//EN_getDataInt ( (pFace) parent, MD_lookupMeshDataId ("mRef"), &parentLevel );
//		for (int i=0; i<fchildren.size(); i++){
////			int levelF;
////			EN_getDataInt ( (pFace) fchildren[i], MD_lookupMeshDataId ("mRef"), &levelF );
//			int childLevel = pErrorAnalysis->getLevelOfRefinement(fchildren[i]);
//			int diff = childLevel - parentLevel;
//			if ( (childLevel < 0) && (diff!= -1) ){
//				//EN_attachDataInt ( (pFace) fchildren[i], MD_lookupMeshDataId ("mRef"), (parentLevel-1) );
//				pErrorAnalysis->setLevelOfRefinement(fchildren[i],parentLevel - 1);
//				correction = true;
//			}
//		}
//		fchildren.clear();
//	}
//}
//
