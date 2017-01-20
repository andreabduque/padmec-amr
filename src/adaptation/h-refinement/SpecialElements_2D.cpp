///*
// * SpecialElements_2D.cpp
// *
// *  Created on: 14/05/2012
// *      Author: rogsoares
// */
//
//#include "H_Refinement_2D.h"
//
//void H_Refinement_2D::removeSpecialElements(pMesh theMesh, ErrorAnalysis *pErrorAnalysis){
//
//#ifdef __ADAPTATION_DEBUG__
//	if (!SE_list_2D.size()){
//		throw Exception(__LINE__,__FILE__,"List of special elements is empty! Maybe you have forgotten to create it before perform the unitary element uneveness.\n");
//	}
//#endif
//
//	/*
//	 * loop over special elements list to remove their respective children (originated by hang-nodes)
//	 * It's faster than looping over the entire element tree
//	 * Before killing leaves, get their Refinement level from Error Analysis and flag parent with highest of them.
//	 */
//	for (SE_iter_2D = SE_list_2D.begin(); SE_iter_2D != SE_list_2D.end(); SE_iter_2D++){
//		mFace* face = *SE_iter_2D;
//
//#ifdef __ADAPTATION_DEBUG__
//		int old_depth = theMesh->getRefinementDepth(face);
//#endif
//
//		// get special element children
//		mFace *f1 = (mFace*)face->get(2,0);
//		mFace *f2 = (mFace*)face->get(2,1);
//
//		// get max refinement level flag from leaves and pass it to parent
//		int rlevel = std::max(pErrorAnalysis->getLevelOfRefinement(f1),pErrorAnalysis->getLevelOfRefinement(f2));
//		pErrorAnalysis->setLevelOfRefinement(face,rlevel);
//
//		// cut elements connection with their parent
//		f1->deleteParent();
//		f2->deleteParent();
//
//		mEdge *edge = 0;
//
//		// get common edge between children
//		edge = f1->commonEdge(f2);
//
//		if (!edge){
//			throw Exception(__LINE__,__FILE__,"Null edge!");
//		}
//
//		// delete: common edge and element children
//		theMesh->DEL(edge);
//		theMesh->DEL(f1);
//		theMesh->DEL(f2);
//
//		// delete adjacencies among special element and its children
//		face->deleteAdjacencies(2);
//		this->setElementAsSpecial(face,false);
//
////#ifdef __ADAPTATION_DEBUG__
////		// A special element must have depth=1 and depth=0 before and after removing its two leaves
////		int new_depth = theMesh->getRefinementDepth(face);
////		char msg[256];
////		int IDs[3];
////
////		if (!old_depth || !new_depth){
////			getTriVerticesIDs(face,IDs);
////		}
////		if (!old_depth){
////			sprintf(msg,"Special element: %d %d %d. Refinement depth different to 1 before removing its leaves.",IDs[0],IDs[1],IDs[2]);
////			throw Exception(__LINE__,__FILE__,msg);
////		}
////		if (!new_depth){
////			sprintf(msg,"Special element: %d %d %d. Refinement depth different to 0 after removing its leaves.",IDs[0],IDs[1],IDs[2]);
////			throw Exception(__LINE__,__FILE__,msg);
////		}
////#endif
//	}
//	SE_list_2D.clear();
//}
//
//void H_Refinement_2D::createSpecialElementsList(pMesh theMesh, ErrorAnalysis *pErrorAnalysis){
//#ifdef __ADAPTATION_DEBUG__
//	cout << "createSpecialElementsList:   START\n";
//#endif
//	/*
//	 * Loop over all elements (triangles) and find those with hang-nodes.
//	 * For 2-D, elements with hang-nodes have only one edge with two children.
//	 */
//	FIter fit = M_faceIter (theMesh);
//	while (pEntity face = FIter_next(fit)){
//		// special elements have any children (for while)
//		if (!theMesh->getRefinementDepth(face)){
//
//			// A special element is found when the sum of refinement depths of all its three edges is equal 1.
//			int sum = 0;
//			for (int i=0; i<3; i++){
//				sum += theMesh->getRefinementDepth(face->get(1,i));
//			}
//			if (sum==1){
//				SE_list_2D.push_back((mFace*)face);
//				this->setElementAsSpecial(face,true);
//			}
//		}
//	}
//	FIter_delete(fit);
//	if (!SE_list_2D.size()){
//		cout << "WARNING: Any element with hang nodes were found.\n";
//	}
//
//#ifdef __ADAPTATION_DEBUG__
//	cout << "createSpecialElementsList:   END\n";
//#endif
//}
//
//void H_Refinement_2D::eliminateHangNodes(pMesh theMesh){
//	mEdge* edgesChildren[2];
//	// edgeIndices help to get face's edges without children (the third one contains the hang-node).
//	// Ex.: face has three edges: edge[0], edge[1]. edge[2]. If edge[1] has the hangnode, then
//	// edge[0] and edge[2] don't have children, i.e, i = 1 -> {0,2}
//	int edgeIndices[3][2] = {{1,2},{0,2},{0,1}};
//
//	if (!SE_list_2D.size()){
//		throw Exception(__LINE__,__FILE__,"Special element list is empty!");
//	}
//	for (SE_iter_2D = SE_list_2D.begin(); SE_iter_2D != SE_list_2D.end(); SE_iter_2D++){
//		mFace* face = *SE_iter_2D;
//		if (!face->isAdjacencyCreated(2)){
//			mEdge* edges[3] = {(mEdge*)F_edge(face,0), (mEdge*)F_edge(face,1), (mEdge*)F_edge(face,2)};
//			// get the edge with two children. The following loop will run just once
//			for (int i=0; i<3; i++){
//				if (theMesh->getRefinementDepth(edges[i])==1){
//					// get edge's children
//					getEdgesChildren(edges[i],edgesChildren);
//					pVertex hangNode = edgesChildren[0]->commonVertex(edgesChildren[1]);
//					// edge_0 and edge_1 have any children
//					mEdge* edge_0 = (mEdge*)edges[ edgeIndices[i][0] ];
//					mEdge* edge_1 = (mEdge*)edges[ edgeIndices[i][1] ];
//					pVertex v0 = edge_0->commonVertex(edge_1);
//					pVertex v1 = edge_0->commonVertex(edges[i]);
//					pVertex v2 = edge_1->commonVertex(edges[i]);
//
//					// give order to edgesLeaves to make new faces further
//					if ( edge_0->commonVertex(edgesChildren[0]) != v1)
//						std::swap(edgesChildren[0],edgesChildren[1]);
//
//					// create a new edge, splitting the special element into two connecting the hangnode to vertex v.
//					mEdge* newEdge = theMesh->createEdge((mVertex*)hangNode,(mVertex*)v0,face->getClassification());
//
//					// create two new elements face
//					mFace *f1 = theMesh->createFaceWithVertices ((mVertex*)hangNode,(mVertex*)v1,(mVertex*)v0,face->getClassification());
//					mFace *f2 = theMesh->createFaceWithVertices ((mVertex*)hangNode,(mVertex*)v2,(mVertex*)v0,face->getClassification());
//
//					// set face element as f1 and f2 father
//					f1->setParent(face);
//					f2->setParent(face);
//
//					// tell f1 and f2 who are their edges
//					f1->add(newEdge);
//					f1->add(edge_0);
//					f1->add(edgesChildren[0]);
//					f2->add(newEdge);
//					f2->add(edge_1);
//					f2->add(edgesChildren[1]);
//
//					// tell face that it has now two children
//					face->add(f2);
//					face->add(f1);
//					f1->setCommonBdry(face->getCommonBdry());
//					f2->setCommonBdry(face->getCommonBdry());
//				}
//			}
//		}
//	}
//}
