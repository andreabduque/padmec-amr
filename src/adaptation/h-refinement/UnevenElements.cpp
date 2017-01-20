///*
// * UnevenElements.cpp
// *
// *  Created on: 16/05/2012
// *      Author: rogsoares
// */
//
//#include "H_Refinement_2D.h"
//
///*
// * Strategy: Degree of refinement (dor) difference must be equal to 1 for any two adjacent mesh element.
// * It's must be imposed to guarantee smoother mesh.
// *
// *  unevenElements_part1: for an element flagged with the highest dor, elements around it will be flagged with dor-1
// *  unevenElements_part2: correct A's element dor which is closed by two adjacent elements with dor+1. In this case, A element
// *  					  will be set to dor+1
// *
// *  After unevenElements_part2 has been called, inconsistencies on 'dor' uneven can be still presented on mesh. So parts 1 and 2
// *  must be called more times until it has disappeared.
// */
//void H_Refinement_2D::unevenElements(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, std::list<pEntity>& faceList, int max, int min){
////	list<pEntity>::iterator iter;
////	list<mEntity*>::iterator sit;
////
////	int n;
////	std::vector<pEntity> neighboursFace(3);
////
//////	int count = 0;
//////	bool HangNodeEdge;
////
////	FIter fit = M_faceIter (theMesh);
////	while (pEntity face = FIter_next (fit) ){
////		// get elements without children: Depth = 0
////		if (!theMesh->getRefinementDepth(face)){
////			int markFace, levelFace;
////			markFace = pErrorAnalysis->getLevelOfRefinement(face);
////			levelFace = theMesh->getRefinementLevel(face);
////
////			//funcao recursiva: elementos marcados com maximo nivel ate min nivel
////			if (markFace == max){
////				pEntity root = face->root();
////				if ( root ){
////					list<mEntity*> subTree;
////					root->getAllSubTree(subTree);
////					for (sit = subTree.begin();sit!=subTree.end();sit++){
////						faceList.push_back(*sit);
////					}
////					subTree.clear();
////				}
////				//std::vector <pFace> fvector;
////				cout << __LINE__ << endl;
////				//neighboursFace(theMesh,face, neighboursFace);
////				//printFaceIDs(face);
////				getNeighboursFace(face,neighboursFace,n);
////				cout << __LINE__ << endl;
////
////				for (int i=0; i<n; i++){
////					pFace viz = neighboursFace[i];
////					int markViz = pErrorAnalysis->getLevelOfRefinement(viz);
////					int levelViz = theMesh->getRefinementLevel(viz);
////					//cout << "markViz " << pErrorAnalysis->getLevelOfRefinement(viz)<< endl;
////
////					// Variavel que verifica alteracao na marcacao do elemento e o inclui na fila
////					bool addList = false;
////					if (levelFace == levelViz && markViz < (markFace-1)){
////						pErrorAnalysis->setLevelOfRefinement(viz,markFace-1);
////						addList = true;
////					}
////					else if ( (levelFace == levelViz + 1) && (markViz < markFace) ) {
////						pErrorAnalysis->setLevelOfRefinement(viz,markFace);
////						addList = true;
////					}
////					else if (markViz < (markFace-1)){
////						pErrorAnalysis->setLevelOfRefinement(viz,markFace-1);
////						addList = true;
////					}
////
////					//Compare level between element and neighbors
////					if (addList){
////						// force neighbor's face refinement level flag be one unity lesser than face
////
////						//element with children: add subtree (include element itself)
////						pEntity root = face->root();
////
////						if ( root ){
////							list<mEntity*> subTree;
////							root->getAllSubTree(subTree);
////							for (sit = subTree.begin();sit!=subTree.end();sit++){
////								faceList.push_back(*sit);
////							}
////							subTree.clear();
////						}
////						//element without children: add element itself
////						else{
////							faceList.push_back(face);
////						}
////					}
////				}
////				//fvector.clear();
////			}
////		}
////	}
////	FIter_delete (fit);
////	if ((max - min) > 1){
////		unevenElements(pErrorAnalysis,theMesh,faceList,max-1,min);
////	}
////
////	meshRegularization(pErrorAnalysis,theMesh,min,faceList);
//}
//
////void H_Refinement_2D::unevenElements_part1(ErrorAnalysis* pErrorAnalysis, pMesh theMesh){
////	int dor; 					// degree of refinement
////	int neighbor_dor;			// face's neighbor degree of refinement
////	int numFaces;				// Number of face's neighbors
////	pFace faces[3];				// Faces may have one, two or three neighbor(s).
////	for (dor=4; dor>=0; dor--){
////		FIter fit = M_faceIter(theMesh);
////		while (pEntity face = FIter_next(fit)){
////			// get final leaves only
////			if (!theMesh->getRefinementDepth(face) && pErrorAnalysis->getLevelOfRefinement(face)==dor){
////				// get faces around face
////				getFacesAroundFace(face,faces,numFaces);
////				for (int i=0; i<numFaces; i++){
////					neighbor_dor = pErrorAnalysis->getLevelOfRefinement(faces[i]);
////					// guarantee unitary unevenness
////					int diff = fabs(neighbor_dor - dor);
////					if ( diff > 1){
////						int setlevel = fabs(dor-1); // modulo deve ser retirado
////						pErrorAnalysis->setLevelOfRefinement(faces[i],setlevel);
////					}
////				}
////			}
////		}
////		FIter_delete(fit);
////	}
////}
////
////void H_Refinement_2D::unevenElements_part2(ErrorAnalysis* pErrorAnalysis, pMesh theMesh){
////	int dor; 					// degree of refinement
////	int neighbor_dor;			// face's neighbor degree of refinement
////	int numFaces;				// Number of face's neighbors
////	pFace faces[3];				// Faces may have one, two or three neighbor(s).
////	for (dor=3; dor>=0; dor--){
////		FIter fit = M_faceIter(theMesh);
////		while (pEntity face = FIter_next(fit)){
////			// get final leaves only
////			if (!theMesh->getRefinementDepth(face) && pErrorAnalysis->getLevelOfRefinement(face)==dor){
////				// get faces around face
////				getFacesAroundFace(face,faces,numFaces);
////				// count how many faces around face whose refinement degree flag is dor+1
////				int count = 0;
////				for (int i=0; i<numFaces; i++){
////					neighbor_dor = pErrorAnalysis->getLevelOfRefinement(faces[i]);
////					if (neighbor_dor == dor+1) count++;
////				}
////				if (count>=2){
////					int setlevel = dor+1;
////					pErrorAnalysis->setLevelOfRefinement(face,setlevel);
////				}
////			}
////		}
////		FIter_delete(fit);
////	}
////}
//
//int H_Refinement_2D::checkMaxDORdifference(ErrorAnalysis* pErrorAnalysis, pMesh theMesh){
//	// max_diff = max degree of freedom difference between two adjacent elements
////	int max_diff = -1000;
////	int numFaces;				// Number of face's neighbors
////	pFace faces[3];
////
////	FIter fit = M_faceIter(theMesh);
////	while (pEntity face = FIter_next(fit)){
////		// get final leaves only
////		if (!theMesh->getRefinementDepth(face)){
////			int dor = pErrorAnalysis->getLevelOfRefinement(face);
////			// get faces around face
////			getFacesAroundFace(face,faces,numFaces);
////			for (int i=0; i<numFaces; i++){
////				int neighbor_dor = pErrorAnalysis->getLevelOfRefinement(faces[i]);
////				int diff = fabs(neighbor_dor - dor);
////				if ( diff > max_diff )
////					max_diff = diff;
////			}
////		}
////	}
////	FIter_delete(fit);
////	return max_diff;
//}
//
//// Essa funcao nao está geral. Nao funcionaria para 3D
////void neighboursFace(pMesh theMesh, pEntity face, std::vector <pFace> &fvector){
////	int ID1, ID2;
////	int ID1viz, ID2viz,ID3viz;
////	std::vector <pFace> neighboursFace;
////	std::vector<pVertex> vertices;
////
////	// Percorrendo arestas da face
////	for (int i=0; i<3; i++){
////		cout << __LINE__ << endl;
////		pEdge edge = F_edge(face,i);
////		getEdgesNodeIDs(edge,ID1,ID2);
////		cout << __LINE__ << endl;
////		if (isNotEdgeOnBdry(edge)){
////
////			// Checar se adjacencias da aresta foram criadas
////			if (!edge->isAdjacencyCreated(2)){
////				throw Exception(__LINE__,__FILE__,Exception::UNKNOWN_EXIT);
////			}
////
////			//Sempre vai retornar duas faces
////			cout << __LINE__ << endl;
////			neighboursFaceofEdge(theMesh, edge, neighboursFace);
////			cout << __LINE__ << endl;
////
////			for (int j=0; j<2; j++){
////				cout << __LINE__ << endl;
////				if (neighboursFace[j]!=face){
////					cout << __LINE__ << endl;
////				// Considerar vizinhos apenas das arestas que nao foram refinadas
////					if (theMesh->getRefinementDepth(edge)==0){
////						cout << __LINE__ << endl;
////						cout << "neighboursFace[j]: " << neighboursFace[j] << endl;
////						M_GetVertices(neighboursFace[j],vertices);
////						cout << __LINE__ << endl;
////						ID1viz = EN_id(vertices[0]);
////						cout << __LINE__ << endl;
////						ID2viz = EN_id(vertices[1]);
////						cout << __LINE__ << endl;
////						ID3viz = EN_id(vertices[2]);
////						cout << __LINE__ << endl;
////						cout << __LINE__ << endl;
////						cout << "IDSvizinhos1  " << ID1viz << " " << ID2viz << " "<< ID3viz <<endl;
////						vertices.clear();
////						neighboursFace.clear();
////						cout << __LINE__ << endl;
////						fvector.push_back(neighboursFace[j]);
////						cout << __LINE__ << endl;
////					}
////
////				// Para arestas que possuam filhos
////					else if (theMesh->getRefinementDepth(edge) == 1){
////						cout << __LINE__ << endl;
////						std::vector<pFace> echildren;
////						getChildren (theMesh, edge,echildren);
////						cout << "criancas da aresta: "<< echildren.size() << endl;
////
////						for (int k=0; k<2; k++){
////							neighboursFaceofEdge(theMesh, echildren[k], neighboursFace);
////
////							M_GetVertices(echildren[k],vertices);
////							ID1viz = EN_id(vertices[0]);
////							ID2viz = EN_id(vertices[1]);
////							cout << "********************NUM FACE " << E_numFaces(echildren[k]) << endl;
////							cout << "********************IDSarestas  " << ID1viz << " " << ID2viz << endl;
////							vertices.clear();
////							neighboursFace.clear();
////
////							for (int m=0; m<2; m++){
////								if (neighboursFace[m] !=face){
////									M_GetVertices(neighboursFace[m],vertices);
////									ID1viz = EN_id(vertices[0]);
////									ID2viz = EN_id(vertices[1]);
////									ID3viz = EN_id(vertices[2]);
////									cout << "IDSvizinhos2  " << ID1viz << " " << ID2viz << " "<< ID3viz <<endl;
////									vertices.clear();
////									neighboursFace.clear();
////									fvector.push_back(neighboursFace[m]);
////								}
////							}
////						}
////						echildren.clear();
////					}
////				}
////			}
////		}
////
////	}
////}
//
//// Retorna faces vizinhas de uma aresta
//void neighboursFaceofEdge(pMesh theMesh, pEntity edge, std::vector <pFace> &fvector){
//	// Verificando as duas faces vizinhas da aresta
//	//if (isNotEdgeOnBdry(edge)){
//	// CUIDADO COM A ARESTA DE CONTORNO!!!!!!!!
////	fvector.clear();
////		for (int i=0; i< 2; i++){
////			pFace viz = E_face (edge,i);
////			fvector.push_back(viz);
////		}
//	//}
//}
//
////Versao Erika
///*void neighboursFace(pMesh theMesh, pEntity face, std::vector <pFace> &fvector){
//	int ID1, ID2,ID3;
//	pEntity vizRef;
//	pEntity parentViz;
//	for (int i=0; i<3; i++){
//		//cout << __LINE__ << endl;
//		pEdge edge = F_edge(face,i);
//		//cout << __LINE__ << endl;
//		getEdgesNodeIDs(edge,ID1,ID2);
//		if (isNotEdgeOnBdry(edge)){
//			if (!edge->isAdjacencyCreated(2)){
//				// Tratamento
//				//cout << __LINE__ << endl;
//				mEdge* edgesChildren[2];
//				getEdgesChildren(edge,edgesChildren);
//				pFace viz_1 = E_face (edgesChildren[1],0);
//				pFace viz_2 = E_face (edgesChildren[1],1);
//
//				if (viz_1!= face)
//					vizRef = viz_1;
//				else if(viz_2!= face)
//					vizRef = viz_2;
//
//				parentViz = vizRef->parent();
//
//				fvector.push_back(parentViz);
//				std::vector<pVertex> vertices;
//				M_GetVertices(parentViz,vertices);
//				int ID1viz, ID2viz,ID3viz;
//				ID1viz = EN_id(vertices[0]);
//				ID2viz = EN_id(vertices[1]);
//				ID3viz = EN_id(vertices[2]);
//				//cout << "IDSvizinhos  " << ID1viz << " " << ID2viz << " "<< ID3viz<<endl;
//				vertices.clear();
//
//				//throw Exception(__LINE__,__FILE__,Exception::UNKNOWN_EXIT);
//			}
//			else {
//				//cout << "Antes vizinhos" << endl;
//				pFace viz_1 = E_face (edge,0);
//				//cout << "Depois vizinhos" << endl;
//				pFace viz_2 = E_face (edge,1);
//				//cout << "Depois vizinhos2" << endl;
//
//				if (viz_1!= face){
//					fvector.push_back(viz_1);
//
//					//cout << "Depois vizinhos" << endl;
//
//					std::vector<pVertex> vertices;
//					M_GetVertices(viz_1,vertices);
//					int ID1viz, ID2viz,ID3viz;
//					ID1viz = EN_id(vertices[0]);
//					ID2viz = EN_id(vertices[1]);
//					ID3viz = EN_id(vertices[2]);
//					//cout << "IDSvizinhos  " << ID1viz << " " << ID2viz << " "<< ID3viz<<endl;
//					vertices.clear();
//				}
//				else if(viz_2!= face){
//					fvector.push_back(viz_2);
//					//cout << "Depois vizinhos" << endl;
//
//					std::vector<pVertex> vertices;
//					M_GetVertices(viz_2,vertices);
//
//					int ID1viz, ID2viz,ID3viz;
//					ID1viz = EN_id(vertices[0]);
//					ID2viz = EN_id(vertices[1]);
//					ID3viz = EN_id(vertices[2]);
//					cout << "IDSvizinhos  " << ID1viz << " " << ID2viz << " "<< ID3viz<<endl;
//					vertices.clear();
//
//				}
//			}
//		}
//	}
//}*/
//
//void getEdgesNodeIDs(pEdge edge, int &ID1, int &ID2){
//	std::vector<pVertex> vertices;
//	M_GetVertices(edge,vertices);
//	ID1 = EN_id(vertices[0]);
//	ID2 = EN_id(vertices[1]);
//	vertices.clear();
//}
//
//bool isNotEdgeOnBdry(pEdge edge){
//	return ( E_numFaces(edge) != 1);
//}
//
//void getChildren (pMesh theMesh, pEntity parent, std::vector<pFace> &fchildren){
//	list<mEntity*>::iterator lit;
//	list<mEntity*> leaves;
//	//parent was refined only once
//	if (theMesh->getRefinementDepth(parent) == 1){
//		parent->getLeaves(leaves);
//		for (lit = leaves.begin(); lit!=leaves.end(); lit++)
//			fchildren.push_back(*lit);
//	}
//	//parent was refined more than once
//	else if (theMesh->getRefinementDepth(parent) > 1){
//		parent->getAllSubTree(leaves);
//		for (lit = leaves.begin(); lit != leaves.end(); lit++){
//			//get only the parent children (difference between level is 1)
//			int levelDiff = theMesh->getRefinementLevel(*lit) - theMesh->getRefinementLevel(parent);
//			if ( levelDiff == 1 ){
//				fchildren.push_back(*lit);
//			}
//		}
//	}
//	leaves.clear();
//}
