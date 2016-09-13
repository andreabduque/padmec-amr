/*
 * GeomData_initialize.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: rogerio
 */

#include "GeomData.h"

namespace PRS{

	/*
	 * Based on pre-processor data, matrix structure is created and data transferred
	 */
	void GeomData::initilize(pMesh theMesh, const std::set<int> &setOfDomain, int FVM){
		// number of domains stay the same no matter how many adaptations happen
		int i = 0;
		int ndom = (int)setOfDomain.size();
		int *domainList = new int[ndom];
		for(set<int>::iterator iter = setOfDomain.begin(); iter!=setOfDomain.end(); iter++){
			domainList[i++] =  *iter;
		}
		setNumDomains(ndom);
		setDomainList(domainList);
		delete[] domainList; domainList = 0;

//		if (theMesh->getDim()==2){
//			numElem = M_numFaces(theMesh);
//			elemtype = 3;
//		}
//		else{
//			numElem = M_numRegions(theMesh);
//			elemtype = 4;
//		}

		// FVM:
		// CASE 1: classical edge based Finite Volume Formulation
		// CASE 2: modified  CASE 1 (highly heterogeneous porous media)
		if (FVM==1){
			initilize(theMesh);				//every new mesh adaptation, mesh data structure changes
		}
	}

	void GeomData::initilize(pMesh theMesh){
		if (theMesh->getDim()==2){
			numElem = M_numFaces(theMesh);
			elemtype = 3;
		}
		else{
			numElem = M_numRegions(theMesh);
			elemtype = 4;
		}
		calculateNumEdges(theMesh);					// calculate number of data to be stored
		calculateNumElements(theMesh);
		calculateNumBDRYEdges(theMesh);
		calculateNumBDRYFaces(theMesh);
		calculateNumNodes(theMesh);
		calculateNumBdryNodes(theMesh);
		allocatePointers(M_numVertices(theMesh),theMesh->getDim());	// allocate storage
		calculateEdgeProperties(theMesh);							// fill storage
		dataTransfer(theMesh);
		calculate_extFaceVersor(theMesh);
		mappingNodesIds(theMesh);									// map data to find them quickly
	}
}
