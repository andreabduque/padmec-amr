/*
 * GeomData_initialize.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: rogerio
 */

#include "GeomData.h"

namespace PRS{
	void GeomData::calculateNumElements(pMesh theMesh){
		int dim = theMesh->getDim();
		int ndom = getNumDomains();

		numDomElem = new int[ndom];
		pEntity face, tetra;
		if (dim==2){
			for (int i=0; i<ndom; i++){
				int counter = 0;
				FIter fit = M_faceIter(theMesh);
				while ( (face = FIter_next(fit)) ){
					int flag1 = getFaceFlag(face);
					int flag2 = domainList[i];
					if (flag1==flag2){
						counter++;
					}
				}
				FIter_delete(fit);
				numDomElem[i] = counter;
				if (!numDomElem[i]){
					throw Exception(__LINE__,__FILE__,"Number of element: 0.");
				}
			}
		}
		else{
			for (int i=0; i<ndom; i++){
				int counter = 0;
				RIter rit = M_regionIter(theMesh);
				while ( (tetra = RIter_next(rit)) ){
					int flag1 = getTetraFlag(tetra);
					int flag2 = domainList[i];
					if (flag1==flag2){
						counter++;
					}
				}
				RIter_delete(rit);
				numDomElem[i] = counter;
				if (!numDomElem[i]){
					throw Exception(__LINE__,__FILE__,"Number of element: 0.");
				}
			}
		}
	}

	void GeomData::calculateNumNodes(pMesh theMesh){
		int ndom = getNumDomains();
		int dim = theMesh->getDim();
		setMeshNodes(M_numVertices(theMesh));
		numNodesPerDomain = new int[ndom];
		pEntity face, tetra;
		std::set<pEntity> nodeSet;
		for (int i=0; i<ndom; i++){
			if (dim==2){
				FIter fit = M_faceIter(theMesh);
				while ( (face = FIter_next(fit)) ){
					if (getFaceFlag(face)==domainList[i]){
						for (int j = 0; j<3; j++){
							nodeSet.insert(face->get(0,j));
						}
					}
				}
				FIter_delete(fit);
			}
			else{
				RIter rit = M_regionIter(theMesh);
				while ( (tetra = RIter_next(rit)) ){
					int tetra_flag = getTetraFlag(tetra);
					if (tetra_flag==domainList[i]){
						for (int j = 0; j<4; j++){
							nodeSet.insert(tetra->get(0,j));
						}
					}
				}
				RIter_delete(rit);
			}
			numNodesPerDomain[i] = (int)nodeSet.size();
			if (!numNodesPerDomain[i]){
				throw Exception(__LINE__,__FILE__,"Number of nodes: 0.");
			}
			nodeSet.clear();
		}
	}

	void GeomData::calculateNumBdryNodes(pMesh theMesh){
		int ndom = getNumDomains();
		int dim = theMesh->getDim();
		numBdryNodesPerDomain = new int[ndom];
		pEntity face, tetra;
		std::set<pEntity> nodeSet;
		for (int i=0; i<ndom; i++){
			if (dim==2){
				FIter fit = M_faceIter(theMesh);
				while ( (face = FIter_next(fit)) ){
					int faceflag = getFaceFlag(face);
					if (faceflag==domainList[i]){
						for (int j = 0; j<3; j++){
							if (getVertexFlag(face->get(0,j))!=faceflag){
								nodeSet.insert(face->get(0,j));
							}
						}
					}
				}
				FIter_delete(fit);
			}
			else{
				RIter rit = M_regionIter(theMesh);
				while ( (tetra = RIter_next(rit)) ){
					int tetra_flag = getTetraFlag(tetra);
					if (tetra_flag==domainList[i]){
						for (int j = 0; j<4; j++){
							if (getVertexFlag(tetra->get(0,j))!=tetra_flag){
								nodeSet.insert(tetra->get(0,j));
							}
						}
					}
				}
				RIter_delete(rit);
			}
			numBdryNodesPerDomain[i] = (int)nodeSet.size();
			if (!numBdryNodesPerDomain[i]){
				throw Exception(__LINE__,__FILE__,"Number of boundary nodes 0!");
			}
			nodeSet.clear();
		}
	}

	void GeomData::calculateNumEdges(pMesh theMesh){
		int dim = theMesh->getDim();
		int ndom = getNumDomains();
		numDomEdges = new int[ndom];
		//this->domainList = new int[ndom];
		//		for (int i=0; i<ndom; i++){
		//			this->domainList[i] = domainList[i];
		//		}
		pEntity face, tetra;
		std::set<pEntity> edgeList;
		for (int i=0; i<ndom; i++){
			if (dim==2){
				FIter fit = M_faceIter(theMesh);
				while ( (face = FIter_next(fit)) ){
					if (getFaceFlag(face)==domainList[i]){
						for (int j = 0; j<3; j++){
							edgeList.insert(face->get(1,j));
						}
					}
				}
				FIter_delete(fit);
			}
			else{
				RIter rit = M_regionIter(theMesh);
				while ( (tetra = RIter_next(rit)) ){
					int tetra_flag = getTetraFlag(tetra);
					if (tetra_flag==domainList[i]){
						for (int j = 0; j<6; j++){
							edgeList.insert(tetra->get(1,j));
						}
					}
				}
				RIter_delete(rit);
			}
			numDomEdges[i] = (int)edgeList.size();
			if (!numDomEdges[i]){
				throw Exception(__LINE__,__FILE__,"Number of edges: 0. Check if elements flags in mesh and physical.dat files match.");
			}
			edgeList.clear();
		}
	}

	void GeomData::calculateNumBDRYEdges(pMesh theMesh){
		if (theMesh->getDim()==3){
			return;
		}

		int ndom = getNumDomains();
		numDomBDRYEdges = new int[ndom];
		pEntity edge, face;
		std::set<pEntity> edgeList;
		for (int i=0; i<ndom; i++){
			FIter fit = M_faceIter(theMesh);
			while ( (face = FIter_next(fit)) ){
				if (getFaceFlag(face)==domainList[i]){
					for (int j = 0; j<3; j++){
						edge = (pEdge)face->get(1,j);
						//cout << "edge flag: " << getEdgeFlag(edge) << endl;
						if ( getEdgeFlag(edge)!= domainList[i]){
							edgeList.insert(edge);
						}
					}
				}
			}
			FIter_delete(fit);
			numDomBDRYEdges[i] = (int)edgeList.size();
			if (!numDomBDRYEdges[i]){
				throw Exception(__LINE__,__FILE__,"Number of boundary edges: 0. You must set the correct flag to ");
			}
			edgeList.clear();
		}
		numExtBdryEdges = 0;
		EIter eit = M_edgeIter(theMesh);
		while ( (edge=EIter_next(eit)) ){
			if (E_numFaces(edge)==1){
				numExtBdryEdges++;
			}
		}
		EIter_delete(eit);
	}

	void GeomData::calculateNumBDRYFaces(pMesh theMesh){
		if (theMesh->getDim()==2){
			return;	// nothing to do here
		}

		int ndom = getNumDomains();
		numDomBDRYFaces = new int[ndom];
		pEntity tetra, face;
		std::set<pEntity> faceList;
		for (int i=0; i<ndom; i++){
			RIter rit = M_regionIter(theMesh);
			while ( (tetra = RIter_next(rit)) ){
				if (getTetraFlag(tetra)==domainList[i]){
					for (int j = 0; j<4; j++){
						face = (pEntity)tetra->get(2,j);
						int faceflag = getFaceFlag(face);
						if ( faceflag!= domainList[i]){
							faceList.insert(face);
						}
					}
				}
			}
			RIter_delete(rit);
			numDomBDRYFaces[i] = (int)faceList.size();
			//cout << "numDomBDRYFaces = " << numDomBDRYFaces[i] << endl;
			if (!numDomBDRYFaces[i]){
				throw Exception(__LINE__,__FILE__,"Number of boundary faces: 0.");
			}
			faceList.clear();
		}
		//throw 1;

		// Calculate number of external boundary faces (triangles)
		numExtBdryFaces = 0;
		FIter fit = M_faceIter(theMesh);
		while ( (face=FIter_next(fit)) ){
			if (F_numRegions(face)==1){
				numExtBdryFaces++;
			}
		}
		FIter_delete(fit);
		if (!numExtBdryFaces){
			throw Exception(__LINE__,__FILE__,"Any external face detected!");
		}
	}

	void GeomData::calculateEdgeProperties(pMesh theMesh){
		int ndom = getNumDomains();
		double vec[3];
		double elength, Lij[3];
		double Icoord[3], Jcoord[3];
		double delta_x = 1e30;
		int i,j;
		int dim = theMesh->getDim();
		pEntity edge;

		//cout << "Num. Edges: " << M_numEdges(theMesh) << endl;
		for (int dom = 0; dom<ndom; dom++){
			int row = 0;
			int dom_flag = domainList[dom];
			EIter eit = M_edgeIter(theMesh);
			while ( (edge = EIter_next(eit)) ){
				if (edgeBelongToDomain(edge,dom_flag)){
					E_getVerticesCoord(edge,Icoord,Jcoord);
					for (j=0; j<dim; j++){
						Lij[j] = .0;
					}
					makeVector(Jcoord,Icoord,Lij);
					for (i=0; i<dim; i++){
						vec[i] = Lij[i];
					}
					elength = .0;
					for (i=0; i<dim; i++){
						elength += vec[i]*vec[i];
					}
					elength = sqrt(elength);
					if (elength == .0){
						char msg[256]; sprintf(msg,"Edge [%d %d] has null length!",EN_id(edge->get(0,0)),EN_id(edge->get(0,1)));
						throw Exception(__LINE__,__FILE__,msg);
					}
					edge_length[dom].setValue(row,elength);
					for (i=0; i<dim; i++){
						vec[i] /= elength;
					}
					edge_versor[dom].setValue(row,0,vec[0]);
					edge_versor[dom].setValue(row,1,vec[1]);
					edge_versor[dom].setValue(row,2,vec[2]);
					if (delta_x > elength){
						delta_x = elength;
					}
					row++;
				}
			}
			EIter_delete(eit);
		}
		setSmallestEdgeLength(P_getMinDbl(delta_x));

		// Calculate unit vector for external boundary edges
		int tnedges = 0;
		int row = 0;

		if (dim==2){
			EIter eit = M_edgeIter(theMesh);
			while ( (edge = EIter_next(eit)) ){
				tnedges++;
				if (E_numFaces(edge)==1){
					E_getVerticesCoord(edge,Icoord,Jcoord);
					for (j=0; j<dim; j++){
						Lij[j] = .0;
					}
					makeVector(Jcoord,Icoord,Lij);
					for (i=0; i<dim; i++){
						vec[i] = Lij[i];
					}
					elength = .0;
					for (i=0; i<dim; i++){
						elength += vec[i]*vec[i];
					}
					elength = sqrt(elength);

					for (i=0; i<dim; i++){
						vec[i] /= elength;
					}
					versor_ExtBdryElem[0].setValue(row,0,vec[0]);
					versor_ExtBdryElem[0].setValue(row,1,vec[1]);
					versor_ExtBdryElem[0].setValue(row,2,vec[2]);
					row++;
				}
			}
			EIter_delete(eit);
			this->setTotalNumberOfEdges(tnedges);
		}
	}

	// calculate a unit vector normal to each external face
	void GeomData::calculate_extFaceVersor(pMesh theMesh){
		int row = 0;
		double dij[3], norma;
		int ndom = getNumDomains();
		for (int i=0; i<ndom; i++){
			FIter fit = M_faceIter(theMesh);
			while ( pFace face = FIter_next(fit) ){

				// only external boudanry faces
				if (F_numRegions(face)==1){
					if (faceBelongToDomain(face,domainList[i])){
						getDij(face,domainList[i],dij);
						norma = sqrt( inner_product(dij,dij,3) );
						versor_ExtBdryElem[0].setValue(row,0,dij[0]/norma);
						versor_ExtBdryElem[0].setValue(row,1,dij[1]/norma);
						versor_ExtBdryElem[0].setValue(row,2,dij[2]/norma);
						row++;
					}
				}
			}
			FIter_delete(fit);
		}
	}

	void GeomData::calculate_NodeAverage_CDL(){
		int i,dom, row;
		const int* indices = NULL;

		// initialize CDL per node
		for(i=0; i<numNodes; i++){
			node_CDL[i] = 0;
		}

		// loop over elements
		int k = 0;
		for (dom=0; dom<_ndom; dom++){
			int nelements = numDomElem[dom];

			// loop over domains' elements
			for (row=0; row<nelements; row++){
				getElement(dom,row,indices);
				node_CDL[indices[3]] += elem_CDL[k];
				node_CDL[indices[4]] += elem_CDL[k];
				node_CDL[indices[5]] += elem_CDL[k];
			}
		}

		// average node CDL by number of sharing faces
		for(i=0; i<numNodes; i++){
			node_CDL[i] /= numElemSharingVertex[i];
		}
	}

	void GeomData::calculate_CDL(){
		int i,dom, row;
		const int* indices = NULL;

		// elements' coordinates
		const double* coords[4] = {NULL,NULL,NULL,NULL};

		// vector created over elements' edges
		int num_edges = 3*(dim-1);
		double* vectors[6];
		for (i=0;i<num_edges;i++){
			vectors[i] = new double[dim];
		}

		// indices to automate edge vectors creation
		int idx[6][2] = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};

		// Auxiliary variable
		int pos = dim+1;

		int k = 0;
		for (dom=0; dom<_ndom; dom++){
			int nelements = numDomElem[dom];

			// loop over domains' elements
			for (row=0; row<nelements; row++){
				getElement(dom,row,indices);

				// get vertices coordinates
				for (i=0; i<pos; i++){
					getCoordinates(indices[i+pos],coords[i]);
				}

				// create edge vectors
				for (i=0; i<num_edges; i++){
					const double* vtx1 = coords[ idx[i][0] ];
					const double* vtx2 = coords[ idx[i][1] ];
					for (int j=0; j<dim; j++){
						(vectors[i])[j] = vtx1[j] - vtx2[j];
					}
				}

				// vectors lengths summation
				for (i=0; i<num_edges; i++){
					elem_CDL[k] = sqrt( inner_product(vectors[i],vectors[i],dim) );
				}

				// take average lenght
				elem_CDL[k] /= pos;
				k++;
			}
		}
		for (i=0;i<num_edges;i++){
			delete[] vectors[i]; vectors[i] = 0;
		}
	}

	void GeomData::calculateNumElemSharingVertex(){
		int i,j,dom;
		const int *idx = NULL;

		// initialize
		for (i=0; i<numNodes; i++){
			numElemSharingVertex[i] = 0;
		}

		/* ----------------------------------------------------------
					 2-D (triangle)			   3-D (tetrahedra)
				idx[0]	: local index		idx[0]	: local index
				idx[1]	: local index		idx[1]	: local index
				idx[2]	: local index		idx[2]	: local index
											idx[3]	: global index

				idx[3]	: global index		idx[4]	: global index
				idx[4]	: global index		idx[5]	: global index
				idx[5]	: global index		idx[6]	: global index
											idx[7]	: global index
				---------------------------------------------------- */
		int pos1 = dim+1;
		int pos2 = pos1+dim;
		for (dom=0; dom<_ndom; dom++){
			for (i=0; i<numDomElem[dom]; i++){
				idx = elem[dom].getrowconst(i);
				for (j=pos1; j<=pos2; j++){
					numElemSharingVertex[ idx[j] ]++;
				}
			}
		}
	}
}
