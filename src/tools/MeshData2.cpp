/*
 * MeshData2.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: rogerio
 */

#include "MeshData.h"

namespace PRS{

	void MeshData::createGlobalNodeIDMapping(){
		std::set<pEntity> dirichletnodes_set;
		int freeRows_size, dirichletCols_size;
		int* freeRows_array;
		int* dirichletCols_array;

		mapNodeID();
		setGlobalIndices();
		getArrays_free(freeRows_size,&freeRows_array);
		getArrays_dirichlet(dirichletCols_size,&dirichletCols_array,dirichletnodes_set);
		initDirichletVectors(dirichletnodes_set);

		ISCreateGeneral(MPI_COMM_WORLD,freeRows_size,freeRows_array,PETSC_COPY_VALUES,&IS_freeRows);
		ISCreateGeneral(MPI_COMM_WORLD,dirichletCols_size,dirichletCols_array,PETSC_COPY_VALUES,&IS_dirichletCols);
		ISDuplicate(IS_freeRows,&IS_freeCols);
		ISDuplicate(IS_freeRows,&IS_dirichletRows);

		delete[] dirichletCols_array; dirichletCols_array=0;
		delete[] freeRows_array; freeRows_array=0;
	}

	void MeshData::mapNodeID(){
		pEntity node;
		int i,ID;

		i = 0;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			ID = EN_id(node);
			nodeID_map[ID] = i++;
		}
		VIter_delete(vit);
		msize = M_numVertices(theMesh);
	}

	void MeshData::getArrays_free(int& freeRows_size, int** freeRows_array){
		int flag, i;
		std::set<int> freenodes_set;
		pEntity node;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			flag = GEN_tag( node->getClassification() );
			if (pSimPar->isNodeFree(flag)){
				freenodes_set.insert(EN_id(node));
			}
		}
		VIter_delete(vit);

		freeRows_size = (int)freenodes_set.size();
		*freeRows_array = new int[freeRows_size];

		i = 0;
		std::set<int>::iterator iter_set;
		for(iter_set=freenodes_set.begin(); iter_set!=freenodes_set.end(); iter_set++){
			*freeRows_array[i++] = *iter_set;
		}
	}

	void MeshData::getArrays_dirichlet(int& dirichletCols_size, int** dirichletCols_array, std::set<pEntity>& dirichletnodes_set){
		int flag, i;
		pEntity node;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			flag = GEN_tag( node->getClassification() );
			if (!pSimPar->isNodeFree(flag)){
				dirichletnodes_set.insert(node);
			}
		}
		VIter_delete(vit);

		dirichletCols_size = (int)dirichletnodes_set.size();
		*dirichletCols_array = new int[dirichletCols_size];

		i = 0;
		std::set<pEntity>::iterator iter_set;
		for(iter_set=dirichletnodes_set.begin(); iter_set!=dirichletnodes_set.end(); iter_set++){
			*dirichletCols_array[i++] = EN_id(*iter_set);
		}
	}

	void MeshData::initDirichletVectors(std::set<pEntity>& dirichletnodes_set){
		// defines variables and allocate memory for vectors
		int flag, i;
		pEntity node;
		nnfreen = (int)dirichletnodes_set.size();
		nfreen = msize - nnfreen;
		dirichlet_idx = new int[nnfreen];
		dirichlet_data = new double[nnfreen];

		i = 0;
		std::set<pEntity>::iterator iter_set;
		for(iter_set=dirichletnodes_set.begin(); iter_set!=dirichletnodes_set.end(); iter_set++){
			node = (pEntity)*iter_set;
			flag = GEN_tag( node->getClassification() );
			dirichlet_idx[i] = i;
			dirichlet_data[i] = pSimPar->getBC_Value(flag);
			i++;
		}
	}

	void MeshData::setGlobalIndices(){
		pEntity elem;
		int ID[3];							// faces vertices ID
		int nelem = M_numFaces(theMesh);	// number of elements
		int row, i;
		MatIndices.allocateMemory(nelem,3);
		MatIndices.initialize(0);

		row = 0;
		FIter fit = M_faceIter(theMesh);
		while (( elem=FIter_next(fit) )){
			getTriVerticesIDs(elem,ID);
			for (i=0;i<3;i++){
				ID[i] = nodeID_map[ ID[i] ];
				MatIndices.setValue(row,i,ID[i]);
			}
			row++;
		}
		FIter_delete(fit);
	}

	void MeshData::getGlobalIndices(int ith_elem, int ith_edge, int* idxm, int* idxn){
		int row = ith_elem;
		int col = ith_edge;
		idxm[0] = MatIndices.getValue(row,5*col);
		idxm[1] = MatIndices.getValue(row,5*col+1);
		idxn[0] = MatIndices.getValue(row,5*col+2);
		idxn[1] = MatIndices.getValue(row,5*col+3);
		idxn[2] = MatIndices.getValue(row,5*col+4);
	}

	IS MeshData::IS_getFreeRows() const{
		return IS_freeRows;
	}

	IS MeshData::IS_getFreeCols() const{
		return IS_freeCols;
	}

	IS MeshData::IS_getDirichletRows() const{
		return IS_dirichletRows;
	}

	IS MeshData::IS_getDirichletCols() const{
		return IS_dirichletCols;
	}
	int MeshData::nrows() const{
		return msize;
	}

	int MeshData::ncols() const{
		return msize;
	}

	int MeshData::numFreeNodes() const{
		return nfreen;
	}

	const int* MeshData::getDirichlet_idx() const{
		return dirichlet_idx;
	}

	const double* MeshData::getDirichlet_data() const{
		return dirichlet_data;
	}
}
