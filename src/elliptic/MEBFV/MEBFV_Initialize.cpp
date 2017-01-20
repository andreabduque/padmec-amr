/*
 * MEBFV_Initialize.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: rogerio
 */

#include "MEBFV_elliptic.h"

namespace PRS{
	void MEBFV_elliptic::Initialize(pMesh theMesh){

		if (initialize){
			return;
		}
		// initialize matrices/vectors only once
		initialize = true;

		int nrows = pMData->nrows();					// number of mesh nodes
		int ncols = pMData->ncols();					// number of mesh nodes
		int numFreeNodes = pMData->numFreeNodes(); 		// number of mesh free nodes
		int numDirichletNodes = nrows - numFreeNodes;	// number of dirichlet nodes (prescribed)

		// generate vertices mapping for matrix/vector assemblies
		pMData->createGlobalNodeIDMapping();

		// define A matrix
		MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,nrows,ncols,80,PETSC_NULL,80,PETSC_NULL,&A);

		// define b vector
		VecCreate(PETSC_COMM_WORLD,&b);
		VecSetSizes(b,PETSC_DECIDE,numFreeNodes);
		VecSetFromOptions(b);

		// defines solution vector
		VecDuplicate(b,&x);

		// defines SST (Source/Sink term) vector
		VecDuplicate(b,&SST);
		VecZeroEntries(SST);
		pPPData->calculateNodalFlowRate(theMesh,pGCData);
		setSST();

		// setup rhs vector: dirichlet
		VecCreate(PETSC_COMM_WORLD,&dirichletVec);
		VecSetSizes(dirichletVec,PETSC_DECIDE,numDirichletNodes);
		VecSetFromOptions(dirichletVec);
		VecZeroEntries(dirichletVec);
		VecSetValues(dirichletVec,numDirichletNodes,pMData->getDirichlet_idx(),pMData->getDirichlet_data(),INSERT_VALUES);
		VecAssemblyBegin(dirichletVec);
		VecAssemblyEnd(dirichletVec);
	}

	void MEBFV_elliptic::setSST(){
		int nwells, nnodes, i, j, row;
		double Qi;
		nwells = pPPData->getNumWells();			// number of wells prescribed (Neumann boundary value)
		for (i=0; i<nwells; i++){
			nnodes = pPPData->getNumWellNodes(i);	// number of nodes located on well
			for (j=0; j<nnodes; j++){
				pPPData->getFlowRate(i,j,row,Qi);
				VecSetValue(SST,row,Qi,INSERT_VALUES);
			}
		}
		VecAssemblyBegin(SST);
		VecAssemblyEnd(SST);
	}

	void MEBFV_elliptic::freeMemory(){
		MatDestroy(&A);
		MatDestroy(&dirichletMat);
		VecDestroy(&b);
		VecDestroy(&x);
		VecDestroy(&SST);
		VecDestroy(&dirichletVec);
	}
}
