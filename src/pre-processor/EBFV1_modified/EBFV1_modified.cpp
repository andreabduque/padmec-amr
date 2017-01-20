/*
 * EBFV1_modified.cpp
 *
 *  Created on: Oct 17, 2014
 *      Author: rogerio
 */

#include "EBFV1_modified.h"

// for debug purposes
void printMatrix(string name, double* Mat, int nrows, int ncols);

void EBFV1_modified_preprocessor_2D(pMesh theMesh, GeomData* pGCData){

	// loop over triangles elements. For each one calculate:
	// 		Geometric coeficients: Cij, Dij, and nodal volumes;
	// 		matrices E, F and G;

//	pEntity elem;
//	int numElements, dim;
//	double K[4];
//	double Dij[6], Cij[6], vol;
//	double E_ij[8], E_jk[8], E_ik[8];	// E: 2x4
//	double F_ij[8], F_jk[8], F_ik[8];	// F: 4x2
//	double G_ij[4], G_jk[4], G_ik[4];	// G: 2x2
//	double A_ij[6], A_jk[6], A_ik[6];	// A: 2x3
//
//	pGCData->setMeshDim(2);
//
//	// function pointer which returns the absolute permeability tensor
//	RockPropFuncPointer pFUNC_RockProp;
//
//	dim = theMesh->getDim();
//	theMesh->modifyState(2,1);	// create edge data structure
//	theMesh->modifyState(0,2);	// create adjacency around nodes
//	theMesh->modifyState(1,2);	// create faces around edge
//	numElements = M_numFaces(theMesh);
//	pGCData->initializeElementMatrix(numElements);
//	calculateEdgeLength(theMesh,pGCData);
//	pFUNC_RockProp = getRockPropertyFuncPointer();					// based on user inputs, define a function pointer to what rock property
//
//	int k = 0;
//	FIter fit = M_faceIter(theMesh);
//	while ( (elem = FIter_next(fit)) ){
//		pFUNC_RockProp(elem,K);
//		calculateGeomCoefficients(elem,Cij,Dij,vol);				// calculates Cij and Dij for all three element edges and volume as well.
//		calculateMatrix_E(elem,pGCData,K,Cij,E_ij,E_jk,E_ik);		// calculates edge matrix for all three edges
//		calculateMatrix_F(elem,Cij,Dij,vol,F_ij,F_jk,F_ik);			// calculates edge matrix for all three edges
//		calculateMatrix_G(elem,pGCData,K,Cij,G_ij,G_jk,G_ik);		// calculates edge matrix for all three edges
//		assemblyMatrix_A(E_ij,E_jk,E_ik,F_ij,F_jk,F_ik,G_ij,G_jk,G_ik,A_ij,A_jk,A_ik);
//		pGCData->setElementMatrices(k,A_ij,A_jk,A_ik);
//		k++;
//
//
////		printMatrix("Eij",E_ij,2,4);
////		printMatrix("Ejk",E_jk,2,4);
////		printMatrix("Eik",E_ik,2,4);
////
////		printMatrix("Fij",F_ij,4,2);
////		printMatrix("Fjk",F_jk,4,2);
////		printMatrix("Fik",F_ik,4,2);
////
////		printMatrix("Gij",G_ij,2,2);
////		printMatrix("Gjk",G_jk,2,2);
////		printMatrix("Gik",G_ik,2,2);
////
////		printMatrix("Aij",A_ij,3,3);
////		printMatrix("Ajk",A_jk,3,3);
////		printMatrix("Aik",A_ik,3,3);
//	}
//	FIter_delete(fit);
	//STOP();
}

void printMatrix(string name, double* Mat, int nrows, int ncols){
	cout << name << endl;
	cout << setprecision(5) << scientific;
	int k =0;
	for (int i=0; i< nrows; i++){
		for (int j=0; j<ncols; j++){
			cout << Mat[k++] << "  ";
		}
		cout << endl;
	}
}























