/*
 * EBFV1__MatrixFreeSolver.cpp
 *
 *  Created on: 25/08/2012
 *      Author: rogsoares
 */

#include <time.h>
#include "EBFV/EBFV1_elliptic.h"

namespace PRS{           // PRS: Petroleum Reservoir Simulator
	double EBFV1_elliptic::setMatrixFreeOperation(pMesh theMesh){
		CPU_Profile::Start();

		PetscInt m, n;
		VecGetOwnershipRange(output,&m,&n);
		PetscInt nLRows = n - m;
		PetscInt numGF = pMData->getNum_GF_Nodes();
		MatCreateShell(PETSC_COMM_WORLD,nLRows,nLRows,numGF,numGF,matvec_struct,&matrix);
		MatSetFromOptions(matrix);
		MatShellSetOperation(matrix, MATOP_MULT,(void(*)(void))&EBFV1_elliptic::MatMultUser);

		/*
		 * Solve system of equation using a matrix-free approach
		 */

		// create auxiliary vector: used during matrix-vector operations
		VecCreate(PETSC_COMM_WORLD,&matvec_struct->z);
		VecSetSizes(matvec_struct->z,PETSC_DECIDE,matvec_struct->F_nrows);
		VecSetFromOptions(matvec_struct->z);

		// use last pressure solution as guess solution for iterative solver and then reduce number of iterations
		static bool guess_sol = false;
		PetscBool guessNonZero = PETSC_FALSE;
		if (guess_sol){
			int row = 0;
			double p;
			guessNonZero = PETSC_TRUE;
			// get solution from mesh nodes and set output vector for guess solution
			int idx = 0;
			VIter vit = M_vertexIter(theMesh);
			while (pEntity node = VIter_next(vit)){
				if ( pSimPar->isNodeFree( GEN_tag( node->getClassification() )) ){
					pPPData->getPressure(idx,p);
					VecSetValue(output,row,p,INSERT_VALUES);
					//cout << "ID = " << EN_id(node) << "\tidx: " << idx << "\tsize:" << pMData->getNum_GF_Nodes() << "\tflag: " << GEN_tag(node->getClassification()) << endl;
					idx++;
					row++;
				}
			}
			VIter_delete(vit);
		}
		guess_sol = true;

		// call PETSc to solver system of equation
		MatGetSize(matvec_struct->G,&m,&n);
		VecGetSize(matvec_struct->RHS,&m);
		VecGetSize(output,&m);
		KSP_solver(matrix,matvec_struct->G,matvec_struct->RHS,output,pSimPar,guessNonZero,KSPBCGS,PCASM,_KSPiter);
//		printVectorToFile(output,"output.txt");
//		printMatrixToFile(matvec_struct->G,"matvec_structG.txt");

		if (_KSPiter >= 10000){
			throw Exception(__LINE__,__FILE__,"Number of KSP iterations reached the maximum number allowed: 10.000 iterations."
					" It means that the system of equations didn't converged. Please check boundary conditions.\n");
		}
		else{
			cout << "KSP_Solver: Num. iter =  " << _KSPiter << endl;
		}

		CPU_Profile::End("Solver");
		return 0;
	}
}
