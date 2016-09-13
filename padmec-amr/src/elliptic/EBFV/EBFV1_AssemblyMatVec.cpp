/*
 * EBFV1_AssemblyMatVec.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: rogerio
 */


#include "EBFV/EBFV1_elliptic.h"

namespace PRS{

	double EBFV1_elliptic::assembly_EFG_RHS(pMesh mesh){

		// compute assembly matrices timing
		CPU_Profile::Start();

		int np = pMData->getNum_GNodes();
		int numGF = pMData->getNum_GF_Nodes();
		int dim = pGCData->getMeshDim();
		int ndom = pSimPar->getNumDomains();

		matvec_struct->ndom = ndom;
		matvec_struct->F_nrows = np*dim;
		matvec_struct->F_ncols = np;

		const double* Cij = NULL;
		const int* indices = NULL;
		int i, nedges, dom, id0, id1, dom_flag;

		// Matrices assembly based only on geometric data. Allocating memory
		if (Perform_Assembling || pSimPar->adaptation_ocurred()){
			E = new Mat[ndom];
			F = new Mat[ndom];
			MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,np,np,100,PETSC_NULL,100,PETSC_NULL,&G_tmp);
			for (dom=0; dom<ndom; dom++){
				MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,np,np*dim,100,PETSC_NULL,100,PETSC_NULL,&E[dom]);
				MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,np*dim,np,100,PETSC_NULL,100,PETSC_NULL,&F[dom]);
			}
		}

		int counter = 0;
		/*
		 * If adaptation is not present, matrices E, F ang G are assembly just once.
		 * Matrices E ang G are geometric matrices multiplied by the mobility
		 * This multiplication will be done are new time step, but E ang G will remain the same for
		 * the whole simulation
		 *
		 * If adaptation is required, then all matrices MUST be computed again!
		 */
		if (pSimPar->adaptation_ocurred()){
			pMData->rowsToImport(mesh,matvec_struct->nrows,matvec_struct->rows);
			throw_exception(!matvec_struct->nrows,"Number of rows NULL",__LINE__,__FILE__);
			initialize_MAS();
		}

		if (Perform_Assembling || pSimPar->adaptation_ocurred()){
			for (dom=0; dom<ndom; dom++){
				nedges = pGCData->getNumEdgesPerDomain(dom);
				dom_flag = pGCData->getDomFlag(dom);
				for (i=0; i<nedges; i++){
					pGCData->getCij(dom,i,Cij);
					pGCData->getEdge(dom,i,indices);
					// idx0        = indices[0]
					// idx1        = indices[1]
					// idx0_global = indices[2]
					// idx1_global = indices[3]

					pGCData->getID(dom,indices[0],indices[1],id0,id1);
					divergence_E(E[dom],Cij,i,dom,dom_flag,indices[2],indices[3],id0,id1,dim,counter);
					divergence_G(G_tmp,Cij,i,dom,dom_flag,indices[2],indices[3],id0,id1,dim,counter);
					gradient_F_edges(F[dom],Cij,dom,indices[0],indices[1],id0,id1,dim);
					setMASindices(counter,id0,id1);
					counter++;
				}
				gradient_F_bdry(F[dom],dom);
				assemblyMatrix(F[dom]);
			}
			Perform_Assembling = false;
		}
		CPU_Profile::End("MatricesAssembly-01");

		// this multiplication will be performed every new time step.
		multiplyMatricesbyMobility();

		// every time G and E are multiplied by the mobility they MUST be assembled
		G_assembly(G_tmp);
		counter = 0;
		for (dom=0; dom<ndom; dom++){
			E_assembly(E[dom],dom,counter);
		}

		CPU_Profile::Start();
		// Get from matrix G its contribution to RHS. matvec_struct->G correspond to all free nodes
		set_SOE(mesh,G_tmp,matvec_struct->G,true,matvec_struct->RHS,true,true);

		// Get from matrices E and F their contribution to RHS. E and F must be handled domain by domain.
		Mat EF, tmp;
		Vec EF_rhs;
		for (i=0; i<ndom; i++){
			MatMatMult(E[i],F[i],MAT_INITIAL_MATRIX,1.0,&EF);
			set_SOE(mesh,EF,tmp,false,EF_rhs,true,false);		// EF matrix is destroyed inside set_SOE function
			VecAXPY(matvec_struct->RHS,1.0,EF_rhs);
			VecDestroy(&EF_rhs);
			MatDestroy(&EF);
		}

		// Create E and F matrices related to free nodes only. Note that they were been created using all mesh nodes.
		/* TRANSFER VALUES FROM TEMPORARY E MATRIX TO matvec_struct->E */
		matvec_struct->F = new Mat[ndom];
		matvec_struct->E = new Mat[ndom];

		// Gets all E[i] columns. Same for all processors
		int nrows = matvec_struct->nrows;
		int *rows = matvec_struct->rows;

		static bool cvfm = true;
		if (!pSimPar->userRequiresAdaptation()){
			if (cvfm){
				pMData->createVectorsForMatrixF(F[0]);
				cvfm = false;
			}
		}
		else{
			pMData->createVectorsForMatrixF(F[0]);
		}
		CPU_Profile::End("MatricesAssembly-02");

		CPU_Profile::Start();
		IS rows_F, cols_F, rows_E, cols_E;
		ISCreateGeneral(PETSC_COMM_WORLD,pMData->get_F_nrows(),pMData->get_F_rows_ptr(),PETSC_COPY_VALUES,&rows_F);
		ISCreateGeneral(PETSC_COMM_WORLD,numGF,pMData->get_F_cols_ptr(),PETSC_COPY_VALUES,&cols_F);
		ISCreateGeneral(PETSC_COMM_WORLD,nrows,rows,PETSC_COPY_VALUES,&rows_E);
		ISCreateGeneral(PETSC_COMM_WORLD,np*dim,pMData->get_pos_ptr(),PETSC_COPY_VALUES,&cols_E);

		for (i=0; i<ndom; i++){
			MatGetSubMatrix(F[i],rows_F,cols_F,MAT_INITIAL_MATRIX,&matvec_struct->F[i]);
			MatGetSubMatrix(E[i],rows_E,cols_E,MAT_INITIAL_MATRIX,&matvec_struct->E[i]);
		}

		ISDestroy(&rows_F);
		ISDestroy(&cols_F);
		ISDestroy(&rows_E);
		ISDestroy(&cols_E);

		// if the mesh has been modified, then matrices must be deleted and new ones will be created
		if (pSimPar->adaptation_ocurred()){
			for (i=0; i<ndom; i++){
				MatDestroy(&E[i]);
				MatDestroy(&F[i]);
			}
			MatDestroy(&G_tmp);
			finalize_MAS();
		}
		// otherwise, reuse the original matrices and save time allocating memory.
		else{
			MatZeroEntries(G_tmp);
			for (i=0; i<ndom; i++){
				MatZeroEntries(E[i]);
			}
		}

		// Create the output vector
		VecCreate(PETSC_COMM_WORLD,&output);
		VecSetSizes(output,PETSC_DECIDE,numGF);
		VecSetFromOptions(output);
		CPU_Profile::End("MatricesAssembly-03");
		return 0;
	}

	int EBFV1_elliptic::assemblyMatrix(Mat A){
		MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
		return 0;
	}

	int EBFV1_elliptic::set_SOE(pMesh mesh, Mat A, Mat &LHS, bool assemblyLHS, Vec &RHS,bool assemblyRHS, bool includeWell){
		int i=0, j=0;
		int numGF = pMData->getNum_GF_Nodes();	// set global free nodes
		int numGP = pMData->getNum_GP_Nodes();	// set global prescribed (dirichlet) nodes

		// -------------------------------------------------------------------------
		// 					ASSEMBLY MATRIX LHS (STIFFNESS MATRIX)
		// -------------------------------------------------------------------------
		int nrows = matvec_struct->nrows;
		int *rows = matvec_struct->rows;

		//cout << "n_rows: " << nrows << endl;

		IS rowsToImport;
		ISCreateGeneral(PETSC_COMM_WORLD,nrows,rows,PETSC_COPY_VALUES,&rowsToImport);

		IS colsToImport;
		ISCreateGeneral(PETSC_COMM_WORLD,numGF,pMData->get_idxFreecols_ptr(),PETSC_COPY_VALUES,&colsToImport);

#ifdef _SEEKFORBUGS_
		if (!nrows){
			throw Exception(__LINE__,__FILE__,"nrows NULL!");
		}
		if (!rows || !pMData->get_idxFreecols_ptr()){
			throw Exception(__LINE__,__FILE__,"rows NULL!");
		}
#endif

		if (assemblyLHS){
			MatGetSubMatrix(A,rowsToImport,colsToImport,MAT_INITIAL_MATRIX,&LHS);
		}

		// --------------------------- ASSEMBLY RHS VECTOR -------------------------
		// matrix formed by columns associated to prescribed nodes.
		// -------------------------------------------------------------------------
		if (assemblyRHS){
			Mat rhs;
			int m,n;

			IS colsToImport_dirichlet;
			ISCreateGeneral(PETSC_COMM_WORLD,numGP,pMData->get_idxn_ptr(),PETSC_COPY_VALUES,&colsToImport_dirichlet);
			MatGetSubMatrix(A,rowsToImport,colsToImport_dirichlet,MAT_INITIAL_MATRIX,&rhs);
			//MatDestroy(&A);
			ISDestroy(&colsToImport_dirichlet);

			VecCreate(PETSC_COMM_WORLD,&RHS);
			VecSetSizes(RHS,PETSC_DECIDE,numGF);
			VecSetFromOptions(RHS);

			// fill RHS vector with values from rhs matrix. entries from this matrix are multiplied by their corresponding prescribed values
			//printMatrixToFile(rhs,"rhs_nondirichlet.txt");
			MatGetOwnershipRange(rhs,&m,&n);
			for (i=m; i<n; i++){
				double sum = .0;
				MIter prescribedIter = pMData->dirichletBegin();
				for (j=0; j<numGP; j++, prescribedIter++){
					int col = j;
					double val=.0;
					MatGetValues(rhs,1,&i,1,&col,&val);
					sum += -val*prescribedIter->second;
				}
				VecSetValue(RHS,i,sum,ADD_VALUES);

			}

			if ( includeWell ){
				wells_RHS_Assembly__Wells(mesh,RHS);
			}
			VecAssemblyBegin(RHS);
			VecAssemblyEnd(RHS);
			MatDestroy(&rhs);
		}
		ISDestroy(&rowsToImport);
		ISDestroy(&colsToImport);
		return 0;
	}

	void EBFV1_elliptic::multiplyMatricesbyMobility(){
		int i, idx_I, idx_J;
		double Sw_0, Sw_1;

		int nedges = 0;
		int ndom = pGCData->getNumDomains();
		for(i=0; i<ndom; i++){
			nedges += pGCData->getNumEdgesPerDomain(i);
		}

		for (i=0; i<nedges; i++){
			idx_I = pMAS->indices[i][0] - 1;					// it return global vertex ID
			idx_J = pMAS->indices[i][1] - 1;					// it return global vertex ID
			pPPData->getSaturation(idx_I,Sw_0);
			pPPData->getSaturation(idx_J,Sw_1);
			pMAS->edge_lambda[i] = 0.5*(pPPData->getTotalMobility(Sw_0) + pPPData->getTotalMobility(Sw_1));
		}
	}

	int EBFV1_elliptic::G_assembly(Mat G){
		int i, j, k;
		double Gij[4];
		int indices[2];
		int ndom = pGCData->getNumDomains();
		int counter = 0;
		for(i=0; i<ndom; i++){
			int nedges = pGCData->getNumEdgesPerDomain(i);
			for(j=0; j<nedges; j++){
				indices[0] = pMAS->indices[counter][0] - 1;
				indices[1] = pMAS->indices[counter][1] - 1;

				for(k=0; k<4; k++){
					Gij[k] = pMAS->Gij[counter][k]*pMAS->edge_lambda[counter];
				}
				MatSetValues(G,2,indices,2,indices,Gij,ADD_VALUES);
				counter++;
			}
		}
		assemblyMatrix(G);
		return 0;
	}

	int EBFV1_elliptic::E_assembly(Mat E, int dom, int &counter){
		int j,k;
		int dim = pGCData->getMeshDim();
		double Eij[4*dim];
		int idxn[2*dim];
		int idxm[2];
		int pos1, pos2;

		int nrows = 2;
		int ncols = 2*dim;

		int nedges = pGCData->getNumEdgesPerDomain(dom);
		for(j=0; j<nedges; j++){
			// define where to assembly sub-matrix
			idxm[0] = pMAS->indices[counter][0] - 1;
			idxm[1] = pMAS->indices[counter][1] - 1;
			pos1 = dim*idxm[0];
			pos2 = dim*idxm[1];
			for (k=0; k<dim; k++){
				idxn[k] = pos1+k;
				idxn[dim+k] = pos2+k;
			}

			// assembly sub-matrix
			for(k=0; k<4*dim; k++){
				Eij[k] = pMAS->Eij[counter][k]*pMAS->edge_lambda[counter];
			}
			MatSetValues(E,nrows,idxm,ncols,idxn,Eij,ADD_VALUES);
			counter++;
		}
		assemblyMatrix(E);
		return 0;
	}

	void EBFV1_elliptic::setMASindices(int count, int id0, int id1){
		if (id0 > id1){
			std::swap(id0,id1);
		}
		pMAS->indices[count][0] = id0;
		pMAS->indices[count][1] = id1;
	}
}



