/*
 * IAR_gradients.cpp
 *
 * IAR_gradients (Interpolation for Adaptative Remeshing)
 *
 *  Created on: 15/02/2013
 *      Author: rogsoares
 */


#include "interpolation.h"

double calculate_Gradients(InterpolationDataStruct* pIData, int field){
	int NumOfNodes = M_numVertices(pIData->m2);

	//PETSc variables
	Vec Fx, Fy, gradx, grady;
	Mat M;

	//create vectors Fx, Fy, gradx and grady
	VecCreate(PETSC_COMM_WORLD,&Fx);
	VecSetSizes(Fx,PETSC_DECIDE,NumOfNodes);
	VecSetFromOptions(Fx);
	VecDuplicate(Fx,&Fy);
	VecDuplicate(Fx,&gradx);
	VecDuplicate(Fx,&grady);

	//create an sparse matrix
	MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,NumOfNodes,NumOfNodes,80,PETSC_NULL,80,PETSC_NULL,&M);

	//Assembly
	Assembly_Mat_Vec(pIData,field,M,Fx,Fy);


	//Solve the linear system
	KSP_solver(M, M, Fx, gradx);
	KSP_solver(M, M, Fy, grady);

	//VecView(Fy,PETSC_VIEWER_STDOUT_SELF);


	//store first order derivatives
	store_FirstOrderDerivatives(pIData,gradx,grady);

	//Second order derivatives
	calculate_SecondOrderDerivatives(pIData);

	//destroy vectors, matrix and KSP's
	VecDestroy(&Fx);
	VecDestroy(&Fy);
	VecDestroy(&gradx);
	VecDestroy(&grady);
	MatDestroy(&M);
	return 0;
}

double Assembly_Mat_Vec(InterpolationDataStruct* pIData, int field, Mat M, Vec Fx, Vec Fy) {
	double coord1[3],coord2[3],coord3[3], value1, value2, value3, Fx_values[3],Fy_values[3];

	//Assembly
	FIter fit = M_faceIter(pIData->m2);
	while (pEntity entity = FIter_next(fit)){

		//get the coordinates and the id's from the element's vertices
		pEntity vertices[3] = {entity->get(0,0), entity->get(0,1), entity->get(0,2)};
		int IDs[3] = {EN_id(vertices[0]), EN_id(vertices[1]), EN_id(vertices[2])};
		int position[3] = {IDs[0]-1, IDs[1]-1, IDs[2]-1};
		V_coord(vertices[0],coord1);
		V_coord(vertices[1],coord2);
		V_coord(vertices[2],coord3);

		//calculate area
		double At = F_area(coord1,coord2,coord3);
		double At_by_6 = (double)(At/6.);
		double At_by_12 = (double)(At/12.);
		double M_values[9] = {At_by_6, At_by_12, At_by_12, At_by_12, At_by_6, At_by_12, At_by_12, At_by_12, At_by_6};

		//get node values
		pIData->pGetDblFunctions[field](IDs[0]-1,value1);
		pIData->pGetDblFunctions[field](IDs[1]-1,value2);
		pIData->pGetDblFunctions[field](IDs[2]-1,value3);
		double CV[3] = {value1/6., value2/6., value3/6.};

		//calculate the contributions
		for(int i=0;i<3;i++){
			Fx_values[i] = (coord2[1]-coord3[1])*CV[0] + (coord3[1]-coord1[1])*CV[1] + (coord1[1]-coord2[1])*CV[2];
			Fy_values[i] = (coord3[0]-coord2[0])*CV[0] + (coord1[0]-coord3[0])*CV[1] + (coord2[0]-coord1[0])*CV[2];
		}
		VecSetValues(Fx,3,position,Fx_values, ADD_VALUES);
		VecSetValues(Fy,3,position,Fy_values, ADD_VALUES);
		MatSetValues(M,3,position,3,position,M_values, ADD_VALUES);
	}
	FIter_delete(fit);


	MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(Fx);
	VecAssemblyEnd(Fx);
	VecAssemblyBegin(Fy);
	VecAssemblyEnd(Fy);
	return 0;
}

// solve system of equations: Ax=y
double KSP_solver(Mat A, Mat pcMatrix, Vec y, Vec x){
	double startt = MPI_Wtime();
	KSP ksp;
	PC preconditioner;

	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetOperators(ksp,A,pcMatrix);
	KSPSetType(ksp,KSPCG);
	KSPGetPC(ksp,&preconditioner);
	PCSetType(preconditioner,PCJACOBI);
	KSPSetInitialGuessNonzero(ksp,PETSC_FALSE);
	KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
	KSPSetFromOptions(ksp);
	KSPSolve(ksp,y,x);
	KSPDestroy(&ksp);
	double endt = MPI_Wtime();
	return endt-startt;
}

double store_FirstOrderDerivatives(InterpolationDataStruct* pIData, Vec gradx, Vec grady){
	double *dx, *dy;
	VecGetArray(gradx,&dx);
	VecGetArray(grady,&dy);

	//loop on vertices to calculate the secord order derivatives
	VIter vit = M_vertexIter(pIData->m2);
	while ( pEntity vertex = VIter_next(vit) ){
		double coord[3];
		int id = EN_id(vertex);
		V_coord(vertex,coord);

		//store the first order derivatives
		pIData->pGrad(id,0) = dx[id-1];
		pIData->pGrad(id,1) = dy[id-1];
	}
	VIter_delete(vit);
	VecRestoreArray(gradx,&dx);
	VecRestoreArray(grady,&dy);
	return 0;
}

double calculate_SecondOrderDerivatives(InterpolationDataStruct* pIData) {
	//loop on vertices to calculate the secord order derivatives
	VIter vit = M_vertexIter(pIData->m2);
	while ( pEntity vertex = VIter_next(vit) ){
		//find out number of elements made with vertex
		int NumOfNeighbours = V_numFaces(vertex);
		if (!NumOfNeighbours){
			pIData->m2->modifyState(0,2);
			pIData->m2->modifyState(2,0);
			NumOfNeighbours = V_numFaces(vertex);
			if (!NumOfNeighbours){
				throw Exception(__LINE__,__FILE__,"Any neighbors found!\n");
			}
		}
		pEntity vertice1, vertice2, vertice3;
		int id = EN_id(vertex);
		double coord[3];
		V_coord(vertex,coord);

		//calculate second order derivatives on all elements around a certain vertex and compute a media
		double dx2 = 0;
		double dy2 = 0;
		double dxdy = 0;

		//find the elements that have the vertex and calculate their contribution
		for(int i=0; i<NumOfNeighbours; i++){
			double ptn1[3],ptn2[3],ptn3[3];

			//get the coordinates and id's
			pEntity face = vertex->get(2,i);
			if (!face){
				throw Exception(__LINE__,__FILE__,"Null face.\n");
			}

			//get coordinates and id from vertice 1
			vertice1 = face->get(0,0);
			V_coord(vertice1,ptn1);
			int id1 = EN_id(vertice1);
			double Dx1 = pIData->pGrad(id1,0);
			double Dy1 = pIData->pGrad(id1,1);

			//get coordinates and id from vertice 2
			vertice2 = face->get(0,1);
			V_coord(vertice2,ptn2);
			int id2 = EN_id(vertice2);
			double Dx2 = pIData->pGrad(id2,0);
			double Dy2 = pIData->pGrad(id2,1);

			//get coordinates and id from vertice 3
			vertice3 = face->get(0,2);
			V_coord(vertice3,ptn3);
			int id3 = EN_id(vertice3);
			double Dx3 = pIData->pGrad(id3,0);
			double Dy3 = pIData->pGrad(id3,1);

			//calculate area
			double Area = F_area(ptn1,ptn2,ptn3);

			//calculate the contributions
			dx2 = dx2 + 0.5*(ptn2[1]-ptn3[1])*Dx1/Area + 0.5*(ptn3[1]-ptn1[1])*Dx2/Area + 0.5*(ptn1[1]-ptn2[1])*Dx3/Area;
			dy2 = dy2 + 0.5*(ptn3[0]-ptn2[0])*Dy1/Area + 0.5*(ptn1[0]-ptn3[0])*Dy2/Area + 0.5*(ptn2[0]-ptn1[0])*Dy3/Area;
			dxdy = dxdy + 0.5*(ptn2[1]-ptn3[1])*Dy1/Area + 0.5*(ptn3[1]-ptn1[1])*Dy2/Area + 0.5*(ptn1[1]-ptn2[1])*Dy3/Area;
		}

		//calculate the second order derivatives as a media
		pIData->pGrad(id,2) = dx2/NumOfNeighbours;
		pIData->pGrad(id,3) = dy2/NumOfNeighbours;
		pIData->pGrad(id,4) = dxdy/NumOfNeighbours;
	}
	VIter_delete(vit);
	return 0;
}
