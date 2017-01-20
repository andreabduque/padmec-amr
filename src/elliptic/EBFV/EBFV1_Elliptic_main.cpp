#include <time.h>
#include "EBFV/EBFV1_elliptic.h"

namespace PRS{
	EBFV1_elliptic::EBFV1_elliptic(){
	}
	
	EBFV1_elliptic::EBFV1_elliptic(pMesh mesh, PhysicPropData *ppd, SimulatorParameters *sp, GeomData *gcd, MeshData *md){
		matvec_struct = new Data_struct;
		pPPData = ppd;
		pGCData = gcd;
		pSimPar = sp;
		pMData = md;
		theMesh = mesh;
		one_eighth = (double)(1.0/8.0);
		firstVTK = true;

		// Matrix Assembly Support - MAS
		initialize_MAS();

		// Auxiliary vector to assembly distributed matrix system of equations
		pMData->rowsToImport(mesh,matvec_struct->nrows,matvec_struct->rows);
		throw_exception(!matvec_struct->nrows,"Number of rows NULL",__LINE__,__FILE__);
	}
	
	EBFV1_elliptic::~EBFV1_elliptic(){
		if (matvec_struct){
			delete[] matvec_struct->rows;
			matvec_struct->rows = 0;
			delete matvec_struct;
			matvec_struct = 0;
		}
	}
	
	// solves system of equation for pressure field
	double EBFV1_elliptic::solver(pMesh theMesh){
		if (pSimPar->userRequiresAdaptation()){
			if (pSimPar->adaptation_ocurred()){
				matvec_struct = new Data_struct;
			}
		}
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: pressure solver\tIN\n";
		#endif

		if ( !pSimPar->exactSolutionExist() ){
			assembly_EFG_RHS(theMesh);
			setMatrixFreeOperation(theMesh);
		}
		updatePressure(theMesh);
		calculatePressureGradient();
		freeMemory();

		// print VTK only for: Sw = S_initial and p = p(Sw_initial). VTK zero
		if (firstVTK){
			if (!pSimPar->timeToPrintVTK()){
				pSimPar->setTimetoPrintVTK();	// forces VTK to be printed
				//pSimPar->printOutVTK(theMesh,pPPData,0,pSimPar,pGCData,exportSolutionToVTK);
				firstVTK = false;
			}
		}

		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: pressure solver\tOUT\n";
		#endif
		return 0;
	}
	
	double EBFV1_elliptic::updatePressure(pMesh theMesh){
		CPU_Profile::Start();

		PetscScalar val;
		if ( !pSimPar->exactSolutionExist() ){
			int nLIDs, *IDs_ptr;
			pMData->getRemoteIDs(nLIDs,&IDs_ptr);
			for(int i=0; i<nLIDs;i++){
				int ID = pMData->get_PETScToApp_Ordering(IDs_ptr[i]+1);
				VecGetValues(output,1,&i,&val);
				//cout << "p = " << val << endl;
				pPPData->setPressure(ID-1,val);
			}
		}

		static bool key = true;
		if (key){
			int nnodes = M_numVertices(theMesh);
			if ( !pSimPar->exactSolutionExist() ){
				for(int i=0;i<nnodes;i++){
					int ID = pMData->get_AppToPETSc_Ordering(i+1);
					if ( pMData->getDirichletValue(ID,&val) ){
						pPPData->setPressure(i,val);
					}
				}
				key = false;
			}
			else{
				double coords[3], x, y, z;
				for(int i=0;i<nnodes;i++){
					pGCData->getCoordinates(i,coords);
					x = coords[0]; y = coords[1]; z = coords[2];
					pPPData->setPressure(i,pSimPar->exact_solution(x,y,z));
				}
			}
		}

		CPU_Profile::End("updatePressure");
		return 0;
	}
	
	double EBFV1_elliptic::freeMemory(){
		CPU_Profile::Start();

		if ( !pSimPar->exactSolutionExist() ){
			MatDestroy(&matrix);
			MatDestroy(&matvec_struct->G);

			for(int i=0; i<pSimPar->getNumDomains(); i++){
				MatDestroy(&matvec_struct->E[i]);
				MatDestroy(&matvec_struct->F[i]);
			}
			delete[] matvec_struct->F;
			delete[] matvec_struct->E;
			matvec_struct->F = 0;
			matvec_struct->E = 0;

			VecDestroy(&matvec_struct->RHS);
			VecDestroy(&matvec_struct->z);
			VecDestroy(&output);
		}
		
		if (pSimPar->userRequiresAdaptation()){
			if (pSimPar->adaptation_ocurred()){
				// todo: ta dando erro aqui na hora de liberar memoria
				delete[] matvec_struct->rows; matvec_struct->rows = 0;
				delete matvec_struct; matvec_struct = 0;
			}
		}

		CPU_Profile::End("freeMemory");
		return 0;
	}

	void EBFV1_elliptic::initialize_MAS(){
		int dim = pGCData->getMeshDim();
		int nedges = 0;
		for(int i=0;i<pGCData->getNumDomains(); i++){
			nedges += pGCData->getNumEdgesPerDomain(i);
		}
		pMAS = new MAS;
		pMAS->Eij = new double*[nedges];
		pMAS->Gij = new double*[nedges];
		pMAS->indices = new int*[nedges];
		pMAS->edge_lambda = new double[nedges];

		for(int i=0;i<nedges;i++){
			pMAS->Eij[i] = new double[4*dim];
			pMAS->Gij[i] = new double[4];
			pMAS->indices[i] = new int[2];
		}
		Perform_Assembling = true;
	}

	void EBFV1_elliptic::finalize_MAS(){
		int nedges = 0;
		for(int i=0;i<pGCData->getNumDomains(); i++){
			nedges += pGCData->getNumEdgesPerDomain(i);
		}
		for(int i=0;i<nedges;i++){
			delete[] pMAS->Eij[i]; pMAS->Eij[i] = 0;
			delete[] pMAS->Gij[i]; pMAS->Gij[i] = 0;
			delete[] pMAS->indices[i]; pMAS->indices[i] = 0;
		}
		delete[] pMAS->Eij; pMAS->Eij = 0;
		delete[] pMAS->Gij; pMAS->Gij = 0;
		delete[] pMAS->indices; pMAS->indices = 0;
		delete pMAS;
	}
}
