#include "MeshData.h"
#include "GeomData.h"

namespace PRS{

	MeshData::MeshData(){
	}

	MeshData::~MeshData(){
		destroyPointers();
	}

	void MeshData::destroyPointers(){
		if (P_pid()>=1 && rowToImport){
			MatDestroy(&joinNodes);
			MatDestroy(&updateValues);
		}
		if (pMS)         { delete pMS; pMS = 0; }
		if (rowToImport){ delete[] rowToImport; rowToImport = 0; }
		if (F_rows)      { delete[] F_rows; F_rows=0; }
		if (F_cols)      { delete[] F_cols; F_cols=0; }
		if (pos)         { delete[] pos; pos=0; }
		if (idxn)        { delete[] idxn; idxn=0; }
		if (idxFreecols) { delete[] idxFreecols; idxFreecols=0; }
		if (FP_Array)    { delete[] FP_Array; FP_Array=0; }
		if (localIDs)    { delete[] localIDs; localIDs=0; }
		dirichlet.clear();
		setOfDomains.clear();
		localIDNumbering.clear();
		mapPB_nodes.clear();
		AODestroy(&ao);
	}

	MeshData::MeshData(SimulatorParameters *sp, pMesh mesh){
		pSimPar = sp;
		theMesh = mesh;
		structsCreation = true;
		pMS = new UVMN_Struct;
	}

	void MeshData::initialize(pMesh theMesh, GeomData* pGCData){
		reorderVerticesIds(theMesh,0);
		settingFreeAndPrescribedNodes(theMesh);
		createVectorsForRHS(theMesh,theMesh->getDim());
	}

	void MeshData::deallocateData(){
		destroyPointers();
	}

	int MeshData::get_AppToPETSc_Ordering(int n) const {
		AOApplicationToPetsc(ao,1,&n);
		return n;
	}

	int MeshData::get_PETScToApp_Ordering(int n) const{
		AOPetscToApplication(ao,1,&n);
		return n;
	}

	void MeshData::reorderVerticesIds(pMesh theMesh, int (*pFunc_numRemoteCopies)(pEntity)){
		// rank p must take all remote nodes
		set<int> remoteNodesSet;
		set<int>::iterator Iter;
		int numUniqueNodes = 0;

		int i = 0;
		pEntity node;
		VIter vit = M_vertexIter(theMesh);
		while ( (node=VIter_next(vit)) ){
//			if ( pFunc_numRemoteCopies(node) ){
//				remoteNodesSet.insert(EN_id(node));
//			}
//			else{
				numUniqueNodes++;
//			}
		}

		i = 0;
		int numRemoteNodes = remoteNodesSet.size();
		int *remoteNodes = new int[numRemoteNodes];

		// transfer nodes from set to an array
		for (Iter=remoteNodesSet.begin(); Iter != remoteNodesSet.end(); Iter++) remoteNodes[i++] = *Iter;

		// rank must know how many nodes will receive from all other ranks
		int *numAllRemoteNodes = new int[P_size()];
		MPI_Allgather(&numRemoteNodes,1,MPI_INT,numAllRemoteNodes,1,MPI_INT,MPI_COMM_WORLD);

		int total = 0;
		for (i=0; i<P_size(); i++) total += numAllRemoteNodes[i];
		int *allRemoteNodes = new int[total];

		int *displs = new int[P_size()];
		displs[0] = 0;
		for (i=1; i<P_size(); i++) displs[i] = displs[i-1] + numAllRemoteNodes[i-1];

		// receive remote nodes from all other ranks
		MPI_Allgatherv (remoteNodes,numRemoteNodes,MPI_INT,allRemoteNodes,numAllRemoteNodes,displs,MPI_INT,MPI_COMM_WORLD);

		delete[] displs;
		remoteNodesSet.clear();

		// how many nodes from allRemoteNodes will be checked by rank p
		total = 0;
		for (i=0;i<P_pid(); i++) total += numAllRemoteNodes[i];
		delete[] numAllRemoteNodes;

		// transfer to a set container only those remote nodes from allRemoteNodes from ranks
		// lower than rank p
		for (i=0;i<total;i++) remoteNodesSet.insert(allRemoteNodes[i]);
		delete[] allRemoteNodes;

		// duplicated numbers in different ranks must be avoided to no disturb PETSc!
		// each rank will check which remote nodes from rank p belong also to ranks lower than it
		// store on another set container those nodes that do not belong to ranks lower that rank p
		set<int> remoteNodesLowSet;
		for (i=0;i<numRemoteNodes;i++){
			Iter = remoteNodesSet.find(remoteNodes[i]);
			if ( Iter == remoteNodesSet.end()) remoteNodesLowSet.insert(remoteNodes[i]);
		}
		remoteNodesSet.clear();
		delete[] remoteNodes;

		// N means how many nodes rank p should map.
		int N = numUniqueNodes + remoteNodesLowSet.size();
		// rank 0: mapping = 1,2,3,...,N0
		// rank 1: mapping = N0+1,N0+2,N0+3,...,N0+N1
		// rank 2: mapping = N0+N1+1,NO+N1+2,NO+N1+3,...,N0+N1+N2

		this->numGN = P_getSumInt(N);

		// each rank should know all Ni's
		int *recvNs = new int[P_size()];
		MPI_Allgather(&N,1,MPI_INT,recvNs,1,MPI_INT,MPI_COMM_WORLD);

		// start mapping
		int *apOrdering = new int[N];
		int *petscOrdering = new int[N];

		int from=1, ID, j=0;
		for (i=0; i<P_pid(); i++) from += recvNs[i];

		i = from;
		vit = M_vertexIter(theMesh);
		while ( (node=VIter_next(vit)) ){
			ID = EN_id(node);
//			if ( pFunc_numRemoteCopies(node) ){
//				Iter = remoteNodesLowSet.find( ID );
//				if ( Iter != remoteNodesLowSet.end() ){
//					apOrdering[j] = ID;
//					petscOrdering[j] = i++;
//					j++;
//				}
//			}
//			else{
				apOrdering[j] = ID;
				petscOrdering[j] = i++;
				j++;
//			}
		}
		remoteNodesLowSet.clear();
		delete[] recvNs;

		// Petsc will be used to make the parallel job
		AOCreateMapping(PETSC_COMM_WORLD,N,apOrdering,petscOrdering,&ao);
		delete[] apOrdering;
		delete[] petscOrdering;
	}

	void MeshData::settingFreeAndPrescribedNodes(pMesh theMesh){
		if (!numGN)
			throw Exception(__LINE__,__FILE__,"Number of global nodes is unknown. Did you call reorderVerticesIds before?\n");
		FreePrescribedNodes(theMesh);
		mappingUnknowns();
	}

	int MeshData::FreePrescribedNodes(pMesh theMesh){
		// all processors must hold all prescribed nodes
		// (that's ONLY true for EBFV1 elliptic equation!!!)
		getNodesWithKnownValues(theMesh);

		numGP = dirichlet.size();
		numGF = numGN - numGP;
		if (!P_pid()){
			cout << "Number of processes required: " << P_size() << endl;
			cout << "Number of global nodes      : " << numGN << endl;
			cout << "Number of dirichlet nodes   : " << numGP << endl;
			cout << "Number of free nodes        : " << numGF << endl;
			cout << "Number of edges             : " << M_numEdges(theMesh) << endl;
			if (theMesh->getDim()==2){
				cout << "Number of triangles (2D)    : " << M_numFaces(theMesh) << endl;
			}
			else{
				cout << "Number of tetrahedra (3D)   : " << M_numRegions(theMesh) << endl;
			}
			cout << "-------------------------------------------------------------------------------------------------------------------------\n\n";
		}

		// a set container for all unknowns. But first, we set for all nodes and
		// then subtract the prescribed ids. it will be used to assembly the LHS
		// matrix
		//		int i;
		//		set<int> setGFIDs;
		//		for (i=1; i<=numGN; i++) setGFIDs.insert(i);
		//		for (MIter mit=dirichlet.begin(); mit!=dirichlet.end(); mit++) setGFIDs.erase((*mit).first);
		//		setGFIDs.clear();
		MPI_Barrier(MPI_COMM_WORLD);
		return 0;
	}

	int MeshData::getNodesWithKnownValues(pMesh theMesh){
		pEntity node, edge, face;
		int i,ID;

		// search for flagged nodes
		VIter vit = M_vertexIter( theMesh );
		while ( (node = VIter_next(vit)) ){
			int flag = (!node->getClassification())?0:GEN_tag(node->getClassification());
			if ( !pSimPar->isNodeFree(flag) ){
				ID = get_AppToPETSc_Ordering(EN_id(node));
				dirichlet[ID] = pSimPar->getBC_Value(flag);
			}
		}
		VIter_delete(vit);

//		theMesh->modifyState(3,1);
//		theMesh->modifyState(1,3);

		// search for flagged edges
		EIter eit = M_edgeIter( theMesh );
		while ( (edge = EIter_next(eit)) ){
			if (!theMesh->getRefinementDepth(edge)){
				int flag = EN_getFlag(edge);
				#ifdef CRUMPTON_EXAMPLE
				double coord1[3],coord2[3],x1,y1,x2,y2;
				V_coord(edge->get(0,0),coord1); x1 = coord1[0]; y1 = coord1[1];
				V_coord(edge->get(0,1),coord2); x2 = coord2[0]; y2 = coord2[1];
				if (flag >= 2000 && flag <3000 &&  flag != 2005 ){	// take only external boundary edges
					// boundary conditions:
					// u(x,y) = [2*sin(y)+cos(y)]alpha*x + sin(y),	x <= 0
					//			exp(x)cos(y),						x > 0

					ID = get_AppToPETSc_Ordering(EN_id(edge->get(0,0)));
					dirichlet[ID] = (x1 <= .0)?((2.*sin(y1)+cos(y1))*ALPHA*x1 + sin(y1)):exp(x1)*sin(y1);
					ID = get_AppToPETSc_Ordering(EN_id(edge->get(0,1)));
					dirichlet[ID] = (x2 <= .0)?((2.*sin(y2)+cos(y2))*ALPHA*x2 + sin(y2)):exp(x2)*sin(y2);
				}
				#else
				if ( !pSimPar->isNodeFree(flag) )
					for (i=0; i<2; i++){
						ID = get_AppToPETSc_Ordering(EN_id(edge->get(0,i)));
						dirichlet[ID] = pSimPar->getBC_Value(flag);
					}
				#endif
			}
		}
		EIter_delete(eit);

		// search for flagged faces (on boundary only)
		// -------------------------------------------------------
		if (theMesh->getDim()==3){

			// External boundary condition defined
			// -------------------------------------------------------
			if (pSimPar->SimulationHas_BC_ExternalDefinition()){
				double coords[3], x, y, z;
				FIter fit = M_faceIter( theMesh );
				while ( (face = FIter_next(fit)) ){
					if (!theMesh->getRefinementDepth(face)){
						int flag = EN_getFlag(face);
						if ( !pSimPar->isNodeFree(flag) ){
							for (i=0; i<3; i++){
								ID = EN_id(face->get(0,i));
								pVertex v = (pVertex)theMesh->getVertex(ID);
								ID = get_AppToPETSc_Ordering(ID);

								V_coord(v,coords);
								x = coords[0];
								y = coords[1];
								z = coords[2];

								dirichlet[ID] = pSimPar->exact_solution(x,y,z);
							}
						}
					}
				}
				FIter_delete(fit);
			}
			else{
				// conventional (Dirichlet) boundary condition: specified in numeric.dat
				// --------------------------------------------------------------------
				FIter fit = M_faceIter( theMesh );
				while ( (face = FIter_next(fit)) ){
					if (!theMesh->getRefinementDepth(face)){
						int flag = EN_getFlag(face);
						if ( !pSimPar->isNodeFree(flag) ){
							for (i=0; i<3; i++){
								ID = get_AppToPETSc_Ordering(EN_id(face->get(0,i)));
								dirichlet[ID] = pSimPar->getBC_Value(flag);
							}
						}
					}
				}
				FIter_delete(fit);
			}
		}

		for(MIter mit = dirichlet.begin(); mit != dirichlet.end(); mit++){
			//cout << "Node [" << mit->first << "]:\t " << mit->second << endl;
		}

		string msg("Dirichlet nodes were not found. Check if Phisycal command was used properly in .geo file\n");
		throw_exception(!P_getSumInt(dirichlet.size()),msg,__LINE__,__FILE__);


		// go ahead only if parallel
		if (P_size()==1) return 0;
		// now, all partitions must know all prescribed nodes and their flags
		// first of all, let partitions aware of how many prescribed nodes exist on each one
		// if processor p does not have any prescribed node let nLPN=1 because p cannot send 0 element
		int nLPN = dirichlet.size();
		int *recvLP = new int[P_size()];
		MPI_Allgather ( &nLPN, 1, MPI_INT, recvLP, 1, MPI_INT, MPI_COMM_WORLD );
		// number of global prescribed nodes
		// Note that nGPN value is not necessary the real global prescribed nodes
		// Nodes on partition boundary can be counted twice or more
		int nGPN=0;
		for (i=0; i<P_size(); i++) nGPN += recvLP[i];
		// sPIDs = send prescribed IDs    sPFlags = send prescribed flags

		i=0;
		int *sPIDs = new int[nLPN];
		double *sPFlags = new double[nLPN];
		for(MIter mit = dirichlet.begin(); mit != dirichlet.end(); mit++){
			sPIDs[i] = mit->first;
			sPFlags[i] = mit->second;
			i++;
		}

		// rcount says how many values each rank will send
		int *rcounts = recvLP;
		// displs says where to start to read in recv buffer
		int *displs = new int[P_size()];
		displs[0] = 0;
		for (i=1; i<P_size(); i++) displs[i] = displs[i-1]+recvLP[i-1];
		int *rPIDs = new int[nGPN];
		// get all prescribed nodes
		MPI_Allgatherv(sPIDs,nLPN,MPI_INT,rPIDs,rcounts,displs,MPI_INT,MPI_COMM_WORLD);
		double *rPFlags = new double[nGPN];
		// get flags from all prescribed nodes
		MPI_Allgatherv(sPFlags,nLPN,MPI_DOUBLE,rPFlags,rcounts,displs,MPI_DOUBLE,MPI_COMM_WORLD);
		for (i=0; i<nGPN; i++)  dirichlet.insert( pair<int,double>(rPIDs[i],rPFlags[i]) );
		delete[] sPIDs; sPIDs=0;
		delete[] sPFlags; sPFlags=0;
		delete[] rPIDs; rPIDs=0;
		delete[] rcounts; rcounts=0;
		delete[] displs; displs=0;
		return 0;
	}

	bool MeshData::getDirichletValue(int ID, double *val=0){
		bool key = false;
		MIter mit = dirichlet.find( ID );
		if ( mit != dirichlet.end() ){
			if (val) *val = mit->second;
			key = true;
		}
		return key;
	}


	/*
	 * It's supposed the global matrix is assembled for all (free and dirichlet)
	 * nodes. For the current FVM used, where system of equation is assembled in
	 * a sub-domain by sub-domain approach, work becomes easier using all nodes.
	 * After that, it's desired to solve the system of equation only for
	 * the free ones as expected. Each rank needs to import some rows from the
	 * assembled matrix to a new matrix. The function below does this job.
	 *
	 * - nrows: number of rows will be imported and will be local to that rank.
	 * - rows: array with number of rows indices that correspond to a free node.	 *
	 * */
	int MeshData::rowsToImport(pMesh theMesh, int &nrows, int *&rows){
		Mat temp;
		PetscErrorCode ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,numGF,numGF,0,PETSC_NULL,0,PETSC_NULL,&temp);

		// get range of local owned rows of each process
		const PetscInt *ranges[P_size()];
		ierr = MatGetOwnershipRanges(temp,ranges);
		//printf("[%d] - ranges[0] %d, ranges[1] %d, ranges[2] %d, ranges[3] %d\n",P_pid(),ranges[0][0],ranges[0][1],ranges[0][2],ranges[0][3]);

		// each rank should take from matrix A row indices associated to free nodes
		// that fill exactly the number of local rows of LHS, i.e,
		// - rank 0 takes the first n1 rows from A, even if it is not local to rank 0.
		// - rank 1 takes from the first n1 rows to n1+n2
		// - rank 2 takes from the first n1+n2 rows to n1+n2+n3
		// - (...)
		// - rank p takes from the first sum(ni),i=1:p rows to sum(ni),i=i:p+1

		int i,k;
		int RANGE[P_size()];
		for (i=0; i<P_size(); i++){
			RANGE[i] = ranges[0][i+1]-ranges[0][i];
		}

		//printf("[%d] - RANGE %d %d %d\n",P_pid(),RANGE[0],RANGE[1],RANGE[2]);

		const int from = ranges[0][P_pid()];
		const int to = ranges[0][P_pid()+1];
		//printf("[%d] - from %d  to %d\n",P_pid(),from,to); exit(1);
		ierr = MatDestroy(&temp); CHKERRQ(ierr);

		
		// Inform which ROWS from A, on processor p, must be copied to assembly LHS matrix (free nodes). For each rank, the number of 
		// rows associated to free nodes cannot be forecast so a list is used.
		list<int> freerowsList;
		k = 0;
		if (P_pid()==0){
			for (i=1; i<=numGN; i++){
				if ( !getDirichletValue(i,0) && k<to ){
					freerowsList.push_back(i-1);
					k++;
				}
			}
		}
		else{
			for (i=1; i<=numGN; i++){
				if ( !getDirichletValue(i,0)  ){
					if ( k>=from && k<to )
						freerowsList.push_back(i-1);
					k++;
				}
			}
		}

		// transfer row indices from list to array
		nrows = freerowsList.size();
		if ( !nrows ){
			throw Exception(__LINE__,__FILE__,"numLocalRows = 0\n");
		}
		rows = new int[nrows];

		i=0;
		for(list<int>::iterator lit=freerowsList.begin(); lit!=freerowsList.end(); lit++){
			rows[i++] = *lit;
			//	if (P_pid()==1) printf("rows[%d] = %d\n",i-1,rows[i-1]);
		}
		freerowsList.clear();
		//printf("[%d]     rows: %d\n",P_pid(),nrows);

		// Rogerio, try to think in something more clear in the future!
		/*
		 * Create a local array consisting of only local free IDs
		 * */
		set<int> setLFNodes;
		VIter vit = M_vertexIter(theMesh);
		while (pEntity node = VIter_next(vit)){
			int ID = get_AppToPETSc_Ordering(EN_id(node));
			if ( !getDirichletValue(ID,0)  ) {
				setLFNodes.insert(ID);
			}
		}
		if (setLFNodes.size()==0){
			throw Exception(__LINE__,__FILE__,"No local free nodes.\n");
		}
		set<int>::iterator setLFN_Iter = setLFNodes.begin();
		numLocalIDs = setLFNodes.size();
		localIDs = new int[numLocalIDs];
		for (i=0; i<numLocalIDs; i++){
			localIDs[i] = *setLFN_Iter-1;
			setLFN_Iter++;
		}
		setLFNodes.clear();
		return 0;
	}

	void MeshData::getRemoteIDs(int &nLIDs, int** IDs_ptr){
		nLIDs = numLocalIDs;
		*IDs_ptr = localIDs;
	}

	void MeshData::mappingUnknowns(){
		int i, k = 0, j = 0;
		FP_Array = new int[numGN];
		for (i = 0; i<numGN; i++){
			FP_Array[i] = ( !getDirichletValue(i+1,0) )?k++:j++;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	void MeshData::createVectorsForRHS(pMesh theMesh, int dim){
		/*
		 * Allocate memory for auxiliary vectors.
		 * */
		try{
			// inform which COLUMNS from A must be copied to assembly LHS matrix (free nodes)
			// this vector is the same for all processors.
			idxFreecols = new PetscInt[getNum_GF_Nodes()];
			// inform which columns from A must be copied to assembly RHS vector (prescribed nodes)
			// this vector is the same for all processors.
			idxn = new PetscInt[getNum_GP_Nodes()];
			pos = new PetscInt[getNum_GNodes()*dim];
		}
		catch(const std::exception &error){
			std::cerr << "An exception has been caught: " << error.what() << std::endl;
		}

		int i,j=0,k=0;
		int np = getNum_GNodes();
		for (i=1; i<=np; i++){
			if ( !getDirichletValue(i,0) )
				idxFreecols[k++] = i-1;
			else
				idxn[j++] = i-1;
		}
		for (i=0; i<np*dim; i++) pos[i] = i;
	}

	int MeshData::createVectorsForMatrixF(Mat &mat){
		// Gets all F[i] local rows. Varies for all processors
		int i, k, F_m, F_n;
		PetscErrorCode ierr = MatGetOwnershipRange(mat,&F_m,&F_n);CHKERRQ(ierr);
		set_F_nrows(F_n - F_m);
		/*
		 * Allocate memory for auxiliary vectors.
		 * */
		try{
			F_rows = new PetscInt[get_F_nrows()];
			F_cols = new PetscInt[getNum_GF_Nodes()];
		}
		catch(const std::exception &error){
			std::cerr << "An exception has been caught: " << error.what() << std::endl;
		}
		for (i=F_m,k=0; i<F_n; i++) F_rows[k++] = i;
		// Gets F[i] columns related to free nodes
		int np = getNum_GNodes();
		for (i=1,k=0; i<=np; i++){
			if ( !getDirichletValue(i,0) ){
				F_cols[k++] = i-1;
			}
		}
		return 0;
	}

	/*
	 * When user order to unify vectors (gradients) on partition baoudaries the function pointer passed as argument returns a vector.
	 * All coordenates (x,y,z) are unified one-by-one as scalars.
	 * A internal function getScalar/setScalar take the function pointer as argument to return one coordinate at time.
	 */
	double getScalar(int row, UVMN_Struct *pMS);
	void setScalar(int row, double val, UVMN_Struct *pMS);

	double getScalar(int row, UVMN_Struct *pMS){
		double vec[3];
		pMS->pFunc_getVector(pMS->dom,row,vec);
		return vec[pMS->coord_xyz];
	}

	void setScalar(int row, double val, UVMN_Struct *pMS){
		double vec[3];
		pMS->pFunc_getVector(pMS->dom,row,vec);
		vec[pMS->coord_xyz] = val;
		pMS->pFunc_setVector(pMS->dom,row,vec);
	}

	int MeshData::unifyVectorsOnMeshNodes(void (*pFunc_getVector)(int,int,double*),
			                              void (*pFunc_setVector)(int,int,const double*),
			                              GeomData* pGCData, int dom,
			                              bool onlyRemoteNodesOnBoundaries){

//		if (P_size()==1) return 0;
//
//		char tag[8]; sprintf(tag,"%d",dom);
//		pMS->tag = tag;
//		pMS->onlyRNOB = onlyRemoteNodesOnBoundaries;
//		pMS->dom = dom;
//		pMS->dim = pGCData->getMeshDim();
//		pMS->pFunc_getVector = pFunc_getVector;
//		pMS->pFunc_setVector = pFunc_setVector;
//
//		// unify for each vector coordinate
//		for (int i = 0; i < pMS->dim; i++){
//			pMS->coord_xyz = i;
//			unifyScalarsOnMeshNodes(0,0,pGCData,pMS);
//		}
		return 0;
	}

	int MeshData::unifyScalarsOnMeshNodes(double(*pFunc_getScalar)(pEntity),
			                               void (*pFunc_setScalar)(pEntity,double),
			                               GeomData* pGCData, void *ptr){
//		if (P_size()==1) return 0;
//
//		// STEP 1
//		// use a map to store only nodes with remote copies (this number is unknown)
//		map<int,double>::iterator mit;
//		pEntity node,face;
//		int ID,i,m,n,k,row;
//
//		if (structsCreation){
//
//			/*
//			 * pMS->onlyRNOB means one desires to unify vectors only nodes on boundaries
//			 */
//			if (ptr && pMS->onlyRNOB){
//				FIter fit = M_faceIter(theMesh);
//				while ( (face = FIter_next(fit)) ){
//					for (i=0;i<3;i++){
//						node = (pEntity)face->get(0,i);
//						if (0){
//							ID = get_AppToPETSc_Ordering(EN_id(node));					// node ID
//							//mapPB_nodes[ID] = (!ptr)?pFunc_getScalar(node):getScalar(node,pMS);
//							row = pSimPar->getLocalNodeIDNumbering(node,pMS->tag);
//							mapPB_nodes[ID] = (!ptr)?pFunc_getScalar(node):getScalar(row,pMS);
//
//
//						}
//					}
//				}
//				FIter_delete(fit);
//			}
//			else{
//				VIter vit = M_vertexIter(theMesh);
//				while ( (node = VIter_next(vit)) ){
//					if (0){
//						ID = get_AppToPETSc_Ordering(EN_id(node));					// node ID
//						//mapPB_nodes[ID] = (!ptr)?pFunc_getScalar(node):getScalar(node,pMS);
//						row = pSimPar->getLocalNodeIDNumbering(node,pMS->tag);
//						mapPB_nodes[ID] = (!ptr)?pFunc_getScalar(node):getScalar(row,pMS);
//					}
//				}
//				VIter_delete(vit);
//			}
//		}
//		else{
//			for(mit = mapPB_nodes.begin(); mit != mapPB_nodes.end(); mit++){
//				ID = get_PETScToApp_Ordering(mit->first);
//				node = (pEntity)theMesh->getVertex(ID);
//				if (!node) throw Exception(__LINE__,__FILE__,"Null vertex!\n");
//				//mit->second = (!ptr)?pFunc_getScalar(node):getScalar(node,pMS);
//				row = pSimPar->getLocalNodeIDNumbering(node,pMS->tag);
//				mit->second = (!ptr)?pFunc_getScalar(node):getScalar(row,pMS);
//			}
//		}
//
//		// STEP 2
//		// number of nodes on partition bdry
//		int numPB_Nodes = mapPB_nodes.size();
//		// nodes on partition bdry are now known. Let's transfer them to a PETSc column
//		// matrix to sum the contribution from all processor that share the same node.
//
//		int np = getNum_GNodes();
//		if (structsCreation){
//			ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,
//					PETSC_DECIDE,np,1,0,PETSC_NULL,0,PETSC_NULL,&joinNodes);CHKERRQ(ierr);
//		}
//		else{
//			ierr = MatZeroEntries(joinNodes);CHKERRQ(ierr);
//		}
//
//		int col = 0;
//		double data;
//		for(mit = mapPB_nodes.begin(); mit != mapPB_nodes.end(); mit++){
//			int row = mit->first-1;				// -1 to satisfy C/C++ index style
//			data = mit->second;
//			ierr = MatSetValues(joinNodes,1,&row,1,&col,&data,ADD_VALUES); CHKERRQ(ierr);
//		}
//		ierr = MatAssemblyBegin(joinNodes,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//		ierr = MatAssemblyEnd(joinNodes,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//
//
//		// STEP 3
//		// which rows from global matrix each rank must import
//		i = 0;
//		if (structsCreation){
//			rowToImport = new int[numPB_Nodes];
//			for (mit=mapPB_nodes.begin();mit!=mapPB_nodes.end();mit++)
//				rowToImport[i++] = mit->first-1;
//		}
//
//		if (structsCreation){
//			ierr = MatGetSubMatrixRaw(joinNodes,numPB_Nodes,rowToImport,1,&col,
//					PETSC_DECIDE,MAT_INITIAL_MATRIX,&updateValues); CHKERRQ(ierr);
//		}
//		else{
//			ierr = MatGetSubMatrixRaw(joinNodes,numPB_Nodes,rowToImport,1,&col,
//					PETSC_DECIDE,MAT_REUSE_MATRIX,&updateValues); CHKERRQ(ierr);
//		}
//
//		ierr = MatGetOwnershipRange(updateValues,&m,&n); CHKERRQ(ierr);
//		k = m;
//
//		for (i=0; i<numPB_Nodes; i++){
//			ierr = MatGetValues(updateValues,1,&k,1,&col,&data); CHKERRQ(ierr);
//			ID = rowToImport[i] + 1;
//			node = theMesh->getVertex( get_PETScToApp_Ordering(ID) );
//			if (!ptr)
//				pFunc_setScalar(node,data);
//			else{
//				row = pSimPar->getLocalNodeIDNumbering(node,pMS->tag);
//				setScalar(row,data,pMS);
//			}
//			k++;
//		}
//		structsCreation = false;
		return 0;
	}
}
