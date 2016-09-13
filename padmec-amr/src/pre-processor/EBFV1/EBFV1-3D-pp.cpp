/**
 * Edge based finite volume pre-processor for tetrahedral meshes. 
 * Implemented by Rogerio Soares. (NW5-2011)
 */

#include "EBFV1__pre-processors.h"

// Calculates Cij, Dij(boundary faces) and nodal volume. It also works for multi domains meshes
int EBFV1_preprocessor_3D(pMesh theMesh, void *pData, int &ndom){
	//PetscPrintf(PETSC_COMM_WORLD,"Starting EBFV1-3D pre-processor...");

	GeomData *pGCData = (GeomData*)pData;
	pGCData->setMeshDim(theMesh->getDim());
	if (theMesh->getDim() != 3){
		throw Exception(__LINE__,__FILE__,"Only 3-D meshes are allowed. Exiting...\n");
	}

	theMesh->modifyState(3,2);
	theMesh->modifyState(2,3);

	pEntity tetra, edge, face;
	int i,j,k,K,m;
	double tCenter[3], eCenter[3], v[3], *vec[2], val, *tetraFCenter[4];
	double normal[3], I[3], J[3], IJ[3], proj[3];

	// allocate vectors
	for (i=0; i<4; i++) tetraFCenter[i] = new double[3];
	for (i=0; i<2; i++) vec[i] = new double[3];

	// Search over all tetrahedra flags set to define domains. Store all flags on a list and initialize Cij vector
	std::set<int> setOfDomain;
	std::vector<double> Cij(3,.0);
	RIter rit = M_regionIter(theMesh);
	while ( (tetra = RIter_next(rit)) ){
		int flag = GEN_tag(tetra->getClassification());
		setOfDomain.insert( flag );
		for (i=0; i<6; i++){
			edge = (pEntity)tetra->get(1,i);
			pGCData->setCij(edge,flag,Cij);
		}
	}
	RIter_delete(rit);

	// mark faces (boundary)
	bool detectBdryFaces = false;
	std::set<int>::iterator iter = setOfDomain.begin();
	for (; iter != setOfDomain.end(); iter++){
		int countBDRYFaces = 0;
		int dom = *iter;
		char dom_string[256]; sprintf(dom_string,"dom_str_%d",dom);
		RIter rit = M_regionIter(theMesh);
		while ( (tetra = RIter_next(rit)) ){
			int tetraflag = GEN_tag(tetra->getClassification());
			if (tetraflag==dom){
				for (i=0;i<4;i++){
					face = (pFace)tetra->get(2,i);
					int faceflag = GEN_tag(face->getClassification());
					// if flags are different, then face in on boundary
					if (faceflag!=tetraflag){
						detectBdryFaces = true;
						countBDRYFaces++;
						// set edge belonging to domain 'dom'
						EN_attachDataInt(face,MD_lookupMeshDataId( dom_string ),1);
					}
				}
			}
		}
		RIter_delete(rit);
	}
	//throw 1;

	if (!detectBdryFaces){
		throw Exception(__LINE__,__FILE__,"Any boundary face (triangles) were detected. Boundary elements MUST have different flag of internal elements.");
	}

	/// initialize weighted volume
	VIter vit = M_vertexIter(theMesh);
	while (pEntity v = VIter_next(vit))
		pGCData->setWeightedVolume(v,0.0);
	VIter_delete(vit);

	vector<pVertex> tetraVertex(4), edgeVertex(2);
	double vt = .0; // total volume

	// loop over elements
	// for each element:	1 - calculate Cij for each edge and store them on the respective edge
	//						2 - calculate element contribution to volume of control volume
//	PetscPrintf(PETSC_COMM_WORLD,"Loop over tetras - start...  ");MPI_Barrier(MPI_COMM_WORLD);
	rit = M_regionIter(theMesh);
	while ( (tetra = RIter_next(rit)) ){
		// get all four tetrahedron's vertices
		for (i=0; i<4; i++){
			tetraVertex[i] = (pEntity)tetra->get(0,i);
		}

		int tetraFaces[4][3] = { {EN_id(tetraVertex[0]),EN_id(tetraVertex[1]),EN_id(tetraVertex[2])},
				{EN_id(tetraVertex[0]),EN_id(tetraVertex[1]),EN_id(tetraVertex[3])},
				{EN_id(tetraVertex[0]),EN_id(tetraVertex[2]),EN_id(tetraVertex[3])},
				{EN_id(tetraVertex[1]),EN_id(tetraVertex[2]),EN_id(tetraVertex[3])}};

		getFCenter(tetraVertex[0], tetraVertex[1], tetraVertex[2], tetraFCenter[0]);
		getFCenter(tetraVertex[0], tetraVertex[1], tetraVertex[3], tetraFCenter[1]);
		getFCenter(tetraVertex[0], tetraVertex[2], tetraVertex[3], tetraFCenter[2]);
		getFCenter(tetraVertex[1], tetraVertex[2], tetraVertex[3], tetraFCenter[3]);

		// get tetrahedron center
		tCenter[0] = tCenter[1] = tCenter[2] = .0;
		R_center(tetra,tCenter);

		double coord1[3]; V_coord(tetraVertex[0],coord1);
		double coord2[3]; V_coord(tetraVertex[1],coord2);
		double coord3[3]; V_coord(tetraVertex[2],coord3);
		double coord4[3]; V_coord(tetraVertex[3],coord4);

		// step #1: Cij
		// identify which domain the current tetrahedron belongs to
		// edges on domain's partition have differents Cij, one for each domain
		const int dom = GEN_tag(tetra->getClassification());
		char dom_string[256]; sprintf(dom_string,"dom_str_%d",dom);
		markTetraBdryFaces(theMesh,tetra,dom);
		setOfDomain.insert(dom); // store every new domain

		for (int pos1=0; pos1<3; pos1++){
			for (int pos2=pos1+1; pos2<4; pos2++){
				edge = (pEdge)theMesh->getEdge((mVertex*)tetraVertex[pos1],(mVertex*)tetraVertex[pos2]);

				// set edge belonging to domain 'dom'
				EN_attachDataInt(edge,MD_lookupMeshDataId( dom_string ),1);


				//M_GetVertices(edge,edgeVertex);
				for (i=0; i<2; i++) edgeVertex[i] = (pEntity)edge->get(0,i);
				// edge's Cij is the sum of two vectors. Both are orthogonals to planes
				// defined by three vectors. They are:
				// v - ec->tc:		edge center to tetra center
				// vec[0] - ec->fcr:	edge center to right face center
				// vec[1] - ec->fcl:	edge center to left face center
				// Cij = (v)x(vec[0]) + (v)x(vec[1]), Cij points to outside of control volume

				// create vector v:
				for (i=0; i<3; i++) I[i] = J[i]= .0;
				V_coord(edgeVertex[0],I);
				V_coord(edgeVertex[1],J);

				double sign = ( EN_id(edgeVertex[0]) > EN_id(edgeVertex[1]) )?-1.0:1.0;
				for (i=0; i<3; i++){
					IJ[i] = sign*(J[i] - I[i]);
				}

				// get edge center
				eCenter[0] = eCenter[1] = eCenter[2] = .0;
				for (i=0; i<3; i++){
					eCenter[i] = .5*(I[i]+J[i]);
				}

				// create vector v: from edge middle point to tetrahedral centroid
				for (i=0; i<3; i++){
					v[i] = tCenter[i] - eCenter[i];
				}

				// search for tetra faces that share the same edge and get their centers. Of course, there are only two faces
				// sharing the same edge. FMDB can do this job easily, but for 3D big meshes a face structure has a high memory cost.

				K = 0;
				// loop over tetra's faces
				for (i=0; i<4; i++){
					// loop over face's vertices
					for (j=0; j<3; j++){
						if (EN_id(edgeVertex[0]) == tetraFaces[i][j]){
							for (k=0; k<3; k++){
								if (EN_id(edgeVertex[1])==tetraFaces[i][k]){
									for (m=0; m<3; m++) vec[K][m] = -eCenter[m] + tetraFCenter[i][m];
									K++;
								}
							}
						}
					}
				}

				// IJ vector is a reference to Cij. IJ points from node I to node J
				// where: I_id < J_id
				// vec projections on edge IJ must have the same orientation as IJ vector
				double n = sqrt( IJ[0]*IJ[0] + IJ[1]*IJ[1] + IJ[2]*IJ[2] ) ;

				// Cij calculation n ;
				for (j=0; j<3; j++){
					Cij[j] = .0;
				}
				pGCData->getCij(edge,dom,Cij);
				double normal1[3], normal2[3];
				cross(v,vec[0],normal1);
				cross(vec[1],v,normal2);

				for (j=0; j<3; j++){
					proj[j] = .0;
				}
				for (j=0; j<3; j++){
					normal[j] = .5*(normal1[j]+normal2[j]);
				}
				val = dot(normal,IJ)/(n*n);
				for (j=0; j<3; j++){
					proj[j] = val*IJ[j];
				}
				double sinal = (dot(proj,IJ)<.0)?-1.:1.;
				for (j=0; j<3; j++){
					Cij[j] += sinal*normal[j];
				}
				pGCData->setCij(edge,dom,Cij);
			}
		}

		// step #2: volume of control volume
		const double porosity = .0;
		const double tetraVolume = R_Volume(tetra);
		const double nodalVolume = .25*tetraVolume;
		vt = +tetraVolume;

		for (i=0; i<4; i++){
			// total nodal volume
			double v1 = pGCData->getVolume(tetraVertex[i],dom);
			v1 += nodalVolume;
			pGCData->setVolume(tetraVertex[i],dom,v1);
			double v2 = pGCData->getWeightedVolume(tetraVertex[i]);
			v2 += nodalVolume*porosity;
			pGCData->setWeightedVolume(tetraVertex[i],v2);
		}
	}
	RIter_delete(rit);	// END TETRAHEDRALS LOOP
	pGCData->setTotalReservoirVolume(vt);
	//PetscPrintf(PETSC_COMM_WORLD,"Finished\n");MPI_Barrier(MPI_COMM_WORLD);

	// deallocate vectors
	for (i=0; i<4; i++) delete[] tetraFCenter[i];
	for (i=0; i<2; i++) delete[] vec[i];

	/*
	 * Dij vector calculation. A loop over all boundary (external/internal) faces is made. Depending on how FMDB is used, faces on all
	 * tetrahedrals may be created and they must be filtered. If a face does not belong to boundary it will assume the tetrahedral domain
	 * flag. This value is greater than 3000 and is used to filter them from those provided by mesh file.
	 */

	dblarray Dij(3,.0);
	//PetscPrintf(PETSC_COMM_WORLD,"Loop over boundary faces - start... ");///MPI_Barrier(MPI_COMM_WORLD);
	if ( M_numFaces(theMesh) != 0 ){
		FIter fit = M_faceIter(theMesh);
		while ( (face = FIter_next(fit)) ){
			int flag = GEN_tag(face->getClassification());
			iter = setOfDomain.find(flag);
			// get only faces built over surfaces defined by user in geometric model
			// it only work if surface flags defined by user are different from those set to tetrahedrals
			if ( iter==setOfDomain.end() ){
				computeDij(theMesh,face,pGCData);
			}
		}
		FIter_delete(fit);
	}
	//PetscPrintf(PETSC_COMM_WORLD,"Finished\n");MPI_Barrier(MPI_COMM_WORLD);
	calculateEdgeLength(theMesh,pGCData);
	calculateCijNorm(theMesh,pGCData,setOfDomain);
	identifyBoundaryElements(theMesh,pGCData,setOfDomain);
	i = 0;
	ndom = (int)setOfDomain.size();
	int *domlist = new int[ndom];
	for(iter = setOfDomain.begin(); iter!=setOfDomain.end(); iter++){
		domlist[i++] =  *iter;
	}

	//PetscPrintf(PETSC_COMM_WORLD," finished.\n");
	//cout << "-------------------------------------------------------------------------------------------------------------------------\n\n";
	return 0;
}

// before use computeDij(), faces on boundaries (external/internal) should be marked with markTetraBdryFaces()
void computeDij(pMesh theMesh, pEntity face, GeomData *pGCData){
	int ID1=0, dom1=0,ID2=0, dom2=0;
	EN_getDataInt((pEntity)face,MD_lookupMeshDataId("ov1"),&ID1);
	EN_getDataInt((pEntity)face,MD_lookupMeshDataId("dom1"),&dom1);
	EN_getDataInt((pEntity)face,MD_lookupMeshDataId("ov2"),&ID2);
	EN_getDataInt((pEntity)face,MD_lookupMeshDataId("dom2"),&dom2);

	pVertex oppositeVertex = theMesh->getVertex(ID1);
	// NOTE: Dij points outside to dom1
	if ( oppositeVertex ){
		double Dij[3] = {.0,.0,.0};
 		DijVector(face,oppositeVertex,Dij);
		pGCData->setDij(face,dom1,dom2,Dij);
	}
	else{
		throw Exception(__LINE__,__FILE__,"Problems to find opposite vertex to face.\n");
	}
}
