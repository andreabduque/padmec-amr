/*
 * EBFV1__pre-processors.cpp
 *
 *  Created on: 17/02/2012
 *      Author: rogsoares
 */

#include "EBFV1__pre-processors.h"

int EBFV1_preprocessor(pMesh theMesh, void *pData){
	cout<< "EBFV1_preprocessor"<<endl;
#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: Preprocessor\tIN\n";
#endif

	int ndom;
	if (theMesh->getDim()==2){
		EBFV1_preprocessor_2D(theMesh,pData,ndom);
	}
	else{
		EBFV1_preprocessor_3D(theMesh,pData,ndom);
	}

#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: Preprocessor\tOUT\n";
#endif
	return 0;
}

//void DijVector(pFace face, pVertex &oppositeVertex, std::vector<double> &Dij){
//	int i;
//	double xyz[3][3], xyzOV[3];
//
//	std::vector<pVertex> v;
//	M_GetVertices(face,v);
//
//	for (i=0 ;i<3; i++) Dij[i] = 0;
//	for (i=0 ;i<3; i++) V_coord( v[i], xyz[i] );
//	V_coord( oppositeVertex, xyzOV );
//
//	double a[3] = { xyz[1][0]-xyz[0][0], xyz[1][1]-xyz[0][1], xyz[1][2]-xyz[0][2] } ;
//	double b[3] = { xyz[2][0]-xyz[0][0], xyz[2][1]-xyz[0][1], xyz[2][2]-xyz[0][2] } ;
//	double c[3] = { xyzOV[0]-xyz[0][0], xyzOV[1]-xyz[0][1], xyzOV[2]-xyz[0][2] } ;
//
//	double _Dij[3]={.0,.0,.0};
//	computeCrossProduct(a,b,_Dij);
//
//	if ( computeDotProduct( c, _Dij ) > .0 ){
//		//	if ( innerProd > .0 ){
//		for (i=0;i<3;i++)
//			_Dij[i] = -_Dij[i];
//	}
//
//	for (i=0;i<3;i++) Dij[i] = _Dij[i]/6.0;
//	v.clear();
//}

// calcuates Dij vector: It points outer of domain
void DijVector(pFace face, pVertex &oppositeVertex, double *Dij){
	int i;
	double xyz[3][3], xyzOV[3];

	// get all vertices coordinates
	for (i=0 ;i<3; i++){
		V_coord( (pVertex)face->get(0,i), xyz[i] );
	}
	// get opposite vertex to boundary face coordinates
	V_coord( oppositeVertex, xyzOV );

	double a[3] = { xyz[1][0]-xyz[0][0], xyz[1][1]-xyz[0][1], xyz[1][2]-xyz[0][2] } ;	// vector over face's edge
	double b[3] = { xyz[2][0]-xyz[0][0], xyz[2][1]-xyz[0][1], xyz[2][2]-xyz[0][2] } ;	// vector over another face's edge
	double c[3] = { xyzOV[0]-xyz[0][0], xyzOV[1]-xyz[0][1], xyzOV[2]-xyz[0][2] } ;		// vector from face's center to opposite vertex

	// Dij: Dij = a x b
	double _Dij[3]={.0,.0,.0};
	computeCrossProduct(a,b,_Dij);

	// if inner product between Dij and c vector is negative, it means Dij points inside domain
	if ( computeDotProduct( c, _Dij ) > .0 ){
		for (i=0;i<3;i++){
			_Dij[i] = -_Dij[i];
		}
	}
	for (i=0;i<3;i++){
		Dij[i] = _Dij[i]/6.0;
	}
}

void getFCenter(pFace face, double *center){
	vector<pVertex> vertex;
	M_GetVertices(face,vertex);
	getFCenter(vertex[0],vertex[1],vertex[2],center);
	vertex.clear();
}

void getFCenter(pVertex v1, pVertex v2, pVertex v3, double *center){
	int i;
	double coords1[3], coords2[3], coords3[3];
	for (i=0; i<3; i++){
		coords1[i] = coords2[i] = coords3[i] = .0;
	}
	V_coord(v1,coords1);
	V_coord(v2,coords2);
	V_coord(v3,coords3);

	for (i=0; i<3; i++) 	center[i] = (coords1[i] + coords2[i] + coords3[i])/3.0;
}

double F_area(pFace face)
{
	int i;
	vector<pVertex> vertex;
	double xyz[3][3], n[3];

	M_GetVertices(face,vertex);
	for (i=0;i<3;i++) V_coord( vertex[i], xyz[i] );

	double a[3] = { xyz[1][0]-xyz[0][0], xyz[1][1]-xyz[0][1], xyz[1][2]-xyz[0][2] } ;
	double b[3] = { xyz[2][0]-xyz[0][0], xyz[2][1]-xyz[0][1], xyz[2][2]-xyz[0][2] } ;

	computeCrossProduct(a,b,n) ;

	return .5*norm(n);
}

void AllgatherDomains(std::set<int> &setOfDomain){
	int i = 0;
	int numLDomains = (int)setOfDomain.size();
	int domainsarray[numLDomains];
	std::set<int>::iterator iter = setOfDomain.begin();
	for (;iter != setOfDomain.end(); iter++) domainsarray[i++] = *iter;

	int numGDomains[P_size()];
	MPI_Gather(&numLDomains,1,MPI_INT,numGDomains,1,MPI_INT,0,MPI_COMM_WORLD);

	//	if (!P_pid()){
	//		for(i=0; i<P_size(); i++) printf("rank %d receives %d domains from rank %d\n",P_pid(),numGDomains[i],i);
	//	}

	// allocate enough space to receive nodes from all processors
	int *recv_buffer2, *displacements;
	int totalDoms = 0;
	if ( !P_pid() ){
		for(i=0; i<P_size(); i++) totalDoms += numGDomains[i];
		recv_buffer2 = new int[totalDoms]; // only root processor allocates memory
		displacements = new int[P_size()];
		displacements[0] = 0;
		for (int i=1; i<P_size(); i++) displacements[i] = displacements[i-1] + numGDomains[i-1];
	}


	// now it's time to send nodes to root processor
	MPI_Gatherv(domainsarray,numLDomains,MPI_INT,
			recv_buffer2,numGDomains,displacements,MPI_INT,
			0,MPI_COMM_WORLD);
	//	if (!P_pid()){
	//		for(i=0; i<totalDoms; i++) printf("rank %d domains %d\n",P_pid(),recv_buffer2[i]);
	//	}

	// let's filter domains flags to avoid repeated values
	setOfDomain.clear();
	if (!P_pid()){
		for(i=0; i<totalDoms; i++) setOfDomain.insert( recv_buffer2[i] );
	}
	//	printf("rank %d setOfDomain.size() = %d\n",P_pid(),setOfDomain.size());

	//	if (!P_pid()){
	//		for (iter = setOfDomain.begin(); iter != setOfDomain.end(); iter++) printf("rank %d domains %d\n",P_pid(),*iter);
	//	}

	// Send these domains flags to all processes

	i = 0;
	int numGDomains2 = (int)setOfDomain.size();
	numGDomains2 = P_getSumInt(numGDomains2);
	int domainsGarray[numGDomains2];

	//	if (!P_pid()){
	for (iter = setOfDomain.begin(); iter != setOfDomain.end(); iter++) domainsGarray[i++] = *iter;
	//	}

	MPI_Bcast(domainsGarray,numGDomains2,MPI_INT,0,MPI_COMM_WORLD);

	for(i=0; i<numGDomains2; i++) setOfDomain.insert( domainsGarray[i] );
	//printf("rank %d numGDomains2 %d\n",P_pid(),numGDomains2);

	//	for(i=0; i<numGDomains2; i++) printf("rank %d domains %d\n",P_pid(),domainsGarray[i]);
}

void generateFilename(string directory, const string &meshFilename, char *filename, string fileextension){
	string str(meshFilename);
	string::size_type start = str.find_last_of("/");
	start = (start == string::npos)?0:start;
	string::size_type end = str.find_last_of(".");

	char buffer[256];
	memset( buffer, '\0', 256 );
	str.copy(buffer,end-start,start);

	sprintf(filename,"%s%s-%d-of-%d%s",directory.c_str(),buffer,P_pid(),P_size(),fileextension.c_str());
}

void generateFilename(const string &meshFilename, char *filename, PetscBool& flg){
	string str(meshFilename);
	string::size_type start = str.find_last_of("/");
	start = (start == string::npos)?0:start;
	string::size_type end = str.find_last_of(".");

	char buffer[256];
	memset( buffer, '\0', 256 );
	str.copy(buffer,end-start,start);

	//	if (flg)
	//		sprintf(filename,"preprocessor-datafiles/%s.dat",buffer);
	//	else
	sprintf(filename,"preprocessor-datafiles/%s-%d-of-%d.dat",buffer,P_pid(),P_size());
}

// the following function exports geometric coefficients based on finite volume
// method described on Darlan thesis only. Coefficients are: Cij , Dij and nodal
// volume
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// IN
// theMesh - pointer to mesh object
// meshDim - mesh dimension 2D/3D
// setOfDomains - set container containing flags associaated to all sub-domains
// OUT
void exportCoefficients(pMesh theMesh, const int &meshDim, set<int> &setOfDomains,
		const char *meshFilename, GeomData *pGCData,PetscBool flg){
	//	PetscPrintf(PETSC_COMM_WORLD,"Exporting coefficients... \n");MPI_Barrier(MPI_COMM_WORLD);
	//	ofstream fid;
	//	fid.open(meshFilename);
	//
	//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Export pre-processor data file as: %s\n",meshFilename);
	//	PetscSynchronizedFlush(PETSC_COMM_WORLD);
	//	MPI_Barrier(MPI_COMM_WORLD);
	//
	//	fid << setprecision(15) << fixed;
	//
	//	// print in the begining of file all flags associated to all domains
	//	// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	//	std::set<int>::iterator dom;
	//	fid << "Mesh dimension:\n" << theMesh->getDim() << endl;
	//	fid << "Number of domains:\n" << setOfDomains.size() << endl;
	//	fid << "Domains' flag list:\n";
	//	for (dom=setOfDomains.begin(); dom!=setOfDomains.end(); dom++) fid << *dom << endl;
	//
	//	printOnFile_NodalVolume(fid,theMesh,meshDim,setOfDomains,pGCData,flg);
	//	printOnFile_Cij(fid,theMesh,meshDim,setOfDomains,pGCData,flg);
	//	printOnFile_Dij(fid,theMesh,meshDim,setOfDomains,pGCData);
	//
	//	// the line below was added here to guarantee compatibility with earlier
	//	// simulator versions which stop reading preprocessor file at the end of Dij
	//	// list coefficients.
	//	fid << "Number of global edges:\n" << getNumGlobalEdges(theMesh) << endl;
	//	PetscPrintf(PETSC_COMM_WORLD,"Exporting coefficients finished.\n\n");MPI_Barrier(MPI_COMM_WORLD);
	//	fid.close();
}


/*
 * Save mesh in gmsh ascii format version 1.0
 * As numeric formulation is based in an edge data structure, element connecti-
 * vities are not necessary after pre-process stage. For the visualization of pa-
 * rallel simulation results it becomes an issue. The following function should
 * be used to generate a mesh file with element (triangles or tetrahedrals) con-
 * nectivities and which will be processed by the vtk file converter program.
 */
void saveMesh_gmsh1(pMesh theMesh, const char* filename){

	ofstream fid;
	fid.open(filename);

	// Printf nodes list: ID - x - y - z
	// =========================================================================
	double coord[3];
	pEntity node, face, tetra;

	fid << "$NOD\n" << M_numVertices(theMesh) << std::endl;
	VIter vit = M_vertexIter(theMesh);
	while ( (node = VIter_next(vit)) ){
		V_coord(node,coord);
		fid << EN_id(node) << " " << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
	}
	fid << "$ENDNOD\n";


	const int dim = theMesh->getDim();

	// Printf nodes list: ID - x - y - z
	// =========================================================================

	// number of element to be printed:
	fid << "$ELM\n";
	int count = 0;
	if (dim==2){
		fid << M_numFaces(theMesh) << std::endl;
		FIter fit = M_faceIter(theMesh);
		while ( (face = FIter_next(fit)) ){
			int flag = GEN_tag(face->getClassification());
			fid << ++count << " 2 " << flag << " 1 3 ";
			for (int i=0; i<3; i++) fid << EN_id(face->get(0,i)) << " ";
			fid << std::endl;
		}
		FIter_delete(fit);
	}
	else if (dim==3){
		fid << M_numRegions(theMesh) << std::endl;
		RIter rit = M_regionIter(theMesh);
		while ( (tetra = VIter_next(rit)) ){
			int flag = GEN_tag(tetra->getClassification());
			fid << ++count << " 4 " << flag << " 1 4 ";
			for (int i=0; i<4; i++) fid << EN_id(tetra->get(0,i)) << " ";
			fid << std::endl;
		}
		RIter_delete(rit);
	}
	else
		printf("dimension %d not allowed. exiting........\n",dim);

	fid << "$ENDELM\n";
	fid.close();
}

double dot(const double *a, const double *b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Cross product of vectors made by p1p2,p1p3
// =========================================================================================
// IN:
// p1,p2,p3 - pointers for node coordenates
// OUT:
// normal - pointer for normal vector
void cross(const double *p1, const double *p2, const double *p3, double *normal){
	double v1[3], v2[3];
	makeVector(p2,p1,v1);
	makeVector(p3,p1,v2);
	return cross(v1,v2,normal);
}

// Cross product of vectors v1 and v2
// =========================================================================================
// IN:
// v1,v2 - pointers for vectors
// OUT:
// normal - pointer for normal vector
void cross(const double *v1, const double *v2, double* normal){
	normal[0] = (v1[1]*v2[2] - v1[2]*v2[1]); // x component
	normal[1] = (v1[2]*v2[0] - v1[0]*v2[2]); // y component
	normal[2] = (v1[0]*v2[1] - v1[1]*v2[0]); // z component
}

// during tetrahedrals loop to compute Cij coefficientts, faces on boundary are
// marked to inform, during Dij computation, between which domains they are.
// For external boundary, faces belong to only one domain. For internal faces,
// there are two domains.
// how face are marked:
//		1 - for some tetrahedral take all its vertices
//		2 - make a loop over all its faces
//		3 - if face is flagged (is on external[flag=1000] or internal boundary
//			[flag=999]) take all its vertices
//		4 - compare face's vertices with tetra' ones to find the opposite node
//			to face (it will be used to give Dij the corrected orientation)
//		5 - associate (mark) to face tetra's domain flag and opposite node ID
void markTetraBdryFaces(pMesh theMesh, pEntity tetra, const int &dom){
	vector<pVertex> tvertices, fvertices;
	M_GetVertices(tetra,tvertices);

	mVertex *v[4];
	for (int i=0; i<4; i++) v[i] = theMesh->getVertex(EN_id(tvertices[i]));

	mFace *face[4];
	face[0] = theMesh->getTri(v[0],v[1],v[2]);
	face[1] = theMesh->getTri(v[0],v[1],v[3]);
	face[2] = theMesh->getTri(v[0],v[2],v[3]);
	face[3] = theMesh->getTri(v[1],v[2],v[3]);

	for (int i=0; i<4; i++){
		if (face[i]){
			int oppositeVertexID = 0;
			M_GetVertices((pEntity)face[i],fvertices);

			// find opposite vertex
			for (int j=0; j<4; j++){
				bool key = false;

				for (int k=0; k<3; k++)
					if ( EN_id(tvertices[j]) == EN_id(fvertices[k]) ){
						key = true;
						break;
					}

				if (!key){
					oppositeVertexID = EN_id(tvertices[j]);
					break;
				}
			}

			fvertices.clear();
			//			printf("[%d] - face (%d,%d,%d) - oppositeVertexID: %d  ",P_pid(),id0,id1,id2,oppositeVertexID);
			//			printf("tetra(%d,%d,%d,%d) - \n",EN_id(tvertices[0]),EN_id(tvertices[1]),EN_id(tvertices[2]));

			// associate to face domain and opposite node ID
			int val = 0;
			EN_getDataInt((pEntity)face[i],MD_lookupMeshDataId("dom1"),&val);
			if (!val){
				EN_attachDataInt((pEntity)face[i],MD_lookupMeshDataId("ov1"),oppositeVertexID);
				EN_attachDataInt((pEntity)face[i],MD_lookupMeshDataId("dom1"),dom);
			}
			else{
				EN_attachDataInt((pEntity)face[i],MD_lookupMeshDataId("ov2"),oppositeVertexID);
				EN_attachDataInt((pEntity)face[i],MD_lookupMeshDataId("dom2"),dom);
			}
		}
	}
	tvertices.clear();
}

//double norm(const double *n){
//	return sqrt( pow(n[0],2) + pow(n[1],2) + pow(n[2],2) );
//}

int getNumGlobalEdges(pMesh theMesh){
	mEdge *edge;
	std::list<mEdge *> remoteEdges;
	std::list<mEdge *>::iterator reiter;

	const int rank = P_pid();
	const int N = P_size() - 1;
	int tag = 1234321;
	int *edgesIDsBuffer; edgesIDsBuffer=0;
	int i,p;
	int ecount = 0;
	int edgesIDsBufferSize;
	MPI_Status status;

	/**
	How things works:

	Edges on partition boundaries have one (2-D) or more (3-D) remote copies.
	A process of rank 'p' is the owner of an edge if any other process of rank
	'q', where p>q, has this edge, i.e, the lowest 'p' process will be always the
	owner of an edge.
	 */

	for (p=0; p<=N-1; p++){

		// receive data
		if ( rank == p+1 ){
			MPI_Recv(&ecount,1,MPI_INT,p,tag,MPI_COMM_WORLD,&status);
			MPI_Recv(&edgesIDsBufferSize,1,MPI_INT,p,tag,MPI_COMM_WORLD,&status);
			edgesIDsBuffer = new int[edgesIDsBufferSize];
			MPI_Recv(edgesIDsBuffer,edgesIDsBufferSize,MPI_INT,p,tag,MPI_COMM_WORLD,&status);
		}

		// send data
		/*
		Each rank counts its edges with remote copies, but counting starts in rank 0.
		The first 'n' edges from rank 0 are numbered from 1 to n.
		Rank 1, counts its edges from n+1 to k, but if an edge already exist in rank 0
		it doesn't count.
		 */
		if ( rank == p ){
			if ( rank == 0){
				// rank 0 gets all edges with remote copies and packs them into a list 'remoteEdges'
				EIter eit = M_edgeIter (theMesh);
				while ( (edge = (mEdge*)EIter_next (eit)) )
					if ( M_numRemoteCopies(theMesh,edge) ) remoteEdges.push_back( edge );
				EIter_delete(eit);

				i = 0;
				ecount = remoteEdges.size();
				edgesIDsBufferSize = 2*ecount;
				edgesIDsBuffer = new int[edgesIDsBufferSize];
				for(reiter = remoteEdges.begin(); reiter != remoteEdges.end(); reiter++){
					edgesIDsBuffer[2*i] = EN_id( (*reiter)->vertex(0) );
					edgesIDsBuffer[2*i+1] = EN_id( (*reiter)->vertex(1) );
					i++;
				}
			}
			else{
				// store array indices from those edges that do not exist the theMesh
				i=0;
				typedef std::pair<int,int> IDpairs;
				std::list< IDpairs > newEdgestoLocalMesh;
				while (i<edgesIDsBufferSize){
					edge = getEdge(theMesh,edgesIDsBuffer[i],edgesIDsBuffer[i+1]);
					if (!edge)
						newEdgestoLocalMesh.push_back( std::pair<int,int>(edgesIDsBuffer[i],edgesIDsBuffer[i+1]) );
					else
						markEdge(edge);
					i += 2;
				}
				delete [] edgesIDsBuffer; edgesIDsBuffer = 0;

				// 2º part
				// all ranks get all edges with remote copies and not marked, then pack them into 'remoteEdges'.
				EIter eit = M_edgeIter (theMesh);
				while ( (edge = (mEdge*)EIter_next (eit)) )
					if ( M_numRemoteCopies(theMesh,edge) && !isEdgeMarked(edge) ) remoteEdges.push_back( edge );
				EIter_delete(eit);

				// two IDs per edge => 2 times number of edges to be sent to next process
				edgesIDsBufferSize = 2*(newEdgestoLocalMesh.size() + remoteEdges.size());
				edgesIDsBuffer = new int[edgesIDsBufferSize];
				// ecount counts number of edges with remote copies without duplication
				ecount += remoteEdges.size();

				// transfer edges IDs from remoteEdges and newEdgestoLocalMesh containers to edgesIDsBuffer array
				i = 0;
				for(reiter = remoteEdges.begin(); reiter != remoteEdges.end(); reiter++){
					edgesIDsBuffer[2*i] = EN_id( (*reiter)->vertex(0) );
					edgesIDsBuffer[2*i+1] = EN_id( (*reiter)->vertex(1) );
					i++;
				}
				remoteEdges.clear();

				std::list<IDpairs>::iterator listIter;
				for(listIter = newEdgestoLocalMesh.begin(); listIter != newEdgestoLocalMesh.end(); listIter++){
					edgesIDsBuffer[2*i] = listIter->first;
					edgesIDsBuffer[2*i+1] = listIter->second;
					i++;
				}
				newEdgestoLocalMesh.clear();
			}

			if (rank <= N-1){
				MPI_Send(&ecount,1,MPI_INT,p+1,tag,MPI_COMM_WORLD);
				MPI_Send(&edgesIDsBufferSize,1,MPI_INT,p+1,tag,MPI_COMM_WORLD);
				MPI_Send(edgesIDsBuffer,edgesIDsBufferSize,MPI_INT,p+1,tag,MPI_COMM_WORLD);
			}
		}
	}
	if (!edgesIDsBuffer){
		delete [] edgesIDsBuffer; edgesIDsBuffer = 0;
	}

	// the round of communication makes the last rank (N-1) calculate the real number
	// of global edges. Then all other rank must take knowledge.
	ecount += getNumGlobalEdgesWithoutRemoteCopies(theMesh);

	int last_rank = P_size() - 1;

	// last rank sends the real number to all other ranks
	MPI_Bcast(&ecount,1,MPI_INT,last_rank,MPI_COMM_WORLD);

	return ecount;
}

mEdge* getEdge(pMesh theMesh, int id0, int id1){
	mVertex *v1 = theMesh->getVertex(id0);
	mVertex *v2 = theMesh->getVertex(id1);
	return (v1 && v2)?theMesh->getEdge(v1,v2):0;
}

int getNumGlobalEdgesWithoutRemoteCopies(pMesh theMesh){
	// each rank counts number of local edges without remote copies
	int egcount,elcount = 0;
	mEdge *edge;
	EIter eit = M_edgeIter (theMesh);
	while ( (edge = (mEdge*)EIter_next(eit)) )
		if ( M_numRemoteCopies(theMesh,edge)==0 ) elcount++;
	EIter_delete(eit);
	MPI_Allreduce(&elcount,&egcount,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	return egcount;
}

void markEdge(mEdge *edge){
	edge->attachInt(MD_lookupMeshDataId("marked"),1);
}

bool isEdgeMarked(mEdge *edge){
	return edge->getAttachedInt(MD_lookupMeshDataId("marked"));
}

void calculateEdgeLength(pMesh theMesh, GeomData *pGCData){
	int dim = theMesh->getDim();
	//double max_edge=.0, avgLength=.0;
	double elength, Lij[dim];
	double Icoord[3], Jcoord[3];
	double delta_x = 1e30; // infinity delta_x
	double versor[3];
	int i,j;
	dblarray vec(3);
	EIter eit = M_edgeIter(theMesh);
	while (pEntity edge = EIter_next(eit)){
		if (!theMesh->getRefinementDepth(edge)){
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
			//cout << "elength = " << elength << endl;
			pGCData->setEdgeLength(edge,elength);
			for (i=0; i<dim; i++){
				vec[i] /= elength;
			}
			//pGCData->setEdgeVec_Unitary(edge,vec);
			versor[0] = vec[0]; versor[1] = vec[1]; versor[2] = vec[2];
			pGCData->setVersor(edge,versor);

			// get the smallest one
			if (delta_x > elength){
				delta_x = elength;
			}
		}
	}
	EIter_delete(eit);//throw 1;
	// if parallel, get the smallest edge from all ranks and then broadcast it.
	pGCData->setSmallestEdgeLength(P_getMinDbl(delta_x));
}

// Calculate Cij Norm
void calculateCijNorm(pMesh theMesh, GeomData *pGCData, std::set<int> setOfDomains){
	int dim = theMesh->getDim();
	double Cij[3];
	std::set<int>::iterator iter = setOfDomains.begin();
	for(;iter!=setOfDomains.end();iter++){
		EIter eit = M_edgeIter(theMesh);
		pEdge edge;
		while ( (edge = EIter_next(eit)) ){
			if (!theMesh->getRefinementDepth(edge)){
				pGCData->getCij(edge,*iter,Cij);
				double Cij_norm = norm_L2(Cij,dim);
				if (Cij_norm > 0.0){
					pGCData->setCij_norm(edge,*iter,Cij_norm);
					int flag = EN_getFlag(edge);
					pGCData->set_belongsToBoundary(edge, flag != *iter );
				}
			}
		}
		EIter_delete(eit);
	}
}

void identifyBoundaryElements(pMesh theMesh, GeomData *pGCData, std::set<int> setOfDomains){
	bool detected=false;
	pEdge edge,face;
	std::set<int>::iterator iter;
	if (theMesh->getDim()==2){
		EIter eit = M_edgeIter(theMesh);
		while ( (edge = EIter_next(eit)) ){
			pGCData->set_belongsToBoundary(edge, false );
			iter = setOfDomains.find( EN_getFlag(edge) );
			if ( iter==setOfDomains.end() ){
				detected = true;
				pGCData->set_belongsToBoundary(edge, true );
			}

		}
		EIter_delete(eit);
	}
	else{
		FIter fit = M_faceIter(theMesh);
		while ( (face = FIter_next(fit)) ){
			pGCData->set_belongsToBoundary(face,false );
			iter = setOfDomains.find( EN_getFlag(face) );
			if ( iter==setOfDomains.end() ){
				detected = true;
				pGCData->set_belongsToBoundary(face, true );
			}
		}
		FIter_delete(fit);
	}

	if (!detected){
		throw Exception(__LINE__,__FILE__,"Any boundary element detected.");
	}
}

