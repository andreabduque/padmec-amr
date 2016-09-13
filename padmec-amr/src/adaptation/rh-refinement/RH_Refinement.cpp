/**
 * * RH_Refinement.cpp
 ** Created: 07-06-2013
 ** Author: Guilherme Caminha ( gpkc @ cin.ufpe.br )
 **/

#define DEBUG_CHECK printf("\nDEBUG(pressione enter para continuar)"); \
getchar()

#include "RH_Refinement.h"

#include "MAdLib.h"
#include "MeshDataBaseInterface.h"
#include "MeshDataBaseIO.h"
#include "IsoMeshSize.h"
#include <unistd.h>

#include "ModelInterface.h"
#include "AdaptInterface.h"

RH_Refinement::RH_Refinement(int argc, char* argv[]) {
	cout << endl << "======================";
	cout << endl << "STARTING MADLIB MODULE";
	cout << endl << "======================" << endl;
	
	//Initialize MAdLib
	MAdLibInitialize(&argc,&argv);
	
	MAdMesh = NULL;
	MAdModel = NULL;
	sizeField = NULL;
	adapter = NULL;
	GM_create(&MAdModel, "theModel");
	
	displayMemoryUsage("MAdLib initialized");
}


RH_Refinement::~RH_Refinement() {
	if(MAdMesh) {
		M_delete(MAdMesh);
	}
	if(MAdModel) {
		GM_delete(MAdModel);
	}
}

//Run the adaptation procedure
//void RH_Refinement::rodar(ErrorAnalysis *pErrorAnalysis, pMesh & theMesh) {
void RH_Refinement::run(pMesh theMesh, std::list<pEntity>& elementList, std::set<pEntity>& elementSet){
	cout << endl << "=========================================";
	cout << endl << "STARTING MADLIB MESH ADAPTATION PROCEDURE";
	cout << endl << "=========================================" << endl;
	
	t0 = ResMngr.getTime();
	
	if (!theMesh){
		throw Exception(__LINE__,__FILE__,"NULL Mesh!");
	}
	
	// NAO CHAMAR ESTA FUNCAO DAQUI. JA VEM DE FORA
//	std::list<pEntity> elementList;
//	std::set<pEntity> elementSet;
	//pErrorAnalysis->getRefUnrefElementsList(theMesh,elementList,elementSet);

	cout << "Number of elements to be (un)refined: " << elementList.size() << endl;
	
	//Copy the mesh from a FMDB mesh to a MADMESH
	if(!MAdMesh)
		importMesh(theMesh);
	importData(theMesh);
	
	//Create the SizeField on the mesh
	makeSizeField(elementSet);
	
	M_writeMsh(MAdMesh, "MAD1.msh");
	
	//Build the adaptation tool
	adapter = new MAd::MeshAdapter(MAdMesh, sizeField);
	
	for(std::list<MAd::pEntity>::iterator it = listaver.begin(); it != listaver.end(); it++) {
		adapter->setConstraint(*it);
	}
	
	adapter->setCollapseOnBoundary(false,0.);
	adapter->setSwapOnBoundary(false, 0.);
	
	displayMemoryUsage("Mesh adapter built");
	
	adapter->run();
	
	adapter->printStatistics(std::cout);
	M_writeMsh(MAdMesh,"MAD2.msh", 1);
	
	delete adapter;
	delete sizeField;
	M_delete(MAdMesh);
	GM_delete(MAdModel);
	MAdModel = NULL;
	MAdMesh = NULL;
	sizeField = NULL;
	adapter = NULL;
	listaver.clear();
	GM_create(&MAdModel, "theModel");
	/*
	 * GM_create(&MAdModel, "theModel");
	 * MAdMesh = M_new(MAdModel);
	 * MAd::M_load(MAdMesh, "MAD1.msh");
	 * 
	 * cout << "numFaces: " <<  M_numFaces(MAdMesh) << endl;
	 */
	
	theMesh = MS_newMesh(0);
	M_load(theMesh,"MAD2.msh");
	
	cout << "Adaptation procedure ocurred in " << ResMngr.getTime() - t0 << " seconds";
	displayMemoryUsage("After adaptation");
	
	//STOP();
}

bool RH_Refinement::importMesh(pMesh & theMesh) {
	cout << "Loading the mesh..." << endl;
	double cpu_mesh_0 = ResMngr.getTime();
	
	//Used for vertex coordinates
	double xyz[3];
	
	//Clear id-to-id relation map
	MAdToMeshIds.clear();
	MeshToMAdIds.clear();
	
	//Start looping through vertices
	int nVerts = M_numVertices(theMesh);
	
	if(!nVerts) return false;
	
	//Create the mesh
	if(!MAdMesh)
		MAdMesh = M_new(MAdModel);
	
	//Build the vertices
	VIter vit = M_vertexIter(theMesh);
	pEntity ent;
	
	for(int ID = 1; (ID <= nVerts) && (ent = VIter_next(vit)); ID++) {
		xyz[0] = 0.;
		xyz[1] = 0.;
		xyz[2] = 0.;
		
		V_coord(ent, xyz);
		
		/*MAd::pVertex pv = */MAd::M_createV2(MAdMesh, xyz, ID);
		
		MeshToMAdIds[EN_id(ent)] = ID;
		MAdToMeshIds[ID] = EN_id(ent);
	}
	
	VIter_delete(vit);
	
	//Build the Elements
	dim = theMesh->getDim();
	
	if(dim==3) {
		//3D
		
	}
	else if (dim==2) {
		//2D
		FIter fit = M_faceIter(theMesh);
		vector<pVertex> v;
		int ids[3];
		
		while(pEntity ent = FIter_next(fit)) {
			MAd::pGFace geom = MAd::GM_faceByTag(M_model(MAdMesh), GEN_tag(EN_whatIn(ent)));
			ids[0] = MeshToMAdIds[ EN_id(F_vertex(ent, 0)) ];
			ids[1] = MeshToMAdIds[ EN_id(F_vertex(ent, 1)) ];
			ids[2] = MeshToMAdIds[ EN_id(F_vertex(ent, 2)) ];
			//cout << ids[0] << " " << ids[1] << " " << ids[2] << endl;
			/*MAd::pFace pf = */MAd::M_createF(MAdMesh, 3, ids, (MAd::pGEntity)geom);
			
			
			
			
			
			//cout << MAd::EN_id((MAd::pEntity)F_vertex(pf, 0)) << " " << MAd::EN_id((MAd::pEntity)F_vertex(pf, 1)) << " " << MAd::EN_id((MAd::pEntity)F_vertex(pf, 2)) << endl;
			
			v.clear();
		}
		FIter_delete(fit);
		
		MAd::EIter eit = MAd::M_edgeIter(MAdMesh);
		while(MAd::pEdge ent = MAd::EIter_next(eit)) {
			if( MAd::E_numFaces(ent) < 2)  {
				listaver.push_back((MAd::pEntity)ent);
				//MAd::E_vertex(ent, 0);
				//MAd::E_vertex(ent, 1);
			}
		}
	}
	
	
	
	M_classifyEntities(MAdMesh);
	
	double cpu_mesh_tot = ResMngr.getTime() - cpu_mesh_0;
	cout << "Loaded the mesh in " << cpu_mesh_tot << " seconds\n";
	
	displayMemoryUsage("Mesh loaded");
	
	return true;
}


void RH_Refinement::displayMemoryUsage(std::string step)
{
	#if LINUX
	double ram, virt;
	
	process_mem_usage(virt,ram);
	
	#ifdef PARALLEL
	double send[2] = {virt,ram};
	double recv[2];
	MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	virt  = recv[0];
	ram = recv[1];
	#endif
	
	std::cout << endl << "Memory usage at step \'"<<step<<"\': " 
	<< ram  << " Mb (resident), "
	<< virt << " Mb (virtual)\n" << endl;
	#endif
}

#if LINUX
void RH_Refinement::process_mem_usage(double& vm_usage, double& resident_set)
{
	using std::ios_base;
	using std::ifstream;
	using std::string;
	
	vm_usage     = 0.0;
	resident_set = 0.0;
	
	// 'file' stat seems to give the most reliable results
	//
	ifstream stat_stream("/proc/self/stat",ios_base::in);
	
	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;
	
	// the two fields we want
	//
	unsigned long vsize;
	long rss;
	
	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
	>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
	>> utime >> stime >> cutime >> cstime >> priority >> nice
	>> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
	
	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage     = vsize / 1024.0 / 1024;
	resident_set = rss * page_size_kb / 1024;
}
#endif

bool RH_Refinement::importMesh(std::string filename) {
	
}

void RH_Refinement::makeSizeField(std::set<pEntity> elementList) {
	cout << "Starting SizeField Creation" << endl;
	sizeField = new MAd::PWLSField(MAdMesh);
	
	sizeField->setCurrentSize();
	
	for(std::set<pEntity>::iterator it = elementList.begin(); it != elementList.end(); ++it) {
		double height;
		EN_getDataDbl(*it, MD_lookupMeshDataId( "elem_height" ), &height);
		
		MAd::pVertex pv = M_vertex(MAdMesh, MeshToMAdIds[ EN_id(*it) ]);
		
		MAd::IsoMeshSize* ps = (MAd::IsoMeshSize*)sizeField->findSize(pv);
		
		void * temp1 = 0;
		if(ps->size() < height) {
			// 			if(dim == 2) {
			// 				MAd::pPList pl = V_edges(pv);
			// 				cout << V_numEdges(pv) << endl;
			// 				while(MAd::pEdge pe = (MAd::pEdge)PList_next(pl, &temp1)) {
			// 					cout << EN_id((pEntity)pe) << endl;
			// 					if(E_numFaces(pe) ==1) continue;
			// 				}
			// 				
			// 				cout << ps->size() << " " << height << endl;
			// 			}
			//continue;
		}
		
		sizeField->setSize((MAd::pEntity)pv, height);
	}
}

void RH_Refinement::importData(pMesh & theMesh) {
	
	VIter vit = M_vertexIter(theMesh);
	
	// 	while(pEntity ent = VIter_next(vit)) {
	// 		double height;
	// 		EN_getDataDbl(*it, MD_lookupMeshDataId( "elem_height" ), &height);
	// 		
	// 		MAd::pVertex pv = M_vertex(MAdMesh, MeshToMAdIds[ EN_id(*it) ]);
	// 		
	// 		MAd::IsoMeshSize* ps = (MAd::IsoMeshSize*)sizeField->findSize(pv);
	// 		
	// 		if(ps->size() < height) {
	// 			cout << ps->size() << " " << height << endl;
	// 			continue;
	// 		}
	// 		
	// 		sizeField->setSize((MAd::pEntity)pv, height);
	// 	}
	
}

void RH_Refinement::clear() {
	
}
