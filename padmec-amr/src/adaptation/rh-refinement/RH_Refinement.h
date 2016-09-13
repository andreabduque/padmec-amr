/**
 * * RH_Refinement.h
 ** Created: 07-06-2013
 ** Author: Guilherme Caminha ( gpkc @ cin.ufpe.br )
 ** 
 ** RH-Refinement module using MAdLib (Mesh Adaption Library), which performs global node repositioning
 ** and mesh adaptation.
 **/


#ifndef RH_Refinement_Include
#define RH_Refinement_Include

#include "includes.h"

#include "AMR.h"

#ifndef NOADAPTATION
//----------------------------------------------------------------
#include "ModelInterface.h"
#include "MeshDataBaseInterface.h"
#include "AdaptInterface.h"
#include "PWLinearSField.h"
#include "MAdResourceManager.h"
//----------------------------------------------------------------
#endif

#define MEM 1

#define LINUX (defined(__linux) || defined(linux)) && MEM

class RH_Refinement : public AMR {
public:
	RH_Refinement(int argc, char* argv[]);
	~RH_Refinement();
	
	//void rodar(ErrorAnalysis *pErrorAnalysis, pMesh & theMesh);
	void run(pMesh theMesh, std::list<pEntity>&, std::set<pEntity>&);
	
#ifndef NOADAPTATION
private:
	std::list<MAd::pEntity> listaver;
	
	//Dimension
	int dim;
	
	//Mesh and Model pointers
	MAd::pMesh MAdMesh;
	MAd::pGModel MAdModel;
	
	//Resource Manager
	MAd::MAdResourceManager ResMngr;
	double t0;
#if LINUX
	void displayMemoryUsage(std::string step="");
	void process_mem_usage(double& vm_usage, double& resident_set);
#endif
	
	//Mesh import functions (from FMDB mesh or .msh)
	bool importMesh(pMesh & theMesh);
	void importData(pMesh & theMesh);
	bool importMesh(std::string filename);
	
	//Mesh export functions (to FMDB mesh or .msh)
	void exportMesh(pMesh & theMesh);
	bool exportMesh(std::string filename);
	
	//Ligação entre ID's do FMDB e do MAdLib
	std::map<int,int> MAdToMeshIds;
	std::map<int,int> MeshToMAdIds;
	
	//CallBack function for data interpolation
	void PADMEC_CBFunction(MAd::pPList before, MAd::pPList after, void *data,
			 MAd::operationType type, MAd::pEntity ppp);
	
	
	//The SizeField used in the current operation
	MAd::PWLSField * sizeField;
	
	//Set the mesh size field
	void makeSizeField(std::set<pEntity>);
	
	//The Mesh Adapter
	MAd::MeshAdapter* adapter;
#endif
public:
	//Clear all data
	void clear();
};

#endif
