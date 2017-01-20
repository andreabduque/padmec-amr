#include "exportVTK.h"
#include "CPU_Profiling.h"

using namespace PRS;

void exportSolutionToVTK(pMesh theMesh, void *pData1, void *pData2, void *pData3, void *pData4, string filename){
	CPU_Profile::Start();

	// open file for output
	// -------------------------------------------------------------------------
	ofstream fid;
	open_file(fid,filename,__LINE__,__FILE__);

	// initialize pointers
	// -------------------------------------------------------------------------
	PhysicPropData *pPPData  = (PhysicPropData*)pData1;
	throw_exception(!pPPData,"PhysicPropData *pPPData = NULL!",__LINE__,__FILE__);

	ErrorAnalysis *pEA  = (ErrorAnalysis*)pData2;
	throw_exception(!pEA,"ErrorAnalysis* pEA = NULL!",__LINE__,__FILE__);

	GeomData *pGCData = (GeomData*)pData4;
	throw_exception(!pGCData,"GeomData *pGCData = NULL!",__LINE__,__FILE__);
	
	// define constants
	// -------------------------------------------------------------------------
	int dim = pGCData->getMeshDim();
	int numElements = pGCData->getNumElements();
	int nnodes = pGCData->getNumNodes();

	// start printing
	// -------------------------------------------------------------------------
	print_headers(fid,nnodes);
	printVerticesCoordenates(fid,theMesh);
	printElementConnectivities(fid,theMesh,dim,numElements);
	printCellTypeList(fid,dim,numElements);

	// LIST HERE ALL NODAL FIELDS
	// -------------------------------------------------------------------
	fid << "\nPOINT_DATA "<< M_numVertices(theMesh) << endl;
	printPressure(fid,pGCData,pPPData);
	printSaturation(fid,pGCData,pPPData);
	//printPressureGradient(fid,pGCData,pPPData);
	//printSaturationGradient(fid,pGCData,pPPData);
	//printWeightedHeight(fid,pGCData,pEA);

	// LIST HERE ALL ELEMENT FIELDS
	// -------------------------------------------------------------------
//	fid << "\nCELL_DATA "<< pGCData->getNumElements() << endl;
//	printElementError(fid,pGCData,pEA);
//	print_h_ratio(fid,pGCData,pEA);
//	print_singular_regions(fid,pGCData,pEA);
//	print_elements_to_remove(fid,pGCData,pEA);
	
	fid.close();
	CPU_Profile::End("VTK");
}

void print_headers(ofstream &fid, int numNodes){
	fid << "# vtk DataFile Version 2.0\n";
	fid << "Two phases flow simulation\n";
	fid << "ASCII\nDATASET UNSTRUCTURED_GRID\n";
	fid << "POINTS " << numNodes << " float\n";
}

void printVerticesCoordenates(ofstream &fid, pMesh theMesh){
	pEntity e;
	int count = 0;
	VIter vit = M_vertexIter(theMesh);
	while( (e = VIter_next(vit)) ) {
		double coord[3] = {.0, .0, .0};
		V_coord(e,coord);
		for (int i=0; i<3; i++){
			fid << coord[i] << " ";
		}
		fid << endl;
		EN_attachDataInt(e,MD_lookupMeshDataId("mLN"),count++);
	}
	VIter_delete(vit);
}

void printElementConnectivities(ofstream &fid, pMesh theMesh, int dim, int numElements){
	fid << "\nCELLS " << numElements << " " << (dim+2)*numElements << endl;
	pEntity elem;
	if (dim==2){
		FIter fit = M_faceIter(theMesh);
		while( (elem = FIter_next(fit)) )
			if ( !theMesh->getRefinementDepth(elem) )
				printElementConnectivities(fid,elem,dim);
		FIter_delete(fit);
	}
	else{
		RIter rit = M_regionIter(theMesh);
		while ( (elem = RIter_next(rit)) )
			if ( !theMesh->getRefinementDepth(elem) )
				printElementConnectivities(fid,elem,dim);
		RIter_delete(rit);
	}
}

void printElementConnectivities(ofstream &fid, pEntity elem, int dim){
	int mappedLNodes;
	fid << dim+1 << " ";
	for(int i=0; i<dim+1; i++){
		EN_getDataInt((pVertex)elem->get(0,i),MD_lookupMeshDataId("mLN"),&mappedLNodes);
		fid << mappedLNodes << " ";
	}
	fid << endl;
}

void printCellTypeList(ofstream &fid, int dim, int numElements){
	fid << "\nCELL_TYPES " << numElements << endl;
	int type = (dim==2)?5:10;
	for(int i=0; i<numElements; i++){
		fid << type << endl;
	}
}

void printPressure(ofstream &fid, GeomData* pGCData, PhysicPropData *pPPData){
	fid << "SCALARS Pressure float 1\n";
	fid << "LOOKUP_TABLE default\n";
	fid << setprecision(8) << fixed;
	double p;
	for(int i=0; i<pGCData->getNumNodes(); i++){
		pPPData->getPressure(i,p);
		fid << p << endl;
	}
}

void printSaturation(ofstream &fid, GeomData* pGCData, PhysicPropData* pPPData){
	fid << "SCALARS Saturation float 1\n";
	fid << "LOOKUP_TABLE default\n";
	fid << setprecision(8) << fixed;
	double Sw;
	for(int i=0; i<pGCData->getNumNodes(); i++){
		pPPData->getSaturation(i,Sw);
		fid << Sw << endl;
	}
}

void printWeightedHeight(ofstream &fid, GeomData* pGCData, ErrorAnalysis *pEA){
	fid << "SCALARS WeightedHeight float 1\n";
	fid << "LOOKUP_TABLE default\n";
	fid << setprecision(8) << fixed;
	for(int i=0; i<pGCData->getNumNodes(); i++){
		fid << pEA->getWE_node(i) << endl;
	}
}

void printElementError(ofstream &fid, GeomData* pGCData, ErrorAnalysis *pEA){
	fid << "SCALARS Element_Error float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	int k = 0;
	for (int dom=0; dom<pGCData->getNumDomains(); dom++){
		for (int row=0; row<pGCData->getNumElemPerDomain(dom); row++){
			fid << pEA->getElementError(k++) << endl;
		}
	}
}

void printPressureGradient(ofstream& fid, GeomData* pGCData, PhysicPropData *pPPData){
	fid << "VECTORS p_grad float\n";
	double* p_grad = NULL;
	for (int dom=0; dom<pGCData->getNumDomains(); dom++){
		for (int node=0; node<pGCData->getNumNodesPerDomain(dom); node++){
			pPPData->get_pw_Grad(dom,node,p_grad);
			fid << p_grad[0] << " " << p_grad[1] << " " << p_grad[2] << endl;
		}
	}
}

void printSaturationGradient(ofstream& fid, GeomData* pGCData, PhysicPropData *pPPData){
	fid << "VECTORS Sw_grad float\n";
	double* Sw_grad = NULL;
	for (int node=0; node<pGCData->getNumNodes(); node++){
		pPPData->get_Sw_Grad(node,Sw_grad);
		fid << Sw_grad[0] << " " << Sw_grad[1] << " " << Sw_grad[2] << endl;
	}
}

void print_h_ratio(ofstream &fid, GeomData* pGCData, ErrorAnalysis *pEA){
	fid << "SCALARS hnew_by_hold float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	int k = 0;
	for (int dom=0; dom<pGCData->getNumDomains(); dom++){
		for (int row=0; row<pGCData->getNumElemPerDomain(dom); row++){
			fid << pGCData->getElem_HR(k++) << endl;
		}
	}
}

void print_singular_regions(ofstream &fid, GeomData* pGCData, ErrorAnalysis *pEA){
	fid << "SCALARS singular float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	int k = 0;
	for (int dom=0; dom<pGCData->getNumDomains(); dom++){
		for (int row=0; row<pGCData->getNumElemPerDomain(dom); row++){
			fid << pEA->isSingular(k++) << endl;
		}
	}
}

void print_elements_to_remove(ofstream &fid, GeomData* pGCData, ErrorAnalysis *pEA){
	fid << "SCALARS singular float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	int k = 0;
	for (int dom=0; dom<pGCData->getNumDomains(); dom++){
		for (int row=0; row<pGCData->getNumElemPerDomain(dom); row++){
			fid << pEA->isToRemove(k++) << endl;
		}
	}
}
