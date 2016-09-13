/*
 * export-VTK2.cpp
 *
 *  Created on: Sep 23, 2014
 *      Author: Rogerio
 */

#include "exportVTK.h"

void exportSolutionToVTK(string& filename, const PhysicPropData* pPPData, const GeomData* pGCData){
	double t1 = MPI_Wtime();
	ofstream fid;
	fid.open(filename.c_str());
	checkFileOpened(fid);
	print_headers(fid,pPPData->getNumNodes());
	print_coordinates(fid,pGCData);
	print_connectivities(fid,pGCData);
	print_celltypes(fid,pGCData->getDim(),pGCData->getNumElements());
	print_saturation(fid,pPPData);
	print_pressure(fid,pPPData);
	fid.close();
	cout << setprecision(2) << scientific << "Time to print VTK file: " << MPI_Wtime - t1 << "[s]\n";
}



void print_Coordinates(ofstream &fid, const GeomData* pGCData){
	//double coords[3];
	const double* coords = NULL;
	int nnodes;
	pGCData->getMeshNodes(nnodes);
	for(int i=0; i<nnodes; i++){
		pGCData->getCoordinates(i,coords);
		fid << coord[0] << " " << coord[1] << " " << coord[2] << endl
	}
}

void print_connectivities(ofstream &fid, const GeomData* pGCData){
	int connectivities[4];
	int dim = pGCData->getDim();
	int numElem = pGCData->getNumElements();
	fid << "\nCELLS " << numElem << " " << (dim+2)*numElem << endl;
	for (int i=0; i<numElem; i++){
		pGCData->getConnectivities(i,connectivities);
		for (int j=0; j<dim; j++){
			fid << connectivities[i] << " ";
		}
		fid << endl;
	}
}

void print_celltype(ofstream& fid, int dim, int numElem){
	fid << "\nCELL_TYPES " << numElem << endl;
	int type = (dim==2)?5:10;				// 5 - triangles;	10 - tetrahedra
	for(int i=0; i<numElem; i++){
		fid << type << endl;
	}
}

void print_saturation(ofstream &fid, const PhysicPropData* pPPData){
	double Sw;
	fid << "\nPOINT_DATA "<< M_numVertices(theMesh) << endl;
	fid << "SCALARS Saturation float 1\n";
	fid << "LOOKUP_TABLE default\n";
	for (int idx=0; idx<pPPData->getNumNodes(); idx++){
		pPPData->getSaturation(idx,Sw);
		fid << Sw << endl;
	}
}

void print_pressure(ofstream &fid, const PhysicPropData* pPPData){
	double p;
	fid << "SCALARS Pressure float 1\n";
	fid << "LOOKUP_TABLE default\n";
	for (int idx=0; idx<pPPData->getNumNodes(); idx++){
		pPPData->getPressure(idx,p);
		fid << p << endl;
	}
}

void checkFileOpen(ofstream& fid){
	if ( !fid.is_open() ){
		char msg[512];
		sprintf(msg,"File '%s' could not be opened or it does not exist.\nCheck if directory was typed correctly.\n",filename.c_str());
		throw Exception(__LINE__,__FILE__,msg);
	}
}
