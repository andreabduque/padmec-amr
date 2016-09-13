/*
 * GeomData2.cpp
 *
 *  Created on: Oct 26, 2014
 *      Author: rogerio
 */

#include "GeomData.h"

namespace PRS{

	double GeomData::getElem_CDL(int i) const{
		return elem_CDL[i];
	}

	int GeomData::getNumNodes() const{
		return numNodes;
	}

	int GeomData::getNumElemPerDomain(int i) const{
		return numDomElem[i];
	}

	int GeomData::getNumEdgesPerDomain(int i) const{
		return numDomEdges[i];
	}

	int GeomData::getNumBDRYEdgesPerDomain(int i) const{
		return numDomBDRYEdges[i];
	}

	int GeomData::getNumBdryFacesPerDomain(int i) const{
		return numDomBDRYFaces[i];
	}

	int GeomData::getNumNodesPerDomain(int i) const{
		return numNodesPerDomain[i];
	}

	void GeomData::getEdge(int dom, int row, const int* &ptr){
		ptr = edges[dom].getrowconst(row);
	}

	void GeomData::getEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global){
		idx_0 = edges[dom].getValue(row,0);
		idx_1 = edges[dom].getValue(row,1);
		idx0_global = edges[dom].getValue(row,2);
		idx1_global = edges[dom].getValue(row,3);
	}

	void GeomData::getEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global, int &flag1, int &flag2){
		idx_0 = edges[dom].getValue(row,0);
		idx_1 = edges[dom].getValue(row,1);
		idx0_global = edges[dom].getValue(row,2);
		idx1_global = edges[dom].getValue(row,3);
		flag1 = edges[dom].getValue(row,4);
		flag2 = edges[dom].getValue(row,5);
	}

	void GeomData::getNodeIdx_Global(int dom, int i, int &idx){
		idx = nodes[dom].getValue(i);
	}

	void GeomData::setVolume(int idx, double v){
		volume_global[0].setValue(idx,v);
	}

	void GeomData::getVolume(int idx,double &v){
		v = volume_global[0].getValue(idx);
	}

	void GeomData::getVolume(int dom, int idx, double& vol){
		vol = volume[dom].getValue(idx);
	}

	void GeomData::getVolume(int dom, int idx_0, int idx_1, double& volumeI, double& volumeJ){
		volumeI = volume[dom].getValue(idx_0);
		volumeJ = volume[dom].getValue(idx_1);
	}

	int GeomData::getNodeID(int i_th) const{
		return pNodeID[i_th];
	}

	void GeomData::getID(int dom, int idx_0, int idx_1, int& id0, int &id1){
		id0 = ID[dom].getValue(idx_0);
		id1 = ID[dom].getValue(idx_1);
	}

	void GeomData::getID(int dom, int idx, int& id){
		id = ID[dom].getValue(idx);
	}

	void GeomData::getBdryEdge(int dom, int row, int &idx_0, int &idx_1){
		idx_0 = edges_bdry[dom].getValue(row,0);		// index number for boundary vertex ID for domain k
		idx_1 = edges_bdry[dom].getValue(row,1);		// index number for boundary vertex ID for domain k
	}

	void GeomData::getBdryEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global){
		idx_0 = edges_bdry[dom].getValue(row,2);		// index number for vertex ID for domain k
		idx_1 = edges_bdry[dom].getValue(row,3);		// index number for vertex ID for domain k
		idx0_global = edges_bdry[dom].getValue(row,4);	// global index number for vertex ID
		idx1_global = edges_bdry[dom].getValue(row,5);	// global index number for vertex ID
	}

	// return
	void GeomData::getBdryFace(int dom, int row, int &idx_0, int &idx_1, int &idx_2, int &idx0_global, int &idx1_global, int &idx2_global){
		idx_0 = faces_bdry[dom].getValue(row,3);		// index number for vertex ID for domain k
		idx_1 = faces_bdry[dom].getValue(row,4);		// idem
		idx_2 = faces_bdry[dom].getValue(row,5);		// idem
		idx0_global = faces_bdry[dom].getValue(row,6);	// global index number for vertex ID
		idx1_global = faces_bdry[dom].getValue(row,7);	// idem
		idx2_global = faces_bdry[dom].getValue(row,8);	// idem
	}

	void GeomData::getElement(int dom, int row, const int* &idx){
		idx = elem[dom].getrowconst(row);
	}

	void GeomData::getElement(int dom, int row, int* idx){
		for(int i=0;i<dim+1;i++){
			idx[i] = elem[dom].getValue(row,i);
			idx[i+dim+1] = elem[dom].getValue(row,i+dim+1);
		}
	}

	// get control volume for one vertex
	void GeomData::getBdryVolume(int dom, int idx, double& vol){
		vol = volume_bdry[dom].getValue(idx);
	}

	// get control volume for trhee vertices
	void GeomData::getBdryVolume(int dom, int idx0, int idx1, int idx2, double *vol){
		vol[0] = volume_bdry[dom].getValue(idx0);
		vol[1] = volume_bdry[dom].getValue(idx1);
		vol[2] = volume_bdry[dom].getValue(idx2);
	}

	// get control volume for three vertices
	void GeomData::getBdryVolume(int dom, int idx0, int idx1, int idx2, double& volumeI, double& volumeJ, double& volumeK){
		volumeI = volume_bdry[dom].getValue(idx0);
		volumeJ = volume_bdry[dom].getValue(idx1);
		volumeK = volume_bdry[dom].getValue(idx2);
	}

	// get control volume for two vertices
	void GeomData::getBdryVolume(int dom, int idx_0, int idx_1, double& volumeI, double& volumeJ){
		volumeI = volume_bdry[dom].getValue(idx_0);
		volumeJ = volume_bdry[dom].getValue(idx_1);
	}

	// get ID for two vertices on boundary
	void GeomData::getBdryID(int dom, int idx_0, int idx_1, int& id0, int &id1){
		id0 = ID_bdry[dom].getValue(idx_0);
		id1 = ID_bdry[dom].getValue(idx_1);
	}

	// get ID for three vertices on boundary
	void GeomData::getBdryID(int dom, int idx_0, int idx_1, int idx_2, int& id0, int &id1, int& id2){
		id0 = ID_bdry[dom].getValue(idx_0);
		id1 = ID_bdry[dom].getValue(idx_1);
		id2 = ID_bdry[dom].getValue(idx_2);
	}

	// get ID for vertex on boundary
	void GeomData::getBdryID(int dom, int idx, int& id){
		id = ID_bdry[dom].getValue(idx);
	}

	// get flag number defined by user on .geo file for domain i (i=0,1,2,...,n-1), where n is the number of domains
	int GeomData::getDomFlag(int i) const{
		return domainList[i];
	}

	// return edge length
	void GeomData::getLength(int dom, int idx, double &length) const{
		length = edge_length[dom].getValue(idx);
	}

	// return versor built over edge pointing from node I to node J, where Node ID I is ALWAYS less than node ID J.
	void GeomData::getVersor(int dom, int idx, double *v) const{
		v[0] = edge_versor[dom].getValue(idx,0);
		v[1] = edge_versor[dom].getValue(idx,1);
		v[2] = edge_versor[dom].getValue(idx,2);
	}

	void GeomData::setCij_norm(int dom, int idx, double val){
		Cij_norm[dom].setValue(idx,val);
	}

	void GeomData::getCij_norm(int dom, int idx, double &val){
		val = Cij_norm[dom].getValue(idx);
	}

	void GeomData::getVersor_ExternalBdryElement(int idx, double* versor){
		versor[0] = versor_ExtBdryElem[0].getValue(idx,0);
		versor[1] = versor_ExtBdryElem[0].getValue(idx,1);
		versor[2] = versor_ExtBdryElem[0].getValue(idx,2);
	}

	// Used by Saturation Gradient
	void GeomData::getExternalBdryEdges(int idx, int &idx0_global, int &idx1_global, int &flag1, int &flag2){
		idx0_global = external_bdry_elem[0].getValue(idx,0);
		idx1_global = external_bdry_elem[0].getValue(idx,1);
		flag1 = external_bdry_elem[0].getValue(idx,2);
		flag2 = external_bdry_elem[0].getValue(idx,3);
	}

	void GeomData::getExternalBdryFaces(int idx, int &idx0_global, int &idx1_global, int &idx2_global, int &flag1, int &flag2, int &flag3){
		idx0_global = external_bdry_elem[0].getValue(idx,0);
		idx1_global = external_bdry_elem[0].getValue(idx,1);
		idx2_global = external_bdry_elem[0].getValue(idx,2);
		flag1 = external_bdry_elem[0].getValue(idx,3);
		flag2 = external_bdry_elem[0].getValue(idx,4);
		flag3 = external_bdry_elem[0].getValue(idx,5);
	}

	// number of external boundary edges
	int GeomData::getNumExternalBdryEdges() const{
		return numExtBdryEdges;
	}

	// number of external boundary edges
	int GeomData::getNumExternalBdryFaces() const{
		return numExtBdryFaces;
	}

	void GeomData::setTotalNumberOfEdges(int n){
		numGEdges = n;
	}

	void GeomData::getTotalNumberOfEdges(int &n) const{
		n = numGEdges;
	}

	void GeomData::setMeshNodes(int n){
		numNodes = n;
	}

	void GeomData::getMeshNodes(int &n) const{
		n = numNodes;
	}

	void GeomData::setNumDomains(int n){
		_ndom = n;
	}

	int GeomData::getNumDomains() const{
		return _ndom;
	}

	const int* GeomData::getDomainList() const{
		return domainList;
	}

	int GeomData::getNumElements() const{
		return numElem;
	}

	void GeomData::getConnectivities(int row, int *connectivities){
		for(int i=0; i<elemtype; i++){
			connectivities[i] = pConnectivities->getValue(row,i);
		}
	}

	void GeomData::getCoordinates(int row, double *coords){
		for(int i=0; i<3; i++){
			coords[i] = pCoords->getValue(row,i);
		}
	}

	void GeomData::getCoordinates(int row, const double* &coords){
		coords = pCoords->getrowconst(row);
	}

	void GeomData::getCij(int dom, int row, double* cij){
		cij[0] = Cij[dom].getValue(row,0);
		cij[1] = Cij[dom].getValue(row,1);
		cij[2] = Cij[dom].getValue(row,2);
	}

	void GeomData::getCij(int dom, int row, const double* &cij){
		cij = Cij[dom].getrowconst(row);
	}

	void GeomData::setCij(int dom, int row, double* cij){
		Cij[dom].setValue(row,0,cij[0]);
		Cij[dom].setValue(row,1,cij[1]);
		Cij[dom].setValue(row,2,cij[2]);
	}

	void GeomData::getDij(int dom, int row, double* dij){
		dij[0] = Dij[dom].getValue(row,0);
		dij[1] = Dij[dom].getValue(row,1);
		dij[2] = Dij[dom].getValue(row,2);
	}

	void GeomData::getDij(int dom, int row, const double* &dij){
		dij = Dij[dom].getrowconst(row);
	}

	void GeomData::setDij(int dom, int row, double* dij){
		Dij[dom].setValue(row,0,dij[0]);
		Dij[dom].setValue(row,1,dij[1]);
		Dij[dom].setValue(row,2,dij[2]);
	}

	int GeomData::getNumFacesSharingVertex(int node) const{
		return numElemSharingVertex[node];
	}
}
