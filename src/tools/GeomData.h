#ifndef GEOMETRICCOEFFICIENTSDATA_H_
#define GEOMETRICCOEFFICIENTSDATA_H_

#include "AttachData.h"
#include "Matrix.h"


namespace PRS
{
	typedef std::vector<double> DataArray;
	typedef map<int,DataArray> MapDataArray;
	typedef map<int,double> MapData;
	typedef map<int,double>::iterator MData_Iter;

	/*! \struct Coefficients GeomData.h
	 *  \brief Set of data that is associated to some specific mesh entity.
	 */
	struct Coefficients{
		MapDataArray Cij;
		DataArray Dij;
		MapData volume;
		double weightedvolume;

		DataArray sumIJ;

		// domains flag for Dij
		int dom1;
		int dom2;

		DataArray edgeVector;
		DataArray edgeVectorUnitary;
		double length;
		double Cij_norm;

		int numRemoteCopies;
		std::map<int,int> numRC; // number of remote copies associated to a domain
		int flag;	// may have several meanings
		bool isBoundary;
	};


	/*! \class GeomData GeomData.h
	 *  \brief GeomData is designed to set/get set of data associated to mesh entities.
	 *  The set of data is defined by the structure Coefficients above.
	 */
	class GeomData : public AttachData
	{
	public:

		GeomData();
		~GeomData();

		void initilize(pMesh);
		void initilize(pMesh, const std::set<int>&);
		void dataTransfer(pMesh);
		void transferCijData(pMesh);
		void transferDijData(pMesh);
		void transferVolData(pMesh);

		// get Dij for 3-D meshes (the old one generate memory linking)
		bool getDij(pEntity face, int dom, double *Dij);
		void setDij(pEntity face, int dom1, int dom2, double *Dij);
		void getCij(pEntity , const int &, DataArray &);
		void getCij(pEntity edge, const int &dom, double *Cij);
		bool getDij(pEntity , const int &, DataArray &);
		double getVolume(pEntity , const int &);
		double getWeightedVolume(pEntity);
		void getVolume(pEntity node, int dom, double &vol);
		int getFlag(pEntity);
		void setCij(pEntity , const int &, DataArray );
		void setCij_norm(pEntity , const int &, double );
		double getCij_norm(pEntity , const int &);
		void setDij(pEntity , const int &, const int &, const DataArray &);
		void setVolume(pEntity , const int &, const  double &);
		void setWeightedVolume(pEntity, const  double &);
		void setFlag(pEntity, int);
		void setEdgeLength(pEntity , double );
		double getEdgeLength(pEntity );
		void setEdgeVector(pEntity, dblarray);
		void getEdgeVector(pEntity, dblarray&);
		void setEdgeVec_Unitary(pEntity, dblarray);
		void getEdgeVec_Unitary(pEntity, dblarray&);

		// check if an edge belongs to domain
		bool nodeBelongToDomain(pEntity , const int &);
		bool edgeBelongToDomain(pEntity , const int &);
		int getDomainFlag(pEntity node);

		// check if an edge belongs to boundary domain (any domain)
		void set_belongsToBoundary(pEdge,bool);
		bool belongsToBoundary(pEdge);
		void setTotalReservoirVolume(double );
		double getReservoirVolume() const;

		static int getNumRemoteCopies(pEntity node){
//			GeomData pGCData;
//			Coefficients* pCoeffnt = pGCData.getAttachedData_pointer<Coefficients>(node);
//			return pCoeffnt->numRemoteCopies;
			return 0;
		}

		void setNumRemoteCopies(pEntity node, int nrc);

		// set the number of copies of an entity associated to a domain
		// An edge can have several remote copies but associated to a unique domain
		// this number can be small.
		void setNumRC(pEntity, int, int);

		static int getNumRC(pMesh theMesh, pEntity ent){
			return 0;//M_numRemoteCopies(theMesh,ent);
		}

		void setMeshDim(int mdim) { dim = mdim; }
		int getMeshDim() const { return dim; }

		void setNumGEdges(int nge) { numGEdges = nge; };
		int getNumGEdges() const { return numGEdges; }

		void setSmallestEdgeLength(double sel) { smallestEdge = sel; }
		double getSmallestEdgeLength() const { return smallestEdge; }

		double getReservoirHeight() const { return reservoirHeight; }
		void setReservoirHeight(double h) { reservoirHeight = h; }

		// used to verify coefficients summation
		void setSumIJ(pEntity,DataArray);

		// used to verify coefficients summation
		void getSumIJ(pEntity,DataArray&);

		void getCij(int dom, int row, double* cij);
		void setCij(int dom, int row, double* cij);
		void getDij(int dom, int row, double* dij);
		void setDij(int dom, int row, double* dij);
		void calculateNumEdges(pMesh);
		void calculateNumFaces(pMesh);
		void calculateNumFacesTmp(pMesh);
		void calculateNumBDRYEdges(pMesh);
		void calculateNumNodes(pMesh);
		void calculateNumBdryNodes(pMesh );
		void calculateEdgeProperties(pMesh theMesh);
		void allocatePointers(int);
		void deallocatePointers();

		int getNumEdgesPerDomain(int i) const{
			return numDomEdges[i];
		}

		int getNumBDRYEdgesPerDomain(int i) const{
			return numDomBDRYEdges[i];
		}

		int getNumNodesPerDomain(int i) const{
			return numNodesPerDomain[i];
		}

		void cleanData(pMesh);

		void mappingNodesIds(pMesh theMesh);

		// creates a mapping for saturation and pressure solution based in a new mesh (the adapted mesh)
		void mappingNodesIds_Tmp(pMesh theMesh);

		void getFace(int dom,int row,int& idx0,int& idx1,int& idx2,int& idx0_global,int& idx1_global,int& idx2_global){
			idx0 = faces[dom].getValue(row,0);
			idx1 = faces[dom].getValue(row,1);
			idx2 = faces[dom].getValue(row,2);
			idx0_global = faces[dom].getValue(row,3);
			idx1_global = faces[dom].getValue(row,4);
			idx2_global = faces[dom].getValue(row,5);
		}

		void getEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global){
			idx_0 = edges[dom].getValue(row,0);
			idx_1 = edges[dom].getValue(row,1);
			idx0_global = edges[dom].getValue(row,2);
			idx1_global = edges[dom].getValue(row,3);
		}

		void getEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global, int &flag1, int &flag2){
			idx_0 = edges[dom].getValue(row,0);
			idx_1 = edges[dom].getValue(row,1);
			idx0_global = edges[dom].getValue(row,2);
			idx1_global = edges[dom].getValue(row,3);
			flag1 = edges[dom].getValue(row,4);
			flag2 = edges[dom].getValue(row,5);
		}

		void getNodeIdx_Global(int dom, int i, int &idx){
			idx = nodes[dom].getValue(i);
		}

		void setVolume(int idx, double v){
			volume_global[0].setValue(idx,v);
		}

		void getVolume(int idx,double &v){
			v = volume_global[0].getValue(idx);
		}

		void getVolume(int dom, int idx, double& vol){
			vol = volume[dom].getValue(idx);
		}

		void getVolume(int dom, int idx_0, int idx_1, double& volumeI, double& volumeJ){
			volumeI = volume[dom].getValue(idx_0);
			volumeJ = volume[dom].getValue(idx_1);
		}

		void getID(int dom, int idx_0, int idx_1, int& id0, int &id1){
			id0 = ID[dom].getValue(idx_0);
			id1 = ID[dom].getValue(idx_1);
		}

		void getID(int dom, int idx, int& id){
			id = ID[dom].getValue(idx);
		}

		void getBdryEdge(int dom, int row, int &idx_0, int &idx_1){
			idx_0 = edges_bdry[dom].getValue(row,0);
			idx_1 = edges_bdry[dom].getValue(row,1);
		}

		void getBdryEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global){
			idx_0 = edges_bdry[dom].getValue(row,2);
			idx_1 = edges_bdry[dom].getValue(row,3);
			idx0_global = edges_bdry[dom].getValue(row,4);
			idx1_global = edges_bdry[dom].getValue(row,5);
		}

		void getBdryVolume(int dom, int idx, double& vol){
			vol = volume_bdry[dom].getValue(idx);
		}

		void getBdryVolume(int dom, int idx0, int idx1, int idx2, double *vol){
			vol[0] = volume_bdry[dom].getValue(idx0);
			vol[1] = volume_bdry[dom].getValue(idx1);
			vol[2] = volume_bdry[dom].getValue(idx2);
		}

		void getBdryVolume(int dom, int idx_0, int idx_1, double& volumeI, double& volumeJ){
			volumeI = volume_bdry[dom].getValue(idx_0);
			volumeJ = volume_bdry[dom].getValue(idx_1);
		}

		void getBdryID(int dom, int idx_0, int idx_1, int& id0, int &id1){
			id0 = ID_bdry[dom].getValue(idx_0);
			id1 = ID_bdry[dom].getValue(idx_1);
		}

		void getBdryID(int dom, int idx, int& id){
			id = ID_bdry[dom].getValue(idx);
		}

		int getDomFlag(int i) const{
			return domainList[i];
		}

		void getLength(int dom, int idx, double &length) const{
			length = edge_length[dom].getValue(idx);
		}

		void getVersor(int dom, int idx, double *v) const{
			v[0] = edge_versor[dom].getValue(idx,0);
			v[1] = edge_versor[dom].getValue(idx,1);
			v[2] = edge_versor[dom].getValue(idx,2);
		}

		void setCij_norm(int dom, int idx, double val){
			Cij_norm[dom].setValue(idx,val);
		}

		void getCij_norm(int dom, int idx, double &val){
			val = Cij_norm[dom].getValue(idx);
		}

		void getEBE(int idx, double* versor){
			versor[0] = EBE_1[0].getValue(idx,0);
			versor[1] = EBE_1[0].getValue(idx,1);
			versor[2] = EBE_1[0].getValue(idx,2);
		}

		void getEBE(int idx, int &idx0_global, int &idx1_global, int &flag1, int &flag2){
			idx0_global = EBE_2[0].getValue(idx,0);
			idx1_global = EBE_2[0].getValue(idx,1);
			flag1 = EBE_2[0].getValue(idx,2);
			flag2 = EBE_2[0].getValue(idx,3);
		}

		int getNumEBE() const{
			return numExtBdryEdges;
		}

		void setTotalNumberOfEdges(int n){
			numGEdges = n;
		}

		void getTotalNumberOfEdges(int &n) const{
			n = numGEdges;
		}

		void setMeshNodes(int n){
			numNodes = n;
		}

		void getMeshNodes(int &n) const{
			n = numNodes;
		}

		void setNumDomains(int n){
			_ndom = n;
		}

		void setDomainList(const int* domlist);

		int getNumDomains() const{
			return _ndom;
		}

		const int* getDomainList() const{
			return domainList;
		}

	private:
		int _ndom;
		int* domainList;
		double reservoirHeight;
		double reservoirVolume;
		int dim;
		double smallestEdge;
		int numGEdges;					// number or global edges
		int* numDomEdges;				// number of edges per domain
		int* numDomFaces;				// number of face per domain
		int* numDomFaces_tmp;			// number of face per domain adapted mesh
		int* numDomBDRYEdges;			// number of edges per domain
		int* numDomBDRYFaces;			// number of edges per domain
		int* numNodesPerDomain;			// number of nodes per domain
		int* numBdryNodesPerDomain;		// number of nodes per domain
		int numExtBdryEdges;
		int numNodes;

		Matrix<int> *ID;				// node ID per domain
		Matrix<int> *edges;				// edges id0-id1 per domain, where id0 and id1 are array indices and not node IDs.
		Matrix<int> *faces;				// stores local (per domain) and global indices for face node's IDs
		Matrix<int> *faces_tmp;			// (for interpolation):stores local (per domain) and global indices for face node's IDs
		Matrix<int> *nodes;				// same as edges
		Matrix<double>* volume;			//node volumes per domain
		Matrix<double>* volume_global;

		Matrix<int> *ID_bdry;			// node ID per domain
		Matrix<int> *edges_bdry;		// edges id0-id1 per domain, where id0 and id1 are array indices and not node IDs.
		Matrix<int> *faces_bdry;
		Matrix<double>* volume_bdry;	//node volumes per domain

		Matrix<double>* EBE_1;
		Matrix<int>* EBE_2;

		Matrix<double>* Cij;			// Cij vector
		Matrix<double>* Dij;

		Matrix<double>* edge_versor;
		Matrix<double>* edge_length;
		Matrix<double>* Cij_norm;
	};
}

#endif /*GEOMETRICCOEFFICIENTSDATA_H_*/

