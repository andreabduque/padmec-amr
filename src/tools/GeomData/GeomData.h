#ifndef GEOMETRICCOEFFICIENTSDATA_H_
#define GEOMETRICCOEFFICIENTSDATA_H_


#include "AttachData.h"
#include "Matrix.h"

namespace PRS{

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
	class GeomData: public AttachData{
	public:

		GeomData();
		~GeomData();

		void initilize(pMesh);
		void initilize(pMesh, const std::set<int>&, int FVM);
		void dataTransfer(pMesh);
		void transferCijData(pMesh);
		void transferDijData(pMesh);
		void transferVolData(pMesh);
		void transferMesh(pMesh);

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

		void getEdgeLength(pEntity, double&);
		void getEdgeVec_Unitary(pEntity, double*);

		// check if an edge belongs to domain
		bool nodeBelongToDomain(pEntity , const int &);
		bool edgeBelongToDomain(pEntity , const int &);
		bool faceBelongToDomain(pEntity , const int &);
		int getDomainFlag(pEntity node);

		// check if an edge belongs to boundary domain (any domain)
		void set_belongsToBoundary(pEdge,bool);
		bool belongsToBoundary(pEdge);
		void setTotalReservoirVolume(double );
		double getReservoirVolume() const;


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
		void getCij(int dom, int row, const double* &cij);
		void setCij(int dom, int row, double* cij);
		void getDij(int dom, int row, double* dij);
		void getDij(int dom, int row, const double* &dij);
		void setDij(int dom, int row, double* dij);
		void calculateNumEdges(pMesh);
		void calculateNumElements(pMesh);
		void calculateNumFacesTmp(pMesh);
		void calculateNumBDRYEdges(pMesh);
		void calculateNumBDRYFaces(pMesh);
		void calculateNumNodes(pMesh);
		void calculateNumBdryNodes(pMesh );

		// calculate: edges length, edges versor (domain and boundary)
		void calculateEdgeProperties(pMesh theMesh);
		void allocatePointers(int,int);
		void deallocatePointers(int);
		int getNumEdgesPerDomain(int i) const;
		int getNumElemPerDomain(int i) const;
		int getNumBDRYEdgesPerDomain(int i) const;
		int getNumBdryFacesPerDomain(int i) const;
		int getNumNodesPerDomain(int i) const;
		int getNumNodes() const;
		void cleanData(pMesh);
		// needed to create an indexed data structured and though get fast access (direct access, O(1)) to any data.
		void mappingNodesIds(pMesh theMesh);
		// creates a mapping for saturation and pressure solution based in a new mesh (the adapted mesh)
		void mappingNodesIds_Tmp(pMesh theMesh);

		void getEdge(int dom, int row, const int* &ptr);
		void getEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global);
		void getEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global, int &flag1, int &flag2);

		void getNodeIdx_Global(int dom, int i, int &idx);
		void setVolume(int idx, double v);
		void getVolume(int idx,double &v);
		void getVolume(int dom, int idx, double& vol);
		void getVolume(int dom, int idx_0, int idx_1, double& volumeI, double& volumeJ);

		// return the i_th global node ID
		int getNodeID(int i_th) const;

		void getID(int dom, int idx_0, int idx_1, int& id0, int &id1);
		void getID(int dom, int idx, int& id);
		void getBdryEdge(int dom, int row, int &idx_0, int &idx_1);
		void getBdryEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global);
		// return
		void getBdryFace(int dom, int row, int &idx_0, int &idx_1, int &idx_2, int &idx0_global, int &idx1_global, int &idx2_global);
		// get control volume for one vertex
		void getBdryVolume(int dom, int idx, double& vol);
		// get control volume for trhee vertices
		void getBdryVolume(int dom, int idx0, int idx1, int idx2, double *vol);
		// get control volume for three vertices
		void getBdryVolume(int dom, int idx0, int idx1, int idx2, double& volumeI, double& volumeJ, double& volumeK);
		// get control volume for two vertices
		void getBdryVolume(int dom, int idx_0, int idx_1, double& volumeI, double& volumeJ);
		// get ID for two vertices on boundary
		void getBdryID(int dom, int idx_0, int idx_1, int& id0, int &id1);
		// get ID for three vertices on boundary
		void getBdryID(int dom, int idx_0, int idx_1, int idx_2, int& id0, int &id1, int& id2);
		// get ID for vertex on boundary
		void getBdryID(int dom, int idx, int& id);
		// get flag number defined by user on .geo file for domain i (i=0,1,2,...,n-1), where n is the number of domains
		int getDomFlag(int i) const;
		// return edge length
		void getLength(int dom, int idx, double &length) const;
		// return versor built over edge pointing from node I to node J, where Node ID I is ALWAYS less than node ID J.
		void getVersor(int dom, int idx, double *v) const;
		void setCij_norm(int dom, int idx, double val);
		void getCij_norm(int dom, int idx, double &val);
		void getVersor_ExternalBdryElement(int idx, double* versor);

		// Used by Saturation Gradient
		void getExternalBdryEdges(int idx, int &idx0_global, int &idx1_global, int &flag1, int &flag2);
		void getExternalBdryFaces(int idx, int &idx0_global, int &idx1_global, int &idx2_global, int &flag1, int &flag2, int &flag3);


		int getNumExternalBdryEdges() const;
		int getNumExternalBdryFaces() const;

		void calculate_extFaceVersor(pMesh theMesh);

		void setTotalNumberOfEdges(int n);
		void getTotalNumberOfEdges(int &n) const;
		void setMeshNodes(int n);
		void getMeshNodes(int &n) const;
		void setNumDomains(int n);
		void setDomainList(const int* domlist);
		int getNumDomains() const;
		const int* getDomainList() const;
		int getNumElements() const;
		void getConnectivities(int row, int *connectivities);
		void getCoordinates(int row, double *coords);
		void getCoordinates(int row, const double* &coords);

		/*MEBFV: function for the Modified EBFV*/
		void initializeElementMatrix(int numElements);
		void setElementMatrices(int row, const double* Aij, const double* Ajk, const double* Aik);
		void getElementMatrices(int row, double* Aij, double* Ajk, double* Aik);
		void getVolume_MEBFV(int,double&);

		void setVersor(pEntity edge, double* versor);
		void getVersor(pEntity edge, double* versor) const;

		// idx is an array of size 2*(dim+1) which stores local and global indices for an element
		// if triangle: 3 locals and 3 globals: total = 2*(2+1) = 6 indices
		// if tetra: 4 locals and 4 globals:    total = 2*(3+1) = 8 indices
		void getElement(int dom, int row, int* idx);

		void getElement(int dom, int row, const int* &idx);

		// calculate number of elements sharing the same vertex for all mesh vertices
		void calculateNumElemSharingVertex();

		// return number of elements sharing a vertex
		int getNumFacesSharingVertex(int elem) const;


		// calculate characteristic dimesion lentgh (CDL) of an element
		// it simply consists to compute the arithmetic mean of all element's edges lengths
		void calculate_CDL();
		double getElem_CDL(int i) const;

		// after calculate_CDL has been called, calculate an average CDL per node
		void calculate_NodeAverage_CDL();

		// return element height: h_new/h_old
		double getElem_HR(int i) const{
			return elem_HR[i];
		}

		// element height ratio: h_new/h_old
		void setElem_HR(int i, double hr){
			elem_HR[i] = hr;
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
		int* numDomElem;				// number of elements per domain
		int* numDomFaces_tmp;			// number of face per domain adapted mesh
		int* numDomBDRYEdges;			// number of edges per domain
		int* numDomBDRYFaces;			// number of edges per domain
		int* numNodesPerDomain;			// number of nodes per domain
		int* numBdryNodesPerDomain;		// number of nodes per domain
		int* numElemSharingVertex;		// number of elements sharing a vertex for all mesh vertices
		int numExtBdryEdges;
		int numExtBdryFaces;
		int numNodes;					// total number of mesh nodes
		double* elem_CDL;				// characteristic dimension length of an element
		double* node_CDL;
		double* elem_HR;				// store element height ration: h_new/h_old

		int* pNodeID;

		int elemtype;					// elemtype = 3 (2-D triangle: 3 nodes), elemtype = 4 (3-D tetrahedron: 4 nodes)
		int numElem;					// number of mesh elements (2D/3D)
		Matrix<int> *pConnectivities;	// matrix for element connectivities
		Matrix<double>* pCoords;		// Mesh node's coordinates: x, y, z

		Matrix<int> *ID;				// node ID per domain
		Matrix<int> *edges;				// edges id0-id1 per domain, where id0 and id1 are array indices and not node IDs.
		Matrix<int> *elem;				// stores local (per domain) and global indices for face node's IDs
		Matrix<int> *faces_tmp;			// (for interpolation):stores local (per domain) and global indices for face node's IDs
		Matrix<int> *nodes;				// same as edges
		Matrix<double>* volume;			//node volumes per domain
		Matrix<double>* volume_global;

		Matrix<int> *ID_bdry;			// node ID per domain
		Matrix<int> *edges_bdry;		// edges id0-id1 per domain, where id0 and id1 are array indices and not node IDs.
		Matrix<int> *faces_bdry;
		Matrix<double>* volume_bdry;	//node volumes per domain

		Matrix<double>* versor_ExtBdryElem;	// versor for external elements (2-D: edges, 3-D: triangles)
		Matrix<int>* external_bdry_elem;	// external_bdry_elem: idx0_global, idx1_global, idx2_global, flag1, flag2, flag3

		Matrix<double>* Cij;			// Cij vector
		Matrix<double>* Dij;

		Matrix<double>* edge_versor;
		Matrix<double>* edge_length;
		Matrix<double>* Cij_norm;

		/*MEBFV*/
		// a matrix to store all edge matrices (Aij, Ajk, and Aik) associated to element edges
		// for triangle elements, we have:
		// geoElementMat[i,:] = {Aij, Ajk, Aik};
		// where: Aij = {Aij_00, Aij_01, Aij_02, Aij_10, Aij_11, Aij__12}
		//        Ajk = {Ajk_00, Ajk_01, Ajk_02, Ajk_10, Ajk_11, Ajk__12}
		//        Aik = {Aik_00, Aik_01, Aik_02, Aik_10, Aik_11, Aik__12}
		//
		// geoElementMat size: m=number of elements, n = 6x3 = 18
		Matrix<double> geoElementMat;
		double* volume_MEBFV;			// vector for all volume control volume
	};
}

#endif /*GEOMETRICCOEFFICIENTSDATA_H_*/

