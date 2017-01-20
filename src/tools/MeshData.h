#ifndef MESHDATA_H_
#define MESHDATA_H_

#include "SimulatorParameters.h"

namespace PRS
{

	/**
	 * MeshData provides methods to get informations about free and prescribed
	 * nodes. It's very important to use it to assembly system of equations
	 * specially in parallel simulations. User does not need to worry about
	 * where a specific node equation should be assembled in global matrix.
	 * MeshData makes use of Petsc AO (Application Ordering) to make parallel
	 * matrix assembly with the minimum remote communication.
	 **/

//	typedef map<int,double> Map;
	typedef map<int,double>::const_iterator MIter;
	typedef set<int>::iterator SIter;

	// _Unify _Vectors On _Mesh _Nodes Struct
	struct UVMN_Struct{
		bool onlyRNOB;	// onlyRemoteNodesOnBoundaries
		int dim;		// dimension
		int dom;		// domain
		char* tag;		//tag for domain
		int coord_xyz;	// vector position: x, y or z;
		void (*pFunc_getVector)(int,int,double*);
		void (*pFunc_setVector)(int,int,const double*);
	};

	class MeshData
	{
	public:

		MeshData();
		MeshData(SimulatorParameters *,pMesh);
		~MeshData();

		/**
		 * Initialize operations:
		 * 		reorderVerticesIds
		 * 		settingFreeAndPrescribedNodes
		 * 		createVectorsForRHS
		 */
		void initialize(pMesh theMesh, GeomData *);

		void deallocateData();
		void destroyPointers();

		// Seek on mesh entities nodes with prescribed values (dirichlet condition)
		int getNodesWithKnownValues(pMesh theMesh);

		// If node is set as a dirichlet it returns true, otherwise returns false
		bool getDirichletValue(int ID, double *val);

		inline int getNum_GF_Nodes() const { return numGF; }	/// global free nodes
		inline int getNum_GP_Nodes() const { return numGP; }	/// global dirichlet nodes
		inline int getNum_LF_Nodes() const { return numLF; }	/// local free nodes
		inline int getNum_GNodes() const { return numGN; }		/// global nodes

		// Get ordering number to be used in matrix assembling
		int get_AppToPETSc_Ordering(int n) const;

		// Get previous ordering number to be used in matrix assembling
		int get_PETScToApp_Ordering(int n) const;

		// Generates node ID ordering
		void reorderVerticesIds(pMesh theMesh, int(*)(pEntity));
		
		// returns a map iterator to access all dirichlet nodes and their prescribed values
		inline MIter dirichletBegin() const { return dirichlet.begin(); }
		inline MIter dirichletEnd() const { return dirichlet.end(); }
		void findDirichletNodes();

		// return mapped position for all nodes
		inline int FPArray(int pos) const { return FP_Array[pos]; }
		inline int get_idxFreecols(int pos) const { return idxFreecols[pos]; }
		inline int get_idxn(int pos) const { return idxn[pos]; }

		// Return pointer to assemble matrices
		inline const PetscInt* get_F_cols_ptr() { return F_cols; }
		inline const PetscInt* get_F_rows_ptr() { return F_rows; }
		inline const PetscInt* get_idxFreecols_ptr() { return idxFreecols; }
		inline const PetscInt* get_pos_ptr() { return pos; }
		inline const PetscInt* get_idxn_ptr() { return idxn; }
		inline PetscInt get_F_nrows() const { return F_nrows; }
		inline void set_F_nrows(int val) { F_nrows = val; }

		// Inform which COLUMNS from global matrix A must be copied to assembly LHS matrix (free nodes). 
		// this vector is the same for all processors.
		PetscInt *idxFreecols;
		
		// inform which columns from A must be copied to assembly RHS vector (prescribed nodes) this vector is the same for all processors.
		PetscInt *idxn;
		PetscInt *idxm;
		PetscInt *pos;
		PetscInt *F_rows;
		PetscInt *F_cols;

		void createVectorsForRHS(pMesh theMesh, int);
		int createVectorsForMatrixF(Mat &);

		int rowsToImport(pMesh theMesh, int &, int *&);
		void settingFreeAndPrescribedNodes(pMesh theMesh);

		// Rogerio, try to think in something more clear in the future!
		void getRemoteIDs(int &nLIDs, int** IDs_ptr);

		/*
		 * unifyScalarsOnMeshNodes function works like MPI_Allreduce with operator
		 * MPI_SUM for each node on partition boundary. For each node with remote
		 * copy, it sums all node's values and broadcast the summation to all nodes'
		 * copies.
		 *
		 * Ex.: node ID = 1736 has 4 remote copies and each copy has the following
		 * 		values:
		 * 		p1 = 23.34;
		 * 		p2 = 873.88;
		 * 		p3 = -12.87;
		 * 		p4 = 73.09;
		 *
		 * 		sum = p1 + p2 + p3 + p4 = 957.44
		 * 		after unifyScalarsOnMeshNodes has been called:
		 *
		 * 		p1 = 957.44;
		 * 		p2 = 957.44;
		 * 		p3 = 957.44;
		 * 		p4 = 957.44;
		 *
		 * User must provide two pointer functions: one to get and another to set
		 * the scalars values
		 */
		int unifyScalarsOnMeshNodes(double(*pFunc_getScalar)(pEntity),
								    void (*pFunc_setScalar)(pEntity,double),
								    GeomData*,void *ptr=0);
//		int unifyVectorsOnMeshNodes(void (*pFunc_getVector)(pEntity,int,dblarray&),
//				                    void (*pFunc_setVector)(pEntity,int,dblarray),
//				                    GeomData* pgc, int, bool rn=false);

		int unifyVectorsOnMeshNodes(void (*pFunc_getVector)(int,int,double*),
								    void (*pFunc_setVector)(int,int,const double*),
								    GeomData* pgc, int, bool rn=false);


		/***********************************************************************************************/
		// MEBFV:
		/***********************************************************************************************/
		int nrows() const;					// number of mesh nodes
		int ncols() const;					// number of mesh nodes
		int numFreeNodes() const; 			// number of mesh free nodes
		void createGlobalNodeIDMapping();	// create all mapping needed to assembly matrices and vectors
		void mapNodeID();					// gives a sequential numbering for all mesh node IDs. From 0 to N-1,
											// where N is the number of mesh nodes

		// for global matrix assembly
		void getArrays_free(int& freeRows_size, int** freeRows_array);
		void getGlobalIndices(int ith_elem, int ith_edge, int* idxm, int* idxn);
		void setGlobalIndices();

		// for matrix solver assembly: extracts rows and columns from global matrix A related to free vertices
		IS IS_getFreeRows() const;
		IS IS_getFreeCols() const;

		// for RHS assembly only
		void getArrays_dirichlet(int& dirichletCols_size, int** dirichletCols_array, std::set<pEntity>& dirichletnodes_set);
		void initDirichletVectors(std::set<pEntity>& dirichletnodes_set);
		const int* getDirichlet_idx() const;		// return pointer of indices for prescribed mesh nodes
		const double* getDirichlet_data() const;	// return pointer of prescribed values assigned to mesh nodes
		IS IS_getDirichletRows() const;					// from global matrix
		IS IS_getDirichletCols() const;					// from global matrix


	private:
		PetscErrorCode ierr;

		AO ao;			// application mapping
		int numLF;		// number of local prescribed nodes
		int numGF;		// number of global free nodes
		int numGN;
		int numGP;
		int *FP_Array;
		int F_nrows;

		// map local prescribed nodes: ID<int> -> prescribed_value<double>
		map<int,double> dirichlet;
		// container with all domains flags
		set<int> setOfDomains;

		SimulatorParameters *pSimPar;
		pMesh theMesh;

		int FreePrescribedNodes(pMesh theMesh);
		void mappingUnknowns();
		int M_numGVertices();

		// Rogerio, try to think in something more clear in the future!
		// array of vertices IDs used to get remote values from a distributed
		// column matrix. Bring to processor those values that are located in
		// other ones.
		int *localIDs;
		int numLocalIDs;

		/*
		 * Set node's Id a crescent numbering.
		 * IF node IDs are: 12,45,79 -> 0,1,2
		 */
		map<int,int> localIDNumbering;

		// unifyScalarsOnMeshNodes's stuffs
		int *rowToImport;
		Mat joinNodes;
		Mat updateValues;
		map<int,double> mapPB_nodes;	// map partition boundary nodes

		/*
		 * structsCreation allows that data structs like vectors and matrices are
		 * created just once. This saves time because dynamic memory allocation
		 * has a high computational cost.
		 */
		bool structsCreation;

		// unifyVectorsOnMeshNodes's stuffs
		UVMN_Struct *pMS;	// _Unify _Vectors On _Mesh _Nodes Struct

		/***********************************************************************************************/
		// MEBFV:
		/***********************************************************************************************/
		std::map<int,int> nodeID_map;

		int* dirichlet_idx;		// a vector to assembly rhs direchlet part: rhs = sst - dirichlet, sst=source/sink term
								// it's size is the number of free nodes and store: 0, 1, 2, ..., n-1, n: number fo free vertices
		double* dirichlet_data;
		int msize;		// number of mesh nodes
		int nfreen;		// number of free (not prescribed) mesh nodes
		int nnfreen;	// number of not free (prescribed) mesh nodes

		// stores indices for sub-matrices related to element
		//                                0      1      2       3      4
		// 2-D (triangles): MatIndices: idxm_0 idxm_1 idxn_0 idxn_1 idxn_2
		Matrix<int> MatIndices;
		IS IS_freeRows;
		IS IS_freeCols;
		IS IS_dirichletRows;
		IS IS_dirichletCols;
	};
}
#endif /*MESHDATA_H_*/
