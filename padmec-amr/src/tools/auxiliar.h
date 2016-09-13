#ifndef AUXILIAR_H_
#define AUXILIAR_H_

#include "includes.h"

struct PADMEC_mesh{
	int numVertices;
	int numElements;
	int numbeges;
	double* coords;
	int *ID;
	int* bedges;
	int* elements;
	double *field1;
	double *field2;
};

void open_file(ofstream& fid, string filename, int line, const char* sourcefile);

void throw_exception(bool,string, int line, const char* sourcefile);

// allocating memory
// -----------------------------------------------
void alloc_BOOL_vector(int LINE, const char* FILE, bool* &p, int size);
void dealloc_BOOL_vector(bool* &p);

void alloc_INT_vector(int LINE, const char* FILE, int* &p, int size);
void dealloc_INT_vector(int* &p);
void alloc_DOUBLE_vector(int LINE, const char* FILE, double* &p, int size);
void dealloc_DOUBLE_vector(double* &p);

enum FIELD {PRESSURE, SATURATION};

const double pi = 3.14159265359;

double F_area(const double *ptn1, const double *ptn2, const double *ptn3);

double norm(const double *n);

void setBarrier();

double E_length(pEdge edge);

void getEdgeVector(pEdge edge, std::vector<double> &edgVec);

double strToDouble(string &str);

int strToInteger(string &str);

const char* getSubString(string &str);

int getVertexFlag(pVertex vertex);

void makeVector(const double *A, const double *B, double *v);

//double* getIdentityMatrix(const int &dim);

double getSmallestEdgeLength(pMesh theMesh);

int getVertexFlag(pVertex node);
int getEdgeFlag(pEdge edge);
int getFaceFlag(pFace face);
int getTetraFlag(pRegion tetra);
int getEntityFlag(int i, pEntity ent);

void replaceAllOccurencesOnString(string &, string::size_type, string, string);
void getIJnodes(pEdge, std::vector<pEntity>&);

//void F_getEdges(pMesh, pEntity, std::vector<pEntity> &);

const double qsi = 1e-10; // qsi is used to avoid division by zero


/*
 * Get edge vertices. M_GetVertices returns the vertices for every entity (edge,
 * triangle, tetra, etc. E_vertices avoid the 'if' statement used by the former
 * function allowing some CPU time reduction.
 */
void E_vertices(pEntity edge, std::vector<pEntity> & vertex);

/*
 * Get coordenates (x,y,z) for edge's vertices I and J.
 */
void E_getVerticesCoord(pEntity edge, double *I, double *J);

typedef Trellis_Util::mPoint VPoint;

/*
 * Get vertex point from an edge
 */
VPoint E_getVertexPoint(pEntity,int);

void E_IJvector(double*, double*, double*);

void printSimulationHeader();

typedef std::vector<double> dblarray;

// let's use C++ <numeric> and <algorithm> to help us doing cool things like
// calculating dot product and euclidian norm
double inner_product(const dblarray &a1, const dblarray &a2);

double inner_product(const double*, const double*, int);

void unitary_vector(const dblarray &a1, dblarray &a2);

void unitary_vector(dblarray &a1);

double norm2(const dblarray &a1);

void checklinepassing(int, const char*, std::string);

enum LOG_FILES {OPENLG, UPDATELG, CLOSELG};

void LogFiles(double timeStep, double assemblyT, double solverT, double gradT, int KSPiter, double hyperbolicCPU,LOG_FILES LG, string path, bool restart, int last_step, double cumulativeSTime_Restart, double CPUTime_Restart);

void failOpeningFile(string, int, const char *);

/*
 * converts seconds to hours minutes seconds
 */
void convertSecToTime(double t, double *h, double *m, double *s);

void STOP();


int printMatrixToFile(Mat&,const char*);
int printVectorToFile(Vec&,const char*);


/*
 * Define types for arrays of pointer functions (Scalars)
 */
typedef double(*GetPFunction)(pEntity);
typedef void(*SetPFunction)(pEntity,double);

typedef GetPFunction* GetFunctionArray;
typedef SetPFunction* SetFunctionArray;

/*
 * Define types for arrays of pointer functions (scalars/gradients)
 */
typedef double (*GetPFuncScalar)(pEntity);
typedef void (*SetPFuncScalar)(pEntity,double);
typedef void (*GetPFuncGrad)(int, int, double*);
typedef void (*SetPFuncGrad)(int, int, double*);

typedef GetPFuncScalar* GetFuncScalarArray;
typedef SetPFuncScalar* SetFuncScalarArray;
typedef GetPFuncGrad* GetFuncGradArray;
typedef SetPFuncGrad* SetFuncGradArray;

//typedef void(*FuncPointer_GetGradient)(FIELD,int,int,int,double*);

/*! \brief: Get leaves of an edge with refinementDepth equal 1 (ONLY two children)
 * \param pEdge edge: edge where to get children from
 * \param pEdge* edgeChildren: pointer array with two children
 */
void getEdgesChildren(pEdge, mEdge**);

/*! \brief: Get faces around a face. In general, face has three neighbors, but it has boundary edges, it may have two or only one neighbor
 * \param face face to find its neighbors
 * \param faces face's neighbors
 * \param numFaces how many neighbor were found
 */
void getFacesAroundFace(pFace face, pFace *faces, int &numFaces);


double dot(const double *u, const double *v, int dim);
double norm_L2(const double *p, int dim);

int EN_getFlag(pEntity);

void getTriVerticesIDs(pFace face, int* IDs);

void printFaceIDs(pEntity);

//void neighboursFace(pMesh theMesh, pEntity face, std::vector <pFace> &fvector);
//void neighboursFaceofEdge(pMesh theMesh, pEntity edge, std::vector <pFace> &fvector);
void getNeighboursFace(pEntity face, std::vector<pEntity> &neighbousFace, int &);
void getChildren (pMesh theMesh, pEntity parent, std::vector<pFace> &fchildren);
void getEdgesNodeIDs(pEdge edge, int &ID1, int &ID2);
bool isNotEdgeOnBdry(pEdge edge);

/*! \brief: Makes a copy from mesh m1 to mesh m2
 * \param pMesh m1 from
 * \param pMesh m2 tp
 */
void makeMeshCopy(pMesh m1, pMesh m2, 
				  void(*pSetPressure)(pEntity,double), double(*pGetPressure)(pEntity),
				  void(*pSetSaturation)(pEntity,double), double(*pGetSaturation)(pEntity));

void makeMeshCopy2(pMesh, PADMEC_mesh*, void(*pGetPressure)(int,double&), void(*pGetSaturation)(int,double&));
void makeMeshCopy2(PADMEC_mesh*,pMesh , void(*pSetPressure)(int,double), void(*pSetSaturation)(int,double));

/*! \brief: delete elements mesh
 * \param pMesh m
 */
void deleteMesh(pMesh m);
void deleteMesh(PADMEC_mesh* pm);
void checkNumEdges(pMesh, int, char*);
void checkMesh(pMesh m);
void calculateNumFacesAroundVertices(pMesh m, std::map<int,int> &facesAroundVertices);
void readmesh(pMesh m, const char* filename);

// reads a mesh file from gmsh (file version 1.0) with parametric coordinates
struct parametric_coords{
	int npc;	// number of parametric coordinates
	double pcoord1;
	double pcoord2;
	double pcoord3;
};
void readmesh_parametric(pMesh m, const char* filename, std::list<parametric_coords> &parametric_list);

void getdomains(pMesh m, int &n, int* domlist);

#endif /*AUXILIAR_H_*/
