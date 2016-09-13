#ifndef EBFV1PREPROCESSORS_H
#define EBFV1PREPROCESSORS_H

#include "EdgeInfo.h"
#include "GeomData.h"

using namespace PRS;

struct PP_Parameters{
	pMesh theMesh;
	set<int> setOfDomain;
	const char *meshFilename;
	const char* exportFilename;
	GeomData *pGCData;
	PetscBool flg;
	int dim;
};

void validete_coefficients(pMesh theMesh, std::set<int>& setOfDomain, GeomData* pGCData);

void computeDij(pMesh, pFace, GeomData *);

int unifyCijAmongProcessors(pMesh, const set<int>&, GeomData*);

void identifyBoundaryElements(pMesh theMesh, GeomData *pGCData, std::set<int> setOfDomains);

int unifyVolumesAmongProcessors(pMesh, const set<int>&, const string&, GeomData*);

void getFCenter(pFace face, double *);

void getFCenter(pVertex, pVertex, pVertex, double *);

double F_area(pFace);

void AllgatherDomains(std::set<int> &);

void generateFilename(string,const string&,char*,string);

void generateFilename(const string&,char*,PetscBool&);

void exportCoefficients(pMesh, const int&, set<int>&, const char*, GeomData *,PetscBool);

void saveMesh_gmsh1(pMesh, const char*);

void makeVector(const double *A, const double *B, double *v);

void cross(const double *, const double *, const double *, double *);

void cross(const double *, const double *, double* );

void markTetraBdryFaces(pMesh, pEntity, const int &);

double dot(const double *,const double *);

void DijVector(pFace , pVertex &, dblarray &);

void DijVector(pFace, pVertex &, double *);

void computeDij(pMesh, pFace, GeomData *);

int setCorrectNumberOfRemoteCopies(PP_Parameters* pPPP);

//double norm(const double *n);

// return the number of edges for a serial or partitioned mesh
int getNumGlobalEdges(pMesh);


// returns an edge based on its vetices IDs
mEdge* getEdge(pMesh, int, int);

// returns the sum of all local edges without remote copies
int getNumGlobalEdgesWithoutRemoteCopies(pMesh);

// mark a edge for a specific porpouses
void markEdge(mEdge*);


// check if an edge was marked by 'markEdge' function
bool isEdgeMarked(mEdge*);

int EBFV1_preprocessor(pMesh, void*);
int EBFV1_preprocessor_2D(pMesh, void*, int&);
int EBFV1_preprocessor_3D(pMesh, void*, int&);

void calculateEdgeLength(pMesh, GeomData*);
void calculateCijNorm(pMesh, GeomData*, std::set<int>);


void initializeCoefficients(pMesh theMesh, GeomData *pGCData);

/*! Date: 05/09/2012
 *! -------------------------------------------------------------------------------------------------------
 *! Check geometric coefficients: Summation of all vectors for each control surface must be zero (<10e-15)
 *! -------------------------------------------------------------------------------------------------------
 */
typedef std::map<int, std::vector<double> > mapSUM;
//!brief Main function for geometric coefficients validation. If validation fails, program terminates.
void validate_EBFV1(GeomData *pGCData, pMesh theMesh, std::set<int> &setOfDomain);

//!brief Initialize sum = 0
void initializeCoeffSummation(GeomData *pGCData, pMesh theMesh, mapSUM &mapNodeCVsum);

//! brief Take Cij contribution
void getBoundaryEdgescontribution(GeomData *pGCData, pMesh theMesh, mapSUM &mapNodeCVsum, int flag_domain);

//!brief Take Dij contribution if control volume is on boundary.
void getEdgesContribution(GeomData *pGCData, pMesh theMesh, mapSUM &mapNodeCVsum, int flag_domain);


//! brief Checking validation
void calculateCVSummation(GeomData *pGCData, pMesh theMesh, mapSUM &mapNodeCVsum, int flag_domain);
#endif
