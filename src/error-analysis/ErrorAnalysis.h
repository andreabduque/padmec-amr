/*
 * ErrosAnalysis.h
 *
 *  Created on: 20/02/2012
 *      Author: rogsoares
 */

#ifndef ERRORANALYSIS_H_
#define ERRORANALYSIS_H_

/*
 * Get estimated error from mesh elements.
 */

#include "SimulatorParameters.h"
#include "GeomData.h"
#include "Matrix.h"
#include "CPU_Profiling.h"

using namespace PRS;

class ErrorAnalysis;

bool analyzeField(FIELD, ErrorAnalysis*, SimulatorParameters*, GeomData*, void(*)(FIELD,int,int,int,const double*&));
bool analyzeField2(FIELD, ErrorAnalysis*, SimulatorParameters*, GeomData*, void(*)(FIELD,int,int,int,const double*&));
bool calculate_ErrorAnalysis(ErrorAnalysis*, SimulatorParameters*, GeomData*, void(*)(FIELD,int,int,int,const double*&), std::list<int>&, std::map<int,double>&);
bool calculate_ErrorAnalysis(ErrorAnalysis*, SimulatorParameters*, GeomData*, void(*)(FIELD,int,int,int,const double*&));

class ErrorAnalysis{

public:

	ErrorAnalysis(){
		init = false;
	}
	~ErrorAnalysis(){}

	//calculate functions
	void calculate_ElementsError(SimulatorParameters *pSimPar, GeomData*, void(*)(FIELD,int,int,int,const double*&), FIELD);
	void calculate_SmoothedGradientNorm(SimulatorParameters *pSimPar, GeomData*, void(*)(FIELD,int,int,int, const double*&), FIELD);
	void calculate_SmoothedGradientNorm_Singularity(SimulatorParameters *pSimPar, GeomData*, void(*)(FIELD,int,int,int,const double*&), FIELD);
	void calculate_DegreeOfRefinement(GeomData*, RefinementStrategies, bool, FIELD);

	double getElementError(int i) const{
		return pElemError[i];
	}

	void calculate_GlobalError(GeomData* pGCData);
	void calculate_GlobalError_Singularity(GeomData* pGCData);
	void calculate_AvgError(int,double,bool);
	void initialize(GeomData* pGCData, SimulatorParameters *pSimPar);
	void resetVariables(int);

	// returns a list of connectivities of elements to be removed/modified for mesh adaptation and its size as well.
	void getElementsForAdaptation(double param1, double param2, GeomData* pGCData, std::list<int> &elemList);
	
	// returns a map< vertex ID, weighted element heights> for all mesh nodes
	void getNodesForAdaptation(GeomData*, std::map<int,double>& nodesMap);

	void deletePointers();

	void countElements(int nelem, bool seekforSingular);

	double calculate_ErrorSum(GeomData* pGCData, bool excludingSingularities);

	// set/get functions
	void setNumElements_Singularity(int n) { 
		numElements_Singularity = n;
	}
	
	int getNumElements_Singularity() const{
		return numElements_Singularity;
	}

	void setSmoothedGradNorm(double n){
		SGN = n;
	}

	double getSmoothedGradNorm() const {
		return SGN;
	}

	void setSmoothedGradNorm_Singularity(double n) {
		SGN_sing = n;
	}

	double getSmoothedGradNorm_Singularity() const {
		return SGN_sing;
	}

	void setNumElements(int n) {
		numElements = n;
	}

	int getNumElements() const{
		return numElements;
	}

	int isSingular(int i) const{
		return isElementSingular[i];
	}
	
	void setElementAsSingular(int i){
		isElementSingular[i] = true;
	}

	// flag element as no singular, which means it not belongs to a region with high gradients
	void setAllElementsAsNotSingular(int);

	void setGlobalError(double ge){ 
		globalError = ge;
	}

	double getGlobalError() const{
		return globalError;
	}
	
	void setGlobalError_Singularity(double ge){
		globalError_Singularity = ge;
	}

	double getGlobalError_Singularity() const{
		return globalError_Singularity;
	}

	void setAverageError(double ae){
		averageError = ae;
	}

	double getAverageError() const{
		return averageError;
	}

	void setAverageError_Singularity(double ae){ 
		averageError_Singularity = ae;
	}

	double getAverageError_Singularity() const{
		return averageError_Singularity;
	}

	void setElementError(int i, double error){
		pElemError[i] = error;
	}

	void set_h_min(GeomData* pGCData, SimulatorParameters* pSimPar);

	void set_h_new(double h_new, int i, int j){
		p_h_new->setValue(i,j,h_new);
	}

	double get_h_new(int i, int j) const{
		return p_h_new->getValue(i,j);
	}

	double get_h_min() const{
		return h_min;
	}
	
	bool isToRemove(int k) const{
		return pElmToRemove[k];
	}

	bool adapt;
	std::list<pEntity> singularElemList;

	void calculate_heights(GeomData*, FIELD, bool);
	void calculate_h_new(GeomData*, double, FIELD);
	void identify_singular_regions(GeomData*, FIELD);
	void calculate_h_ratio(GeomData*);

	double getWE_node(int i) const{
		return pWH_node[i];
	}

private:

	int numElements;					// number or mesh's elements
	int numElements_Singularity;		// number or mesh's elements excluding those flagged as singular
	double SGN;							// Smoothed Gradient Norm
	double SGN_sing;					// Smoothed Gradient Norm excluding singular elements
	double globalError;					// global error for all mesh's element
	double globalError_Singularity;		// global error for all mesh's element excluding those flagged as singular
	double averageError;				// average global error for all mesh's element
	double averageError_Singularity;	// average global error for all mesh's element excluding those flagged as singular
	double h_min;						// minimum element height does not allow that elements smaller than h_min be removed from mesh
	double* pElemError;					// where element error are stored
	bool* isElementSingular;			// says if element is singular or not
	Matrix<double>* p_h_new;			// save h_new for all elements for each field
	bool* pElmToRemove;					// says if an element is to be removed or not
	double* pWH_node;					// pointer to array to store Weighted height per mesh node. (for visualization purpose)

	bool init;
};


#endif /* ERROSAnalysis_H_ */
