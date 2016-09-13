#ifndef PHYSICALPROPERTIESDATA_H_
#define PHYSICALPROPERTIESDATA_H_

#include "MeshData.h"
#include "Matrix.h"
#include "GeomData.h"

extern Matrix<double> *pGrad_dom;

/*
 * Saturation gradient is ALWAYS continuous even for heterogeneous problems
 * SwGrad_dom is an auxiliary matrix (DO NOT REMOVE IT FROM THE CODE!!!!!!!!!!)
 * SwGrad is the matrix where saturation must be get from
 */
extern Matrix<double> *SwGrad_dom; // DO NOT REMOVE IT FROM THE CODE!!!!!!!!!!
extern Matrix<double> SwGrad;

extern Matrix<double> pressure;
extern Matrix<double> Sw;
extern Matrix<double> Sw_old;

// temporary matrices (nx1) to store interpolated values from old to new mesh
extern Matrix<double> Sw_tmp;
extern Matrix<double> p_tmp;

namespace PRS{

/**
 * PhysicPropData class provides methods to set/get physical properties such as pressure, saturation and their respective gradients
 * as well as rock-fluid properties (fractional flow and total mobilities). Fluid velocity is also managed by this class. All proper-
 * ties are associated to a mesh entity (node or edge) through the class AttachData that is in charge of store all data (physical,
 * numeric, geometric,...) to mesh entities.
 */
class PhysicPropData{
public:

	PhysicPropData();
	~PhysicPropData();

	void initialize(MeshData *, SimulatorParameters *, pMesh, bool, GeomData* pGCData);
	void allocateTemporaryData(int);
	void allocateData(SimulatorParameters*, GeomData*, int);
	void deallocateData(SimulatorParameters*);

	// transfer data from Sw_tmp and p_tmp to Sw and pressure matrices
	void transferTmpData();

	// update Sw and pressure matrices (nx1) to receive interpolated values from old mesh
	void update_Sw_p(int);

	/*
	 *  set/get gradient pressure in new format. As pGrad_matrix is declared outside of class definition, it will not give errors in a static
	 *  member function. If pGrad_matrix was declared as a member of class, a class object will be required do be declared inside function, which
	 *  will call the constructor everytime set/get gradient was called.
	 *
	 *  dom - domain flag. As computation is performed by domain, it follows crescent order (0,1,2,3,...). Then, domain's flag is not necessary
	 *  to know, only the order (if it is the first, second, ....)
	 *  row - local node numbering. Do not use global ID or it will not work!!!!
	 *  grad - pointer to an array of three positions. Ex.: double p[3];
	 */
	static void set_pw_Grad(int dom, int row, const double* grad){
		pGrad_dom[dom].setValue(row,0,grad[0]);
		pGrad_dom[dom].setValue(row,1,grad[1]);
		pGrad_dom[dom].setValue(row,2,grad[2]);
	}

	static void get_pw_Grad_const(int dom, int row, const double* &grad){
		grad = pGrad_dom[dom].getrowconst(row);
	}

	void get_pw_Grad(int dom, int row, double* &grad){
		grad = pGrad_dom[dom].getrow(row);
	}

	static void setPressure(int idx, double p){
		pressure.setValue(idx,p);
	}

	static void getPressure(int idx, double& p){
		p = pressure.getValue(idx);
	}

	void get_Sw_Grad(int dom, int row, double* &grad){
			grad = SwGrad_dom[dom].getrow(row);
	}

	void get_Sw_Grad(int idx_global, double* &grad){
		grad = SwGrad.getrow(idx_global);
	}

	static void setPressure_NM(int idx, double v){
		p_tmp.setValue(idx,v);
	}

	static void setSaturation_NM(int idx, double v){
		Sw_tmp.setValue(idx,v);
	}

	// that's for error analysis
	static void getGradient(FIELD field, int dom, int idx, int idx_global, const double* &grad){
		switch (field){
		case PRESSURE:
			grad = pGrad_dom[dom].getrowconst(idx);
			break;
		case SATURATION:
			grad = SwGrad.getrowconst(idx_global);
			break;
		default:
			throw Exception(__LINE__,__FILE__,"Unknown field.");
		}
	}

	static void setSaturation(int idx, double v){
		Sw.setValue(idx,v);
	}

	static void getSaturation(int idx, double& sw){
		sw = Sw.getValue(idx);
	}

	void setSw_old(int idx, double v){
		Sw_old.setValue(idx,v);
	}

	static void getSw_old(int idx, double& sw){
		sw = Sw_old.getValue(idx);
	}

	void setInitialSaturation(pMesh, SimulatorParameters*);
	void setSw_max(pEntity,double);
	void setSw_min(pEntity,double);
	void setS_Limit(pEntity,double);

	void setVelocity_new(int dom, int row, double* vel){
		velocity[dom].setValue(row,0,vel[0]);
		velocity[dom].setValue(row,1,vel[1]);
		velocity[dom].setValue(row,2,vel[2]);
	}
	void setVelocity_old(int dom, int row, double* vel){
		velocity[dom].setValue(row,3,vel[0]);
		velocity[dom].setValue(row,4,vel[1]);
		velocity[dom].setValue(row,5,vel[2]);
	}

	void getVelocity_new(int dom, int row, double* vel){
		vel[0] = velocity[dom].getValue(row,0);
		vel[1] = velocity[dom].getValue(row,1);
		vel[2] = velocity[dom].getValue(row,2);
	}

	void getVelocity_old(int dom, int row, double* vel){
		vel[0] = velocity[dom].getValue(row,3);
		vel[1] = velocity[dom].getValue(row,4);
		vel[2] = velocity[dom].getValue(row,5);
	}

	// GET Volume/mobility,fractionalflux
	//double getVolume(pEntity, const int&);
	//double getWeightedVolume(pEntity);
	double getTotalMobility(pEntity);
	double getTotalMobility(double Sw);
	double getFractionalFlux(const double&);
	double getOilFractionalFlux(const double &);
	double get_ksw(const double&);
	double get_kso(const double&);

	void retrieveSwField(pMesh);
	void storeSwField(pMesh);


//	// get set/get pointers to arrays of pointer functions
//	FuncPointer_GetGradient get_FuncPointer_GetGradient(){
//		return getGradient;
//	}

	/*! brief For steady-state simulations, total mobility must be equal 1, otherwise it must be calculated.
	 * \param state if true, lambda_total = 1.0;
	 */
	void setSimulationState(bool state){
		steady_state = state;
	}

	int getNumNodes() const{
		return nnodes;
	}

	int getNumNodeFreeWells() const{
		return nfree;
	}

	int getNumNodesWells() const{
		return nneumann;
	}

	void getNeumannIndex(int i, int &idx) const{
		idx = pWellsNeumann_index[i];
	}

	void getFreeIndex(int i, int &idx) const{
		idx = pWellsFree_index[i];
	}

	void getNonvisc(int idx, double &val){
		val = nonvisc.getValue(idx);
	}

	void setNonvisc(int idx, double val){
		nonvisc.setValue(idx,val);
	}

	void resetNonvisc(double&);

	void initializeNonvisc(){
		nonvisc.initialize(.0);
	}

	bool isInjectionWell(int idx){
		return injectionWell.getValue(idx);
	}

	void setProjectedSwgrad(int k, bool val){
		projectedSw_grad.setValue(k,val);
	}

	bool getProjectedSwgrad(int k){
		return projectedSw_grad.getValue(k);
	}

	/***********************************************************************************************/
	// MEBFV

	int getNumWells() const;						// number of wells prescribed (Neumann boundary value)
	int getNumWellNodes(int) const;					// number of nodes located on well
	void getFlowRate(int,int,int&,double&) const;	// get nodal weighted flow rate (Qi = Qt*(Vi/Vt), where:
													// 						Qt: total well flow rate;
													//						Vt: total well volume
													//						Vi: nodal volume (control volume)
	void calculateNodalFlowRate(pMesh,GeomData*);	// calculate Nodal Flow Rate
	void getTotalWellFlowRate(int flag, double& Qt);

private:
	bool steady_state;	// check if steady-state or transient simulations
	double Swr;			// Irreducible water saturation
	double Sor;			// Residual oil saturation
	double mi_w;		// water viscosity
	double mi_o;		// oil viscosity
	int ksModel;		// rel. permeability model flag

	Matrix<double>* velocity;
	Matrix<double> nonvisc;
	//Matrix<double> *SwGrad_dom;
	Matrix<bool> projectedSw_grad;
	Matrix<bool> injectionWell; // idx = 0 (not injectio/ could be production), idx = 1 free well
	int nfree;
	int nneumann;
	int* pWellsFree_index;
	int* pWellsNeumann_index;

	int nnodes;
	bool Allocated;				// says if Sw and pressure must be allocated


	/******************************************************************************************************************************/
	// MEBFV
	int numWells;						// number of wells on reservoir
	Matrix<int> numWellNodes;
	int** idxWellNodeMat;			//store node indices for well p
	double** flowrateWellNodeMat;	//store weighted flow rate (Qi=Qt*(Vi/Vt)) for each node on well p.
	std::map<int,double> totalFlowRate_map;
};
}
#endif /*PHYSICALPROPERTIESDATA_H_*/
