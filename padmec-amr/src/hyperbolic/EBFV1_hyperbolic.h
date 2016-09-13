/*
 * EBFV1-hyperbolic.h
 *
 *  Created on: 09/01/2009
 *      Author: rogerio
 */

#ifndef EBFV1HYPERBOLI_H_
#define EBFV1HYPERBOLI_H_

#include "Hyperbolic_equation.h"
#include "CPU_Profiling.h"

namespace PRS{

	class EBFV1_hyperbolic : public Hyperbolic_equation{
	public:

		EBFV1_hyperbolic();
		EBFV1_hyperbolic(pMesh, PhysicPropData*, SimulatorParameters*,GeomData*, MeshData*, OilProductionManagement*, ErrorAnalysis*);
		~EBFV1_hyperbolic();
		virtual double solver(pMesh, double&);

	protected:

		// Main functions related to the advective equation presented in order that they must be called.
		double calculateVelocityField(int,int);
		double calculateIntegralAdvectiveTerm(int, double&);
		void calculateExplicitAdvanceInTime(double);

		void nodeWithOut_Wells(double);
		void nodeWith_Wells(double);
		
		/// For mesh adaptation, saturation field interpolation between old and new mesh over Sw_t, and not Sw_t+1
		void saveSwField();

		void setRecoveredOil(double val){
			oilRecovered = val;
		}

		double getCumulativeOil() const{
			return _cumulativeOil;
		}

		void setCumulativeOil(double co){
			_cumulativeOil = co;
		}

		double getRecoveredOil() const{
			return oilRecovered;
		}

		double oilRecovered;

		// Numeric parameter to compute nonvisc term
		double alpha_max;

		OilProductionManagement *pOPManager;
		MeshData *pMData;
		GeomData *pGCData;
		PhysicPropData *pPPData;
		SimulatorParameters *pSimPar;
		ErrorAnalysis *pEA;

		// Saturation gradient is calculated for all domains at once. Nodes on boundary domains contains one gradient vector for each domain.
		void calculateSaturationGradient();

	private:
		// every new time-step, nodal gradient must be set to zero and start a new calculation
		void resetSaturationGradient();
		
		// loop over all edges (omega domain
		void calc_Sw_grad_1(int, int);
		
		// loop over all boundary edges (2-D) or all external faces (3-D)
		void calc_Sw_grad_2(int, int);
		
		// Averaging by Total Volume in 3-D (area in 2-D problems)
		void calc_Sw_grad_3(int);

		// Imposition of Homogeneus Neumman Boundary Conditions
		void calc_Sw_grad_4(int);
		
		// If restart is required, then oil production file data is used to update the cumulated oil value
		void setCumulativeOilProd();

		double _cumulativeOil;

		PetscErrorCode ierr;

		bool PRINT_DEBUG;

		ofstream fid_PRINT_DEBUG;
	};
}


#endif /* EBFV1HYPERBOLI_H_ */
