/*
 * EBFV1_hyperbolic_MIMPES_Adapt.h
 *
 *  Created on: 05/05/2009
 *      Author: rogerio
 */

#ifndef EBFV1_HYPERBOLIC_MIMPES_ADAPT_H_
#define EBFV1_HYPERBOLIC_MIMPES_ADAPT_H_

#include "EBFV1_hyperbolic.h"

namespace PRS{

/*
-------------------------------------------------------------------------------------------------------------------------------------------
class: EBFV1_hyperbolic_MIMPES_Adapt
-------------------------------------------------------------------------------------------------------------------------------------------
This class implements a modified IMPES (Implicit Pressure Explicit Saturation) formulation for the two-phase oil-water flow in porous media.
As pressure filed varies slowly during entire simulation, it can be kept constant for some number of time-steps. At this moment, pressure
gradient field is also hold constant and only the velocity and saturation fields are updated. The frequency of updates is controlled by an
empirical factor called DVTOL which is based in the velocity norm variation. It's expected to reduce computational CPU time significantly
where it can reduce by a factor of 10 when compared to the classical IMPES formulation. For more details about this strategy see:

Numerical Simulation of Oil-Water Displacements Using a Higher Order Control Volume Formulation in Parallel Computers with Distributed Memory
Proceedings of COBEM 2011, 21st, Brazilian Congress of Mechanical Engineering
Copyright © 2011 by ABCM
October 24-28, 2011, Natal, RN, Brazil

-------------------------------------------------------------------------------------------------------------------------------------------
*/

	class EBFV1_hyperbolic_MIMPES_Adapt : public EBFV1_hyperbolic {
	public:

		EBFV1_hyperbolic_MIMPES_Adapt();
		EBFV1_hyperbolic_MIMPES_Adapt(pMesh, PhysicPropData *, SimulatorParameters *,GeomData *, MeshData *, OilProductionManagement *, ErrorAnalysis*);
		~EBFV1_hyperbolic_MIMPES_Adapt();
		double solver(pMesh, double&);

	private:

		int sumNum_p_DT;						// summation of number of pressure timesteps
		int sumNum_Sw_DT;						//  summation of number of Sw timesteps
		double p_timestepOld;
		double cumulative_p_timestep;
		int vel_ts_counter;						// counts number of implicit time-steps
		int sat_ts_counter;						// counts number of explicit time-steps
		double DVTOL;							// velocity variation norm tolerance
		double vel_factor;						// turn velocity vector dimensionless
		double time_factor;						// turn sat_ts dimensionless
		double sat_ts;							// explicit time advance
		double vel_ts;							// implicit time advance
		double vel_ts_old;						// last implicit time advance
		bool initial_vel_ts_old;				// set vel_ts_old as sat_ts at the beginning of simulation
		bool   allowAdaptativeTimeAdvance;		// controls saturations field time advance
		bool   allowVelocityCalculation;		// make velocities filed be calculated just once per explicit (sat_ts) time-step
		double cumulativeExplicitTS;			// cumulativeExplicitTS = sum(sat_ts)
		double dv_norm;							// implicitTS returns a dimensionsless delta T

		/**
		 * Print on file adaptative data:
		 * 	- time (PVI)
		 *  - explicit_counter/implicit_counter ratio
		 *  - velocity norm
		 */
		ofstream fid;
		bool print_ATE_Data;
		void MIMPES_output(double Sw_timestep_sum, double p_timestep, int numCFL_Steps, int sumNum_p_DT, int sumNum_Sw_DT, double DV_norm);

		// Compute a new implicit time-step (vel_ts) if necessary. This happens either at the beginning of simulation or when sum of sat_ts be greater than DT.
		void calculateImplicitTS(double &p_timestep, double timeStep, double &DV_norm, bool &go_ImplicitTS, bool &go_MIMPES);

		// Compute a new implicit time-step as function of the last implicit time-step, velocity norm and DVTOL. DVTOL comes from the input data file 'numeric.dat'.
		void calculateVelocityVariationNorm(double &DV_norm, int dim);

		// Reset logical variables which control programming flux.
		void resetVariables();

		// Time_factor and vel_factor are used to turn implicit time-step and velocity vector dimensionless variables
		void createDimensionLessFactors();

		// Verify if the explicit time advance (summation of sat_ts) do not exceeded the limit imposed by the implicit time advance (vel_ts)
		void correct_p_TS(double &p_timestep, bool &go_MIMPES);
		void correct_Sw_TS(double &timeStep, bool &go_MIMPES);

		void setCumulative_p_TS(double ts){
			cumulative_p_timestep = ts;
		}

		double getCumulative_p_TS() const{
			return cumulative_p_timestep;
		}
	};
}
#endif /* EBFV1_HYPERBOLIC_MIMPES_ADAPT_H_ */
