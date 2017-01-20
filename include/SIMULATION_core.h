#ifndef _IMPES_FORMULATION_H_
#define _IMPES_FORMULATION_H_

#include "solverEquation_methods.h"
#include "EBFV1__pre-processors.h"

#include "RH_Refinement.h"
#include "H_Refinement_2D.h"
#include "AdaptiveRemeshing.h"
#include "ErrorAnalysis.h"
#include "interpolation.h"


//!PRS: Petroleum Reservoir Simulator
namespace PRS{

void initializeParameters(pMesh theMesh, ErrorAnalysis* pErrorAnalysis, AMR* pMeshAdapt);

	class SIMULATION_core{
	public:

		SIMULATION_core();
		~SIMULATION_core();

		/*! brief Load data from files, initialize pointers and structures
		 * \param argc stantard C/C++ input argument
		 * \param argv stantard C/C++ input argument
		 */
		int initialize(int argc, char **argv);

		void initialize_adaptation(int argc, char **argv);

		/*! Call it when performing mesh adaptation. After mesh has been modified, many informations (e.g. number of degrees of freedom) must be recalculated.
		 */
		void updatePointersData(pMesh);

		/*!Initialize a pointer to call the solver defined by user for the elliptic term in IMPES method
		 */
		Elliptic_equation* init_EllipticSolverPointer(int );

		/*! Initialize a pointer to call the solver defined by user for the hiperbolic term in IMPES method
		 */
		Hyperbolic_equation* init_HyperbolicSolverPointer(int );

		/*! \brief: Decides which sort of simulation must be performed: steady state or transient
		 */
		int solver();

		/*! \brief: performs a transient simulation
		 */
		int transient();
		
		/*! \brief: performs a steady state simulation
		 */
		int steadyState();

		bool adaptation();

		/*! Free pointers memory, output files
		 */
		int finalize();

	private:

		/// pointer to solver pressure and saturation fields
		Elliptic_equation* pElliptic_eq;
		Hyperbolic_equation* pHyperbolic_eq;

		/**
		 *  Auxiliary pointers to access all data needed to simulation
		 */

		/// set/get pressure, saturation, velocity, flux, gradients
		PhysicPropData *pPPData;
		/// set/get parameters: CFL, filenames, methods
		SimulatorParameters *pSimPar;
		/// set/get Cij, Dij, volume, reservoir dimensions
		GeomData *pGCData;
		/// parallel stuff
		MeshData *pMData;
		/// calculate oil production
		OilProductionManagement *pOilProduction;
		/// handles mesh entities
		pMesh theMesh;

		// start mesh adaptation procedure
		AMR* pMeshAdapt;

		ErrorAnalysis* pErrorAnalysis;

		InterpolationDataStruct* pIData;

		// Function Pointer for Interpolation functions
		void (*pInterpolateData)(InterpolationDataStruct*);

		int simFlag;
		enum SIMULATION_States{STEADY_STATE, TRANSIENT, MIMPES_ADAPT};
	};
}
#endif
