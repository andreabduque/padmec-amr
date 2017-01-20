#include "SIMULATION_core.h"
#include "CPU_Profiling.h"

double CPU_Profile::tic;
std::map< string, std::list<double> > CPU_Profile::cpumap;

namespace PRS {

	SIMULATION_core::SIMULATION_core(){
		pElliptic_eq = 0;
		pHyperbolic_eq = 0;
		theMesh = MS_newMesh(0);
	}

	SIMULATION_core::~SIMULATION_core(){
	}

	int SIMULATION_core::initialize(int argc, char **argv){

		// Initialize parallel libraries
		ParUtil::Instance()->init(argc,argv);
		PetscInitialize(&argc,&argv,(char *)0,(char *)0);

		if ( argc<2 ){
			char msg[256];
			sprintf(msg,"You MUST type: ./PADAMEC_AMR.exe a b, where:\na = 0 or 1 (Steady State or Transient\n");
			throw Exception(__LINE__,__FILE__, msg );
		}

		simFlag = atoi(argv[1]);
		printSimulationHeader();

		// Initialize simulation pointers
		pPPData = new PhysicPropData;
		pGCData = new GeomData;
		pSimPar = new SimulatorParameters(theMesh);
		pMData = new MeshData(pSimPar,theMesh);

		// Load data from files
		pSimPar->setCommandLineArguments(argv,argc);
		pSimPar->inputDataFromFiles(pGCData);

		/*
		 * If adaptation required, initialize:
		 * 		1- Adaptation Pointer (h-Refinement or adaptive remeshing)
		 * 		2- Error Analysis Pointer
		 * 		3- Interpolation Function
		 */
		initialize_adaptation(argc,argv);

		/*
		 *  Initialization procedure based on previous loaded data from file:
		 *  	- well flow
		 *  	- physical properties
		 *  	- boundary and initial conditions, mapping
		 */
		pGCData->initilize(theMesh,pSimPar->setOfDomains,pSimPar->getEllipticSolver());
		pSimPar->initialize(pGCData,theMesh);
		pPPData->initialize(pMData,pSimPar,theMesh,false,pGCData);
		if (pSimPar->getEllipticSolver()==1){
			pMData->initialize(theMesh,pGCData);
		}

		// Oil production output
		pOilProduction = new OilProductionManagement(pSimPar->getOutputPathName(),pSimPar->getInitialOilVolume(),pSimPar->getTotalInjectionFlowRate(),pSimPar->useRestart());

		// Initialize elliptic and hyperbolic solver pointers
		pElliptic_eq = init_EllipticSolverPointer( pSimPar->getEllipticSolver() );
		pHyperbolic_eq = init_HyperbolicSolverPointer( pSimPar->getHyperbolicSolver() );
		return 0;
	}

	Elliptic_equation* SIMULATION_core::init_EllipticSolverPointer(int elliptic_method){
		switch ( elliptic_method ){
		// EBFV1: A vertex centered edge based finite volume formulation
		case 1:
			return new EBFV1_elliptic(theMesh,pPPData,pSimPar,pGCData,pMData);
		// MEBFV: A modified vertex centered element based finite volume formulation for highly heterogeneous porous media
		case 2:
			return new MEBFV_elliptic(theMesh,pPPData,pSimPar,pGCData,pMData);
		default:
			throw Exception(__LINE__,__FILE__,"Could not initialize a pointer to pElliptic_eq. Unknown method.\n");
		}
	}

	Hyperbolic_equation* SIMULATION_core::init_HyperbolicSolverPointer(int hyperbolic_method){
		switch ( hyperbolic_method ){
		case 1:
			return new EBFV1_hyperbolic(theMesh,pPPData,pSimPar,pGCData,pMData,pOilProduction,pErrorAnalysis);
		case 2:
			return new EBFV1_hyperbolic_MIMPES(theMesh,pPPData,pSimPar,pGCData,pMData,pOilProduction,pErrorAnalysis);
		case 3:
			return new EBFV1_hyperbolic_MIMPES_Adapt(theMesh,pPPData,pSimPar,pGCData,pMData,pOilProduction,pErrorAnalysis);
		default:
			throw Exception(__LINE__,__FILE__,"Could not initialize a pointer to pHiperbolic_eq. Unknown method.\n");
		}
	}

	void SIMULATION_core::updatePointersData(pMesh theMesh){
		cout<< "UPDATEPOINTERS"<<endl;
	#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: updating Pointers\tIN\n";
	#endif
		// starting deallocating data related to simulation pointers
		pSimPar->deallocateData();
		pPPData->deallocateData(pSimPar);
		pMData->deallocateData();


		// initialize them once more
		pSimPar->initialize(pGCData,theMesh);
		pPPData->initialize(pMData,pSimPar,theMesh,true,pGCData);
		pMData->initialize(theMesh,pGCData);
	#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: updating Pointers\tOUT\n";
	#endif
	}

	void SIMULATION_core::initialize_adaptation(int argc, char **argv){
		pErrorAnalysis = new ErrorAnalysis;

		//#ifndef NOADAPTATION
			pIData = new InterpolationDataStruct;
			//pIData->getLevelOfRefinement = ErrorAnalysis::getLevelOfRefinement;
			pIData->numFields = 2;
			pIData->pGetDblFunctions = new GetDblFunction[2];
			pIData->pSetDblFunctions = new SetDblFunction[2];

			// get data from old mesh
			pIData->pGetDblFunctions[0] = pPPData->getPressure;
			pIData->pGetDblFunctions[1] = pPPData->getSaturation;

			// set data (interpolated) to new mesh
			pIData->pSetDblFunctions[0] = pPPData->setPressure_NM;		// set pressure for New Mesh
			pIData->pSetDblFunctions[1] = pPPData->setSaturation_NM;	// set saturarion for New Mesh

//			switch( pSimPar->getRefStrategy() ){
//			case H_REFINEMENT:
//				//pMeshAdapt = new H_Refinement_2D;
//				//pIData->isElementSpecial = H_Refinement_2D::isElementSpecial;
//				break;
//			case ADAPTIVE_REMESHING:
//				//pMeshAdapt = new AdaptiveRemeshing(argc, argv);
//				break;
//			case RH_REFINEMENT:
//				//pMeshAdapt = new RH_Refinement(argc, argv);
//				break;
//			default:
//				throw Exception(__LINE__,__FILE__,"Unknown adaptation strategy.");
//			}
	//	#endif
	}

	int SIMULATION_core::finalize(){
		// Write to file oil production output. Only rank 0 is in charge of it.
		string path = pSimPar->getOutputPathName();

		//char tmp[256]; sprintf(tmp,"%s_PETSc_summary_nproc%d.log",path.c_str(),P_size());
		//PetscErrorCode ierr = PetscLogPrintSummary(MPI_COMM_WORLD,tmp); CHKERRQ(ierr);

		// free memory
		delete pElliptic_eq;
		delete pHyperbolic_eq;
		delete pPPData;
		delete pSimPar;
		delete pGCData;
		delete pMData;

		// finalize MPI
		ParUtil::Instance()->Finalize();
		return 0;
	}
}
