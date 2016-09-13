/*
 * EBFV1-hyperbolic.cpp
 *
 *  Created on: 09/01/2009
 *      Author: rogerio
 */

#include "EBFV1_hyperbolic.h"

namespace PRS{

	EBFV1_hyperbolic::EBFV1_hyperbolic(){
	}

	EBFV1_hyperbolic::EBFV1_hyperbolic(pMesh mesh, PhysicPropData *ppd,SimulatorParameters *sp, GeomData *gcd,MeshData *md, OilProductionManagement *popm, ErrorAnalysis *pEAna){
		pGCData = gcd;
		pOPManager = popm;
		pPPData = ppd;
		pSimPar = sp;
		pMData = md;
		_cumulativeOil = .0;
		pEA = pEAna;
		setCumulativeOilProd();
	}

	EBFV1_hyperbolic::~EBFV1_hyperbolic(){
	}

	double EBFV1_hyperbolic::solver(pMesh theMesh, double &timeStep){
	#ifdef _SEEKFORBUGS_
		if (!P_pid()) std::cout << "Hyperbolic solver...\n";
		if (!M_numEdges(theMesh)){
			throw Exception(__LINE__,__FILE__,"Number of edges: 0!\n");
		}
	#endif
		static int timestep_counter = 0;			// counts number of time steps every new VTK
		int dim = theMesh->getDim();

		// initialize time step with a very high number
		timeStep = 1.0e+10;
		int ndom = (int)pSimPar->setOfDomains.size();
		for (int dom=0; dom<ndom; dom++){
			calculateVelocityField(dom,dim);
		}

		// calculate saturation gradient if adaptation or high order approximation were required
		if (  pSimPar->userRequiresAdaptation() || pSimPar->useHOApproximation()){
			calculateSaturationGradient();
		}

		pPPData->resetNonvisc(alpha_max);
		for (int dom=0; dom<ndom; dom++){
			calculateIntegralAdvectiveTerm(dom,timeStep);
		}

		pSimPar->correctTimeStep(timeStep);					// correct time-step value to print out the desired simulation moment
		pSimPar->saveCurrentSimulationTimes();				// it must be called before cumulative simulation time
		pSimPar->setCumulativeSimulationTime(timeStep); 	// AccSimTime = AccSimTime + timeStep
		calculateExplicitAdvanceInTime(timeStep);			// Calculate saturation field: Sw(n+1)

	#ifdef _SEEKFORBUGS_
		std::cout << "<_SEEKFORBUGS_>  timeStep = : " << timeStep << endl;
		throw_exception(timeStep < 1e-10,"Time step NULL!",__LINE__,__FILE__);
	#endif

		timestep_counter++;
		// oil production output
		if (pSimPar->timeToPrintVTK()){
			pOPManager->printOilProduction(timeStep,pSimPar->getCumulativeSimulationTime(),pSimPar->getSimTime(),getRecoveredOil(),getCumulativeOil(),timestep_counter);
			timestep_counter = 0;
		}

		pSimPar->printOutVTK(theMesh,pPPData,pEA,pSimPar,pGCData,exportSolutionToVTK);
		if (!P_pid()) std::cerr << "done.\n\n";
		return 0;
	}

	void EBFV1_hyperbolic::setCumulativeOilProd(){
		// set cumulative oil production
		if (pSimPar->useRestart()){
			string expofn;
			pSimPar->getExportFileName(expofn);
			char str[512]; sprintf(str,"%s_oil-production-%d.csv",expofn.c_str(),P_size());
			ifstream fid;
			fid.open(str);

			string strline;
			for (int i=0; i<=pSimPar->getLastPVI(); i++){
				getline(fid,strline);
				//cout << "i = " << i << "  " << strline << endl;
			}
			string data[5];
			fid >> data[0] >> data[1] >> data[2] >> data[3] >> _cumulativeOil >> data[4];
			cout << "\t" << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << _cumulativeOil << "\t" << data[4];
			_cumulativeOil *= pOPManager->getInitialOilVolume();
			cout << "\n_cumulativeOil = " << _cumulativeOil << endl;
			cout << "pSimPar->getLastPVI() = " << pSimPar->getLastPVI() << endl;
			fid.close();
			//exit(1);
		}
	}
}
