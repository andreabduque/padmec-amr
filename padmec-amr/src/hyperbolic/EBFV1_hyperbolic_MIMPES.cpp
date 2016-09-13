/*
 * EBFV1-hyperbolic.cpp
 *
 *  Created on: 09/01/2009
 *      Author: rogerio
 */

#include "EBFV1_hyperbolic_MIMPES.h"
#include "exportVTK.h"

namespace PRS{

	EBFV1_hyperbolic_MIMPES::EBFV1_hyperbolic_MIMPES(){
	}

	EBFV1_hyperbolic_MIMPES::EBFV1_hyperbolic_MIMPES(pMesh mesh,PhysicPropData *ppd, SimulatorParameters *sp,GeomData *gcd,MeshData *md, OilProductionManagement *popm, ErrorAnalysis *pEAna){
		pMData = md;
		pOPManager = popm;
		pPPData = ppd;
		pGCData = gcd;
		pSimPar = sp;
		pEA = pEAna;
		DVTOL = pSimPar->getDVTOL();
		createDimensionLessFactors();

		string filename = pSimPar->getOutputPathName() + "_MIMPES_monitor.csv";
		fid.open(filename.c_str());
		if (!fid.is_open()) {
			char msg[256]; sprintf(msg,"File '%s' could not be opened or it does not exist.\n",filename.c_str());
			throw Exception(__LINE__,__FILE__,msg);
		}
		fid << "PVI CFL_dt/imp_dt sum_numSw_dt/sum_nump_dt vel_norm" << endl;
		sumNum_p_DT = 0;
		sumNum_Sw_DT = 0;
		p_timestepOld = -1.0;
		cumulative_p_timestep = 0;
	}

	EBFV1_hyperbolic_MIMPES::~EBFV1_hyperbolic_MIMPES(){
		fid.close();
	}

	double EBFV1_hyperbolic_MIMPES::solver(pMesh theMesh, double &timeStep){
		if (!P_pid()) std::cout << "Adaptative Hyperbolic solver...\n";

		bool go_MIMPES = true;                                                  // says to keep Sw advance while velocity is hold constant
		bool go_ImplicitTS = true;                                              // calculate implicit time step once while Sw advances
		static int timestep_counter = 0;                                // counts number of time steps every new VTK
		int dim = theMesh->getDim();
		timeStep = 1.0e+10;                                                             // initialize time step with a very high number
		int ndom = (int)pSimPar->setOfDomains.size();

		// calculate velocity field
		for (int dom=0; dom<ndom; dom++){
			calculateVelocityField(dom,dim);
		}
		double DV_norm,  p_timestep, Sw_timestep_sum = .0;
		calculateVelocityVariationNorm(DV_norm,dim);            // implicitTS returns a dimensionless delta T
		// calculate how long saturation equation will be solved while pressure field is kept constant
		// p_DT = sum (Sw_DT)
		int numCFL_Steps = 0;
		sumNum_p_DT++;

		do{     // keep velocity field constant if possible
			// calculate saturation gradient if adaptation or high order approximation were required
			if (  pSimPar->userRequiresAdaptation() || pSimPar->useHOApproximation()){
				calculateSaturationGradient();
			}

			timeStep = 1.0e+10;
			pPPData->resetNonvisc(alpha_max);
			for (int dom=0; dom<ndom; dom++){
				calculateIntegralAdvectiveTerm(dom,timeStep);
			}

			calculateImplicitTS(p_timestep,timeStep,DV_norm,go_ImplicitTS,go_MIMPES);
			correct_Sw_TS(timeStep,go_MIMPES);                                                              // do not allow Sw_timestep_sum be greater than p_timestep
			pSimPar->setCumulativeSimulationTime(timeStep);                                 // AccSimTime = AccSimTime + timeStep
			calculateExplicitAdvanceInTime(timeStep);                                               // Calculate saturation field: Sw(n+1)
			timestep_counter++;
			sumNum_Sw_DT++;                                                                                                 // for whole simulation
			numCFL_Steps++;                                                                                                 // for implicit time steps only
			Sw_timestep_sum += timeStep;                                                                    // CFL timestep summation

			// oil production output
			if (pSimPar->timeToPrintVTK()){
				pOPManager->printOilProduction(timeStep,pSimPar->getCumulativeSimulationTime(),pSimPar->getSimTime(),getRecoveredOil(),getCumulativeOil(),timestep_counter);
				timestep_counter = 0;
			}
			pSimPar->printOutVTK(theMesh,pPPData,pEA,pSimPar,pGCData,exportSolutionToVTK);
		}while ( go_MIMPES );
		MIMPES_output(Sw_timestep_sum,p_timestep,numCFL_Steps,sumNum_p_DT,sumNum_Sw_DT,DV_norm);
		if (!P_pid()) cout << "done\n";
		return 0;
	}

	// Compute a new implicit time-step (p_timestep) if necessary.
	// This happens either at the beginning of simulation or when sum of timeStep is greater than DT.
	void EBFV1_hyperbolic_MIMPES::calculateImplicitTS(double &p_timestep, double timeStep, double &DV_norm, bool &go_ImplicitTS, bool &go_MIMPES){
		if (go_ImplicitTS){
			CPU_Profile::Start();
			if (p_timestepOld<.0){							// initialize p_timestepOld for the beginning of simulation (t=0)
				p_timestepOld = timeStep;
			}
			p_timestepOld *= time_factor;					// make implicit p_timestepOld variable dimensionless
			p_timestep = (DVTOL/DV_norm)*p_timestepOld;		// Update new implicit time step: p_timestep
			double RT = (double)(p_timestep/p_timestepOld);	// dimensionless DT ratio
			if (RT>1.25){
				p_timestep = 1.25*p_timestepOld;
			}
			else if ( RT<0.75 ){
				p_timestep = 0.75*p_timestepOld;
			}
			p_timestep /= time_factor;								// gives implicit p_timestep variable dimension of time
			if (p_timestep < timeStep || (RT>0.75 && RT<1.25)){     // It doesn't make sense p_timestep be less than timeStep
				p_timestep = timeStep;
			}
			correct_p_TS(p_timestep,go_MIMPES);						// do not allow p_timestep go beyond total simualtion time
			go_ImplicitTS = false;
			p_timestepOld = p_timestep;
			CPU_Profile::End("ImplicitTS");
		}
	}

	void EBFV1_hyperbolic_MIMPES::correct_p_TS(double &p_timestep, bool &go_MIMPES){
		double p = getCumulative_p_TS() + p_timestep;
		if ( p > pSimPar->getSimTime() ){
			p_timestep = pSimPar->getSimTime()-getCumulative_p_TS();
			p = pSimPar->getSimTime();
			pSimPar->stopSimulation();
		}
		setCumulative_p_TS(p);
	}

	void EBFV1_hyperbolic_MIMPES::correct_Sw_TS(double &timeStep, bool &go_MIMPES){
		// CFL time step must be corrected to satisfy VTK print out or MIMPES procedure
		double cum_VTK = pSimPar->getPrintOutVTKFrequency();
		double cum_p_TS = getCumulative_p_TS();
		double cum_ST = pSimPar->getCumulativeSimulationTime();

		// fit timeStep to MIMPES timestep
		if (cum_ST+timeStep > cum_p_TS){
			timeStep = cum_p_TS - cum_ST;
			if (timeStep<0){throw Exception(__LINE__,__FILE__,"time step negative!");}
			// fit timeStep to PVI time for right VTK print out
			if (cum_ST+timeStep > cum_VTK){
				timeStep = cum_VTK - cum_ST;
				if (timeStep<0){throw Exception(__LINE__,__FILE__,"time step negative!");}
				pSimPar->allowPrintVTK();
			}
			else{
				go_MIMPES = false;
			}
		}
		if (cum_ST+timeStep > cum_VTK){
			timeStep = cum_VTK - cum_ST;
			pSimPar->allowPrintVTK();
		}
		if (cum_ST+timeStep >= pSimPar->getSimTime()){
			timeStep = pSimPar->getSimTime() - cum_ST;
			go_MIMPES = false;
			pSimPar->allowPrintVTK();
		}
	}

	void EBFV1_hyperbolic_MIMPES::calculateVelocityVariationNorm(double &DV_norm, int dim){
		CPU_Profile::Start();

		int i, dom, ndom, edge, nedges, numGEdges;
		double Cij[3], v_new[3], v_old[3], Dv[3], DV_domNorm[4], Cij_norm, aux, dot, DV_norm_sum;

		pGCData->getTotalNumberOfEdges(numGEdges);
		ndom = pSimPar->getNumDomains();
		for (dom=0; dom<ndom; dom++){
			nedges = pGCData->getNumEdgesPerDomain(dom);
			DV_domNorm[dom] = .0;
			for (edge=0; edge<nedges; edge++){
				pGCData->getCij(dom,edge,Cij);
				pGCData->getCij_norm(dom,edge,Cij_norm);
				pPPData->getVelocity_new(dom,edge,v_new);
				pPPData->getVelocity_old(dom,edge,v_old);
				for (i=0; i<dim; i++){
					Dv[i] = vel_factor*(v_new[i]-v_old[i]);
				}
				dot = inner_product(Dv,Cij,dim);
				aux = dot/Cij_norm;
				DV_domNorm[dom] += aux*aux;
			}
			DV_domNorm[dom] /= numGEdges;
		}

		DV_norm_sum = 0;
		for (dom=0; dom<ndom; dom++){
			DV_norm_sum += DV_domNorm[dom];
		}
		aux = (double)(DV_norm_sum/ndom);
		DV_norm = sqrt(aux);

		CPU_Profile::End("VelocityVariationNorm");
	}

	void EBFV1_hyperbolic_MIMPES::MIMPES_output(double Sw_timestep_sum, double p_timestep, int numCFL_Steps, int sumNum_p_DT, int sumNum_Sw_DT, double DV_norm){
		fid << setprecision(6) << fixed
				<< (double)(pSimPar->getCumulativeSimulationTime()/pSimPar->getSimTime()) << " "
				<< (double)(Sw_timestep_sum/numCFL_Steps)/p_timestep << " "
				<< (double)(sumNum_Sw_DT/sumNum_p_DT) << " "
				<< DV_norm << endl;
	}

	void EBFV1_hyperbolic_MIMPES::createDimensionLessFactors(){
		double L, H, phi;
		phi = pSimPar->getPorosity(3300);
		pSimPar->getReservoirGeometricDimensions(L,H);
		double Q = pSimPar->getTotalInjectionFlowRate();
		time_factor = Q/(phi*(L*L)*H);
		vel_factor = L*H/Q;
	}
}

