/*
 * SIMULATION_core__solvers.cpp
 *
 *  Created on: 27/08/2012
 *      Author: rogsoares
 */

#include "SIMULATION_core.h"

namespace PRS{

	int SIMULATION_core::solver(){
		double timeStep;
		if (simFlag==STEADY_STATE){
			bool adapt = false;	// ............ If adaptation is not required, get out while loop.
			int count = 0;
			do{
				pElliptic_eq->solver(theMesh);
				pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,pGCData,exportSolutionToVTK);
				#ifndef NOADAPTATION
				adapt = adaptation();
				#endif
				count++;
			}while (adapt);
			if (adapt){
				cout<< "Loops :"<< count <<endl;
			}
		}
		else if (simFlag==TRANSIENT){
			while ( !pSimPar->finishSimulation() ){
				pElliptic_eq->solver(theMesh);
				pHyperbolic_eq->solver(theMesh,timeStep);
				//#ifndef NOADAPTATION
				adaptation();
				//#endif
			}
		}
		else if (simFlag==MIMPES_ADAPT){
			while ( !pSimPar->finishSimulation() ){
				pElliptic_eq->solver(theMesh);
				pHyperbolic_eq->solver(theMesh,timeStep);
				adaptation();
			}
		}
		return 0;
	}
}
