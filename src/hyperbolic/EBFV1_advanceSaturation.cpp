#include "EBFV1_hyperbolic.h"

namespace PRS{

	// Euler Foward - time advance:  (compute saturation explicitly). For each new time step compute a new saturation.
	void EBFV1_hyperbolic::calculateExplicitAdvanceInTime(double delta_T){
		CPU_Profile::Start();

		saveSwField();					// save Sw field (Sw_t) before calculate Sw_t+1
		nodeWithOut_Wells(delta_T);		// Advance time for all free node not located in wells
		nodeWith_Wells(delta_T);		// Advance time for nodes in production well (fractional-steps)

		CPU_Profile::End("ExplicitAdvanceInTime");
	}

	void EBFV1_hyperbolic::nodeWithOut_Wells(double delta_T){
		double Sw,Sw_old,nonvisc,volume;
		int idx;
		int nnode = pPPData->getNumNodeFreeWells();
		for(int i=0; i<nnode; i++){
			pPPData->getFreeIndex(i,idx);
			pPPData->getSaturation(idx,Sw_old);
			pPPData->getNonvisc(idx,nonvisc);
			pGCData->getVolume(idx,volume);

			// alterar para usar valor de porosidade lido de arquivo
			Sw = Sw_old - (delta_T/(0.2*volume))*nonvisc;

			if ( Sw<1.0e-6 ) Sw = .0;
			pPPData->setSaturation(idx,Sw);
	//#ifdef _SEEKFORBUGS_
			if (Sw > 1.01 || Sw < .0){
				char msg[256]; sprintf(msg,"Sw[%d/%d] = %.4f. 0.0 <= Sw <-1.0 \n",i,nnode,Sw);
				throw Exception(__LINE__,__FILE__,msg);
			}
	//#endif
		}
	}

	void EBFV1_hyperbolic::nodeWith_Wells(double delta_T){
		const int N = pSimPar->getWellTimeDiscretion();		// get number of delta subdivisions
		double dt_well = (double)delta_T/N;					// time step for well nodes saturation

		int i,j, well_idx;
		double Sw,Sw0,Sw_old,Vt,Qi,Vi,Qwi,Qt,fw,fo,nonvisc,wp,cml_oil,Qo,Qw;

		// TODO: nao usar constantes para identificar pocos!
		Qt = pSimPar->getFlowrateValue(51);					// source/sink term
		Vt = pSimPar->getWellVolume(51);					// for node i, Qi is a fraction of total well flow rate (Qt)
		cml_oil = .0;
		Qo = .0;
		Qw = .0;
		int nnodes = pPPData->getNumNodesWells();
		for (i=0; i<nnodes; i++){
			pPPData->getNeumannIndex(i,well_idx);
			pPPData->getNonvisc(well_idx,nonvisc);
			pPPData->getSaturation(well_idx,Sw_old);
			pGCData->getVolume(well_idx,Vi);
			Qi = Qt*(Vi/Vt);								// Fluid (water+oil) flow rate through node i
			wp = 0.2*Vi;
			for (j=0; j<N; j++){
				Sw0 = Sw_old - (dt_well/wp)*(nonvisc);
				fw = pPPData->getFractionalFlux(Sw0);
				Qwi = fabs(fw*Qi);
				Sw = Sw0 - dt_well*(Qwi/wp);
				Sw_old = Sw;
			}
	#ifdef _SEEKFORBUGS_
			if (Sw > 1.01 || Sw < .0){
				char msg[256]; sprintf(msg,"Water saturation is out-of-bound [0 1]. Sw = %.4f\n",Sw);
			}
	#endif

			pPPData->setSaturation(well_idx,Sw);
			fw = pPPData->getFractionalFlux(Sw);	    // oil fractional flux
			fo = pPPData->getOilFractionalFlux(Sw);		// oil fractional flux
			Qo += fabs(Qi*fo);
			Qw += fabs(Qi*fw);
			cml_oil += Qo;
		}

		setRecoveredOil(Qo/(Qo+Qw));
		cml_oil = cml_oil*delta_T + getCumulativeOil();
		setCumulativeOil(cml_oil);
	}

	void EBFV1_hyperbolic::saveSwField(){
		double Sw;
		int i, nnodes;
		pGCData->getMeshNodes(nnodes);
		for(i=0; i<nnodes; i++){
			pPPData->getSaturation(i,Sw);
			pPPData->setSw_old(i,Sw);
		}
	}
}
