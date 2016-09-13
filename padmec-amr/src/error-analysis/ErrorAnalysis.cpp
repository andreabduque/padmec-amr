#include "ErrorAnalysis.h"

using namespace PRS;

bool calculate_ErrorAnalysis(ErrorAnalysis *pEA, SimulatorParameters *pSimPar, GeomData* pGCData, void(*pFunc_getGrad)(FIELD,int,int,int,const double*&), std::list<int>& elemList, std::map<int,double>& nodeMap){
	CPU_Profile::Start();
#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: Error Analysis\tIN\n";
#endif

	pEA->initialize(pGCData,pSimPar);

	bool Sw_adapt = analyzeField(SATURATION,pEA,pSimPar,pGCData,pFunc_getGrad);
	bool p_adapt = analyzeField(PRESSURE,pEA,pSimPar,pGCData,pFunc_getGrad);

	if (Sw_adapt || p_adapt){
		double param1 = pSimPar->Remeshing_param1();
		double param2 = pSimPar->Remeshing_param2();
		pEA->calculate_h_ratio(pGCData);
		pEA->getElementsForAdaptation(param1,param2,pGCData,elemList);
		pEA->getNodesForAdaptation(pGCData,nodeMap);
	}


#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: Error Analysis\tOUT\n";
#endif

	CPU_Profile::End("Error Analysis");
	return Sw_adapt || p_adapt;
}

bool analyzeField(FIELD field, ErrorAnalysis *pEA, SimulatorParameters *pSimPar,GeomData* pGCData, void(*pFunc_getGrad)(FIELD,int,int,int,const double*&)){
	double tol1, tol2;
	bool adapt = false;

	const bool NO_SINGULARITY = false;
	const bool WITH_SINGULARITY = true;

	switch (field){
	case PRESSURE:
		cout << "Pressure:\n";
		cout << "====================================================================================\n";
		tol1 = pSimPar->getp_Tol1();
		tol2 = pSimPar->getp_Tol2();
		break;
	case SATURATION:
		cout << "Saturation:\n";
		cout << "====================================================================================\n";
		tol1 = pSimPar->getSw_Tol1();
		tol2 = pSimPar->getSw_Tol2();
		break;
	}

	// reset global error, element error, etc.
	pEA->resetVariables(pGCData->getNumElements());

	// calculate errors for all mesh elements
	pEA->calculate_ElementsError(pSimPar,pGCData,pFunc_getGrad,field);
	pEA->calculate_SmoothedGradientNorm(pSimPar,pGCData,pFunc_getGrad,field);
	pEA->calculate_GlobalError(pGCData);

//	if (pEA->getGlobalError() > tol1){
		pEA->calculate_AvgError(pGCData->getNumElements(),tol1,NO_SINGULARITY);
		pEA->calculate_h_new(pGCData, pEA->getAverageError(),field);
		pEA->identify_singular_regions(pGCData,field);

		// calculate errors for all mesh elements excluding those flagged as singular
		if (pEA->getNumElements_Singularity()) {

			// Repeat calculations for all mesh elements excluding those flagged as singular
			pEA->calculate_SmoothedGradientNorm_Singularity(pSimPar,pGCData,pFunc_getGrad,field);
			pEA->calculate_GlobalError_Singularity(pGCData);

			//if (pEA->getGlobalError_Singularity() > tol2){
				pEA->calculate_AvgError(pGCData->getNumElements(),tol2,WITH_SINGULARITY);
				pEA->calculate_h_new(pGCData,pEA->getAverageError_Singularity(),field);
			//}
		}
		adapt = true;
//	}
	return adapt;
}
