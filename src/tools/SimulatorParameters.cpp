#include "SimulatorParameters.h"
#include <sstream>

namespace PRS{
	
	SimulatorParameters::SimulatorParameters(){
	}
	
	SimulatorParameters::SimulatorParameters(pMesh mesh):theMesh(mesh){
		TimeStepCounter = 0;
		exportIter = 0;
		cumTS = .0;
		stop_simulation = false;
		useHOApp = false;
		restart = false;
		vtk_step = 0;
		
		PVI_increment = 0.01;	// every 1% of total simulation time a new VTK file will be printed
		PVI_cumulative = .0;	// summation of all PVI_increments
		
		allowPrintingVTK = false;
		pctype = PCNONE;
		firstVTKupdate = true;
		EBFV1_pressureSolver_scheme = false;
		_doAdaptation = false;
		vtk_time_frequency = PVI_increment*getSimTime();
		set_adapt_occur(false);
	}
	
	SimulatorParameters::~SimulatorParameters(){
		MIter_RockProperties mit = mapRockProp.begin();
		for (; mit != mapRockProp.end(); mit++){
			double *K = mit->second->K;
			delete[] K; K = 0;
		}
		mapRockProp.clear();
		mapBC.clear();
	}
	
	void SimulatorParameters::checkPermeabilityTensor(){
		const double* K = getPermeability(3300);//BACALHO
		// Kxx*Kyy >= Kxy*Kyx
		if (K[0]*K[3] < K[1]*K[2]) throw Exception(__LINE__,__FILE__,"Permeability tensor must obey the following relation: K_xx*K_yy >= K_xy*K_yx\n");
		K_Isotropic = (K[1]*K[2] == .0)?true:false;
	}
	
	const double* SimulatorParameters::getPermeability(const int &dom){
		MIter_RockProperties mit = mapRockProp.find(dom);
		if (mit != mapRockProp.end()){
			return mit->second->K;
		}
		cout << "Warning: no permeability tensor associated to domain " << dom << ".\n";
		cout << __FILE__ << "\t at line " << __LINE__ << endl;
		return 0;
	}
	
	void SimulatorParameters::getDomains(){
		for (MIter_RockProperties mit=mapRockProp.begin(); mit!=mapRockProp.end(); mit++)
			setOfDomains.insert(mit->first);
	}
	
	double SimulatorParameters::getBC_Value(const int &flag){
		MapFlagIter mIter = mapBC.find(flag);
		if (mIter != mapBC.end()){
			return mIter->second->val;
		}
		cout << "Warning: getBC_Value() return null value\n";
		cout << __FILE__ << "\t at line " << __LINE__ << endl;
		return 0;
	}
	
	bool SimulatorParameters::hasNodeWell(const int &flag){
		// if node has a prescribed saturation it's supposed that exist a well
		if ( mapSaturation.size() ){
			MapFlagIter MIter = mapSaturation.find(flag);
			return (MIter != mapSaturation.end());
		}
		else{
			isNodeFree(flag);
			return well;
		}
	}
	
	bool SimulatorParameters::isNodeFree(const int &flag){
		MapFlagIter mIter = mapBC.find(flag);
		if (mIter == mapBC.end()) return true;
		else{
			BdryConditionType *bct = mIter->second;
			return ( !bct->type.compare("dirichlet") )?false:true;
		}
	}
	
	bool SimulatorParameters::finishSimulation(){
		std:: cout << setprecision(2) << fixed;
		if (stop_simulation || (getCumulativeSimulationTime() >= getSimTime()) ){
			if (!P_pid()){
				std::cout << "#################################\n";
				std::cout << "Simulation " << (double)100.0*getCumulativeSimulationTime()/getSimTime() << "% concluded.\n";
				std::cout << "End of Simulation\n";
				std::cout << "#################################\n";
			}
			return 1;
		}
		else{
			if (!P_pid()) std::cout << "Simulation " << (double)100.0*getCumulativeSimulationTime()/getSimTime() << "% concluded.\n";
			return 0;
		}
	}
	
	// physical parameters
	double SimulatorParameters::getPorosity(const int &dom){
		MIter_RockProperties mit = mapRockProp.find(dom);
		if (mit != mapRockProp.end()) return mit->second->porosity;
		
		cout << "Warning: porosity() return null value.";
		cout << "Domain "<< dom << " not found.\n";
		cout << __FILE__ << "\t at line " << __LINE__ << endl;
		
		return 0;
	}
	
	// For a specific well, the flow rate for each node will be weighted by nodal volume.
	// Qt = (V1/Vt)Qt + (V2/Vt)Qt + ... + (Vn/Vt)Qt
	// where:
	// 		Qt 				- total flow rate (read from file)
	// 		Q1 = (V1/Vt)Qt 	- flow rate in node 1
	// 		Q2 = (V2/Vt)Qt 	- flow rate in node 2
	// 		Qn = (Vn/Vt)Qt 	- flow rate in node n
	void SimulatorParameters::weightWellFlowRateByVolume(pMesh theMesh, GeomData *pGCData){
		double Vi, Vt;
		pVertex node;
		map<int,set<int> >::iterator miter;	// for each well flag
		SIter_const siter;			// get all flagged node IDs for that well
		
		#ifdef _SEEKFORBUGS_
		if (!M_numVertices(theMesh)){
			throw Exception(__LINE__,__FILE__,"Number of vertices NULL");
		}
		if (!mapNodesOnWell.size()){
			throw Exception(__LINE__,__FILE__,"mapNodesOnWell size NULL");
		}
		#endif 
		
		// FOR EACH PRODUCTION WELL
		for (miter=mapNodesOnWell.begin(); miter!=mapNodesOnWell.end(); miter++){
			int flag = miter->first;
			// FOR EACH NODE ON PRODUCTION WELL
			Vt = .0;
		//	cout << "Number of nodes on well: " << miter->second.size() << endl;
			for (siter = miter->second.begin(); siter!=miter->second.end();siter++){
				int id = *siter;
				// to which domain node belongs get its volume
				for (SIter_const sit=setOfDomains.begin(); sit!=setOfDomains.end(); sit++)	{
					int dom = *sit;	// domains flag
					node = (mEntity*)theMesh->getVertex( id );
					
					#ifdef _SEEKFORBUGS_
					if (!node){
						char msg[256]; sprintf(msg,"NULL node with ID %d",id);
						throw Exception(__LINE__,__FILE__,msg);
					}
					#endif
					pGCData->getVolume(node,dom,Vi);
					Vt += Vi;
				}
			}
			well_info = MWells[flag];
			well_info.wellVolume = Vt;
			MWells[flag] = well_info;
//			cout << setprecision(7);
//			cout << "flag = " << flag << endl;
//			cout << "Vt   = " << Vt << endl;
//			//cout << "ID   = " << id << endl;
		}
		//STOP();
		if (!MWells.size()){
			throw Exception(__LINE__,__FILE__,"MWells size = 0\n");
		}
	}
	
	// get node ids associated to wells
	void SimulatorParameters::getWells(pMesh theMesh, int dim){
		// if user has provided some well, count how many nodes are flagged for each one
		if ( MWells.size() ){
			pEntity node;
			
			// associate a set container to each well flag
			// EX.: two injection and one production wells provided
			// flag		-	nodes IDS associated to that flag
			//  10		-	12,34,76,285,99,43,19,...
			//  11		-	1,3114,56,2285,399,643,189,...
			//  51		-	133,36,7,2855,929,463,149,...
			// initialize map container with a set container for each flag
			//cout << __LINE__ << "\t" << __FILE__ << endl;
			for (MWIter miter = MWells.begin();miter!=MWells.end();miter++){
				setNodesOnWell setNOW;
				mapNodesOnWell[miter->first] = setNOW;
			}

			int flag;
			// loop over nodes searching for wells
			VIter vit = M_vertexIter( theMesh );
			while ( (node = VIter_next(vit)) ){
				// get node's flag
				flag = (!node->getClassification())?0:GEN_tag( node->getClassification() );
				
				// check for INJECTION/PRODUCTION wells
				if ( this->isInjectionWell(flag) || this->isProductionWell(flag) ){
					MNOWIter = mapNodesOnWell.find(flag);
					// insert node ID for flag
					if ( MNOWIter!=mapNodesOnWell.end() )
						MNOWIter->second.insert( EN_id(node) );
				}
			}
			VIter_delete(vit);
		}
		else{
			throw Exception(__LINE__,__FILE__,"MWells size = 0");
		}
	}
	
	bool SimulatorParameters::isInjectionWell(int flag) const{
		MWCIter mit = MWells.find( flag );
		return (mit!=MWells.end())?mit->second.isInjection:false;
	}
	
	double SimulatorParameters::getInitialSaturation(pEntity node){
		int flag = GEN_tag(node->getClassification());
		MapFlagIter MIter = mapSaturation.find(flag);
		return (MIter == mapSaturation.end())?Sw_initial():(1.0 - Sor());
	}
	
	bool SimulatorParameters::isProductionWell(int flag) const{
		MWCIter mit = MWells.find( flag );
		return (mit!=MWells.end())?mit->second.isProduction:false;
	}
	
	double SimulatorParameters::getFlowrateValue(int flag) const{
		MWCIter mit = MWells.find( flag );
		return mit->second.flowRate;
	}
	
	double SimulatorParameters::getTotalInjectionFlowRate() const{
		double Q = 0.0;
		for (MWCIter mwiter=MWells.begin(); mwiter!=MWells.end(); mwiter++){
			//printf("well flag %d. Flow rate: %f.  Volume: %f\n",mwiter->first,mwiter->second.flowRate,mwiter->second.wellVolume);
			if ( isInjectionWell( mwiter->first ) ){
				Q += mwiter->second.flowRate;
			}
		}
		return Q;
	}
	
	double SimulatorParameters::getWellVolume(int well_flag) const{
		MWCIter mit = MWells.find( well_flag );
		return mit->second.wellVolume;
	}

	/*
	 * double GeomData::getReservoirVolume() returns the total reservoir volume.
	 * To compute the initial oil volume it's necessary to take account rock po-
	 * rosity that can vary from domain to domain. So, setInitialOilVolume() function
	 * make a loop over all control volumes from all partitions (if there are two
	 * or more)
	 */
	void SimulatorParameters::setInitialOilVolume(pMesh theMesh, GeomData *pGCData){
		int  dom, ndom, nnodes, node;
		double vol_total = .0, vol;
		ndom = pGCData->getNumDomains();
		for (dom=0; dom<ndom; dom++){
			nnodes = pGCData->getNumNodesPerDomain(dom);
			double vol_domain = .0;
			for (node=0; node<nnodes; node++){
				pGCData->getVolume(dom,node,vol);
				vol_domain += vol;
			}
			vol_total += 0.2*vol_domain;
		}
		_IOV = vol_total*(1.0 - Sw_initial());		// initial oil volume (_IOV)
		if (_IOV < 1e-8){
			throw Exception(__LINE__,__FILE__,"Initial Oil Volume null!");
		}
	}
	
	// Every N time-steps a new VTK file must be printed. It occurs at every 0.1PVI. If the sum of N time-steps exceeds 0.01PVI,
	// then the last one must be corrected to guarantee the exactness of all PVI
	void SimulatorParameters::correctTimeStep(double &timeStep){
		if (firstVTKupdate){
			updatePrintOutVTKFrequency();  
		}
		double timeFrequency = getPrintOutVTKFrequency();
		double cumST = timeStep + getCumulativeSimulationTime();
		if ( cumST > timeFrequency ){
			timeStep = timeFrequency - getCumulativeSimulationTime();
			cumST = timeFrequency;
			allowPrintingVTK = true;
		}
	}
	
	void SimulatorParameters::allowPrintVTK(){
		allowPrintingVTK = true;
	}

	double SimulatorParameters::getPrintOutVTKFrequency(){
		if (firstVTKupdate){
			updatePrintOutVTKFrequency();
		}
		return vtk_time_frequency;
	}

	void SimulatorParameters::printOutVTK(pMesh theMesh, void *pData1, void *pData2, void *pData3, void *pData4, pFunc_PrintVTK printVTK){
		allowPrintingVTK = true;
		if (allowPrintingVTK){
			int theStep = getStepOutputFile();
			incrementeStepOutputFile();
			char fname[256];
			sprintf(fname,"%s__%d-of-%d__step-%d.vtk",expofName.c_str(),P_pid(),P_size(),theStep);
			PetscPrintf(PETSC_COMM_WORLD,"VTK Output: step file #%d\n",theStep);
			printVTK(theMesh,pData1,pData2,pData3,pData4,fname);
			updatePrintOutVTKFrequency();
			allowPrintingVTK = false;
		}
		//throw 1;
	}
	
	void SimulatorParameters::updatePrintOutVTKFrequency(){
		firstVTKupdate = false;
		
		// a new VTK file will be printed at each 0.05 PVI (5% of total simulation time)
		PVI_cumulative += getPVIincrement();
		
		// when next vtk file must be print out
		vtk_time_frequency = PVI_cumulative*getSimTime();
	}
}
