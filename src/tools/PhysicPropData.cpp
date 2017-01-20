#include "PhysicPropData.h"

Matrix<double> *pGrad_dom;		// each pointer is a matrix nx3
Matrix<double> *SwGrad_dom; // DO NOT REMOVE IT FROM THE CODE!!!!!!!!!!
Matrix<double> pressure;
Matrix<double> SwGrad;			// a unique matrix nnodesx3 for whole mesh
Matrix<double> Sw;				// for domain
Matrix<double> Sw_old;			// for domain
Matrix<double> Sw_tmp;
Matrix<double> p_tmp;

namespace PRS{

	PhysicPropData::PhysicPropData(){
		steady_state = false;
		Allocated = false;
	}

	PhysicPropData::~PhysicPropData(){
	}

	void PhysicPropData::update_Sw_p(int n){
		allocateTemporaryData(n);
	}

	void PhysicPropData::allocateTemporaryData(int n){
		Sw_tmp.allocateMemory(n);
		Sw_tmp.initialize(0);
		p_tmp.allocateMemory(n);
		p_tmp.initialize(0);
	}

	void PhysicPropData::transferTmpData(){
		// deallocate Sw and pressure matrices first
		Sw.freeMemory();
		pressure.freeMemory();
		Sw_old.freeMemory();

		// allocate with new dimension
		int rows,cols;
		Sw_tmp.getsize(rows,cols);
		Sw.allocateMemory(rows);
		pressure.allocateMemory(rows);
		Sw_old.allocateMemory(rows);
		Sw_old.initialize(0);

		// transfer data
		for(int i=0;i<rows;i++){
			Sw.setValue(i,Sw_tmp.getValue(i));
			pressure.setValue(i,p_tmp.getValue(i));
		}
		// deallocate temporary data
		Sw_tmp.freeMemory();
		p_tmp.freeMemory();
	}

	void PhysicPropData::allocateData(SimulatorParameters *pSimPar, GeomData* pGCData, int numnodes){
		int ndom = pGCData->getNumDomains();
		pGrad_dom = new Matrix<double>[ndom];
		SwGrad_dom = new Matrix<double>[ndom];
		velocity = new Matrix<double>[ndom];
		for (int k=0; k<ndom; k++){
			int nrows = pGCData->getNumNodesPerDomain(k);
			int nedges = pGCData->getNumEdgesPerDomain(k);
			pGrad_dom[k].allocateMemory(nrows,3);
			pGrad_dom[k].initialize(.0);
			SwGrad_dom[k].allocateMemory(nrows,3);
			SwGrad_dom[k].initialize(.0);
			velocity[k].allocateMemory(nedges,6);
			velocity[k].initialize(.0);
		}
		nnodes = numnodes;
		SwGrad.allocateMemory(nnodes,3);
		SwGrad.initialize(.0);
		injectionWell.allocateMemory(nnodes);
		projectedSw_grad.allocateMemory(nnodes);
		nonvisc.allocateMemory(nnodes);

		// if adaptation used, Sw and pressure will be allocated in transferTmpData. This is for the very first time. Beginning of simulation.
		if (!Allocated){
			Sw.allocateMemory(nnodes);
			Sw.initialize(0);
			Sw_old.allocateMemory(nnodes);
			Sw_old.initialize(0);
			pressure.allocateMemory(nnodes);
			pressure.initialize(0);
			Allocated = true;
		}
	}

	void PhysicPropData::deallocateData(SimulatorParameters *pSimPar){
		int ndom = pSimPar->getNumDomains();
		for (int k=0; k<ndom; k++){
			pGrad_dom[k].freeMemory();
			SwGrad_dom[k].freeMemory();
			velocity[k].freeMemory();
		}
		SwGrad.freeMemory();
		injectionWell.freeMemory();
		projectedSw_grad.freeMemory();
		nonvisc.freeMemory();
		delete[] pGrad_dom; pGrad_dom = 0;
		delete[] SwGrad_dom; SwGrad_dom = 0;
		delete[] velocity; velocity = 0;

		if (!Allocated){
			pressure.freeMemory();
			Sw.freeMemory();
			Sw_old.freeMemory();
		}
	}

	void PhysicPropData::initialize(MeshData *pMData, SimulatorParameters *pSimPar, pMesh theMesh, bool isUpdate, GeomData* pGCData){
		if (pSimPar->getEllipticSolver()==2){
			return;
		}
		allocateData(pSimPar,pGCData,M_numVertices(theMesh));
		Swr = pSimPar->Swr();				// Irreducible water saturation
		Sor = pSimPar->Sor();				// Residual oil saturation
		mi_w = pSimPar->waterViscosity();	// water viscosity
		mi_o = pSimPar->oilViscosity();		// oil viscosity

		// do not set initial saturation field for update (mesh adaptation)
		if (!isUpdate){
			setInitialSaturation(theMesh,pSimPar);
		}

		double Sw;
		int idx = 0;
		int k = 0;
		nfree = 0, nneumann = 0;
		pEntity node;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			int flag = GEN_tag(node->getClassification());
			getSaturation(idx,Sw);
			if ( Sw < 1e-8 ){
				Sw = .0;
			}

			this->setSaturation(idx,Sw);
			// todo: TIRAR ESSE BACALHO DAQUI!!!
			if ( pSimPar->isProductionWell(flag) ){
				nneumann++;
				injectionWell.setValue(k++,false);
			}
			if ( !pSimPar->isInjectionWell(flag) && !pSimPar->isProductionWell(flag) ){
				nfree++;
				injectionWell.setValue(k++,false);
			}
			idx++;
		}
		VIter_delete(vit);

		pWellsFree_index = new int[nfree];
		pWellsNeumann_index = new int[nneumann];
		idx = 0;
		int free_idx = 0;
		int neumann_idx = 0;
		vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			int flag = GEN_tag(node->getClassification());
			if ( pSimPar->isProductionWell(flag) ){
				this->pWellsNeumann_index[neumann_idx] = idx;
				neumann_idx++;
			}
			if ( !pSimPar->isInjectionWell(flag) && !pSimPar->isProductionWell(flag) ){
				this->pWellsFree_index[free_idx] = idx;
				free_idx++;
			}
			idx++;
		}
		VIter_delete(vit);
		ksModel = pSimPar->ksModel();

#ifdef __SEEKFORBUGS__
		if (!numNodes || !numEdges || !ndom)
			throw Exception(__LINE__,__FILE__,"Null value! Exiting....\n");
#endif
	}

	void PhysicPropData::setInitialSaturation(pMesh theMesh, SimulatorParameters *simPar){
		pVertex node;
		double Sw;
		int idx = 0, well = 0;

		VIter vit = M_vertexIter(theMesh);
		// If restart is required, then load saturation field using the last VTK file generated.
		if (simPar->useRestart()){
			ifstream fid;
			string strline;
		//	cout << "restart file name:  " << simPar->getRestartFilename().c_str() << endl;
			fid.open(simPar->getRestartFilename().c_str());

			if (!fid.is_open()){
				throw Exception(__LINE__,__FILE__,"File could not be opened!\n");
			}

			// set position to start reading
			do{
				getline(fid,strline);
			}while( strline.compare("SCALARS Saturation float 1") );
			getline(fid,strline);

			// load data
			//cout << "Loading initial saturation field using the last VTK file generated.\n";
			while ( (node = VIter_next(vit)) ){
				fid >> Sw;
				setSaturation(idx,Sw);
				idx++;
			}
			well = 1;
			fid.close();
		}
		else{
			while ( (node = VIter_next(vit)) ){
				Sw = simPar->getInitialSaturation(node);
				if ( Sw > .0 ){
					//printf("Injection well located in node %d Sw = %f\n",EN_id(node),Sw);
					well = 1;
				}
				setSaturation(idx,Sw);
				idx++;
			}
		}
		VIter_delete(vit);

		// rank with injection well must say to all other ranks that it's OK, a injection well was being informed!
		well = P_getMaxInt(well);

		if (!well){
			throw Exception(__LINE__,__FILE__,"Injection wells are missing!");
		}
	}

	double PhysicPropData::getTotalMobility(double Sw){
		double krw = get_ksw(Sw);
		double kro = get_kso(Sw);
		return  krw/mi_w + kro/mi_o;
	}

	double PhysicPropData::getFractionalFlux(const double &Sw){
		double krw = pow((Sw - Swr)/(1. - Swr - Sor),2);
		double kro = pow((1. - Sw - Swr)/(1. - Swr - Sor),2);
		return krw/( krw + (mi_w/mi_o)*kro );
	}

	double PhysicPropData::getOilFractionalFlux(const double &Sw){
		double krw = pow((Sw - Swr)/(1. - Swr - Sor),2);
		double kro = pow((1. - Sw - Swr)/(1. - Swr - Sor),2);
		return (kro/mi_o)/( krw/mi_w + kro/mi_o );
	}

	double PhysicPropData::get_ksw(const double &Sw){
		switch(ksModel){
		case 1:
			return Sw;
		case 2:
			return pow(Sw,2);
		case 3:
			return pow((Sw - Swr)/(1. - Swr - Sor),2);
		default:
			throw Exception(__LINE__,__FILE__,"Unknown model for relative permeability.\n");
		}
	}

	double PhysicPropData::get_kso(const double &Sw){
		switch(ksModel){
		case 1:
			return 1. - Sw;
		case 2:
			return pow(1. - Sw,2);
		case 3:
			return pow((1. - Sw - Swr)/(1. - Swr - Sor),2);
		default:
			throw Exception(__LINE__,__FILE__,"Unknown model for relative permeability.\n");
		}
	}

	void PhysicPropData::resetNonvisc(double &alpha_max){
		alpha_max = .0;
		for (int i=0; i<nnodes; i++){
			nonvisc.setValue(i,.0);
		}
	}
}
