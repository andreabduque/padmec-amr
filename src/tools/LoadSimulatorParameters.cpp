/*
 * LoadSimulationParameters.cpp
 *
 *  Created on: 19/01/2011
 *      Author: rogerio
 */

#include "SimulatorParameters.h"
#include "EBFV1__pre-processors.h"
#include "EBFV1_modified.h"
#include <sstream>

namespace PRS{

	void SimulatorParameters::inputDataFromFiles(GeomData *pGCData){
		// define how the dirichlet and neumann values will be evaluated (from numeric.dat file or from Boundary_conditions.h)
		defineExactSolution();

		getParametersDirectoryPath();
		loadNumericParameters(pGCData);
		load_preprocessorData(pGCData);
		loadPhysicalParameters(pGCData->getMeshDim());
		loadSolverParameters();
		loadHighOrderParameter();
		loadMeshAdaptationParameters();
	}

	void SimulatorParameters::getParametersDirectoryPath(){

		// read numeric parameters
		ifstream fid;
		fid.open("simulation-parameters/filenameParameters.dat");
		if ( !fid.is_open() ){
			throw Exception(__LINE__,__FILE__,"filenameParameters.dat could not be opened or it does not exist\n");
		}

		setPositionToRead(fid,"path:");
		fid >> parametersPath;

		// check if numeric.dat physical.dat slope_limiters.dat and solve.dat files are presented
		string checkForNumeric(parametersPath + "/numeric.dat");
		string checkForPhysical(parametersPath + "physical.dat");
		string checkForSlope_limiters(parametersPath + "/slope_limiters.dat");
		string checkForSolve(parametersPath + "/solve.dat");

		ifstream fid1, fid2, fid3, fid4;
		fid1.open(checkForNumeric.c_str());
		fid2.open(checkForPhysical.c_str());
		fid3.open(checkForSlope_limiters.c_str());
		fid4.open(checkForSolve.c_str());

		string msg1("could not be opened or it does not exist\n");
		char msg[256];

		if ( !fid1.is_open() ) sprintf(msg,"numeric.dat %s",msg1.c_str());
		if ( !fid2.is_open() ) sprintf(msg,"physical.dat %s",msg1.c_str());
		if ( !fid3.is_open() ) sprintf(msg,"slope-limiters.dat %s",msg1.c_str());
		if ( !fid4.is_open() ) sprintf(msg,"solve.dat %s",msg1.c_str());

		string failMsg(msg);

		if (!failMsg.length()){
			throw Exception(__LINE__,__FILE__,msg);
		}
		else{
			fid1.close();
			fid2.close();
			fid3.close();
			fid4.close();
		}
	}

	void SimulatorParameters::initialize(GeomData *pGCData, pMesh theMesh){
#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: SimulatorParameters::initialize\tIN\n";
#endif
		if (!M_numFaces(theMesh)){
			throw Exception(__LINE__,__FILE__,"Number of mesh elements null!");
		}

		if (!this->useRestart()){
			setStepOutputFile(0);
		}
		getDomains();
		getWells(theMesh,pGCData->getMeshDim());
		weightWellFlowRateByVolume(theMesh,pGCData);
		if (getEllipticSolver()==1){
			setInitialOilVolume(theMesh,pGCData);
		}

#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: SimulatorParameters::initialize\tOUT\n";
#endif
	}

	void SimulatorParameters::deallocateData(){
		mapNodesOnWell.clear();
		_mapNodesDomain.clear();
		delete[] numNodesDom; numNodesDom = 0;
	}

	void SimulatorParameters::load_preprocessorData(void *pData){
		int ndom;
		ifstream fid;
		fid.open(prepFilename().c_str());
		if ( !fid.is_open() ){
			char msg[256]; sprintf(msg,"File: %s \ncould not be opened or it doesn't exist. Verify if path is correct.\n",prepFilename().c_str());
			throw Exception(__LINE__,__FILE__,msg);
		}
		fid.close();

//		std::list<parametric_coords> parametric_list;
//		readmesh_parametric(theMesh,prepFilename().c_str(),parametric_list);
//		cout << "Size parametric list: " << parametric_list.size() << endl;
//		STOP();

		M_load(theMesh,prepFilename().c_str());		// load mesh using FMDB
		int dim = theMesh->getDim();

		if (dim<2 || dim>3){
			throw Exception(__LINE__,__FILE__,"Mesh dimension unknown.");
		}

		cout << "Eliptic solver: " << getEllipticSolver() << endl;
		switch ( getEllipticSolver() ){
		case 1:	//.............................................. EBFV1
			if (dim==2){
				EBFV1_preprocessor_2D(theMesh,pData,ndom);
			}
			else{
				EBFV1_preprocessor_3D(theMesh,pData,ndom);
			}
			break;
		case 2:	//.............................................. EBFV1 MODIFIED
			if (dim==3){
				throw Exception(__LINE__,__FILE__,"MEBFV1 method is available for 2-D domains only.");
			}
			EBFV1_modified_preprocessor_2D(theMesh,(GeomData*)pData);
			break;
		default:
			throw Exception(__LINE__,__FILE__,"Available methods: 1 - EBFV1 and 2 - MEBFV1.");
		}
	}

	void SimulatorParameters::loadNumericParameters(GeomData *pGCData){
		string str;
		ifstream fid;

		// read numeric parameters
		string numeric(parametersPath + "/numeric.dat");
		fid.open(numeric.c_str());
		if ( !fid.is_open() ){
			char msg[256]; sprintf(msg,"File: %s \ncould not be opened or it doesn't exist\n",numeric.c_str());
			throw Exception(__LINE__,__FILE__,msg);
		}


		// start reading file
		// -------------------------------------------------------------------------
		setPositionToRead(fid,"type method ID right below:");
		fid >> ellipticSolver;

		setPositionToRead(fid,"type method ID right below:");
		fid >> hyperbolicSolver;

		setPositionToRead(fid,"Use high order approximation:");
		fid >> str;

		if (!str.compare("yes")){
			setUseHOApproximation();
		}

		setPositionToRead(fid,"C.F.L (Courant-Friedilich-Lendroff):");
		fid >> _CFL;

		setPositionToRead(fid,"Discretization time on production well:");
		fid >> _timeStepWell;

		if (!str.compare("yes") && _CFL > .5)
			throw Exception(__LINE__,__FILE__,"High order approximation cannot use CFL condition number greater than 0.5!\n");

		setPositionToRead(fid,"mesh filename:");
		fid >> prepfName;

		setPositionToRead(fid,"export filename (to VTK format without extension .vtk):");
		fid >> expofName;

		string strYes;
		setPositionToRead(fid,"Start simulation from <restart> file?");
		getline(fid,strYes);

		setPositionToRead(fid,"Dirichlet (Specify: flags and respective pressure values):");
		mapFlag(fid,"dirichlet");

		setPositionToRead(fid,"Injection Wells: 10 - 50");
		getWells(fid);

		setPositionToRead(fid,"Production Wells: 51 - 90");
		getWells(fid);

		setPositionToRead(fid,"Specify water saturation in injection wells:");
		mapFlag(fid,"saturation");

		setPositionToRead(fid,"water saturation at t=0:");
		fid >> _Sw_initial;

		setPositionToRead(fid,"simulation-time:");
		fid >> _ST;

		setPositionToRead(fid,"dvtol:");
		fid >> _dvtol;

		setPositionToRead(fid,"Reservoir Dimensions: L - H:");
		fid >> reservoir_Length >> reservoir_Height;

		pGCData->setReservoirHeight(reservoir_Height);

		fid.close();

		// RESTART: check if simulation starts from the beginning or from where it stopped.
		string::size_type pos = strYes.find("yes",0);
		if (pos != string::npos){
			restart = true;

			ifstream rfid;
			// seek for the last VTK file generated
			for (int i=0; i<=100; i++){
				char strtmp[512]; sprintf(strtmp,"%s__0-of-1__step-%d.vtk",expofName.c_str(),i);
				rfid.open(strtmp);
				if (rfid.is_open()){
					rfid.close();
				}
				else{
					cout << "\nRestart:\n--------------------------------------------------------------------------------------\n";
					if (!i){
						throw Exception(__LINE__,__FILE__,"Any VTK file has been generated using the provided case file. Please, check carefully for mistypes\n");
					}
					lastpvi = i-1; // lastpvi=i is the VTK file that should be created by simulation was terminated. So the last is i-1.
					PVI_cumulative = double(lastpvi/100.0);
					char strtmp2[512]; sprintf(strtmp2,"%s__0-of-1__step-%d.vtk",expofName.c_str(),lastpvi);
					restartFilename = strtmp2;
					rfid.close();

					cout << "\nSimulation will be resumed from VTK file:\n" << restartFilename << endl;

					// set cumulative simulation time variable: cumTS = summation of all time steps (from the beginning to PVI i-1
					string strline;
					char strtmp4[512]; sprintf(strtmp4,"%s_simulation-monitor-%d.csv",expofName.c_str(),P_size());
					rfid.open(strtmp4);
					if (!rfid.is_open()){
						throw Exception(__LINE__,__FILE__,"Simulation monitor file could not be opened.\n");
					}
					break;
				}
			}
			this->setStepOutputFile(lastpvi+1);
			cumTS = PVI_cumulative*this->getSimTime();
			cout << "Simulation Time           : " << getSimTime();
			cout << "\nCumulative Simulation Time: " << cumTS << "  [" << setprecision(2) << double(100.*cumTS/getSimTime()) << "% concluded]" << endl;
			cout << "lastpvi= " << lastpvi << endl;
		}

		if (!P_pid()){
			cout << "SIMULATION INFO:\n";
			cout << "-------------------------------------------------------------------------------------------------------------------------";
			cout << "\nPath and file name:\n";
			cout << prepfName << "\n";
			cout << "Output path and file name:\n";
			cout << expofName << endl;
		}

	}

	void SimulatorParameters::loadPhysicalParameters(int dim){
		ifstream fid;
		string str;
		// read physical parameters
		string physical(parametersPath + "/physical.dat");
		fid.open(physical.c_str());
		if ( !fid.is_open() ) throw Exception(__LINE__,__FILE__,"physical.dat could not be opened or it doesn't exist.\n");

		// start reading file
		setPositionToRead(fid,"Oil density");
		fid >> str; _oil_density = strtod(str.c_str(), 0);

		setPositionToRead(fid,"Water density:");
		fid >> str; _water_density = strtod(str.c_str(), 0);

		setPositionToRead(fid,"Oil viscosity:");
		fid >> str; _oil_viscosity = strtod(str.c_str(), 0);

		setPositionToRead(fid,"Water viscosity:");
		fid >> str; _water_viscosity = strtod(str.c_str(), 0);

		setPositionToRead(fid,"# dom, K and phi:");
		readRockProperties(fid,dim);

		setPositionToRead(fid,"Residual Oil Satutarion");
		fid >> str; _Sor = strtod(str.c_str(), 0);

		setPositionToRead(fid,"Irreducible Water Satutarion");
		fid >> str; _Swr = strtod(str.c_str(), 0);

		setPositionToRead(fid,"type model ID right above:");
		fid >> _ksModel;

		//		printf("_ksModel = %d\n",_ksModel); throw 1;

		fid.close();
		getDomains();
	}

	void SimulatorParameters::loadSolverParameters(){
		ifstream fid;
		string str;
		// read solver parameters
		string solver(parametersPath + "/solver.dat");
		fid.open(solver.c_str());
		if ( !fid.is_open() ) throw Exception(__LINE__,__FILE__,"solver.dat could not be opened or it doesn't exist.\n");

		// start reading file
		setPositionToRead(fid,"absolute convergence tolerance");
		fid >> _abstol;

		setPositionToRead(fid,"relative convergence tolerance for KSP solver");
		fid >> _rtol;

		setPositionToRead(fid,"relative convergence tolerance for external iteration");
		fid >> _rtol2;

		setPositionToRead(fid,"convergence tolerance in terms of the norm of the change in the solution between steps");
		fid >> _dtol;

		setPositionToRead(fid,"maximum number of iterations");
		fid >> _maxit;

		setPositionToRead(fid,"Sets number of iterations at which GMRES, FGMRES and LGMRES restarts:");
		fid >> _Krylov_restart;


		// Load a 'database' to set a preconditioner chosen by user. <private>
		setPreconditionerDataBase();

		// Read solver.dat and check which preconditioner user chose.
		string::size_type pos;
		int numpc = (int)mapPC.size(); // number of preconditioners presented on data base.
		setPositionToRead(fid,"Define a preconditioner (Default: none)");
		for (int i=0; i<numpc; i++) {
			getline(fid,str,'\n');

			// search for chosen option
			pos = str.find("(x)",0);

			// if marked with 'x', search in 'str' for pre-conditioner available
			if (pos != string::npos){
				std::map<std::string,PCType>::iterator miter = mapPC.begin();
				while (miter != mapPC.end()){
					pos = str.find(miter->first,0);
					if (pos != string::npos){				// miter->first: string name
						pctype = miter->second;				// miter->second preconditioner
					}
					miter++;
				}
			}
		}

		// read pressure solver scheme for EBFV1 formulation
		str.clear();
		setPositionToRead(fid,"Use <defect-correction> scheme for EBFV1 pressure solver?: (default: matrix-free scheme)");
		getline(fid,str);
		pos = str.find("(x)",0);
		if (pos != string::npos) {
			setUseDefectCorrection();
			if (!P_pid()) printf("Using defect-correction scheme for EBFV1 pressure solver.\n");
		}

		fid.close();
	}

	void SimulatorParameters::loadHighOrderParameter(){
		ifstream fid;
		string str;
		// read solver parameters
		string slimiters(parametersPath + "/slope_limiters.dat");
		fid.open(slimiters.c_str());
		if ( !fid.is_open() ) throw Exception(__LINE__,__FILE__,"solve-limiters.dat could not be opened or it doesn't exist.\n");

		// start reading file
		const int HOM = 5;
		std::string theString[HOM];

		// Read from file all slope limiter function options
		setPositionToRead(fid,"Nodal sloper limiter functions (default: MUSCL):");
		for (int i=0; i<HOM; i++) getline(fid,theString[i],'\n');

		// Define a "data bank" for node slope limiters function
		mapSL_node["(x) - Superbee"] = node_Superbee;
		mapSL_node["(x) - Minmod"] = node_Minmod;
		mapSL_node["(x) - Van_Albada"] = node_Van_Albada;
		mapSL_node["(x) - Osher"] = node_Osher;
		mapSL_node["(x) - WoodField"] = node_WoodField;

		std::map<string,NSLF>::iterator mit;

		// define node SL default as MUSCL
		slf_method_Node = node_MUSCL;

		// Find which option user has chosen
		for (int i=0; i<HOM; i++){
			mit = mapSL_node.find(theString[i]);
			if ( mit != mapSL_node.end() ) slf_method_Node = mit->second;
			theString[i].clear();
		}

		// Read from file all slope limiter function options
		setPositionToRead(fid,"Edge sloper limiter functions (default: MUSCL):");
		for (int i=0; i<HOM; i++) getline(fid,theString[i],'\n');

		// Define a "data bank" for edge slope limiters function
		mapSL_edge["(x) - Superbee"] = SUPERBEE;
		mapSL_edge["(x) - Minmod"] = MINMOD;
		mapSL_edge["(x) - Van Albada"] = VAN_ALBADA;
		mapSL_edge["(x) - Van Albada"] = OSHER;
		mapSL_edge["(x) - WoodField"] = WOODFIELD;

		// define edge SL default as MUSCL
		slf_method_Edge = MUSCL;

		// Find which option user has chosen
		std::map<string,ESLF>::iterator mite;
		for (int i=0; i<HOM; i++){
			mite = mapSL_edge.find(theString[i]);
			if ( mite != mapSL_edge.end() ) slf_method_Edge = mite->second;
		}

		// Modified Taylor expansion series extrapolating nodal saturation value on volume control interfaces. Mark an 'x' (without quotes) inside
		// the parentheses to define the coefficient 'k' for the saturation extrapolation.

		// define koef as 0 <default>
		koef = .0;
		map_koef["(x) -  -1,   método de ponderação a montante de segunda ordem"] = -1.;
		map_koef["(x) - 1/3, método de ponderação a montante de terceira ordem"] = 1./3.;
		map_koef["(x) -   1,   método de diferenças centradas de três pontos"] = 1.;

		// Read from file the options
		setPositionToRead(fid,"Taylor expansion coefficient - (default: koef = 0, método de Fromm):");
		for (int i=0; i<4; i++) getline(fid,theString[i],'\n');

		// Find which option the use has chosen
		std::map<string,double>::iterator mit_koef;
		for (int i=0; i<HOM; i++){
			mit_koef = map_koef.find(theString[i]);
			if ( mit_koef != map_koef.end() ){
				koef = mit_koef->second;
			}
		}

		// Delta value for woodfield
		_WFdelta = 0.2;
		//setPositionToRead(fid,"Delta value for Woodfield:");
		//fid >> _WFdelta;
		fid.close();
	}

	void SimulatorParameters::setPositionToRead(ifstream &fid, string str2){
		int count = 0;
		string str1;
		do{
			getline(fid,str1);
			if (++count > 100){
				if (!P_pid()){
					cout << "\n\nWARNING: Could not find <" << str2
							<<">.\nProbably you are using an out of date input file."
							" Default value will be used instead.\n\n";
				}break;
			}
		}while( str1.compare(str2) );
	}

	void SimulatorParameters::readRockProperties(ifstream &fid, int dim){
		RockProperties *rockprop;
		double num;
		string str;
		istringstream stream1;
		std::list<double> K_tmp;

		// This is a way to use the simulator to evaluate elliptic equation without screw-up the input data procedure.
	#ifdef CRUMPTON_EXAMPLE
		RockProperties *rockprop1 = new RockProperties;
		RockProperties *rockprop2 = new RockProperties;

		rockprop1->K = new double[4];	// get permeability tensor
		rockprop1->K[0] = 1; rockprop1->K[1] = 0;
		rockprop1->K[2] = 0; rockprop1->K[3] = 1;
		mapRockProp[3300] = rockprop1;

		rockprop2->K = new double[4];	// get permeability tensor
		rockprop2->K[0] = (double)2*ALPHA; rockprop2->K[1] = (double)1*ALPHA;
		rockprop2->K[2] = (double)1*ALPHA; rockprop2->K[3] = (double)2*ALPHA;
		mapRockProp[3301] = rockprop2;

	#else
		while( true ){
			getline(fid,str);
			if (!str.compare("end")) break;
			rockprop = new RockProperties;
			const int dom = atoi(str.c_str());	// get domain flag
			rockprop->K = new double[dim*dim];	// get permeability tensor

			// read permeability tensor
			getline(fid,str);
			stream1.str(str);
			int count =0;
			// check if file data size match to mesh dimensions
			// 2-D -> K[2][2] or 3-D -> K[3][3]
			while( stream1 >> num ){
				K_tmp.push_back(num);
				count++;
			}
			if ( count != dim*dim ){
				throw Exception(__LINE__,__FILE__,"Mesh dimension and permeability tensor size do not match!\n");
			}

			std::copy(K_tmp.begin(),K_tmp.end(),rockprop->K);

			// get porosity
			getline(fid,str);
			rockprop->porosity = strtod(str.c_str(), 0);
			mapRockProp[dom] = rockprop;
			K_tmp.clear();
			stream1.clear();
		}
		if (mapRockProp.size()==0)
			throw Exception(__LINE__,__FILE__,"You must provide rock properties\n");
	#endif
	}

	void SimulatorParameters::mapFlag(ifstream &fid, string whatmap){
		string str;
		BdryConditionType *bct;

		while( 1 ){
			fid >> str;
			if (!str.compare("end")) break;
			int flag = atoi(str.c_str());
			fid >> str;
			double val = strtod(str.c_str(), 0);

			// for every flag number create a new BdryConditionType pointer and add inform type of condition and its value.
			bct = new BdryConditionType;
			bct->type = whatmap;
			bct->val = val;

			// two map container are used: one for saturation and other to dirichlet/neumann boundary conditions
			if (!whatmap.compare("saturation")){
				mapSaturation[flag] = bct;
			}
			else{
				mapBC[flag] = bct;
			}
		}
	}

	void SimulatorParameters::getWells(ifstream &fid){
		string str;
		while( true ){
			fid >> str;
			if (!str.compare("end")) break;
			int flag = atoi(str.c_str());

			fid >> str;
			double val = strtod(str.c_str(), 0);
			// associate a vector of size 2 to each well flag
			// data[0] = Q, total flow rate
			// data[1] = V, total well volume (to be calculated later)
			//			dblarray data(2);
			//			data[0] = val;
			//			MWells[flag] = data;

			well_info.isInjection = false;
			well_info.isProduction = true;
			well_info.flowRate = val;
			if (val>.0){
				well_info.isInjection = true;
				well_info.isProduction = false;
			}
			MWells[flag] = well_info;
		}
	}

	string SimulatorParameters::getFilename(string meshName){
		char s1[256], s4[256];
		generateFilename(meshName,s1);
		sprintf(s4,"%s_CFL%.2f-%d-of-%d",s1,CFL(),P_pid(),P_size());
		string s5(s4);
		return s5;
	}

	void SimulatorParameters::generateFilename(const string &meshFilename, char *filename){
		string str(meshFilename);
		string::size_type start = str.find_last_of("/");
		start = (start == string::npos)?0:start;
		string::size_type end = str.find_last_of(".");
		char buffer[256];
		memset( buffer, '\0', 256 );
		str.copy(buffer,end-start,start);
		sprintf(filename,"%s",buffer);
	}

	void SimulatorParameters::setPreconditionerDataBase(){
		mapPC["(x) JACOBI"]    = PCJACOBI;
		mapPC["(x) SOR"]       = PCSOR;
		mapPC["(x) LU"]        = PCLU;
		mapPC["(x) SHELL"]     = PCSHELL;
		mapPC["(x) BJACOBI"]   = PCBJACOBI;
		mapPC["(x) MG"]        = PCMG;
		mapPC["(x) EISENSTAT"] = PCEISENSTAT;
		mapPC["(x) ILU"]       = PCILU;
		mapPC["(x) ICC"]       = PCICC;
		mapPC["(x) ASM"]       = PCASM;
	}
}
