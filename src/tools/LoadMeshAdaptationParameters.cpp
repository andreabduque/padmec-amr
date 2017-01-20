/*
 * LoadMeshAdaptationParameters.cpp
 *
 *  Created on: 16/02/2012
 *      Author: rogsoares
 */

#include "SimulatorParameters.h"
#include <sstream>

namespace PRS{

	void SimulatorParameters::loadMeshAdaptationParameters(){
		//if (!P_pid()) std::cout << "loading mesh adaptation parameters... ";
		ifstream fid;
		string str;
		// read madapt parameters
		string madapt(parametersPath + "/adaptation.dat");
		fid.open(madapt.c_str());
		if ( !fid.is_open() ) throw Exception(__LINE__,__FILE__,"adaptation.dat could not be opened or it doesn't exist.\n");

		// start reading file
		// ---------------------------------------------------------------------
		setPositionToRead(fid,"Use mesh adaptation:");
		fid >> str;
		if (!str.compare("yes")) setAdaptation();

		setPositionToRead(fid,"Refinement Strategy: {h-refinement, remeshing, rh-refinement}");
		fid >> str;
		if (!str.compare("h-refinement")){
			refstrategy = H_REFINEMENT;
		}
		else if (!str.compare("remeshing")){
			refstrategy = ADAPTIVE_REMESHING;
		}
		else if (!str.compare("rh-refinement")){
			refstrategy = RH_REFINEMENT;
		}
		
		setPositionToRead(fid,"Error tolerance for all mesh elements:");
		fid >> str;

		double error_1 = strtod(str.c_str(), 0);

		setPositionToRead(fid,"Error tolerance for all mesh elements excluding those on singularities regions:");
		fid >> str;
		double error_2 = strtod(str.c_str(), 0);

		if (error_1 < error_2){
			ofstream fout;
			fout.open("WARNING__Adaptation.txt");
			fout << "========================================================\n"
					" WARNING:\n"
					" Tolerance for all mesh elements SHOULD be GREATER than the"
					" tolerance for mesh elements excluding those on singularity"
					" regions.\n"
					"========================================================\n";
		}

		setPositionToRead(fid,"Maximum number of 2D element subdivisions (recommended: 4)");
		fid >> str;
		int num2D = atoi(str.c_str());

		setPositionToRead(fid,"Maximum number of 3D element subdivisions (recommended: 3)");
		fid >> str;
		int num3D = atoi(str.c_str());

		setPositionToRead(fid,"Maximum number of refinement steps:");
		fid >> str;
		int numSteps = atoi(str.c_str());

		setTol1(error_1);
		setTol2(error_2);
		setMax_2D(num2D);
		setMax_3D(num3D);
		setNumSubdivision_perStep(numSteps);


		// Interpolation Methods
		// ------------------------------------------------------------------------------------------------------------
		string::size_type pos;
		INTERPOLATION_OPTIONS intpmethod;	// initialize it as an unknown method
		setPositionToRead(fid,"Interpolation methods:");
		for (int i=0; i<8; i++) {
			getline(fid,str,'\n');
			pos = str.find("(x)",0);
			// if marked with 'x', search in 'str' for the interpolation method available
			if (pos != string::npos){
				if      (str == "(x) - h_REFINEMENT") intpmethod = h_REFINEMENT;
				else if (str == "(x) - LINEAR") intpmethod =  LINEAR;
				else if (str == "(x) - QUADRATIC") intpmethod =  QUADRATIC;
				else if (str == "(x) - ADAPTATIVE") intpmethod =  ADAPTATIVE;
				else if (str == "(x) - CONSERVATIVE") intpmethod =  CONSERVATIVE;
				else if (str == "(x) - PURE_INJECTION") intpmethod =  PURE_INJECTION;
				else if (str == "(x) - HALF_WEIGHTING") intpmethod =  HALF_WEIGHTING;
				else if (str == "(x) - FULL_WEIGHTING") intpmethod =  FULL_WEIGHTING;
			}
		}
		setInterpolationMethod(intpmethod);

		setPositionToRead(fid,"Parameters to decide if an element must be (un)refined or not.");
		fid >> str;
		remeshing_param1 = strtod(str.c_str(), 0);
		fid >> str;
		remeshing_param2 = strtod(str.c_str(), 0);

		if (remeshing_param1 > remeshing_param2){
			char msg[256]; sprintf(msg,"Remeshing adaptation: param1 (%.5f) must be less than param2 (%.5f)",remeshing_param1,remeshing_param2);
			throw Exception(__LINE__,__FILE__,msg);
		}

		setPositionToRead(fid,"h_min_allowed = h_min/n, where h_min is the minimum element height of the initial mesh.");
		fid >> str; //cout << str << endl;
		remeshing_param3 = strtod(str.c_str(), 0);

		//if (!P_pid()) std::cout << "done.\n";
	}
}

