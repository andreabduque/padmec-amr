/*
 * SimulatorParameters2.cpp
 *
 *  Created on: Jan 29, 2015
 *      Author: rogerio
 */

#include "SimulatorParameters.h"

namespace PRS{

	void SimulatorParameters::defineExactSolution(){
		bc_external_definition = true;

		// put all command line strings into one string object
		string casename;
		for(int i=0; i<__argc; i++){
			casename.append(__argv[i]);
		}

		/*
		 * Run benchmark case using the exact solution as output.
		 */
		exact_sol_exist = ( casename.find("print_exact_solution") != std::string::npos )?true:false;

		cout << "Print exact solution is " << exact_sol_exist << endl;

		if ( casename.find("benchmark3d_case1") != std::string::npos ){
			exact_solution = Benchmark3D_case1__ES;
			ss_term = Benchmark3D_case1__SST;
			case_problem = CASE_1;
		}
		else if ( casename.find("benchmark3d_case5") != std::string::npos ){
			exact_solution = Benchmark3D_case5__ES;
			ss_term = Benchmark3D_case5__SST;
			case_problem = CASE_5;
		}
		else{
			bc_external_definition = false;
		}
	}
}
