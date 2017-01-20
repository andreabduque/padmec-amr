/*
 * OilProductionManagement.cpp
 *
 *  Created on: 14/01/2009
 *      Author: rogerio
 */

#include "OilProductionManagement.h"

namespace PRS{

	OilProductionManagement::OilProductionManagement(){}

	OilProductionManagement::OilProductionManagement(string path, double iov,double TotalInjectionFlowRate, bool restart){

		char tmp[256];
		sprintf(tmp,"%s_oil-production-%d.csv",path.c_str(),P_size());
		string fname(tmp);

		IOV = iov;
		TIFR = TotalInjectionFlowRate;

		if (restart){
			fid.open(fname.c_str(),ios_base::app);
		}
		else{
			fid.open(fname.c_str());
			fid << "PVI time-step cumulative-time Recovered-Oil Cumulative-Oil Num.Time-steps" << std::endl;
		}

		if (!fid.is_open()) {
			throw Exception(__LINE__,__FILE__,"Oil production file could not be opened or it does not exist.\n");
		}
	}

	OilProductionManagement::~OilProductionManagement(){
		fid << "CLOSE" << std::endl;
		fid.close();
	}

	void OilProductionManagement::printOilProduction(double timeStep,
			                                         double cml_time,
			                                         double total_SimTime,
			                                         double rec_oil,
			                                         double cml_oil,
			                                         int timestep_counter){
		fid << fixed << setprecision(2) << (double)(cml_time/total_SimTime) << " " << setprecision(7)
			<< timeStep << " " << cml_time << " " << rec_oil << " " <<  cml_oil/IOV << " " << timestep_counter << endl;
	}
}
