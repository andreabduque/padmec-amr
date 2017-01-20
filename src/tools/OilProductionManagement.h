/*
 * OilProductionManagement.h
 *
 *  Created on: 14/01/2009
 *      Author: rogerio
 */

#ifndef OILPRODUCTIONMANAGEMENT_H_
#define OILPRODUCTIONMANAGEMENT_H_

#include "auxiliar.h"

namespace PRS{

	class OilProductionManagement{
	public:

		OilProductionManagement();
		OilProductionManagement(string, double, double,bool);
		~OilProductionManagement();

		// print a new oil production value frequency controlled by getPrintStep
		void printOilProduction(double timeStep,double cml_time,double total_SimTime,double rec_oil,double cml_oil,int timestep_counter);

		double getInitialOilVolume() const { return IOV; }

	private:
		ofstream fid;			// output stream for oil production
		double IOV;				// Initial Oil Volume
		double TIFR;			// Total Injection Flow Rate
	};
}

#endif /* OILPRODUCTIONMANAGEMENT_H_ */
