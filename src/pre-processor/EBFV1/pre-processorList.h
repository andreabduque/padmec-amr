/*
 * pre-processorList.h
 *
 *  Created on: 28/12/2008
 *      Author: rogerio
 */

#ifndef PREPROCESSORLIST_H_
#define PREPROCESSORLIST_H_

#include "auxiliar.h"

namespace PRS{
	/**
	 *	Every numeric formulation will require a specific way to load the preporcessor
	 *	data and their respective load function must be defined here. All functions
	 *	prototypes must obey the following syntax:
	 *
	 *					void funcName(pMesh,void*,string);
	 *	where:
	 *		pMesh	- is the FMDB pointer function
	 *		void	- is a pointer to a specific object usually GeometricCoefficientData
	 *		string	- is the file name containing pre-processed data.
	 *
	 *	After that, user must edit load_preprocessorData() function in file SimulatorParameters.cpp
	 *	and include the new preprocessor function defining a new flag.
	 */

	void load_EBFV1_preprocessorData(pMesh, void*, string, int&);
}
#endif /* PREPROCESSORLIST_H_ */
