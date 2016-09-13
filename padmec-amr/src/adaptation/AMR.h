/*
 * ADAPTATIVE MESH REFINEMENT CLASS
 *
 * AMR.h
 *
 *  Created on: 24/01/2013
 *      Author: Prof. Rog�rio Soares
 */

#ifndef AMR_H_
#define AMR_H_

#include "ErrorAnalysis.h"
#include "exportVTK.h"


class AMR{
public:
	AMR(){}
	virtual ~AMR(){}

	virtual void run(pMesh theMesh, std::list<pEntity>&, std::set<pEntity>&)=0;
};


#endif /* AMR_H_ */
