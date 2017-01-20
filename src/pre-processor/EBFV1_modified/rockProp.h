/*
 * rockProp.h
 *
 *  Created on: Oct 31, 2014
 *      Author: rogerio
 */

#ifndef ROCKPROP_H_
#define ROCKPROP_H_

#include "auxiliar.h"
//This file is intent to define functions which represent rock properties such as:
//
//		1) absolute permeability tensor
//		2) porosity

#include "auxiliar.h"

// RPFP: Rock property function pointer
typedef void (*RockPropFuncPointer)(pEntity,int,double*);

void getRockPropertyFuncPointer(RockPropFuncPointer*);


/*
 * USERS MUST DEFINE THEIR FUNCTION BELOW
 */

void Homogeneous2D(pEntity e, double* K);

#endif /* ROCKPROP_H_ */
