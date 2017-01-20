/*******************************************************************************
							UNIVERSIDADE FEDERAL DE PERNAMBUCO

Author:		Rogerio Soares da Silva
Date:		Dec-2007
Licence:	This program is free software; you can redistribute it and/or modify it.
Warranty:	None! Use it at your own risk!
Visit:		www.padmec.org

*******************************************************************************/

#ifndef INCLUDES_H
#define INCLUDES_H

// C/C++ includes
// =============================================================================
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <valarray>
#include <list>
#include <numeric>
#include <algorithm>
#include <new>

using namespace std;

// FMDB includes
// =============================================================================
#include "mMesh.h"
#include "mEntity.h"
#include "mVertex.h"
#include "mEdge.h"
#include "mFace.h"
#include "mTet.h"
#include "AOMD.h"
#include "AOMD_Internals.h"
#include "pmModel.h"
#include "AOMD_cint.h"
#include "AOMD_OwnerManager.h"
#include "mMirrorEntity.h"
#include "mAOMD.h"
#include "ParUtil.h"
#include "mExchangeData.h"
#include "mBuildAdj.h"
using namespace AOMD;

// Petsc includes
// =============================================================================
#include "petscksp.h"
#include "petscao.h"

// Project include
// =============================================================================
#include "Exception.h"

#endif
