/*
 * SimParStrucs.h
 *
 *  Created on: 30/07/2010
 *      Author: rogerio
 */

#ifndef SIMPARSTRUCS_H_
#define SIMPARSTRUCS_H_

/*
 * Boundary conditions structure
 */
	struct BdryConditionType{
		string type;		// dirichlet or neumann
		double val;		// pointer to prescribed value
	};

	typedef map<int,BdryConditionType*> MapFlag;
	typedef MapFlag::iterator MapFlagIter;
	typedef set<int>::const_iterator SIter_const;
	typedef void (*pFunc_PrintVTK)(pMesh, void*, void*, void*, void*, string);

/*
 * Rock Properties structures
 */
	struct RockProperties{
		double porosity;
		double *K;
	};

	typedef map<int,RockProperties*> MapRockProperties;
	typedef MapRockProperties::iterator MIter_RockProperties;

	struct WellInfo{
		double flowRate;
		double wellVolume;
		bool isInjection;
		bool isProduction;
	};

	typedef map<int,WellInfo> MapWells;
	typedef MapWells::iterator MWIter;
	typedef MapWells::const_iterator MWCIter;

	/*
	 * Slope limiter functions parameters
	 */
	// Nodal sloper limiter functions:
	enum NSLF {node_Superbee, node_Minmod, node_MUSCL, node_Van_Albada, node_Osher, node_WoodField};

	// Edge sloper limiter functions:
	enum ESLF {SUPERBEE,MINMOD,MUSCL,VAN_ALBADA,OSHER,WOODFIELD};

	// fractional flux implementations (to avoid bad code for specific test cases):
	// Verma1: fw = Swn², where Swn = (Sw-Swr)/(1-Swr-Sor);
	enum FRACTIONALFLUX {Verma1};

#endif /* SIMPARSTRUCS_H_ */
