/*
 * EBFV1_AssemblyMatVec.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: rogerio
 */

#include "EBFV/EBFV1_elliptic.h"

namespace PRS{

	// NOTE: mapNodesOnWells should not be available in the way it appears here.
	// Some function should be implemented inside SimulatorParameters
	// to provide only the data required.
	void EBFV1_elliptic::wells_RHS_Assembly__Wells(pMesh mesh, Vec &RHS){

		if (!pSimPar->SimulationHas_BC_ExternalDefinition()){
			wells_1(mesh,RHS);
		}
		else{
			switch (pSimPar->getCaseProblem()){
			case CASE_1:{
				double x,y,z,volume,coord[3];
				int node_ID, row, i;

				int nrows = matvec_struct->nrows;	// number of free nodes
				int *rows = matvec_struct->rows;	// free nodes indices
				/*
				 *  1 - Take only free nodes
				 *  2 - This a specific problem with only two sub-domains flagged as 3300 and 3301.
				 *  3 - A node located on boundary domains has two control-volumes, each one associated to a sub-domain
				 *      Source/sink term must be applied to both
				 *  4 - When a node's volume is entirely located in a sub-domain, vol1 or vol2 will be equal to zero.
				 */
				for (i=0; i<nrows; i++){
					pGCData->getVolume(0,rows[i],volume);
					pGCData->getCoordinates(rows[i],coord);
					x = coord[0]; y = coord[1]; z = coord[2];
					node_ID = pMData->get_AppToPETSc_Ordering(rows[i]+1);
					row = pMData->FPArray(node_ID-1);
					VecSetValue(RHS,row,volume*pSimPar->ss_term(x,y,z),ADD_VALUES);
				}
			}
			break;
			case CASE_5:{
				double x,y,z,volume,coord[3];
				int node_ID, row, i;

				int nrows = matvec_struct->nrows;	// number of free nodes
				int *rows = matvec_struct->rows;	// free nodes indices
				/*
				 *  1 - Take only free nodes
				 *  2 - This a specific problem with only two sub-domains flagged as 3300 and 3301.
				 *  3 - A node located on boundary domains has two control-volumes, each one associated to a sub-domain
				 *      Source/sink term must be applied to both
				 *  4 - When a node's volume is entirely located in a sub-domain, vol1 or vol2 will be equal to zero.
				 */
				for (i=0; i<nrows; i++){
					cout << "i = " << i << "\tnrows: " << nrows << " ";
					pGCData->getVolume(rows[i],volume); cout << __LINE__ << " ";
					pGCData->getCoordinates(rows[i],coord); cout << __LINE__ << " ";
					x = coord[0]; y = coord[1]; z = coord[2]; cout << __LINE__ << " ";
					node_ID = pMData->get_AppToPETSc_Ordering(rows[i]+1); cout << __LINE__ << " ";
					row = pMData->FPArray(node_ID-1); cout << "row: " << row << "\t" << __LINE__ << " ";
					VecSetValue(RHS,row,volume*pSimPar->ss_term(x,y,z),ADD_VALUES); cout << __LINE__ << "\n";
				}
			}
			break;
			default:
				throw Exception(__LINE__,__FILE__,"Unknown benchmark case. Exiting....");
			}
		}
	}

	int EBFV1_elliptic::wells_1(pMesh mesh, Vec &RHS){
		int node_ID, row;
		double Vt, Vi, Qi, Qt;

		if (!pSimPar->mapNodesOnWell.size()){
			return 0;
			throw Exception(__LINE__,__FILE__,"No wells found!");
		}
		// for each well flag
		map<int,set<int> >::iterator mit = pSimPar->mapNodesOnWell.begin();
		for (; mit!=pSimPar->mapNodesOnWell.end(); mit++){
			int well_flag = mit->first;
			// source/sink term
			Qt = pSimPar->getFlowrateValue(well_flag);
			if ( fabs(Qt)<=1e-7 ){
				//throw Exception(__LINE__,__FILE__,"Flow rate NULL!");
			}

			// get all flagged node IDs for that well
			if (!mit->second.size()){
				throw Exception(__LINE__,__FILE__,"No wells found!");
			}
			SIter sit = mit->second.begin();
			for (; sit!=mit->second.end(); sit++){
				Vt = pSimPar->getWellVolume(well_flag);
	#ifdef _SEEKFORBUGS_
				if ( Vt<1e-12 ){
					char msg[256]; sprintf(msg,"Well with null volume V = %.6f. Vertex (%d)",Vt,node_ID);
					throw Exception(__LINE__,__FILE__,msg);
				}
	#endif //_SEEKFORBUGS_

				Vi = .0;
				node_ID = *sit;
				for (SIter_const dom=pSimPar->setDomain_begin(); dom!=pSimPar->setDomain_end(); dom++){
					pVertex node = (mEntity*)mesh->getVertex( node_ID );
					Vi += pGCData->getVolume(node,*dom);
				}

				// for node i, Q is a fraction of total well flow rate
				Qi = Qt*(Vi/Vt);

				// FPArray ('F'ree 'P'rescribed array) maps node id: node_ID -> row
				// row: position in Petsc GlobalMatrix/RHSVBector where node must be assembled
				node_ID = pMData->get_AppToPETSc_Ordering(node_ID);
				row = pMData->FPArray(node_ID-1);

				// Do not include well flux on nodes with prescribed pressure
				if (pSimPar->isNodeFree(well_flag)){
					VecSetValue(RHS,row,Qi,ADD_VALUES);
				}
			}
		}
		return 0;
	}


	/*
	 * This is a way to use the simulator to evaluate elliptic equation without screw-up
	 * the input data procedure. Treating source/sink terms:
	 */
	int EBFV1_elliptic::wells_2(pMesh mesh, Vec &RHS){
	#ifdef CRUMPTON_EXAMPLE
		double x,y,coord[3];
		double srcsnk, vol1, vol2, f_xy1, f_xy2;

		int nrows = matvec_struct->nrows;	// number of free nodes
		int *rows = matvec_struct->rows;	// free nodes indices
		/*
		 *  1 - Take only free nodes
		 *  2 - This a specific problem with only two sub-domains flagged as 3300 and 3301.
		 *  3 - A node located on boundary domains has two control-volumes, each one associated to a sub-domain
		 *      Source/sink term must be applied to both
		 *  4 - When a node's volume is entirely located in a sub-domain, vol1 or vol2 will be equal to zero.
		 */
		for (int i=0; i<nrows; i++){
			pVertex vertex = mesh->getVertex(rows[i]+1);
			if (vertex){
				V_coord(vertex,coord); x = coord[0]; y = coord[1];
				vol1 = pGCData->getVolume(vertex,3300); 			// it supposes coord_x <= 0
				vol2 = pGCData->getVolume(vertex,3301); 			// it supposes coord_x  > 0
				f_xy1 = -(2.*sin(y) + cos(y))*ALPHA*x - sin(y);		// source/sink for coord_x <= 0
				f_xy2 = 2.*ALPHA*exp(x)*cos(y);						// source/sink for coord_x  > 0
				srcsnk = vol1*f_xy1 + vol2*f_xy2;					// source/sink term
				node_ID = pMData->get_AppToPETSc_Ordering(rows[i]+1);
				row = pMData->FPArray(node_ID-1);
				VecSetValue(RHS,row,-srcsnk,ADD_VALUES);
			}
		}
	#endif
		return 0;
	}
}
