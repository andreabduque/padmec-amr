/*
 * PhysicPropData2.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: rogerio
 */

#include "PhysicPropData.h"


namespace PRS{

	int PhysicPropData::getNumWells() const{
		return numWells;
	}

	int PhysicPropData::getNumWellNodes(int row) const{
		return numWellNodes.getValue(row,0);
	}

	void PhysicPropData::getFlowRate(int well_idx, int ith_node, int& row_idx, double& Qi) const{
		row_idx = idxWellNodeMat[well_idx][ith_node];
		Qi = flowrateWellNodeMat[well_idx][ith_node];
	}

	void PhysicPropData::calculateNodalFlowRate(pMesh theMesh, GeomData* pGCData){
		int flag,n,i,j;
		double Qt, Vi, TWV;
		std::set<pEntity> nodes_set;
		std::map<int,std::set<pEntity> > wellNodes_map;

		// search for flagged nodes: we are supposing that flags belonging to range 1 to 1999
		pEntity node, edge;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			flag = GEN_tag( node->getClassification() );
			if (flag>1 && flag<1999){
				nodes_set = wellNodes_map[flag];
				nodes_set.insert(node);
				wellNodes_map[flag] = nodes_set;
			}
		}
		VIter_delete(vit);

		EIter eit = M_edgeIter(theMesh);
		while ( (edge = EIter_next(eit)) ){
			flag = GEN_tag( edge->getClassification() );
			if (flag>1 && flag<1999){
				nodes_set = wellNodes_map[flag];
				nodes_set.insert((pEntity)edge->get(0,0));
				nodes_set.insert((pEntity)edge->get(0,1));
				wellNodes_map[flag] = nodes_set;
			}
		}
		EIter_delete(eit);

		nodes_set.clear();
		std::set<pEntity>::iterator set_iter;
		std::map<int,std::set<pEntity> >::iterator map_iter = wellNodes_map.begin();
		numWells = (int)wellNodes_map.size();		// defines number of wells
		flowrateWellNodeMat = new double*[numWells];
		idxWellNodeMat = new int*[numWells];
		for (i=0; i<numWells; i++){

			// first: Allocate memory
			flag = map_iter->first;
			n = (int)map_iter->second.size();			// size of node set belonging to i_th well
			flowrateWellNodeMat[i] = new double[n];		// allocate memory
			idxWellNodeMat[i] = new int[n];

			// second: compute total well volume (TWV) (TWV = sum(Vi), Vi = volume of nodal control volume which belongs to well p)
			// loop over set of well p nodes
			TWV = .0;
			for (set_iter = map_iter->second.begin(); set_iter!=map_iter->second.end(); set_iter++){
				node = *set_iter;
				pGCData->getVolume_MEBFV(EN_id(node),Vi);
				TWV += Vi;
			}

			// third: compute weighted flow rate per node for well p.
			set_iter = map_iter->second.begin();
			getTotalWellFlowRate(flag,Qt);
			for (j=0; j<n; j++){
				node = *set_iter;
				pGCData->getVolume_MEBFV(EN_id(node),Vi);
				flowrateWellNodeMat[i][j] = (double)Qt*(Vi/TWV);
				idxWellNodeMat[i][j] = EN_id(node)-1;
				set_iter++;
			}
			map_iter++;
		}
		wellNodes_map.clear();
	}

	void PhysicPropData::getTotalWellFlowRate(int flag, double& Qt){
		Qt = totalFlowRate_map[flag];
	}
}
