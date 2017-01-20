#include "GeomData.h"

namespace PRS{

	void GeomData::dataTransfer(pMesh theMesh){
		transferCijData(theMesh);
		transferDijData(theMesh);
		transferVolData(theMesh);
		transferMesh(theMesh);
	}

	void GeomData::transferCijData(pMesh theMesh){
		dblarray Cij(3,.0);
		double cij[3], norm;
		int ndom = getNumDomains();
		for (int i=0; i<ndom; i++){
			int dom = domainList[i];
			int row = 0;
			EIter eit = M_edgeIter(theMesh);
			while ( pEdge edge = EIter_next(eit) ){
				if ( edgeBelongToDomain(edge,dom) ){
					getCij(edge,dom,Cij);
					norm = getCij_norm(edge,dom);
					cij[0] = Cij[0];
					cij[1] = Cij[1];
					cij[2] = Cij[2];
					setCij(i,row,cij);
					setCij_norm(i,row,norm);
					row++;
				}
			}
			EIter_delete(eit);
		}
	}

	void GeomData::transferDijData(pMesh theMesh){
		dblarray Dij(3,.0);
		double dij[3];
		int ndom = getNumDomains();
		for (int i=0; i<ndom; i++){
			int dom = domainList[i];
			int row = 0;

			if (dim==2){
				EIter eit = M_edgeIter(theMesh);
				while ( pEdge edge = EIter_next(eit) ){
					if ( edgeBelongToDomain(edge,dom) ){
						if (belongsToBoundary(edge)){
							getDij(edge,dom,Dij);
							dij[0] = Dij[0];
							dij[1] = Dij[1];
							dij[2] = Dij[2];
							setDij(i,row,dij);
							row++;
						}
					}
				}
				EIter_delete(eit);
			}
			else{
				FIter fit = M_faceIter(theMesh);
				while ( pFace face = FIter_next(fit) ){
					if ( faceBelongToDomain(face,dom) ){
						if (belongsToBoundary(face)){
							getDij(face,dom,dij);
							setDij(i,row,dij);
							row++;
						}
					}
				}
				FIter_delete(fit);
			}
		}
	}

	void GeomData::transferVolData(pMesh theMesh){
		int idx = 0;
		pEntity node;
		double volume;
		int ndom = getNumDomains();
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			volume = .0;
			for (int i=0; i<ndom; i++){
				volume += getVolume(node,domainList[i]);
			}
			setVolume(idx,volume);
			idx++;
		}
		VIter_delete(vit);
	}

	void GeomData::transferMesh(pMesh theMesh){
		int ID,i;
		pEntity node, face, tetra;
		int row = 0;
		double coords[3];
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			V_coord(node,coords);
			for(i=0; i<3; i++){
				pCoords->setValue(row,i,coords[i]);
			}
			row++;
		}
		VIter_delete(vit);

		row = 0;
		if (theMesh->getDim()==2){
			FIter fit = M_faceIter(theMesh);
			while ( (face = FIter_next(fit)) ){
				for(i=0; i<3; i++){
					node = (pEntity)face->get(0,i);
					ID = EN_id(node);
					pConnectivities->setValue(row,i,ID-1);	// ID-1: VTK purposes
				}
				row++;
			}
			FIter_delete(fit);
		}
		else{
			RIter rit = M_regionIter(theMesh);
			while ( (tetra = RIter_next(rit)) ){
				for(i=0; i<4; i++){
					node = (pEntity)tetra->get(0,i);
					ID = EN_id(node);
					pConnectivities->setValue(row,i,ID-1);	// ID-1: VTK purposes
				}
				row++;
			}
			RIter_delete(rit);
		}
	}
}
