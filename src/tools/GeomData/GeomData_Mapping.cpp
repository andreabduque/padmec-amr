#include "GeomData.h"

namespace PRS{
	void GeomData::mappingNodesIds(pMesh theMesh){
		int ndom = getNumDomains();
		const int* domainList = getDomainList();
		std::map<int,int> mapIDtoIndex, mapBdryIDtoIndex, mapIDtoIndex_global;
		std::set<int> setID, setBdryID;
		int i, j, k, id0, id1, id;
		pEntity node, edge, face, tetra;

		i = 0;
		// set a sequential numbering for each vertex ID: from 0 to n-1, where n is the number of vertices
		alloc_INT_vector(__LINE__,__FILE__,pNodeID,M_numVertices(theMesh));
		VIter vit = M_vertexIter(theMesh);
		while ( (node=VIter_next(vit))){
			id = EN_id(node);
			pNodeID[i] = id;
			mapIDtoIndex_global[id] = i;
			i++;
		}
		VIter_delete(vit);

		pEntity v1,v2,v3;
		bool extElem = true;// PERFORME ONLY ONCE

		// STEP1: get all vertices ID from domain "dom" and its boundary.
		// All IDs are stored into setID and setBdryID set containers.
		for (k=0; k<ndom; k++){
			int dom = domainList[k];

			if (dim==2){
				FIter fit = M_faceIter(theMesh);
				while ( (face=FIter_next(fit) )){
					int faceflag = getFaceFlag(face);
					if (faceflag==domainList[k]){
						for(i=0; i<3; i++){
							setID.insert( EN_id(face->get(0,i)) );				// nodes domain
						}
						for(i=0; i<3; i++){
							if (getVertexFlag(face->get(0,i))!=faceflag){
								setBdryID.insert(EN_id(face->get(0,i)));		// boundary nodes
							}
						}
					}
				}
				FIter_delete(fit);
			}
			else{
				RIter rit = M_regionIter(theMesh);
				while ( (tetra = RIter_next(rit) )){
					int tetraflag = getTetraFlag(tetra);
					if (tetraflag==domainList[k]){
						for(i=0; i<4; i++){
							setID.insert( EN_id(tetra->get(0,i)) );				// nodes domain
						}
						for(i=0; i<4; i++){
							if (getVertexFlag(tetra->get(0,i))!=tetraflag){
								setBdryID.insert(EN_id(tetra->get(0,i)));		// boundary nodes
							}
						}
					}
				}
				RIter_delete(rit);
			}

			// get external boundary edges/faces only: It's for saturation gradient calculation
			if (extElem){
				i = 0;
				if (dim==2){
					EIter eit = M_edgeIter(theMesh);
					while ( (edge=EIter_next(eit)) ){
						/* ONLY EXTERNAL EDGES */
						if (E_numFaces(edge)==1){
							v1 = edge->get(0,0);
							v2 = edge->get(0,1);
							external_bdry_elem[0].setValue(i,0,mapIDtoIndex_global[EN_id(v1)]);
							external_bdry_elem[0].setValue(i,1,mapIDtoIndex_global[EN_id(v2)]);
							external_bdry_elem[0].setValue(i,2,GEN_tag(v1->getClassification()));
							external_bdry_elem[0].setValue(i,3,GEN_tag(v2->getClassification()));
							i++;
						}
					}
					EIter_delete(eit);
				}
				else{
					FIter fit = M_faceIter(theMesh);
					while ( (face = FIter_next(fit)) ){
						/* ONLY EXTERNAL FACES */
						if ( F_numRegions(face)==1 ){
							v1 = face->get(0,0);
							v2 = face->get(0,1);
							v3 = face->get(0,2);
							external_bdry_elem[0].setValue(i,0,mapIDtoIndex_global[EN_id(v1)]);
							external_bdry_elem[0].setValue(i,1,mapIDtoIndex_global[EN_id(v2)]);
							external_bdry_elem[0].setValue(i,2,mapIDtoIndex_global[EN_id(v3)]);
							external_bdry_elem[0].setValue(i,3,GEN_tag(v1->getClassification()));
							external_bdry_elem[0].setValue(i,4,GEN_tag(v2->getClassification()));
							external_bdry_elem[0].setValue(i,5,GEN_tag(v3->getClassification()));
							i++;
						}
					}
					FIter_delete(fit);
				}
				extElem = false;
			}

			// STEP2: map IDs from step1 like this:
			//		id0 = 0;
			//		id1 = 1;
			//		id2 = 2;
			//		...
			//		idn = n;

			if ((int)setID.size() > this->getNumNodesPerDomain(k)){
				char msg[256]; sprintf(msg,"NUmber of elements collected [%d] is greater than the max [%d]",(int)setID.size(),this->getNumNodesPerDomain(k));
				throw Exception(__LINE__,__FILE__,msg);
			}

			// domain
			i = 0;
			pVertex node;
			for (std::set<int>::iterator iter=setID.begin(); iter!= setID.end(); iter++){
				int id = *iter;			// vertex ID
				mapIDtoIndex[id] = i;	// gives a sequential numbering (0,1,2,...) for all vertices ID for domain k
				ID[k].setValue(i,id);
				node = theMesh->getVertex(id);
				volume[k].setValue(i,getVolume(node,dom) );
				nodes[k].setValue(i,mapIDtoIndex_global[id]);
				i++;
			}
			setID.clear();

			// boundary
			i = 0;
			for (std::set<int>::iterator iter=setBdryID.begin(); iter!= setBdryID.end(); iter++){
				int id = *iter;					// boundary vertex ID
				mapBdryIDtoIndex[id] = i;		// gives a sequential numbering (0,1,2,...) for all boundary vertices ID for domain k
				ID_bdry[k].setValue(i,id);
				node = theMesh->getVertex(id);
				volume_bdry[k].setValue(i,getVolume(node,dom) );
				i++;
			}
			setBdryID.clear();

			// STEP3: use mapped ID to initialize edge struct, where:
			//		edge_id0 = index_mapped
			//		edge_id1 = index_mapped

			// ALL EDGES
			i = 0;
			EIter eit = M_edgeIter(theMesh);
			while ( (edge=EIter_next(eit)) ){
				if ( edgeBelongToDomain(edge,domainList[k]) ){
					v1 = edge->get(0,0);
					v2 = edge->get(0,1);
					id0 = EN_id( v1 );
					id1 = EN_id( v2 );
					edges[k].setValue(i,0,mapIDtoIndex[id0]);					// index number for vertex ID for domain k
					edges[k].setValue(i,1,mapIDtoIndex[id1]); 					// index number for vertex ID for domain k
					edges[k].setValue(i,2,mapIDtoIndex_global[id0]);			// global index number for vertex ID
					edges[k].setValue(i,3,mapIDtoIndex_global[id1]);			// global index number for vertex ID
					edges[k].setValue(i,4,GEN_tag(v1->getClassification()));
					edges[k].setValue(i,5,GEN_tag(v2->getClassification()));
					i++;
				}
			}

			if (!i){
				throw Exception(__LINE__,__FILE__,"Any edge detected!");
			}

			if (dim==2){
				i = 0;
				FIter fit = M_faceIter(theMesh);
				while ( (face=FIter_next(fit)) ){
					int flag = GEN_tag(face->getClassification());
					if ( flag==dom ){
						for (j=0; j<3; j++){
							v1 = face->get(0,j);
							id0 = EN_id( v1 );
							elem[k].setValue(i,j,mapIDtoIndex[id0]);					// index number for vertex ID for domain k
							elem[k].setValue(i,j+3,mapIDtoIndex_global[id0]);			// global index number for vertex ID
						}
						i++;
					}
				}
			}
			else{
				i = 0;
				RIter rit = M_regionIter(theMesh);
				while ( (tetra=RIter_next(rit)) ){
					int flag = GEN_tag(tetra->getClassification());
					if ( flag==dom ){
						for (j=0; j<4; j++){
							v1 = tetra->get(0,j);
							id0 = EN_id( v1 );
							elem[k].setValue(i,j,mapIDtoIndex[id0]);					// index number for vertex ID for domain k
							elem[k].setValue(i,j+4,mapIDtoIndex_global[id0]);			// global index number for vertex ID
						}
						i++;
					}
				}
			}

			i = 0;
			bool bdryElem_detected = false;
			if (dim==2){
				// BOUNDARY EDGES PER DOMAIN
				eit = M_edgeIter(theMesh);
				while ( (edge=EIter_next(eit)) ){
					if ( edgeBelongToDomain(edge,domainList[k]) ){
						if (belongsToBoundary(edge)){
							bdryElem_detected = true;
							for (j=0;j<2;j++){
								id = EN_id( edge->get(0,j) );
								edges_bdry[k].setValue(i,j,mapBdryIDtoIndex[id]);		// index number for boundary vertex ID for domain k
								edges_bdry[k].setValue(i,j+2,mapIDtoIndex[id]); 		// index number for vertex ID for domain k
								edges_bdry[k].setValue(i,j+4,mapIDtoIndex_global[id]);	// global index number for vertex ID
							}
							i++;
						}
					}
				}
			}
			else{
				// BOUNDARY FACES PER DOMAIN
				FIter fit = M_faceIter(theMesh);
				while ( (face = FIter_next(fit)) ){
					if ( faceBelongToDomain(face,domainList[k]) ){
						if ( belongsToBoundary(face) ){
							bdryElem_detected = true;
							for (j=0;j<3;j++){
								id = EN_id( face->get(0,j) );
								faces_bdry[k].setValue(i,j,mapBdryIDtoIndex[id]);		// index number for boundary vertex ID for domain k
								faces_bdry[k].setValue(i,j+3,mapIDtoIndex[id]);			// index number for vertex ID for domain k
								faces_bdry[k].setValue(i,j+6,mapIDtoIndex_global[id]);	// global index number for vertex ID
							}
							i++;
						}
					}
				}
				FIter_delete(fit);
				if (!bdryElem_detected){
					throw Exception(__LINE__,__FILE__,"Any boundary element detected! You must use Physical command in .geo Gmsh to define them.");
				}
			}
			mapBdryIDtoIndex.clear();
			mapIDtoIndex.clear();
		}
		mapIDtoIndex_global.clear();
	}
}
