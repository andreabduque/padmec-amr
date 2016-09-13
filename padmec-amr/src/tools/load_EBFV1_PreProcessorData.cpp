/*
 * load_EBFV1_PreProcessorData.cpp
 *
 *  Created on: 29/12/2008
 *      Author: rogerio
 */

#include "GeomData.h"
#include "pre-processorList.h"

namespace PRS{

	// local functions prototypes
	// =============================================================================
	void loadNodalVolumes(ifstream&, pMesh, GeomData*, int *&setDomains, int &);
	void loadCij(ifstream&, pMesh, GeomData *, int *setDomains, int ndom);
	void loadDij(ifstream&, pMesh, GeomData *, int *setDomains);
	// =============================================================================

	void load_EBFV1_preprocessorData(pMesh theMesh, void* p, string fname, int &ndom){
		if (!P_pid()) cout << "loading pre-processor data begin.\n";
		char filename[256];
		sprintf(filename,"%s-%d-of-%d.dat",fname.c_str(),P_pid(),P_size());
		//cout << "pre-processor file: " << filename << endl;
		///PetscSynchronizedPrintf(PETSC_COMM_WORLD,"pre-processor file: %s\n",filename);
		//PetscSynchronizedFlush(PETSC_COMM_WORLD);
		MPI_Barrier(PETSC_COMM_WORLD);



		ifstream fid;
		fid.open(filename);

		if (!fid.is_open()){
			char str[256];
			sprintf(str,"File %s could not be opened or it does not exist\n",filename);
			throw Exception(__LINE__,__FILE__,str);
		}
		if (M_numFaces(theMesh)==0 && theMesh->getDim()==3)
			throw Exception(__LINE__,__FILE__,"Boundary faces are missing. Exinting...\n");

		//int ndom;
		int *setDomains; setDomains=0;
		GeomData *pGCData = (GeomData*)p;
		loadNodalVolumes(fid,theMesh,pGCData,setDomains,ndom);
		loadCij(fid,theMesh,pGCData,setDomains,ndom);
		loadDij(fid,theMesh,pGCData,setDomains);

		// Number of global edges:
		string line;
		getline(fid,line,'\n');
		getline(fid,line,'\n');
		getline(fid,line,'\n');
		fid >> line; pGCData->setNumGEdges( atoi(line.c_str()) );

		fid.close();

		char msg[256]; sprintf(msg,"\t[%d] - Vertices(%d), Edges(%d), Boundary Faces(%d)\n",
				P_pid(),M_numVertices(theMesh),pGCData->getNumGEdges(),M_numFaces(theMesh));

		if (!P_pid()) {
			std::cout << "Number of domains: " << ndom << std::endl;
			std::cout << msg;
			cout << "loading pre-processor data end.\n";
		}

		// calculate edge length
		int i,j,dim = 3;
		double elength, Lij[dim];
		double Icoord[3], Jcoord[3];
		double delta_x = 1e30; // infinity delta_x
		dblarray vec(3);
		EIter eit = M_edgeIter(theMesh);
		while (pEntity edge = EIter_next(eit)){
			E_getVerticesCoord(edge,Icoord,Jcoord);
			for (j=0; j<dim; j++) Lij[j] = .0;
			makeVector(Jcoord,Icoord,Lij);

			for (i=0; i<dim; i++) vec[i] = Lij[i];
//			pGCData->setEdgeVector(edge,vec);

			elength = .0;
			for (i=0; i<dim; i++) elength += vec[i]*vec[i];
			elength = sqrt(elength);

			if (elength == .0){
				char msg[256]; sprintf(msg,"Edge [%d %d] has null length!",EN_id(edge->get(0,0)),EN_id(edge->get(0,1)));
				throw Exception(__LINE__,__FILE__,msg);
			}

			pGCData->setEdgeLength(edge,elength);

			for (i=0; i<dim; i++) vec[i] /= elength;
			pGCData->setEdgeVec_Unitary(edge,vec);

			// get the smallest one
			if (delta_x > elength) delta_x = elength;

			#ifdef _SEEKFORBUGS_
				if (elength == .0) throw Exception(__LINE__,__FILE__,"Edge length is NULL!\n");
			#endif
		}
		EIter_delete(eit);//throw 1;

		// if parallel, get the smallest edge from all ranks and then broadcast it.
		pGCData->setSmallestEdgeLength(P_getMinDbl(delta_x));

		// calculate edge length
		int count;
		pMeshDataId dataAttached_id = MD_lookupMeshDataId("node-connectivities");
		eit = M_edgeIter(theMesh);
		while (pEntity edge = EIter_next(eit)){
			count = 0;
			EN_getDataInt(edge->get(0,0),dataAttached_id,&count);
			EN_attachDataInt(edge->get(0,0),dataAttached_id,++count);

			count = 0;
			EN_getDataInt(edge->get(0,1),dataAttached_id,&count);
			EN_attachDataInt(edge->get(0,1),dataAttached_id,++count);
		}
		EIter_delete(eit);

		int max_connect = 0;
		int min_connect = 1e+8;
		VIter vit = M_vertexIter(theMesh);
		while (pEntity node = VIter_next(vit)){
			EN_getDataInt(node,dataAttached_id,&count);
			if (count > max_connect) max_connect = count;
			if (count < min_connect) min_connect = count;
		}
		VIter_delete(vit);
	}

	void loadNodalVolumes(ifstream &fid, pMesh theMesh, GeomData *pGCData, int *&setDomains, int &ndom){
		if (!P_pid()) printf("\tloading Vol... ");
		int nrc, ID, i, iGrp;
		double x,y,z,vol,wvol;
		double vt = 0.0;
		string str, sdata[10];
		pEntity v;

		getline(fid,str,'\n');			//Mesh dimension:
		getline(fid,str,'\n');
		int dim = atoi(str.c_str()); 	//3
		pGCData->setMeshDim(dim);
		getline(fid,str,'\n');			//Number of domains:
		getline(fid,str,'\n');
		ndom = atoi(str.c_str());		//2
		setDomains = new int[ndom];
		getline(fid,str,'\n');			//Domains' flag list:
		for (i=0; i<ndom; i++){
			getline(fid,str,'\n');
			setDomains[i] = atoi(str.c_str());
		}
		getline(fid,str,'\n');			//For each node: ID - numRemoteCopies - x - y - z - nodalVolume - weightedvolume
		for (int i=0; i<ndom; i++){
			getline(fid,str,'\n');		//Nodes over domain: XXXX
			while (1){
				fid >> sdata[0];
				if ( !sdata[0].compare("end") ) break;
				ID = atoi(sdata[0].c_str());
				if (!ID) throw Exception(__LINE__,__FILE__,"ID zero\n");

				fid >> sdata[0] >> sdata[1] >> sdata[2] >> sdata[3] >> sdata[4] >> sdata[5] >> sdata[6];

				iGrp = atoi(sdata[0].c_str());
				nrc = atoi(sdata[1].c_str());
				x = strtod(sdata[2].c_str(), 0);
				y = strtod(sdata[3].c_str(), 0);
				z = strtod(sdata[4].c_str(), 0);

				/*
				 * For 2-D domains, multiply volume by reservoir height (H) for 2D/3D
				 * simulations (physics occurs only on 2-D but reservoir volume
				 * is considered)
				 */
				double H = (pGCData->getMeshDim()==2)?pGCData->getReservoirHeight():1.0;

				//cout << "pGCData->getReservoirHeight() " << pGCData->getReservoirHeight() << endl; STOP();

				// do not divide vol by number of remote copies!
				vol = H*strtod(sdata[5].c_str(), 0);
				if (vol == .0){
					char msg[256]; sprintf(msg,"Node %d has null volume!",ID);
					throw Exception(__LINE__,__FILE__,msg);
				}

//				printf("rank: %d nó %d vol: %f\n",P_pid(),ID,vol);
				vt += vol/((double)nrc+1.0);
				wvol = H*strtod(sdata[6].c_str(), 0);

#ifdef _SEEKFORBUGS_
    if (vol == .0) throw Exception(__LINE__,__FILE__,"Volume is NULL!\n");
#endif //_SEEKFORBUGS_

				v = (pEntity)theMesh->getVertex(ID);
				if (!v)	v = (mEntity*)theMesh->createVertex(ID,x,y,z,theMesh->getGEntity(iGrp,0));


				pGCData->setVolume(v,setDomains[i],vol);
				pGCData->setWeightedVolume(v,wvol);
				//pGCData->setNumRemoteCopies(v,nrc);

				///printf("dom = %d,  vol[%d] = %f\n",setDomains[i],ID,pGCData->getVolume(v,setDomains[i]));
			}
			getline(fid,str,'\n');		//Nodes over domain: 3300
		}
		if (!P_pid()) printf("done.\n");
		if (!M_numVertices(theMesh)) throw Exception(__LINE__,__FILE__,"Mesh has any nodes.\n");
		pGCData->setTotalReservoirVolume(vt);

		PetscPrintf(PETSC_COMM_WORLD,"Volume total: %f\n",pGCData->getReservoirVolume());
	}

	void loadCij(ifstream &fid, pMesh theMesh, GeomData *pGCData, int *setDomains, int ndom){
		if (!P_pid()) printf("\tloading Cij... ");
		string str[10];
		const int dim = 3;
		int nrc, iGrp, j, i = 0;
		int edgeIDs[2];
		dblarray Cij(dim);
		double max_edge=.0, avgLength=.0;

		bool hasRC=false;

		fid >> str[0];
		getline(fid,str[0],'\n');
		for (i=0; i<ndom; i++){
			getline(fid,str[0],'\n');
			while ( 1 ){
				fid>>str[0];
				if (str[0] == "end") break;
				iGrp = atoi(str[0].c_str());			// flag
				for (j=1; j<=dim+3; j++) fid >> str[j];
				nrc = atoi(str[1].c_str());				// number of remote copies
				edgeIDs[0] = atoi(str[2].c_str());		// ID0
				edgeIDs[1] = atoi(str[3].c_str());		// ID1

				if (nrc>0) hasRC=true;

				if (edgeIDs[0]>edgeIDs[1]) std::swap(edgeIDs[0],edgeIDs[1]);

				mEdge* edge =  theMesh->getEdge(theMesh->getVertex(edgeIDs[0]),theMesh->getVertex(edgeIDs[1]));
				if( !edge )
					edge = theMesh->createEdge(edgeIDs[0],edgeIDs[1],theMesh->getGEntity(iGrp,1));

				/*
				 * For 2-D domains, multiply Cij by reservoir height (H) for 2D/3D
				 * simulations (physics occurs only on 2-D but reservoir volume
				 * is considered)
				 */

				double H = (pGCData->getMeshDim()==2)?pGCData->getReservoirHeight():1.0;

				for (j=0; j<dim; j++) Cij[j] = strtod(str[j+4].c_str(),0);

				if (Cij[0] == .0 && Cij[1] == .0 && Cij[2] == .0){
					char msg[256]; sprintf(msg,"Cij from edge %d-%d is completely null!",EN_id(edge->get(0,0)),EN_id(edge->get(0,1)));
					throw Exception(__LINE__,__FILE__,msg);
				}


				double aux = 0;
				for (j=0; j<dim; j++) aux += Cij[j]*Cij[j];
				aux = sqrt(aux);

				// BACALHO pra multi-dominios
				if (aux > 0.0){
					pGCData->setCij(edge,setDomains[i],Cij);
					pGCData->setCij_norm(edge,setDomains[i],aux);
					//pGCData->setNumRemoteCopies(edge,nrc);
					//pGCData->setNumRC(edge,setDomains[i],nrc);
					pGCData->set_belongsToBoundary(edge, (setDomains[i] != iGrp) );

					dblarray edIJ(dim,.0);
					double coord1[3],coord2[3];
					V_coord(theMesh->getVertex(edgeIDs[0]),coord1);
					V_coord(theMesh->getVertex(edgeIDs[1]),coord2);
					for (j=0; j<dim; j++) edIJ[j]=coord2[j]-coord1[j];
					double length = sqrt( std::inner_product(edIJ.begin(),edIJ.end(),edIJ.begin(),.0) );
					avgLength += length;
					if (max_edge<length){
						max_edge=length;
					}
					edIJ.clear();
				}
			}
			getline(fid,str[0],'\n');
		}

		if (!hasRC && P_size()>1){
			throw Exception(__LINE__,__FILE__,"No remote copies detected. Weird situation for parallel simulation.\n");
		}
		if (!P_pid()) printf("done.\n");
	}

	void loadDij(ifstream &fid, pMesh theMesh, GeomData *pGCData, int *setDomains){
		if (!P_pid()) printf("\tloading Dij... ");

		int IDs[3], doms[2], j, iGrp;
		int dim = pGCData->getMeshDim();
		int pos = 2*dim + 2;

		std::string str[10];
		dblarray Dij(3);
		mEntity *e;
		mEdge *edge;
		mFace *face;
		getline(fid,str[0],'\n');
		getline(fid,str[0],'\n');
		while ( 1 ){
			// flag dom1 dom2 id0 id1 id2      Dijx                Dijy             Dijz
			// 999  3300 4400 161 160 866 0.000000000000000 -0.000173423583642 -0.000000000000000
			//  0     1   2    3   4   5         6                  7                 8
			fid>>str[0];
			if (str[0] == "end") break;
			iGrp = atoi(str[0].c_str());
			for (j=1; j<=pos; j++) fid >> str[j];
			for (j=0; j<2; j++) 	doms[j] = atoi(str[j+1].c_str());
			for (j=0; j<dim; j++) IDs[j] = atoi(str[j+3].c_str());

			/*
			 * For 2-D domains, multiply Dij by reservoir height (H) for 2D/3D
			 * simulations (physics occurs only on 2-D but reservoir volume
			 * is considered)
			 */
			double H = (pGCData->getMeshDim()==2)?pGCData->getReservoirHeight():1.0;
			for (j=0; j<dim; j++) Dij[j] = strtod(str[j+dim+3].c_str(),0);

			if (dim==2) // edge already exist!
				edge = theMesh->getEdge(theMesh->getVertex(IDs[0]),theMesh->getVertex(IDs[1]));
			else{
				// face must be created!
				face = theMesh->getTri(theMesh->getVertex(IDs[0]),theMesh->getVertex(IDs[1]),theMesh->getVertex(IDs[2]));
				if (!face)
					face = theMesh->createFaceWithVertices(IDs[0],IDs[1],IDs[2],theMesh->getGEntity(iGrp,2));
				else
					throw Exception(__LINE__,__FILE__,"Face already exist!\n");
			}

			e = (dim==2)?(mEntity*)edge:(mEntity*)face;
			pGCData->setDij((pEntity)e,doms[0],doms[1],Dij);
			pGCData->setFlag((pEntity)e,iGrp);
		}
		if (!P_pid()) printf("done.\n");
	}
}
