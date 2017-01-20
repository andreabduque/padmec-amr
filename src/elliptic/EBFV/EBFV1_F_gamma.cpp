#include "EBFV/EBFV1_elliptic.h"

namespace PRS{

	int EBFV1_elliptic::gradient_F_bdry(Mat F, int dom){
		double Dij[3], Fij_column1[4], Fij_column2[4], volumeI, volumeJ, volumeK;
		int nedges, nfaces, id0, id1, id2, idx_0, idx_1, idx_2, idx0_global, idx1_global, idx2_global, i, j, pos1, pos2;

		if (pGCData->getMeshDim()==2){
			nedges = pGCData->getNumBDRYEdgesPerDomain(dom);
			for (j = 0; j<nedges; j++){
				pGCData->getBdryEdge(dom,j,idx_0,idx_1);
				pGCData->getBdryID(dom,idx_0,idx_1,id0,id1);
				pGCData->getBdryVolume(dom,idx_0,idx_1,volumeI,volumeJ);
				pGCData->getDij(dom,j,Dij);
				Fij_column1 [0] = 5.*Dij[0]/(6.*volumeI);
				Fij_column1 [1] = 5.*Dij[1]/(6.*volumeI);
				Fij_column1 [2] = Dij[0]/(6.*volumeJ);
				Fij_column1 [3] = Dij[1]/(6.*volumeJ);
				Fij_column2 [0] = Dij[0]/(6.*volumeI);
				Fij_column2 [1] = Dij[1]/(6.*volumeI);
				Fij_column2 [2] = 5.*Dij[0]/(6.*volumeJ);
				Fij_column2 [3] = 5.*Dij[1]/(6.*volumeJ);
				id0 = pMData->get_AppToPETSc_Ordering(id0);
				id1 = pMData->get_AppToPETSc_Ordering(id1);
				pos1 = 2*(id0-1);
				pos2 = 2*(id1-1);
				int idxm[4] = {pos1,pos1+1,pos2,pos2+1};
				int idxn[2] = {id0-1,id1-1};
				MatSetValues(F,4,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);
				MatSetValues(F,4,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);
			}
		}
		else{
			nfaces = pGCData->getNumBdryFacesPerDomain(dom);
			cout << "Faces: " << nfaces << endl;
			for (j = 0; j<nfaces; j++){
				pGCData->getBdryFace(dom,j,idx_0,idx_1,idx_2,idx0_global,idx1_global,idx2_global);
				pGCData->getBdryID(dom,idx_0,idx_1,idx_2,id0,id1,id2);
				pGCData->getBdryVolume(dom,idx_0,idx_1,idx_2,volumeI,volumeJ,volumeK);
				pGCData->getDij(dom,j,Dij);

				double tmp[3] = {1./(8.*volumeI), 1./(8.*volumeJ), 1./(8.*volumeK)};
				double aux[3][3] = {{6.*tmp[0],tmp[0],tmp[0]},{tmp[1],6.*tmp[1],tmp[1]},{tmp[2],tmp[2],6.*tmp[2]}};

				// fill edge matrix
				double Fij_column1[9], Fij_column2[9], Fij_column3[9];
				for (i=0; i<3; i++){
					Fij_column1[3*i] =   aux[i][0]*Dij[0];
					Fij_column1[3*i+1] = aux[i][0]*Dij[1];
					Fij_column1[3*i+2] = aux[i][0]*Dij[2];

					Fij_column2[3*i] = aux[i][1]*Dij[0];
					Fij_column2[3*i+1] = aux[i][1]*Dij[1];
					Fij_column2[3*i+2] = aux[i][1]*Dij[2];

					Fij_column3[3*i] = aux[i][2]*Dij[0];
					Fij_column3[3*i+1] = aux[i][2]*Dij[1];
					Fij_column3[3*i+2] = aux[i][2]*Dij[2];
				}

				// index for global Fg
				// where edge matrix must be assembled.
				id0 = pMData->get_AppToPETSc_Ordering(id0);
				id1 = pMData->get_AppToPETSc_Ordering(id1);
				id2 = pMData->get_AppToPETSc_Ordering(id2);

				int pos1 = 3*(id0-1);
				int pos2 = 3*(id1-1);
				int pos3 = 3*(id2-1);

				//cout << "id: " << id0 << "  " << id1 << "  " << id2 << endl;

				int idxm[9] = {pos1,pos1+1,pos1+2, pos2,pos2+1,pos2+2, pos3,pos3+1,pos3+2};
				int idxn[3] = {id0-1, id1-1, id2-1};

				MatSetValues(F,9,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);
				MatSetValues(F,9,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);
				MatSetValues(F,9,idxm,1,&idxn[2],Fij_column3,ADD_VALUES);
			}
		}
		return 0;
	}
}
