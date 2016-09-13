#include "EBFV/EBFV1_elliptic.h"

namespace PRS{
//	int EBFV1_elliptic::gradient_F_edges(Mat F, pEntity edge, const int &dom, int dim, double *Cij){
int EBFV1_elliptic::gradient_F_edges(Mat F, const double *Cij, int dom, int idx0, int idx1, int id0, int id1, int dim){
		int i;
		double volumeI, volumeJ;

		// get nodes I and J
//		pEntity I = (pVertex)edge->get(0,0);
//		pEntity J = (pVertex)edge->get(0,1);
//		if (EN_id(edge->get(0,0)) > EN_id(edge->get(0,1))){
//			std::swap(I,J);
//		}
		if (id0 > id1){
			std::swap(id0,id1);
			std::swap(idx0,idx1);
		}

		// Where Fij should be assembled
		id0 = pMData->get_AppToPETSc_Ordering(id0);
		id1 = pMData->get_AppToPETSc_Ordering(id1);
		pGCData->getVolume(dom,idx0,idx1,volumeI,volumeJ);


		double aux1 = (1./(2.*volumeI));
		double aux2 = (-1./(2.*volumeJ));
		double Fij_column1[2*dim], Fij_column2[2*dim];
		for (i=0; i<dim; i++){
			Fij_column1[i] = Cij[i]*aux1;
			Fij_column2[i] = Fij_column1[i];
			Fij_column1[i+dim] = Cij[i]*aux2;
			Fij_column2[i+dim] = Fij_column1[i+dim];
		}

		int pos1 = dim*(id0-1);
		int pos2 = dim*(id1-1);
		int idxm[2*dim];				// says on which rows Fij must be assembled into Fg
		for (i=0; i<dim; i++){
			idxm[i] = pos1+i;
			idxm[dim+i] = pos2+i;
		}
		int idxn[2] = {id0-1,id1-1};	// says on which columns Fij must be assembled into Fg
		MatSetValues(F,2*dim,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);
		MatSetValues(F,2*dim,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);
		return 0;
	}
}
