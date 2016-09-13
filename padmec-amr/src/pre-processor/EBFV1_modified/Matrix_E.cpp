#include "EBFV1_modified.h"

void calculateMatrix_E(pEntity face, GeomData* pGCData, const double* K, const double* Cij, double* E_ij, double* E_jk, double* E_ik){
	calculeMatrix_E((pEntity)F_edge(face,0),pGCData,K,Cij,E_ij);
	calculeMatrix_E((pEntity)F_edge(face,1),pGCData,K,Cij,E_jk);
	calculeMatrix_E((pEntity)F_edge(face,2),pGCData,K,Cij,E_ik);
}

void calculeMatrix_E(pEntity edge, GeomData* pGCData, const double* K, const double* Cij, double *E_ij){
	int dim = 2;
	const double Identity[4] = {1.0,.0,.0,1.0};
	int i, j, k, pos;
	double sign = 1.0;
	double versor[3], matLij[4], matSubtrac_ILij[4];

	int id0 = EN_id(edge->get(0,0));		// get node I ID
	int id1 = EN_id(edge->get(0,1));		// get node J ID
	sign = ( id0 > id1 )?-1.0:1.0; 		// vector IJ must point from the smaller vertex ID to the greater

	pGCData->getVersor(edge,versor);
	versor[0] *= sign;
	versor[1] *= sign;
	k = 0;
	for (i=0; i<dim; i++){
		for (j=0; j<dim; j++){
			matLij[k++] = versor[i]*versor[j];
		}
	}

	for (k=0; k<dim*dim; k++){
		matSubtrac_ILij[k] = -0.5*(Identity[k] - matLij[k]);
	}

	k = 0;
	pos = 0;
	double EA[4] = {.0,.0,.0,.0};
	for (i=0; i<dim; i++){
		for (j=0; j<dim; j++){
			for (k=0; k<dim; k++){
				EA[pos] += K[dim*i+k]*matSubtrac_ILij[dim*k+j];
			}
			pos++;
		}
	}

	double Eij_part1[3] = {.0,.0,.0};
	for (i=0; i<dim; i++){
		for (j=0; j<dim; j++){
			Eij_part1[i] += Cij[j]*EA[dim*j+i];
		}
	}

	E_ij[0] = Eij_part1[0];
	E_ij[1] = Eij_part1[1];
	E_ij[2] = E_ij[0];
	E_ij[3] = E_ij[1];

	E_ij[4] = -E_ij[0];
	E_ij[5] = -E_ij[1];
	E_ij[6] = -E_ij[0];
	E_ij[7] = -E_ij[1];
}

