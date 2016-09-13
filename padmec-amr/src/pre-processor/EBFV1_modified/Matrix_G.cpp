#include "EBFV1_modified.h"



void calculateMatrix_G(pEntity face, GeomData* pGCData, const double* K, const double* Cij, double* G_ij, double* G_jk, double* G_ik){
	calculeMatrix_G((pEntity)F_edge(face,0),pGCData,K,Cij,G_ij);
	calculeMatrix_G((pEntity)F_edge(face,1),pGCData,K,Cij,G_jk);
	calculeMatrix_G((pEntity)F_edge(face,2),pGCData,K,Cij,G_ik);
}

void calculeMatrix_G(pEntity edge, GeomData* pGCData, const double* K, const double* Cij, double *G_ij){
	int dim = 2;
	int i, j;
	double versor[3], length, aux, sign = 1.0;

	int id0 = EN_id(edge->get(0,0));		// get node I ID
	int id1 = EN_id(edge->get(0,1));		// get node J ID
	if (id0 > id1){
		std::swap(id0,id1);
		sign = -sign;
	}

	pGCData->getEdgeLength(edge,length);
	pGCData->getVersor(edge,versor);
	for (i=0; i<dim; i++){
		versor[i] *= sign;
	}

	double KL[3] = {.0, .0, .0};

	for (i=0; i<dim; i++){
		for (j=0; j<dim; j++){
			KL[i] += K[dim*i+j]*versor[j];
		}
	}

	aux = -(Cij[0]*KL[0]+Cij[1]*KL[1])/length;
	G_ij[0] = aux;
	G_ij[1] = -aux;
	G_ij[2] = -aux;
	G_ij[3] = aux;
}
