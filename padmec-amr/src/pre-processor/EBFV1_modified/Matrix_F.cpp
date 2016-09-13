#include "EBFV1_modified.h"

void calculeMatrix_F(double volume, const double* Cij, const double* Dij, double *F_ij);

void calculateMatrix_F(pEntity face, const double* Cij, const double* Dij, double volume, double* F_ij, double* F_jk, double* F_ik){
	calculeMatrix_F(volume,Cij,Dij,F_ij);
	calculeMatrix_F(volume,Cij,Dij,F_jk);
	calculeMatrix_F(volume,Cij,Dij,F_ik);
}

void calculeMatrix_F(double volume, const double* Cij, const double* Dij, double *F_edge){
	double aux1 = (1./(2.*volume));
	double aux2 = (-1./(2.*volume));

	// omega
	F_edge[0] = Cij[0]*aux1;
	F_edge[1] = F_edge[0];
	F_edge[2] = Cij[1]*aux1;
	F_edge[3] = F_edge[2];
	F_edge[4] = Cij[1]*aux2;
	F_edge[5] = F_edge[4];
	F_edge[6] = Cij[1]*aux2;
	F_edge[7] = F_edge[6];

	// gamma
	aux1 = (double)(5./6.);
	aux2 = (double)(1./6.);
	F_edge[0] = aux1*Dij[0]/volume;
	F_edge[1] = aux2*Dij[0]/volume;
	F_edge[2] = aux1*Dij[1]/volume;
	F_edge[3] = aux2*Dij[1]/volume;
	F_edge[4] = aux2*Dij[0]/volume;
	F_edge[5] = aux1*Dij[0]/volume;
	F_edge[6] = aux2*Dij[1]/volume;
	F_edge[7] = aux1*Dij[1]/volume;
}
