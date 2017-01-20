/*
 * assemblyMatrix_A.cpp
 *
 *  Created on: Oct 17, 2014
 *      Author: rogerio
 */
#include "EBFV1_modified.h"

/*
     A11 = (E_ij[0,0]+E_ik[0,0])*(F_ij[0,0]+F_ik[0,0])+
          (E_ij[0,1]+E_ik[0,1])*(F_ij[1,0]+F_ik[1,0])+
           E_ij[0,2]*F_ij[2,0]+
           E_ij[0,3]*F_ij[3,0]+
           E_ik[0,2]*F_ik[2,0]+
           E_ik[0,3]*F_ik[3,0]+
           G_ij[0,0]+G_ik[0,0]


    A12 = (E_ij[0,0]+E_ik[0,0])*F_ij[0,1]+
          (E_ij[0,1]+E_ik[0,1])*F_ij[1,1]+
           E_ij[0,2]*(F_ij[2,1]+F_jk[0,0])+
           E_ij[0,3]*(F_ij[3,1]+F_jk[1,0])+
           E_ik[0,2]*F_jk[2,0]+
           E_ik[0,3]*F_jk[3,0]+
           G_ij[0,1]

    A13 = (E_ij[0,0]+E_ik[0,0])*F_ik[0,1]+
          (E_ij[0,1]+E_ik[0,1])*F_ik[1,1]+
           E_ij[0,2]*F_jk[0,1])+
           E_ij[0,3]*F_jk[1,1])+
           E_ik[0,2]*(F_ik[2,1]+F_jk[2,1])+
           E_ik[0,3]*(F_ik[3,1]+F_jk[3,1])+
           G_ik[0,1]

    A21 =  E_ij[1,0]*(F_ij[0,0]+F_ik[0,0])+
           E_ij[1,1]*(F_ij[1,0]+F_ik[1,0])+
           (E_ij[1,2]+E_jk[0,0])*F_ij[2,0]+
           (E_ij[1,3]+E_jk[0,1])*F_ij[3,0]+
           E_jk[0,2]*F_ik[2,0]+
           E_jk[0,3]*F_ik[3,0]+
           G_ij[1,0]

    A22 =  E_ij[1,0]*F_ij[0,1]+
           E_ij[1,1]*F_ij[1,1]+
           (E_ij[1,2]+E_jk[0,0])*(F_ij[2,1]+F_jk[0,0])+
           (E_ij[1,3]+E_jk[0,1])*(F_ij[3,1]+F_jk[1,0])+
           E_jk[0,2]*F_jk[2,0]+
           E_jk[0,3]*F_jk[3,0]+
           G_ij[1,1]+G_jk[0,0]

    A23 =  E_ij[1,0]*F_ik[0,1]+
           E_ij[1,1]*F_ik[1,1]+
           (E_ij[1,2]+E_jk[0,0])*F_jk[0,1]+
           (E_ij[1,3]+E_jk[0,1])*F_jk[1,1]+
           E_jk[0,2]*(F_ik[2,1]+F_jk[2,1])+
           E_jk[0,3]*(F_ik[3,1]+F_jk[3,1])+
           G_jk[0,1]

    A31 = E_ik[1,0]*(F_ij[0,0]+F_ik[0,0])+
          E_ik[1,1]*(F_ij[1,0]+F_ik[1,0])+
          E_jk[1,0]*F_ij[2,0]+
          E_jk[1,1]*F_ij[3,0]+
          (E_ik[1,2]+E_jk[1,2])*F_ik[2,0]+
          (E_ik[1,3]*E_jk[1,3])*F_ik[3,0]+
          G_ik[1,0]

    A32 = E_ik[1,0]*F_ij[0,1]+
          E_ik[1,1]*F_ij[1,1]+
          E_jk[1,0]*(F_ij[2,1]+F_jk[0,0])+
          E_jk[1,1]*(F_ij[3,1]+F_jk[1,0])+
          (E_ik[1,2]+E_jk[1,2])*F_jk[2,0]+
          (E_ik[1,3]*E_jk[1,3])*F_jk[3,0]+
          G_jk[1,0]

    A33 = E_ik[1,0]*F_ik[0,1]+
          E_ik[1,1]*F_ik[1,1]+
          E_jk[1,0]*F_jk[0,1]+
          E_jk[1,1]*F_jk[1,1]+
          (E_ik[1,2]+E_jk[1,2])*(F_ik[2,1]+F_jk[2,1])+
          (E_ik[1,3]*E_jk[1,3])*(F_ik[3,1]+F_jk[3,1])+
          G_ik[1,1]+G_jk[1,1]

 */

void assemblyMatrix_A(const double* E_ij, const double* E_jk, const double* E_ik,
		const double* F_ij, const double* F_jk, const double* F_ik,
		const double* G_ij, const double* G_jk, const double* G_ik,double* A_ij,double* A_jk, double* A_ik){

	A_ij[0] = E_ij[0]*(F_ij[0]+F_ik[0])+E_ij[1]*(F_ij[2]+F_ik[2])+E_ij[2]*F_ij[4]+E_ij[3]*F_ij[6]+G_ij[0];
	A_ij[1] = E_ij[0]*F_ij[1]+E_ij[1]*F_ij[3]+E_ij[2]*(F_ij[5]+F_jk[0])+E_ij[3]*(F_ij[7]+F_jk[2])+G_ij[1];
	A_ij[2] = E_ij[0]*F_ik[1]+E_ij[1]*F_ik[3]+E_ij[2]*F_jk[1]+E_ij[3]*F_jk[3];
	A_ij[3] = E_ij[4]*(F_ij[0]+F_ik[0])+E_ij[5]*(F_ij[2]+F_ik[2])+E_ij[6]*F_ij[4]+E_ij[7]*F_ij[6]+G_ij[2];
	A_ij[4] = E_ij[4]*F_ij[1]+E_ij[5]*F_ij[3]+E_ij[6]*(F_ij[5]+F_jk[0])+E_ij[7]*(F_ij[7]+F_jk[2])+G_ij[3];
	A_ij[5] = E_ij[4]*F_ik[1]+E_ij[5]*F_ik[3]+E_ij[6]*F_jk[1]+E_ij[7]*F_jk[3];

	A_jk[0] = E_jk[0]*F_ij[4]+E_jk[1]*F_ij[6]+E_jk[2]*F_ik[4]+E_jk[3]*F_ik[6];
	A_jk[1] = E_jk[0]*(F_ij[5]+F_jk[0])+E_jk[1]*(F_ij[7]+F_jk[2])+E_jk[2]*F_jk[4]+E_jk[3]*F_jk[6]+G_jk[0];
	A_jk[2] = E_jk[0]*F_jk[1]+E_jk[1]*F_jk[3]+E_jk[2]*(F_ik[5]+F_jk[5])+E_jk[3]*(F_ik[7]+F_jk[7])+G_jk[1];
	A_jk[3] = E_jk[4]*F_ij[4]+E_jk[5]*F_ij[6]+E_jk[6]*F_ik[4]+E_jk[7]*F_ik[6];
	A_jk[4] = E_jk[4]*(F_ij[5]+F_jk[0])+E_jk[5]*(F_ij[7]+F_jk[2])+E_jk[6]*F_jk[4]+E_jk[7]*F_jk[6]+G_jk[2];
	A_jk[5] = E_jk[4]*F_jk[1]+E_jk[5]*F_jk[3]+E_jk[6]*(F_ik[5]+F_jk[5])+E_jk[7]*(F_ik[7]+F_jk[7])+G_jk[3];

	A_ik[0] = E_ik[0]*(F_ij[0]+F_ik[0])+E_ik[1]*(F_ij[2]+F_ik[2])+E_ik[2]*F_ik[4]+E_ik[3]*F_ik[6]+G_ik[0];
	A_ik[1] = E_ik[0]*F_ij[1]+E_ik[1]*F_ij[3]+E_ik[2]*F_jk[4]+E_ik[3]*F_jk[6];
	A_ik[2] = E_ik[0]*F_ik[1]+E_ik[1]*F_ik[3]+E_ik[2]*(F_ik[5]+F_jk[5])+E_ik[3]*(F_ik[7]+F_jk[7])+G_ik[1];
	A_ik[3] = E_ik[4]*(F_ij[0]+F_ik[0])+E_ik[5]*(F_ij[2]+F_ik[2])+E_ik[6]*F_ik[4]+E_ik[7]*F_ik[6]+G_ik[2];
	A_ik[4] = E_ik[4]*F_ij[1]+E_ik[5]*F_ij[3]+E_ik[6]*F_jk[4]+E_ik[7]*F_jk[6];
	A_ik[5] = E_ik[4]*F_ik[1]+E_ik[5]*F_ik[3]+E_ik[6]*(F_ik[5]+F_jk[5])+E_ik[7]*(F_ik[7]+F_jk[7])+G_ik[3];
}
