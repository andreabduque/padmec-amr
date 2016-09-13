#include "EBFV/EBFV1_elliptic.h"

// =============================================================================
// Matrix E is responsible to project the gradient over an orthogonal plane to
// edge. E is assembled with a loop over all edges. In parallel, if an edge is
// on the partition boundary, Eij must be divided by (NRC) Number of Remote
// Copies. Each partition will add Eij/NRC to the global matrix E and the resul-
// ting Eij contribution will be the same as it was assembled in serial.
// =============================================================================


namespace PRS{
	int EBFV1_elliptic::divergence_E(Mat E, const double *Cij, int edge, int dom, int dom_flag, int idx0_global, int idx1_global, int id0, int id1, int dim, int counter){
		const double I2D[4] = {1.0,.0,.0,1.0};
		const double I3D[9] = {1.0,.0,.0,.0,1.0,.0,.0,.0,1.0};
		int i, j, k;
		double sign = 1.0;

		if (id0 > id1){
			std::swap(id0,id1);
			sign = -sign;
		}

		// index for global Fg
		id0 = pMData->get_AppToPETSc_Ordering(id0);
		id1 = pMData->get_AppToPETSc_Ordering(id1);

		// get absolute K tensor: dim x dim
		const double *K = pSimPar->getPermeability(dom_flag);

		// ####################################################################
		// what is calculated below: Eij = -(K*Mob_IJ)/2*|I-versor*versor I-versor*versor|*Cij
		//												 |versor*versor-I versor*versor-I|
		// ####################################################################

		// product between two vectors: versor*versor => matrix 3x3
		double versor[3];
		pGCData->getVersor(dom,edge,versor);
		for (i=0; i<dim; i++){
			versor[i] *= sign;
		}
		double matLij[dim*dim];
		k = 0;
		for (i=0; i<dim; i++){
			for (j=0; j<dim; j++){
				matLij[k++] = versor[i]*versor[j];
			}
		}

		const double *Identity = (dim==2)?I2D:I3D;
		double matSubtrac_ILij[dim*dim];
		for (k=0; k<dim*dim; k++) matSubtrac_ILij[k] = -0.5*(Identity[k] - matLij[k]);

		dblarray EA(dim*dim,.0);
		k = 0;
		int pos = 0;
		if ( pSimPar->is_K_Isotropic() ){
			int pos1 = 0;
			int pos2 = 0;
			for (i=0; i<dim; i++){
				for (j=0; j<dim; j++) EA[k++] = K[pos1]*matSubtrac_ILij[pos2+j];
				pos1 += dim+1;
				pos2 += dim;
			}
		}else{
			for (i=0; i<dim; i++){
				for (j=0; j<dim; j++){
					for (k=0; k<dim; k++)
						EA[pos] += K[dim*i+k]*matSubtrac_ILij[dim*k+j];
					pos++;
				}
			}
		}

		/*
		 * Se o meio é isotropico:
		 *
		 *     |K00*matSubtrac_ILij00   K00*matSubtrac_ILij01|
		 * EA =|                                             |
		 *     |K11*matSubtrac_ILij10   K11*matSubtrac_ILij11|
		 *
		 *
		 * */

		double Eij_part1[3] = {.0,.0,.0};
		for (i=0; i<dim; i++){
			for (j=0; j<dim; j++){
				Eij_part1[i] += Cij[j]*EA[dim*j+i];
			}
		}

		// fill Eij matrix
		double Eij_row1[2*dim], Eij_row2[2*dim];
		for (i=0; i<dim; i++){
			Eij_row1[i] = Eij_part1[i];
			Eij_row1[i+dim] = Eij_row1[i];
			Eij_row2[i] = -Eij_row1[i];
			Eij_row2[i+dim] = -Eij_row1[i];
		}

		K = 0;

		for (i=0;i<2*dim; i++){
			pMAS->Eij[counter][i] = Eij_row1[i];
			pMAS->Eij[counter][i+2*dim] = Eij_row2[i];
		}
		return 0;
	}
}

