/*
 * calculateGeomCoefficients.cpp
 *
 *  Created on: Oct 17, 2014
 *      Author: rogerio
 */

#include "EBFV1_modified.h"

void calculateGeomCoefficients(pEntity face, double *Cij, double *Dij, double &vol){
	// Cij_x = Cij[0];		Dij_x = Dij[0];
	// Cij_y = Cij[1];		Dij_y = Dij[1];
	// Cik_x = Cij[2];		Dij_x = Dij[2];
	// Cik_y = Cij[3];		Dij_y = Dij[3];
	// Cjk_x = Cij[4];		Dij_x = Dij[4];
	// Cjk_y = Cij[5];		Dij_y = Dij[5];

	int i, pos, id0, id1;
	double faceCenter[3], edgeCenter[3], IJ[2], innerprod, sign;
	double I[3], J[3], v[2];
	pEntity edge;

	pos = 0;
	getFCenter(face,faceCenter);
	for (i=0; i<3; i++){					// loop over all three face's edges to calculate Cij and Dij
		edge = F_edge(face, i);
		edgeCenter[0] = .0;
		edgeCenter[1] = .0;
		E_center(edge,edgeCenter);			// get edge middle point
		V_coord(edge->get(0,0),I);			// get node I coordinate
		V_coord(edge->get(0,1),J);			// get node J coordinate
		id0 = EN_id(edge->get(0,0));		// get node I ID
		id1 = EN_id(edge->get(0,1));		// get node J ID
		sign = ( id0 > id1 )?-1.0:1.0; 		// vector IJ must point from the smaller vertex ID to the greater
		IJ[0] = sign*(J[0]-I[0]);			// edge vector, used as a reference vector
		IJ[1] = sign*(J[1]-I[1]);
		v[0] = edgeCenter[0]-faceCenter[0];	// reference vector
		v[1] = edgeCenter[1]-faceCenter[1];
		innerprod = v[1]*IJ[0] + (-v[0])*IJ[1];	// Cij must point as if I inside the CV and J outside
		if ( innerprod <= .0 ){
			v[0] = -v[0];
			v[1] = -v[1];
		}
		Cij[pos] = v[1];
		Cij[pos+1] = -v[0];
		Dij[pos] = -(J[1]-I[1])/2.0;					// Dij vector is orthogonal to edge (it's unknown Dij orientation)
		Dij[pos+1] =  (J[0]-I[0])/2.0;
		v[0] = edgeCenter[0]-faceCenter[0];				// reference vector
		v[1] = edgeCenter[1]-faceCenter[1];

		innerprod = Dij[pos]*v[0] + Dij[pos+1]*v[1];	// Dij must point to outside element
		if (  innerprod <= .0 ){
			Dij[pos] = -Dij[pos];
			Dij[pos+1] = -Dij[pos+1];
		}
		pos += 2;
	}
	vol = F_area(face)/3.0;								// calculate node control volume
}
