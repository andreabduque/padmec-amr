/*
 * CalculateElementsError.cpp
 *
 *  Created on: 24/02/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis.h"

void ErrorAnalysis::calculate_ElementsError(SimulatorParameters *pSimPar, GeomData* pGCData,void(*pFunc_getGrad)(FIELD,int,int,int,const double*&), FIELD field){
	int i, j;

	// mesh dimension
	int dim = pGCData->getMeshDim();

	// arrays of pointers for arrays: gradients and delta gradients
	const double* grad[4] = {NULL, NULL, NULL, NULL};
//	for (i=0; i<dim+1; i++){
//		grad[i] = new double[dim];
//	}

	// vector created over elements' edges
	int num_edges = 3*(dim-1);
	double* vectors[num_edges];
	double* delt_grad[num_edges];
	for (i=0;i<num_edges;i++){
		vectors[i] = new double[dim];
		delt_grad[i] = new double[dim];
	}

	// indices to automate edge vectors creation
	int idx[6][2] = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};

	// automate delta gradient calculation: triangle has 3 grad. diff | tetrahedron has 6 grad. differences
	int dltGrd_idx[6][2] = {{0,1},{1,2},{2,0},{3,0},{3,1},{3,2}};

	// indices: elements connectivities
	const int* indices = NULL;

	// elements' coordinates
	const double* coords[4] = {NULL,NULL,NULL,NULL};

	// number of domains
	int ndom = pSimPar->getNumDomains();

	// element counter
	int ith_face = 0;

	// Auxiliary avriable
	int pos = dim+1;

	cout << setprecision(5);

	// loop over domains
	for (int dom=0; dom<ndom; dom++){

		// loop over domains' elements
		for (int row=0; row<pGCData->getNumElemPerDomain(dom); row++){
			pGCData->getElement(dom,row,indices);

			// get vertices coordinates
			for (i=0; i<pos; i++){
				pGCData->getCoordinates(indices[i+pos],coords[i]);
//				cout << "Node - " << indices[i+pos] << " ";
//				for (j=0;j<3;j++){
//					cout << setprecision(5) << (coords[i])[j] << " ";
//				}
//				cout << endl;
			}

			// create edge vectors
			for (i=0; i<num_edges; i++){
				//double* vec = vectors[i];
				const double* vtx1 = coords[ idx[i][0] ];
				const double* vtx2 = coords[ idx[i][1] ];
				for (j=0; j<dim; j++){
					(vectors[i])[j] = vtx1[j] - vtx2[j];
				}
			}


//			for (i=0; i<num_edges; i++){
//				cout << "<";
//				for (j=0; j<dim; j++){
//					cout << (vectors[i])[j] << ",";
//				}
//				cout << ">\t";
//			}
//			cout << endl;

			// get nodal gradients for one element
			for (i=0; i<pos; i++){
				pFunc_getGrad(field,dom,indices[i],indices[i+pos],grad[i]);
			}

//			for (i=0; i<dim+1; i++){
//				cout << "<";
//				for (j=0; j<dim; j++){
//					cout << (grad[i])[j] << ",";
//				}
//				cout << ">\t";
//			}
//			cout << endl;

			// gradient variation
			for (i=0; i<num_edges; i++){    /*number of element vertices*/
				//double* deltgrd = delt_grad[i];
				const double* grd_A = grad[ dltGrd_idx[i][0] ];
				const double* grd_B = grad[ dltGrd_idx[i][1] ];
				//cout << grd_A << "\t" << grd_B << "\t\t" << dltGrd_idx[i][0] << "," << dltGrd_idx[i][1] << endl;
				for (j=0; j<dim; j++){  /*number of coordinates*/
					(delt_grad[i])[j] = grd_B[j] - grd_A[j];
				}
			}

//			for (i=0; i<dim+1; i++){
//				cout << "<";
//				for (j=0; j<dim; j++){
//					cout << (delt_grad[i])[j] << ",";
//				}
//				cout << ">\t";
//			}
//			cout << endl;

			// inner product
			double e, error = .0;
			for (i=0; i<num_edges; i++){    /*number of element vertices*/
				e = inner_product(delt_grad[i],vectors[i],dim);
				error += fabs(e);
			}

			error = (double)(error/(dim+1.0)); // divide by the number of element vertices
			setElementError(ith_face,sqrt(error));
			ith_face++;
		}
	}

	// free memory
	for (i=0; i<dim+1; i++){
		//delete[] grad[i]; grad[i] = 0;
		delete[] delt_grad[i]; delt_grad[i] = 0;
	}
	for (i=0;i<num_edges;i++){
		delete[] vectors[i]; vectors[i] = 0;
	}
}
