/*
 * Rogerio: 19/05/2015
 */
#include "EBFV1__pre-processors.h"

void validete_coefficients(pMesh theMesh, std::set<int>& setOfDomain, GeomData* pGCData){

#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: validate pre-processor coefficient calculation\tIN\n";
#endif


	cout << setprecision(8) << scientific;
	pEntity edge;
	std::vector<double> Cij(3), Dij(3);
	std::map<int,double*> coeff_sum;

	std::set<int>::iterator iter = setOfDomain.begin();

	for(;iter!=setOfDomain.end(); iter++){
		int dom = *iter;

		// initialize all nodes for domain "dom"
		double* p;
		VIter vit = M_vertexIter(theMesh);
		while (pEntity node = VIter_next(vit)){
			p = new double[3];
			p[0] = p[1] = p[2] = .0;
			coeff_sum[EN_id(node)] = p;
		}
		VIter_delete(vit);


		// calculate summation of coefficients
		EIter eit = M_edgeIter(theMesh);
		while ( (edge = EIter_next(eit)) ){
			int I = EN_id(edge->get(0,0));
			int J = EN_id(edge->get(0,1));
			pGCData->getCij(edge,dom,Cij);
			double sign = (I>J)?-1.0:1.0;


			p = coeff_sum[I];
			p[0] += sign*Cij[0];
			p[1] += sign*Cij[1];
			p[2] += sign*Cij[2];
			if (pGCData->getDij(edge,dom,Dij)){
				p[0] += Dij[0];
				p[1] += Dij[1];
				p[2] += Dij[2];
			}

			p = coeff_sum[J];
			p[0] += -sign*Cij[0];
			p[1] += -sign*Cij[1];
			p[2] += -sign*Cij[2];
			if (pGCData->getDij(edge,dom,Dij)){
				p[0] += Dij[0];
				p[1] += Dij[1];
				p[2] += Dij[2];
			}
		}
		EIter_delete(eit);

		double summation[3] = {.0,.0,.0};
		cout << "Domain: " << dom << endl;
		//cout << "ID:\t" << "Summ:\n";
		std::map<int,double*>::iterator miter = coeff_sum.begin();
		for(;miter!=coeff_sum.end();miter++){
			//int ID = miter->first;
			p = miter->second;
			//cout << ID << "\t" << p[0] << " " << p[1] << " " << p[2] << endl;
			summation[0] += p[0];
			summation[1] += p[1];
			summation[2] += p[2];
			dealloc_DOUBLE_vector(p);
		}
		cout <<  "\n\nSummation of all summations (Cij/Dij):  " << summation[0] << " " << summation[1] << " " << summation[2] << "\n";

		double vol_control = 0;
		vit = M_vertexIter(theMesh);
		while (pEntity node = VIter_next(vit)){
			vol_control += pGCData->getVolume(node,dom);
		}
		VIter_delete(vit);
		cout <<  "Total Volume                         :  " << vol_control << "\n\n\n";

	}
#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: validate pre-processor coefficient calculation\tOUT\n";
#endif

}
