/*
 * gradientSaturation.cpp
 *
 *  Created on: 14/05/2009
 *      Author: rogerio
 */

#include "EBFV/EBFV1_elliptic.h"

namespace PRS{
	void EBFV1_elliptic::calculatePressureGradient(){
		CPU_Profile::Start();

		int dim = pGCData->getMeshDim();
		resetPressureGradient();
		int ndom = (int)pSimPar->setOfDomains.size();
		for (int dom=0; dom<ndom; dom++){
			calc_p_grad_1(dom,dim);
			calc_p_grad_2(dom,dim);
		}
		calc_p_grad_3(dim);

		CPU_Profile::End("calculatePressureGradient");
	}

	void EBFV1_elliptic::calc_p_grad_1(int dom, int dim){
		double p_I, p_J, val;
		const double* Cij = NULL;
		double* p_grad_I = NULL;
		double* p_grad_J = NULL;
		int i, nedges, edge, idx0, idx1, idx0_global, idx1_global, id0, id1;
		nedges = pGCData->getNumEdgesPerDomain(dom);
		for (edge=0; edge<nedges; edge++){
			pGCData->getCij(dom,edge,Cij);
			pGCData->getEdge(dom,edge,idx0,idx1,idx0_global,idx1_global);
			pGCData->getID(dom,idx0,idx1,id0,id1);
			if (id0 > id1){
				std::swap(idx0_global,idx1_global);
				std::swap(idx0,idx1);
			}
			pPPData->getPressure(idx0_global,p_I);
			pPPData->getPressure(idx1_global,p_J);
			pPPData->get_pw_Grad(dom,idx0,p_grad_I);
			pPPData->get_pw_Grad(dom,idx1,p_grad_J);
			val = 0.5*(p_I + p_J);
			for (i=0; i<dim; i++){
				p_grad_I[i] += val*Cij[i];
				p_grad_J[i] += -val*Cij[i];
			}
		}
	}

	void EBFV1_elliptic::calc_p_grad_2(int dom, int dim){
		const double* Dij = NULL;
		double* p_grad_I = NULL;
		double* p_grad_J = NULL;
		double* p_grad_K = NULL;

		double p_I, p_J;
		int nedges, edge, idx0, idx1, idx2, idx0_global, idx1_global, idx2_global;

		if (dim==2){
			nedges = pGCData->getNumBDRYEdgesPerDomain(dom);
			for (edge=0; edge<nedges; edge++){
				pGCData->getDij(dom,edge,Dij);
				pGCData->getBdryEdge(dom,edge,idx0,idx1,idx0_global,idx1_global);
				pPPData->getPressure(idx0_global,p_I);
				pPPData->getPressure(idx1_global,p_J);
				pPPData->get_pw_Grad(dom,idx0,p_grad_I);
				pPPData->get_pw_Grad(dom,idx1,p_grad_J);
				for (int i=0; i<dim; i++){
					p_grad_I[i] += ((5.*p_I + p_J)/6.)*Dij[i];
					p_grad_J[i] += ((p_I + 5.*p_J)/6.)*Dij[i];
				}
			}
		}
		else{
			int nfaces = pGCData->getNumBdryFacesPerDomain(dom);
			double line1[3] = {6., 1., 1.};
			double line2[3] = {1., 6., 1.};
			double line3[3] = {1., 1., 6.};
			double dot1, dot2, dot3, p_vec[3], aux[3];
			for (int face=0; face<nfaces; face++){
				pGCData->getDij(dom,face,Dij);
				pGCData->getBdryFace(dom,face,idx0,idx1,idx2,idx0_global,idx1_global,idx2_global);
				pPPData->getPressure(idx0_global,p_vec[0]);
				pPPData->getPressure(idx1_global,p_vec[1]);
				pPPData->getPressure(idx2_global,p_vec[2]);
				pPPData->get_pw_Grad(dom,idx0,p_grad_I);
				pPPData->get_pw_Grad(dom,idx1,p_grad_J);
				pPPData->get_pw_Grad(dom,idx2,p_grad_K);

				dot1 = inner_product(line1,p_vec,3);
				dot2 = inner_product(line2,p_vec,3);
				dot3 = inner_product(line3,p_vec,3);
				for (int i=0; i<3; i++){
					aux[i] = Dij[i]/8.0;
				}
				for (int i=0; i<3; i++){
					p_grad_I[i] += aux[i]*dot1;
					p_grad_J[i] += aux[i]*dot2;
					p_grad_K[i] += aux[i]*dot3;
				}
			}
		}
	}

	void EBFV1_elliptic::calc_p_grad_3(int dim){
		double* p_grad = NULL;
		double vol;
		int i, dom, ndom, nnodes, node;

		ndom = (int)pSimPar->setOfDomains.size();
		for (dom=0; dom<ndom; dom++){
			nnodes = pGCData->getNumNodesPerDomain(dom);
			for (node=0; node<nnodes; node++){
				pPPData->get_pw_Grad(dom,node,p_grad);
				pGCData->getVolume(dom,node,vol);
				for (i=0; i<dim; i++){
					p_grad[i] /= vol;
				}
			}
		}
	}

	void EBFV1_elliptic::resetPressureGradient(){
		double* p_grad = NULL;
		int dom, ndom, nnodes, node, i;
		ndom = (int)pSimPar->setOfDomains.size();
		for (dom=0; dom<ndom; dom++){
			nnodes = pGCData->getNumNodesPerDomain(dom);
			for (node=0; node<nnodes; node++){
				pPPData->get_pw_Grad(dom,node,p_grad);
				for (i=0; i<3; i++){
					p_grad[i] = 0;
				}
			}
		}
	}
}
