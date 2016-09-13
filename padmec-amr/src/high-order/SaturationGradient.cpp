/*
 * gradientSaturation.cpp
 *
 *  Created on: 14/05/2009
 *      Author: rogerio
 */

#include "EBFV1_hyperbolic.h"

namespace PRS{
	double EBFV1_hyperbolic::calculateSaturationGradient(pMesh theMesh){
		int dim = pGCData->getMeshDim();
		resetSaturationGradient(theMesh);
		int ndom = (int)pSimPar->setOfDomains.size();
		for (int dom=0; dom<ndom; dom++){
			calc_Sw_grad_1(dom,dim);
			calc_Sw_grad_2(dom,dim);
		}
		calc_Sw_grad_3(theMesh,dim);
		calc_Sw_grad_4(dim);
		double Sw_grad[3];
		int nnodes = M_numVertices(theMesh);
		for (int node=0; node<295; node++){
			pPPData->get_Sw_Grad(node,Sw_grad);
			//cout << node+1 << setprecision(6) << fixed << " " << Sw_grad[0] << " " << Sw_grad[1] << endl;
		}
		//STOP();
		return 0;
	}

	void EBFV1_hyperbolic::calc_Sw_grad_1(int dom, int dim){
		double Cij[3], Sw_grad_I[3], Sw_grad_J[3], Sw_I, Sw_J, val;
		int i,j,nedges, edge, idx0, idx1, idx0_global, idx1_global, id0, id1, dom_flag;
		nedges = pGCData->getNumEdgesPerDomain(dom);
		for (edge=0; edge<nedges; edge++){
			pGCData->getCij(dom,edge,Cij);
			pGCData->getEdge(dom,edge,idx0,idx1,idx0_global,idx1_global);
			pGCData->getID(dom,idx0,idx1,id0,id1);
			if (id0 > id1){
				std::swap(idx0_global,idx1_global);
				std::swap(idx0,idx1);
			}
			Sw_I = pPPData->getSaturation(idx0_global);
			Sw_J = pPPData->getSaturation(idx1_global);
			pPPData->get_Sw_Grad(dom,idx0,Sw_grad_I);
			pPPData->get_Sw_Grad(dom,idx1,Sw_grad_J);
			val = 0.5*(Sw_I + Sw_J);
			for (i=0; i<dim; i++){
				Sw_grad_I[i] += val*Cij[i];
				Sw_grad_J[i] += -val*Cij[i];
			}
			pPPData->set_Sw_Grad(dom,idx0,Sw_grad_I);
			pPPData->set_Sw_Grad(dom,idx1,Sw_grad_J);
		}
	}

	void EBFV1_hyperbolic::calc_Sw_grad_2(int dom, int dim){
		double Dij[3], Sw_grad_I[3], Sw_grad_J[3], Sw_I, Sw_J;
		int i,j,nedges, edge, idx0, idx1,idx0_global, idx1_global;
		nedges = pGCData->getNumBDRYEdgesPerDomain(dom);
		for (edge=0; edge<nedges; edge++){
			pGCData->getDij(dom,edge,Dij);
			pGCData->getBdryEdge(dom,edge,idx0,idx1,idx0_global,idx1_global);
			Sw_I = pPPData->getSaturation(idx0_global);
			Sw_J = pPPData->getSaturation(idx1_global);
			pPPData->get_Sw_Grad(dom,idx0,Sw_grad_I);
			pPPData->get_Sw_Grad(dom,idx1,Sw_grad_J);
			for (int i=0; i<dim; i++){
				Sw_grad_I[i] += ((5.*Sw_I + Sw_J)/6.)*Dij[i];
				Sw_grad_J[i] += ((Sw_I + 5.*Sw_J)/6.)*Dij[i];
			}
			pPPData->set_Sw_Grad(dom,idx0,Sw_grad_I);
			pPPData->set_Sw_Grad(dom,idx1,Sw_grad_J);
		}
	}

	void EBFV1_hyperbolic::calc_Sw_grad_3(pMesh theMesh, int dim){
		double Sw_grad_tmp[3], Sw_grad[3], volume;
		int i, dom, ndom, nnodes, node, idx;
		ndom = (int)pSimPar->setOfDomains.size();
		// performe sw_grad accumulation for multidomains
		for (dom=0; dom<ndom; dom++){
			nnodes = pGCData->getNumNodesPerDomain(dom);
			for (node=0; node<nnodes; node++){
				pGCData->getNodeIdx_Global(dom,node,idx);
				pPPData->get_Sw_Grad(idx,Sw_grad);
				pPPData->get_Sw_Grad(dom,node,Sw_grad_tmp);
				//cout << node+1 << " " << Sw_grad_tmp[0] << " " << Sw_grad_tmp[1] << "\t";
				for (i=0; i<dim; i++){
					Sw_grad[i] += Sw_grad_tmp[i];
				}
				pPPData->set_Sw_Grad(idx,Sw_grad);
				//cout << node+1 << " " << Sw_grad[0] << " " << Sw_grad[1] << endl;
			}
		}
		// weight cumulated sw_grad per node volume
		nnodes = M_numVertices(theMesh);
		for (node=0; node<nnodes; node++){
			pPPData->get_Sw_Grad(node,Sw_grad);
			pGCData->getVolume(node,volume);
			for (i=0; i<dim; i++){
				Sw_grad[i] /= volume;
			}
			pPPData->set_Sw_Grad(node,Sw_grad);
			//cout << node+1 << " " << Sw_grad[0] << " " << Sw_grad[1] << endl;
		}
		//STOP();
	}

	void EBFV1_hyperbolic::calc_Sw_grad_4(int dim){
		double Sw_grad_I[3], Sw_grad_J[3], versor[3], innerp1, innerp2;
		int i,nedges, edge, idx0_global, idx1_global, flag1, flag2;

		nedges = pGCData->getNumEBE();
		for (edge=0; edge<nedges; edge++){
			pGCData->getEBE(edge,versor);
			pGCData->getEBE(edge,idx0_global,idx1_global,flag1,flag2);
			pPPData->get_Sw_Grad(idx0_global,Sw_grad_I);
			pPPData->get_Sw_Grad(idx1_global,Sw_grad_J);

			// Parallel Component of The Saturation Gradient
			innerp1 = inner_product(Sw_grad_I,versor,dim);
			innerp2 = inner_product(Sw_grad_J,versor,dim);
//			for (i=0; i<dim; i++){
//				innerp1 += Sw_grad_I[i]*versor[i];
//				innerp2 += Sw_grad_J[i]*versor[i];
//			}
			for (i=0; i<dim; i++){
				Sw_grad_I[i] = innerp1*versor[i];
				Sw_grad_J[i] = innerp2*versor[i];
			}

			// exclude injection wells
			if ( !pSimPar->isInjectionWell(flag1)  && !pPPData->getProjectedSwgrad(idx0_global)){
				pPPData->set_Sw_Grad(idx0_global,Sw_grad_I);
				pPPData->setProjectedSwgrad(idx0_global,true);
			}
			if ( !pSimPar->isInjectionWell(flag2)  && !pPPData->getProjectedSwgrad(idx1_global) ){
				pPPData->set_Sw_Grad(idx1_global,Sw_grad_J);
				pPPData->setProjectedSwgrad(idx1_global,true);
			}
		}
	}

	void EBFV1_hyperbolic::resetSaturationGradient(pMesh theMesh){
		double Sw_grad[3]={.0,.0,.0};
		int dom, ndom, nnodes, node;
		ndom = (int)pSimPar->setOfDomains.size();
		for (dom=0; dom<ndom; dom++){
			nnodes = pGCData->getNumNodesPerDomain(dom);
			for (node=0; node<nnodes; node++){
				pPPData->set_Sw_Grad(dom,node,Sw_grad);
			}
		}
		nnodes = M_numVertices(theMesh);
		for (node=0; node<nnodes; node++){
			pPPData->set_Sw_Grad(node,Sw_grad);
			pPPData->setProjectedSwgrad(node,false);
		}
	}
}
