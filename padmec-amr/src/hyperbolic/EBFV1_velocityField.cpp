/*
 * velocityField.cpp
 *
 *  Created on: 14/01/2009
 *      Author: rogerio
 */

#include "EBFV1_hyperbolic.h"

namespace PRS{

	double EBFV1_hyperbolic::calculateVelocityField(int dom, int dim){

		CPU_Profile::Start();

		const double* Cij = NULL;
		const double* pGrad_I = NULL;
		const double* pGrad_J = NULL;

		double p_I, p_J, Sw_I, Sw_J, MobI, MobJ, MobIJ, length, flag;
		double vel[3], versor[3], pw_grad_IJ[3];
		int i, j, idx0_global, idx1_global, idx0, idx1, id0, id1;

		int nedges = pGCData->getNumEdgesPerDomain(dom);
		for(int edge=0; edge<nedges; edge++){
			pGCData->getCij(dom,edge,Cij);
			pGCData->getEdge(dom,edge,idx0,idx1,idx0_global,idx1_global);
			pGCData->getID(dom,idx0,idx1,id0,id1);
			double sign = 1.0;
			if (id0 > id1){
				std::swap(idx0_global,idx1_global);
				std::swap(idx0,idx1);
				sign = -sign;
			}
			pPPData->getPressure(idx0_global,p_I);
			pPPData->getPressure(idx1_global,p_J);
			pPPData->get_pw_Grad_const(dom,idx0,pGrad_I);
			pPPData->get_pw_Grad_const(dom,idx1,pGrad_J);
			pGCData->getLength(dom,edge,length);
			pGCData->getVersor(dom,edge,versor);

			for (i=0; i<dim; i++){
				versor[i] *= sign;
			}
			for (i=0; i<dim; i++){
				versor[i] = -versor[i];
			}

			for (i=0; i<dim; i++){
				pw_grad_IJ[i] = .5*(pGrad_I[i] + pGrad_J[i]);
			}
			double VIJN[3] = {.0,.0,.0};
			double a = .0;
			for (i=0; i<dim; i++){
				a += versor[i]*pw_grad_IJ[i];
			}
			for (i=0; i<dim; i++){
				VIJN[i] = pw_grad_IJ[i] - a*versor[i];
			}
			a = (p_J - p_I)/length;
			double VIJFD[3] = {versor[0]*a, versor[1]*a, versor[2]*a};
			double aux[3] = {VIJFD[0]+VIJN[0], VIJFD[1]+VIJN[1], VIJFD[2]+VIJN[2]};
			pPPData->getVelocity_new(dom,edge,vel);
			pPPData->setVelocity_old(dom,edge,vel);
			flag = pGCData->getDomFlag(dom);
			const double *K = pSimPar->getPermeability(flag);
			pPPData->getSaturation(idx0_global,Sw_I);
			pPPData->getSaturation(idx1_global,Sw_J);
			MobI = pPPData->getTotalMobility(Sw_I);
			MobJ = pPPData->getTotalMobility(Sw_J);
			MobIJ = 0.5*(MobI + MobJ);

			for (i=0; i<dim; i++){
				vel[i] = 0;
			}
			for (i=0; i<dim; i++){
				for (j=0; j<dim; j++){
					vel[i] += K[dim*i+j]*aux[j];
				}
			}

			for (i=0; i<dim; i++){
				vel[i] *= -MobIJ;
			}
			pPPData->setVelocity_new(dom,edge,vel);
		}

		CPU_Profile::End("VelocityField");
		return 0;
	}
}

