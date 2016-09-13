#include "EBFV1_hyperbolic.h"

namespace PRS{
	double EBFV1_hyperbolic::calculateIntegralAdvectiveTerm(int dom, double &timeStep){
		CPU_Profile::Start();
		int dim = pGCData->getMeshDim();
		const double* Cij = NULL;
		double* Sw_grad_I = NULL;
		double* Sw_grad_J = NULL;
		double Sw_I, Sw_J, fwII, fwJJ, fwIJ, df_dsIJ, n, alpha;
		double  non_visc_fv, non_visc_ad, nonvisc_I, nonvisc_J;
		double courant, phi, length;
		double vel[3], FluxIJ[3], Cij_norm, versor[3], edIJ[3];//, Sw_grad_I[3], Sw_grad_J[3];
		int i,j,idx0_global, idx1_global, idx0, idx1, id0, id1,flag1,flag2;
		double sign, koef = (double)1./3.;

		int nedges = pGCData->getNumEdgesPerDomain(dom);
		for(j=0; j<nedges; j++){
			pGCData->getCij(dom,j,Cij);
			pGCData->getCij_norm(dom,j,Cij_norm);
			pGCData->getEdge(dom,j,idx0,idx1,idx0_global,idx1_global,flag1,flag2);
			pGCData->getID(dom,idx0,idx1,id0,id1);
			// todo: colocar este codigo em algum lugar de forma que so seja feito uma unica vez.
			sign = -1.0;
			if (id0 > id1){
				std::swap(idx0_global,idx1_global);
				sign = -sign;
			}

			// considering up-wind approximation for saturation field
			pPPData->getSaturation(idx0_global,Sw_I);
			pPPData->getSaturation(idx1_global,Sw_J);

			// get high order approximation for saturation
			if ( pSimPar->useHOApproximation() ){
				double SLI, SLJ, ratioI, ratioJ, delta_Sw, slimit_I, slimit_J,SLII,SLJJ,DSwII, DSwJJ;
				slimit_I = slimit_J =1.0;
				pGCData->getLength(dom,j,length);
				pGCData->getVersor(dom,j,versor);
				for (i=0; i<dim; i++){
					edIJ[i] = sign*versor[i]*length;
				}
				pPPData->get_Sw_Grad(idx0_global,Sw_grad_I);
				pPPData->get_Sw_Grad(idx1_global,Sw_grad_J);
				delta_Sw = (Sw_J - Sw_I);
				DSwII = 2.*inner_product(Sw_grad_I,edIJ,dim) - delta_Sw;
				DSwJJ = 2.*inner_product(Sw_grad_J,edIJ,dim) - delta_Sw;
				ratioI = (2.*DSwII*delta_Sw + qsi) / ( DSwII*DSwII + delta_Sw*delta_Sw + qsi);
				ratioJ = (2.*DSwJJ*delta_Sw + qsi) / ( DSwJJ*DSwJJ + delta_Sw*delta_Sw + qsi);
				SLII = (ratioI + fabs(ratioI) + qsi)/(1. + fabs(ratioI) + qsi);
				SLJJ = (ratioJ + fabs(ratioJ) + qsi)/(1. + fabs(ratioJ) + qsi);
				SLI = SLII*slimit_I;
				SLJ = SLJJ*slimit_J;
				if ( !pSimPar->isInjectionWell(flag1) ){
					Sw_I = Sw_I + (SLI/4.)*((1.-koef)*DSwII + (1.+koef)*delta_Sw);
				}
				if ( !pSimPar->isInjectionWell(flag2) ){
					Sw_J = Sw_J - (SLJ/4.)*((1.-koef)*DSwJJ + (1.+koef)*delta_Sw);
				}
			}

			fwII = pPPData->getFractionalFlux(Sw_I);
			fwJJ = pPPData->getFractionalFlux(Sw_J);
			fwIJ = 0.5*(fwII + fwJJ);

			// mid-edge total velocity
			pPPData->getVelocity_new(dom,j,vel);

			// Numerical Flux Function
			for (i=0; i<dim; i++){
				FluxIJ[i] = fwIJ*vel[i];
			}

			// Fractional Flux Flow Function Derivative (with respect to saturation)
			df_dsIJ = ( fabs(Sw_I-Sw_J) > 1.e-12 )?fabs((fwJJ-fwII)/(Sw_J-Sw_I)) : .0;

			// Approximate Eigenvalue (Note that we are using the linearized form of df_dsIJ)
			n = .0;
			for (i=0; i<dim; i++){
				n += vel[i]*vel[i];
			}
			alpha = sqrt(n)*df_dsIJ;

			// get the maximum alpha to compute the time step
			alpha_max = std::max(alpha,alpha_max);

			// Central difference Contribution
			non_visc_fv = .0;
			for (i=0; i<dim; i++){
				non_visc_fv +=  FluxIJ[i]*Cij[i];
			}

			// Numerical Diffusion
			non_visc_ad = 0.5*Cij_norm*alpha*(Sw_J - Sw_I);

			// Computing "Non-Viscous" Terms
			pPPData->getNonvisc(idx0_global,nonvisc_I);
			pPPData->getNonvisc(idx1_global,nonvisc_J);
			nonvisc_I = nonvisc_I + (non_visc_fv - non_visc_ad);
			nonvisc_J = nonvisc_J - (non_visc_fv - non_visc_ad);

			// update nonvisc term. it will be set to 0 at the next time iteration
			pPPData->setNonvisc(idx0_global,nonvisc_I);
			pPPData->setNonvisc(idx1_global,nonvisc_J);
		}

		courant = pSimPar->CFL();
		int flag = pGCData->getDomFlag(dom);
		phi = pSimPar->getPorosity(flag);
		length = pGCData->getSmallestEdgeLength();
		timeStep = std::min(timeStep,(courant*length*phi)/alpha_max);

		CPU_Profile::End("IntegralAdvectiveTerm");
		return 0;
	}
}
