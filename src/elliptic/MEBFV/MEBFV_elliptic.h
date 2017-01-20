/*
 * MEBFV_elliptic.h
 *
 *  Created on: Oct 23, 2014
 *      Author: rogerio
 */

#ifndef MEBFV_ELLIPTIC_H_
#define MEBFV_ELLIPTIC_H_

#include "Elliptic_equation.h"
#include "GeomData.h"
#include "CPU_Profiling.h"

/*
 * class implementation for calculation of elliptic equation:
 * 				div(-K.grad_p) = f
 *
 * 	It's a vertex centered element based finite volume formulation for highly heterogeneous porous media.
 */

namespace PRS{
	class MEBFV_elliptic : public Elliptic_equation{
	public:
		MEBFV_elliptic();
		MEBFV_elliptic(pMesh, PhysicPropData *, SimulatorParameters *, GeomData *, MeshData *);
		~MEBFV_elliptic();

		double solver(pMesh);	// virtual function implementation for solving Ax=b
		void Initialize(pMesh);		// allocate memory for all matrices and vector needed by elliptic term
		void setSST();			// fill SST vector with source/sink terms (Neumann boundary values)
		void setdirichletVec();	// fill vector with prescribed values for corresponding node index
		void Assembly_A();		// assembly global stiffness matrix
		void Assembly_b();		// assembly rhs vector
		double Solve();			// solves Ax=b

		void freeMemory();		// free memory from all matrices and vectors

	private:

		Mat A;				// Stiffness matrix
		Mat A_free;			// Global matrix for free vertices
		Vec b;				// Right hand side (rhs) vector
		Vec x;				// Solution vector

		Mat dirichletMat;	// matrix contribution from prescribed nodes
		Vec dirichletVec;	// vector to be added to rhs vector: rhs = SST - dirichletVec
		Vec SST;			// Source Sink Terms vector for rhs

		bool initialize;	// set variation

		PhysicPropData *pPPData;
		SimulatorParameters *pSimPar;
		GeomData *pGCData;
		MeshData *pMData;
		pMesh theMesh;
	};
}
#endif /* MEBFV_ELLIPTIC_H_ */
