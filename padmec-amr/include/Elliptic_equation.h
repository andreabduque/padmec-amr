#ifndef _ELLIPTIC_EQUATION_H_
#define _ELLIPTIC_EQUATION_H_

#include "exportVTK.h"

/**
 * Base class for elliptic solvers.
 */
namespace PRS           // PRS: Petroleum Reservoir Simulator
{
	class Elliptic_equation{
	public:

		Elliptic_equation(){}
		virtual ~Elliptic_equation(){}
		virtual double solver(pMesh)=0;
		//virtual void getCPUtime(double &assembly, double &solver, double &gradient, int &KSPiter);

	protected:
		Vec output;

		double calculatePressureGradient();

		// solve system of equations: Ax=y
		double KSP_solver(Mat A, Mat pcMatrix, Vec y, Vec x, SimulatorParameters *pSimPar, PetscBool guessNonZero, KSPType ksptype, PCType pctype,PetscInt &its){
			double startt = MPI_Wtime();
			KSP ksp;
			PC preconditioner;
			KSPCreate(PETSC_COMM_WORLD,&ksp);
			KSPSetOperators(ksp,A,pcMatrix);
			KSPSetType(ksp,ksptype);
			KSPGetPC(ksp,&preconditioner);
			PCSetType(preconditioner,pctype);
			KSPSetInitialGuessNonzero(ksp,guessNonZero);
			KSPGMRESSetRestart(ksp,pSimPar->getKrylov_restart());
			KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
			KSPSetFromOptions(ksp);
			KSPSolve(ksp,y,x);
			KSPGetIterationNumber(ksp,&its);
			KSPDestroy(&ksp);
			double endt = MPI_Wtime();
			return endt-startt;
		}
	};
}
#endif
