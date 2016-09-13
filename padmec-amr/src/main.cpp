#include "SIMULATION_core.h"

int main(int argc, char **argv){
	/*@
	 * SIMULATION_core managers all simulation steps.
	 * Simulation is divided into three steps:
	 *   	initialization    : pointer creation, data reading, variables defining
	 *    	time-step solving : call equation functions
	 *    	stop simulation   : write output data to files, free memory
	 *
	 * Note.: During the execution of any of these three steps an error signal
	 * 		  (Exception) can be thrown. If it happens, the exception is caught
	 *        and simulation terminates after showing the cause and where it
	 *        occurred.
	 @*/
	PRS::SIMULATION_core sim;
	try{
		sim.initialize(argc,argv);
		double t1 = MPI_Wtime();
		PetscPrintf(PETSC_COMM_WORLD,"\n\nStart simulation:\n-----------------------------------------------\n");
		sim.solver();
		PetscPrintf(PETSC_COMM_WORLD,"\n\nEnd of simulation:\n----------------------------------------------\n");
		double t2 = MPI_Wtime();

		double h,m,s;
		convertSecToTime(t2-t1,&h,&m,&s);
		cout << setprecision(0) << fixed << "\n\nCPU time elapsed: " << h << "h " << m << "m " << s << "s\n\n";

		CPU_Profile::StatisticOutput(argc,argv);
	}
	catch (Exception excp) {
		excp.showExceptionMessage();
	}

	// If an exception is thrown, finalize all necessary things before terminate running.
	sim.finalize();
	return 0;
}
