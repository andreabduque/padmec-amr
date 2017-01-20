/*
 * CPU_Profiling.h
 *
 *  Created on: Sep 19, 2014
 *      Author: Rogerio Soares
 */

#ifndef CPU_PROFILE_H_
#define CPU_PROFILE_H_

#include "auxiliar.h"

/* This class provides a simple way to get CPU statistics related to a specific code block.
 *
 * See example below to understand how this class works.
 *
 * CPU_Profile::Start();		// Start clock counting
 * foo_1();						// function to be profiled
 * CPU_Profile::End("foo_1");	// Stop clock time and give a name to identify it
 *
 * CPU_Profile::Start();		// Start clock counting
 * foo_2(); 					// function to be profiled
 * CPU_Profile::End("foo_2");	// Stop clock time and give a name to identify it
 *
 * CPU_Profile::Start();		// Start clock counting
 * foo_N(); 					// function to be profiled
 * CPU_Profile::End("foo_n");	// Stop clock time and give a name to identify it
 *
 * // Stop profile procedure. Give a file name to where statistic will be printed
 * CPU_Profile::StatisticsOutPut("Filename.txt");
 *
 * Output:
 *
 * CPU profile: Displays statistics about CPU time
 * ===========================================================================================
 * Function                    CPU[%]                       h/m/s
 *
 * foo_1                        18                          0 53 21
 * foo_2                        38                          1 12 11
 * ...
 * foo_N                         2                          0 1 32
 *
 * A counter can be used to get an average CPU time after a set of function are called several times.
 */

class CPU_Profile{
public:

	CPU_Profile(){
	}
	~CPU_Profile(){}

	static void Start(){
		tic = MPI_Wtime();
	}

	static void End(string whatprofiling){
		CPUtimeList cpulist = cpumap[whatprofiling];
		cpulist.push_back(MPI_Wtime() - tic);
		cpumap[whatprofiling] = cpulist;
	}

	static void End2(string group, string whatprofiling){
		CPUtimeMap getCPUGroup = cpugroup[group];
		CPUtimeList cpulist = getCPUGroup[whatprofiling];

		cpulist.push_back(MPI_Wtime() - tic);
		getCPUGroup[whatprofiling] = cpulist;
		cpugroup[group] = getCPUGroup;
	}

	static void StatisticOutput(int argc, char** argv){
		int sizemap = (int)cpumap.size();
		if (!sizemap){
			cout << "CPU profiling statistics could not be printed.\n";
			return;
		}

		int i;
		double cumulated, h, m, s, total;
		double arrayprofile[sizemap];

		ofstream fid;
		if (argc==2){
			fid.open("StatisticOutput.txt");
		}
		else{
			fid.open(argv[2]);
		}

		i = 0;
		total = .0;
		for (CPUtimeMapIter iter1 = cpumap.begin(); iter1 != cpumap.end(); iter1++){
			CPUtimeList cpulist = (*iter1).second;
			cumulated = .0;
			for (CPUtimeListIter iter2 = cpulist.begin(); iter2 != cpulist.end(); iter2++){
				cumulated += *iter2;
			}
			int sizelist = (int)cpulist.size();
			arrayprofile[i] = (double)cumulated/sizelist; // take average time
			total += (double)arrayprofile[i];
			i++;
		}
		cout << total << endl;

		fid << "CPU profile: Displays statistics about CPU time\n";
		fid << "===========================================================================================\n";
		fid << "Function                    CPU[%]                       h/m/s\n\n";

		const int WIDTH = 33;
		int linewidth;
		i = 0;
		string funcname;
		for (CPUtimeMapIter iter1 = cpumap.begin(); iter1 != cpumap.end(); iter1++){
			convertSecToTime(arrayprofile[i],&h,&m,&s);
			funcname = (*iter1).first;
			linewidth = WIDTH - (int)funcname.length();	// (*iter1).first.size() = size of function name string
			fid << (*iter1).first <<  setw(linewidth) << setprecision(2) << fixed << 100.0*arrayprofile[i]/total << "\t\t\t";
			fid << setprecision(0) << h << " " << m << " " << s << endl;
			i++;
			funcname.clear();
		}
		cpumap.clear();
	}

	static void StatisticOutput2(string group, string filename){
		int sizemap = (int)cpugroup.size();
		if (!sizemap){
			cout << "CPU profiling statistics could not be printed.\n";
			return;
		}


		double cumulated, h, m, s, total;
		int i;

		CPUtimeMap getCPUGroup = cpugroup[group];
		sizemap = (int)getCPUGroup.size();
		double arrayprofile[sizemap];

		ofstream fid;
		fid.open(filename.c_str());

		i = 0;
		total = .0;
		for (CPUtimeMapIter iter1 = getCPUGroup.begin(); iter1 != getCPUGroup.end(); iter1++){
			CPUtimeList cpulist = (*iter1).second;
			cumulated = .0;
			//cout << (*iter1).first << ":\n";
			for (CPUtimeListIter iter2 = cpulist.begin(); iter2 != cpulist.end(); iter2++){
				cumulated += *iter2;
				//cout << cumulated << "\n";
			}
			int sizelist = (int)cpulist.size();
			//cout << "total: " << total << "\tsizelist: " << sizelist << endl;
			arrayprofile[i] = (double)cumulated/sizelist; // take average time
			total += (double)arrayprofile[i];
			//cout << "avg: " << arrayprofile[i] << endl;
			i++;
		}

		fid << "CPU profile: Displays statistics about CPU time\n";
		fid << "===========================================================================================\n";
		fid << "Function                    CPU[%]                       h/m/s\n\n";

		const int WIDTH = 31;
		int linewidth;
		i = 0;
		string funcname;
		for (CPUtimeMapIter iter1 = getCPUGroup.begin(); iter1 != getCPUGroup.end(); iter1++){
			//cout << "arrayprofile[i]: " << arrayprofile[i] << endl;
			convertSecToTime(arrayprofile[i],&h,&m,&s);
			funcname = (*iter1).first;
			linewidth = WIDTH - (int)funcname.length();	// (*iter1).first.size() = size of function name string
			fid << (*iter1).first <<  setw(linewidth) << setprecision(2) << fixed << 100.0*arrayprofile[i]/total << "%\t\t\t";
			fid << setprecision(0) << h << " " << m << " " << s << endl;
			i++;
			funcname.clear();
		}
		//cout << "Statistics were cleaned from memory and saved in file.\n";
		getCPUGroup.clear();
	}

	typedef std::list<double> CPUtimeList;
	typedef std::list<double>::iterator CPUtimeListIter;
	typedef std::map<string,CPUtimeList> CPUtimeMap;
	typedef std::map<string,CPUtimeList>::iterator CPUtimeMapIter;
	typedef std::map<string,CPUtimeMap> MapGroup;

	static MapGroup cpugroup;						// allow multiple profiles
	static CPUtimeMap cpumap;						// store all CPU time for every time _whatprofiling is called
	static double tic;								// initial counting
};


#endif /* CPU_PROFILING_H_ */
