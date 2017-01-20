#include "includes.h"

Exception::Exception(const int i, const char *f, const char *text){
	string str1(f);
	int loc = str1.find_last_of("/");
	string str2 = str1.substr(loc+1, str1.size()-loc);
	sprintf(exceptionMsg,"\n\nAn exception has been thrown by rank %d\nFile:\t%s\n"
						 "Line:\t%d\n"
						 "Cause:\t%s\n",P_pid(),str2.c_str(),i,text);

}

Exception::Exception(int line, const char *f, MSGIDs msgId){
	string str1(f);
	int loc = str1.find_last_of("/");
	string file = str1.substr(loc+1, str1.size()-loc);
	string msg = getMessage(msgId);
	sprintf(exceptionMsg,"\n\nAn exception has been thrown by rank %d\nFile:\t%s\n"
			"Line:\t%d\n"
			"Cause:\t%s\n",P_pid(),file.c_str(),line,msg.c_str());
}

void Exception::showExceptionMessage() const{
	// print message on screen
	std::cout << exceptionMsg << std::endl;
	// print message on file
	ofstream fid;
	fid.open("ExceptionThrown.log");
	fid << "Application has been terminated by an exception throw:\n\n";
	fid << exceptionMsg << std::endl;
	MPI_Abort( MPI_COMM_WORLD, 190 );
}

Exception::~Exception(){
}
