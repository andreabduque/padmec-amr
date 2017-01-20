#ifndef EXCEPTION_H
#define EXCEPTION_H


class Exception
{
private:
	char exceptionMsg[512];

public:

	enum MSGIDs{UNKNOWN_EXIT,INIT_ERROR};


	Exception(){}
	Exception(int, const char*, const char*);
	Exception(int, const char*, MSGIDs);
	void showExceptionMessage() const;
	virtual ~Exception();


	string getMessage(MSGIDs msgId){
		string msg;
		switch (msgId){
		case INIT_ERROR:
			msg = "You must type:\n./SimAdapt_project.exe n\n\nwhere n can be:\n0 - Steady State\n1 - Transient";
			break;
		case UNKNOWN_EXIT:
			msg = "Simulation terminated by user.";
			break;
		default:
			msg = "Unknown Error Message.";
		}
		return msg;
	}
};

#endif
