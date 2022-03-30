#include "msgtool.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "time.h"
#include <unistd.h>
using namespace std;

//usage:
//arg1: input file
//- one line per query


int main(int argc, char *argv[])
{
	//load all queries
	ifstream fin(argv[1], ios::in);
	int MAX_LINE_LENGTH=1000;
	char line[MAX_LINE_LENGTH];
	int type=1;
	vector<string> queries;
	char* pch;
	bool server_exit = false;
	while(fin.getline(line, MAX_LINE_LENGTH))
	{
		string q=line;
		if(q == "server_exit")
		{
			server_exit =true;
			break;
		}
		queries.push_back(q);
	}
	fin.close();
	//do batch processing
	int n=queries.size();
	msg_queue_client client;
	msg_queue_receiver receiver;
	double start_time=get_current_time();
	int query_num=0;
	while(query_num < n)
	{
		if(client.send_msg(type, queries[query_num].c_str()))
			query_num++;
		else
			usleep(1);
	}

	if(server_exit)
		client.send_msg(type, "server_exit");

	//all queries sent, obtain the answers
	while(query_num > 0)
	{
		while(receiver.recv_msg(type) == false);//busy waiting if no notification is received
		do
		{
			//process current notification
			char* notif=receiver.get_msg();
			pch=strtok(notif, "\n");
			cout<<"[INFO] Query "<<pch<<" is done."<<endl;
			//------
			query_num--;
		}
		while(receiver.recv_msg(type));
	}
	cout<<"Total query processing time: "<<(get_current_time()-start_time)<<" seconds"<<endl;
	//----
	return 0;
}
