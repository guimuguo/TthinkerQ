#include "tools/msgtool.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "utils/time.h"
using namespace std;

//usage:
//arg1: input file
//- one line per query
//- format: query_id query_text

//arg2: N
//- number of queries per batch for processing

//arg3: output file (optional)
//- default file: batch_out.txt

struct id_time
{
	int id;
	double time;
	double rate;

	inline bool operator<(const id_time& rhs) const
	{
		return id < rhs.id;
	}

	inline bool operator==(const id_time& rhs) const
	{
		return id == rhs.id;
	}
};

int main(int argc, char *argv[])
{
	//load all queries
	int N=atoi(argv[2]);
	int capacity=0;
	ifstream fin(argv[1], ios::in);
	int LINE_LENGTH=100;
	char line[LINE_LENGTH];
	int type=1;
	vector<string> queries;
	char* pch;
	while(fin.getline(line, LINE_LENGTH))
	{
		string q=line;
		queries.push_back(q);
	}
	fin.close();
	//do batch processing
	int n=queries.size();
	vector<id_time> results;
	msg_queue_client client;
	msg_queue_receiver receiver;
	double start_time=get_current_time();
	int nxtQ=0;
	while(nxtQ < n)
	{
		while(capacity < N && nxtQ < n)
		{
			client.send_msg(type, queries[nxtQ].c_str());
			nxtQ++;
			capacity++;
		}
		while(receiver.recv_msg(type) == false);//busy waiting if no notification is received
		do
		{
			//process current notification
			char* notif=receiver.get_msg();
			id_time entry;
			pch=strtok(notif, " ");
			entry.id=atoi(pch);
			pch=strtok(NULL, " ");
			entry.time=atof(pch);
			pch=strtok(NULL, "\n");
			entry.rate=atof(pch);
			results.push_back(entry);
			//------
			capacity--;
		}
		while(receiver.recv_msg(type));
	}
	//all queries sent, not all answers obtained
	while(capacity > 0)
	{
		while(receiver.recv_msg(type) == false);//busy waiting if no notification is received
		do
		{
			//process current notification
			char* notif=receiver.get_msg();
			id_time entry;
			pch=strtok(notif, " ");
			entry.id=atoi(pch);
			pch=strtok(NULL, " ");
			entry.time=atof(pch);
			pch=strtok(NULL, "\n");
			entry.rate=atof(pch);
			results.push_back(entry);
			//------
			capacity--;
		}
		while(receiver.recv_msg(type));
	}
	cout<<"Total query processing time: "<<(get_current_time()-start_time)<<" seconds"<<endl;
	//output results to a logfile
	char* outfile="batch_out.txt";
	if(argc > 3)
	{
		outfile=argv[3];
	}
	ofstream out(outfile);
	sort(results.begin(), results.end());
	for(int i=0; i<n; i++)
	{
		string q=queries[i];
		id_time en=results[i];
		out<<en.id<<": "<<q<<", response time "<<en.time<<" seconds, access rate = "<<en.rate<<endl;
	}
	out<<"Total query processing time: "<<(get_current_time()-start_time)<<" seconds"<<endl;
	out.close();
	//----
	return 0;
}
