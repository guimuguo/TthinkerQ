#pragma once

#include "../system/workerOL.h"
#include "../system/task.h"
#include <fstream>
#include <assert.h>
#include <climits>
#include <unordered_set>
#include "../system/TaskProgMap.h"

#include "graph.h"

#define RECORD_RESULT

typedef unsigned long long int ULL;

float TIME_THRESHOLD = 0.1;
int BIG_THRESHOLD = 100;
Graph g;

struct ContextValue
{
    int start, end, length;
    std::unordered_set<int> visited; 
    Path path;

    ContextValue() 
    {	
	}
    ~ContextValue() 
    {
	}
};


ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{
	m >> c.start;
	m >> c.end;
	m >> c.length;
    m >> c.path;

    size_t size;
    m >> size;
    for (size_t i = 0; i < size; i++) {
        int key;
        m >> key;
        c.visited.insert(key);
    }

    return m;
}

ifbinstream & operator<<(ifbinstream & m, const ContextValue & c)
{
	m << c.start;
	m << c.end;
	m << c.length;
    m << c.path;

    m << c.visited.size();
    for (auto it = c.visited.begin(); it != c.visited.end(); ++it) {
        m << *it;
    }

    return m;
}

// check if input file is empty
bool empty(ifstream &pFile)
{
    return pFile.peek() == ifstream::traits_type::eof();
}


struct CEQuery
{
    int start, end, length;
    vector<ULL> counters;

    CEQuery()
    {
        counters.assign(32, 0);
    }
};

typedef Task<ContextValue> CETask;

class CEComper: public Comper<CETask, CEQuery>
{
public:
    int counter;

    void dfs(int s, int d, int length, Path &path, std::unordered_set<int> &visited, FILE *gfpout)
    {    
        double drun_time;
        timeb cur_time;
        visited.insert(s);
        path.push_back(s);
        if(s == d) {
            // report simple path

#ifdef RECORD_RESULT
            for(int i=0; i<path.size(); i++)
                fprintf(gfpout, "%d ", path[i]);
            fprintf(gfpout, "\n");
            fflush(gfpout);
#endif

            counter++;
        } else {
            if(path.size() <= length) {
                for(int i=1; i<=g.mppadj_lists_o[s][0]; i++) {
                    int ne = g.mppadj_lists_o[s][i];
                    if(visited.find(ne) == visited.end()) {
                        ftime(&cur_time);
					    drun_time = cur_time.time-g.gtime_start[thread_id].time+(double)(cur_time.millitm-g.gtime_start[thread_id].millitm)/1000;
                        if(drun_time < TIME_THRESHOLD) {
                            dfs(ne, d, length, path, visited, gfpout);
                        } else {
                            // to split
                            CETask *t = new CETask();
                            // outdegree of ne, larger than threshold
                            t->context.start = ne; // next time start with ne
                            t->context.end = d;
                            t->context.length = length;
                            t->context.visited = visited;
                            t->context.path = path;
                            add_task(t);
                        }
                    }
                }
            }
        }
        path.pop_back();
        visited.erase(s);
    }

    virtual bool toQuery(string& line, CEQuery& q)
	{
        istringstream istrStream(line);
        istrStream>>q.start>>q.end>>q.length;
        return true;
	}

    virtual bool task_spawn(CEQuery &q)
	{
		CETask *task = new CETask();
        task->context.start = q.start;
        task->context.end = q.end;
        task->context.length = q.length;

        add_task(task);
        return true;
	}

    virtual void compute(ContextT &context, CEQuery &q)
    {   
        counter = q.counters[thread_id];
        ftime(&g.gtime_start[thread_id]);
        dfs(context.start, context.end, context.length, context.path, context.visited, gfpout);
        q.counters[thread_id] = counter;
    }

    virtual bool is_bigTask(ContextValue &context)
	{   
        if(g.mppadj_lists_o[context.start][0] >= BIG_THRESHOLD) 
            return true;
        else
		    return false;
	}

    virtual bool postprocess(CEQuery &q) 
    {
        cout<<"Query "<<get_queryID()<<" total count: "<<(ULL)accumulate(q.counters.begin(), q.counters.end(), 0)<<endl;
        return false;
    }
};

class CEWorker : public Worker<CEComper>
{
public:
    CEWorker(int num_compers) : Worker(num_compers)
    {  
    }

    ~CEWorker()
    {	
    }

    void load_data(const char* file_path)
    {
        g.loadGraphFromFile(file_path);
    }
};