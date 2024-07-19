#pragma once

#include "../system/workerOL.h"
#include "../system/task.h"
#include "graph.h"

#include <vector>
using namespace std;

#define TIME_THRESHOLD 0.1


Graph data_graph;

typedef unsigned long long int ULL;

struct ContextValue
{   
    vector<ui> R;
    vector<ui> P;
    vector<ui> X;
};

ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{
    m >> c.R;
    m >> c.P;
    m >> c.X;
    return m;
}
ifbinstream & operator<<(ifbinstream & m, const ContextValue & c) 
{
    m << c.R;
    m << c.P;
    m << c.X;
    return m;
}

struct MCQuery
{
    ui src;
    vector<ULL> counters;

    struct timeb start_t;
    struct timeb end_t;

    MCQuery()
    {
        counters.assign(32, 0);
    }
};

typedef Task<ContextValue> MCTask;

void ordered_insert(vector<ui> & newR, ui to_add)
{
	auto it = newR.begin();
	for ( ; it != newR.end(); ++it) {
		if (*it > to_add) {
			newR.insert(it, to_add);
            break;
		}
	}
	if (it == newR.end())
		newR.insert(it, to_add);
}

class MCComper: public Comper<MCTask, MCQuery>
{
public:
    ULL counter;

    MCComper() {}

    virtual bool toQuery(string& line, MCQuery& q)
    {
        q.src = stoi(line);
        ftime(&q.start_t);
        return true;
    }

    virtual bool task_spawn(MCQuery &q)
    {
        // for(ui i=0; i<data_graph.getVerticesCount(); i++)
        // {
            ui i = (ui)q.src;
            vector<ui> R {i};
            vector<ui> P, X;
            ui nbr_count;
            const ui *nbrs = data_graph.getVertexNeighbors(i, nbr_count);
            for (int j=0; j<nbr_count; ++j)
            {
                const ui neighbor = nbrs[j];
                P.push_back(neighbor);
            }
            MCTask *t = new MCTask();
            t->context.R = move(R);
            t->context.P = move(P);
            t->context.X = move(X);
            add_task(t);
            return true;
        // }
    }

    void BK(vector<ui> &R, vector<ui> &P, vector<ui> &X)
    {
        struct timeb cur_time;
		double drun_time;

        if (P.size() == 0 && X.size() == 0)
        {
            // report R;
            // cout << "************** ";
            // for (auto x: R) cout << x << ' ';
            // cout << endl;
            counter++;
            return;
        }
        for (ui i=0; i<P.size(); ++i)
        {
            ui v = P[i];
            vector<ui> newR = R;
            ordered_insert(newR, v);

            ui nbr_count;
            const ui *nbrs = data_graph.getVertexNeighbors(v, nbr_count);
            vector<ui> tmp(nbrs, nbrs + nbr_count);

            vector<ui> newP;
			set_intersection(tmp.begin(), tmp.end(), P.begin() + i, P.end(), back_inserter(newP));

            vector<ui> newX;
			set_intersection(tmp.begin(), tmp.end(), X.begin(), X.end(), back_inserter(newX));

            ftime(&cur_time);
            drun_time = cur_time.time-data_graph.gtime_start[thread_id].time+(double)(cur_time.millitm-data_graph.gtime_start[thread_id].millitm)/1000;

            if(drun_time < TIME_THRESHOLD) {
                BK(newR, newP, newX);
            } else {
                MCTask *t = new MCTask();
                t->context.R = move(newR);
                t->context.P = move(newP);
                t->context.X = move(newX);
                add_task(t);
            }
            X.push_back(v);
        }
    }

    virtual void compute(ContextT &context, MCQuery &q)
    {
        ftime(&data_graph.gtime_start[thread_id]);
        counter = q.counters[thread_id];
        BK(context.R, context.P, context.X);
        q.counters[thread_id] = counter;
    }

    virtual bool postprocess(MCQuery &q)
    {
        
        ftime(&q.end_t);
        double totaltime = q.end_t.time-q.start_t.time+(double)(q.end_t.millitm-q.start_t.millitm)/1000;
        cout<<"Query "<<get_queryID()<<" total time: "<<totaltime<<endl;
        ULL total_results = 0;
        for(ui i=0; i<32; i++)
        {
            total_results += q.counters[i];
        }
        cout<<"Query "<<get_queryID()<<" total count: "<< total_results <<endl;

        return false;
    }
};

class MCWorker : public Worker<MCComper>
{
public:
    MCWorker(ui num_compers) : Worker(num_compers)
    {
    }

    ~MCWorker()
    {
    }

    void load_data(char* file_path)
    {
        std::string fp = std::string(file_path);
        data_graph.loadGraphFromFile(fp);
    }
};