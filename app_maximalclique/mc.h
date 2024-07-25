#pragma once

#include "../system/workerOL.h"
#include "../system/task.h"
#include "graph.h"

#include <vector>
using namespace std;

#define TIME_THRESHOLD 0.1

 using namespace std::chrono;


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

    // struct timeb start_t;
    // struct timeb end_t;

    std::chrono::time_point<std::chrono::steady_clock> start_t, end_t;
    ui max_sz = 0;

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
            cout << "-----------" << i << ", " << nbr_count << endl;
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
            // ftime(&q.start_t);
            q.start_t = std::chrono::steady_clock::now();
            return true;
        // }
    }

    
    void set_intersection(const ui *first1, const ui *last1,
                          const ui *first2, const ui *last2,
                            vector<ui> &out)
    {
        while (first1!=last1 && first2!=last2)
        {
            if (*first1<*first2) ++first1;
            else if (*first2<*first1) ++first2;
            else {
                out.push_back(*first1);
                ++first1; 
                ++first2;
            }
        }
    }

    void BK(vector<ui> &R, vector<ui> &P, vector<ui> &X, MCQuery &q)
    {
        struct timeb cur_time;
		double drun_time;

        if (P.size() == 0 && X.size() == 0)
        {
            // report R;
            // cout << "************** ";
            // for (auto x: R) cout << x << ' ';
            // cout << endl;
            q.max_sz = q.max_sz > R.size() ? q.max_sz : R.size();
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
            // vector<ui> tmp(nbrs, nbrs + nbr_count);

            vector<ui> newP;
			// set_intersection(tmp.begin(), tmp.end(), P.begin() + i, P.end(), back_inserter(newP));
            set_intersection(nbrs, nbrs+nbr_count, P.data() + i, P.data() + P.size(), newP);

            vector<ui> newX;
			// set_intersection(tmp.begin(), tmp.end(), X.begin(), X.end(), back_inserter(newX));
            set_intersection(nbrs, nbrs+nbr_count, X.data(), X.data() + X.size(), newX);

            ftime(&cur_time);
            drun_time = cur_time.time-data_graph.gtime_start[thread_id].time+(double)(cur_time.millitm-data_graph.gtime_start[thread_id].millitm)/1000;

            if(drun_time < TIME_THRESHOLD) {
                BK(newR, newP, newX, q);
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
        BK(context.R, context.P, context.X, q);
        q.counters[thread_id] = counter;
    }

    virtual bool postprocess(MCQuery &q)
    {
        
        // ftime(&q.end_t);
        q.end_t = std::chrono::steady_clock::now();
        // double totaltime = q.end_t.time-q.start_t.time+(double)(q.end_t.millitm-q.start_t.millitm)/1000;
        cout<<"Query "<<get_queryID()<<" total time (us): "<<std::chrono::duration_cast<std::chrono::microseconds>(q.end_t - q.start_t).count()<<endl;
        ULL total_results = 0;
        for(ui i=0; i<32; i++)
        {
            total_results += q.counters[i];
        }
        cout<<"Query "<<get_queryID()<<" total count: "<< total_results << ", max size: " << q.max_sz <<endl;

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