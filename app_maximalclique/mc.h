#pragma once

#include "../system/workerOL.h"
#include "../system/task.h"
#include "graph.h"

#include <vector>
using namespace std;


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
    m >> R;
    m >> P;
    m >> X;
    return m;
}
ifbinstream & operator<<(ifbinstream & m, const ContextValue & c) 
{
    m << R;
    m << P;
    m << X;
    return m;
}

struct MCQuery
{
    ui src;
    vector<ULL> counters;

    GMQuery() 
    {
        counters.assign(32, 0);
    }
};

typedef Task<ContextValue> MCTask;

void ordered_insert(vector<int> & newR, int to_add)
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
    virtual bool toQuery(string& line, MCQuery& q)
    {
        q.src = stoi(line);
        return true;
    }

    virtual bool task_spawn(GMQuery &q)
    {
        // for(ui i=0; i<data_graph.getVerticesCount(); i++)
        // {
            ui i = q.src;
            vector<ui> R {i};
            vector<ui> P, X;
            ui nbr_count;
            const ui *nbrs = data_graph.getVertexNeighbors(i, nbr_count);
            for (int j=0; j<nbr_count; ++j)
            {
                const ui neighbor = nbrs[j];
                if (nbr < i)
                    X.push_back(nbr);
                else 
                    R.push_back(nbr);
            }
            MCTask *t = new MCTask();
            t->context.R = move(R);
            t->context.P = move(P);
            t->context.X = move(X);
            add_task(t);
        // }
    }

    void BK(vector<ui> &R, vector<ui> &P, vector<ui> &X)
    {
        if (P.size() == 0 && X.size() == 0)
        {
            // report R;
            cout << "************** ";
            for (auto x: R) cout << x << ' ';
            cout << endl;
        }
        for (ui i=0; i<P.size(); ++i)
        {
            ui v = P[i];
            vector<ui> newR = R;
            ordered_insert(newR, v);

            ui nbr_count;
            const ui *nbrs = getVertexNeighbors(v, nbr_count);
            vector<int> tmp(nbrs, nbrs + nbr_count);

            vector<int> newP;
			set_intersection(tmp.begin(), tmp.end(), P.begin() + i, P.end(), back_inserter(newP));

            vector<int> newX;
			set_intersection(tmp.begin(), tmp.end(), X.begin(), X.end(), back_inserter(newX));

            BK(newR, newP, newX);

            X.push_back(v);
        }
    }

    virtual void compute(ContextT &context, MCQuery &q)
    {
        BK(context.R, context.P, context.X);
    }
};

class MCWorker : public Worker<MCComper>
{
public:
    GMWorker(ui num_compers) : Worker(num_compers)
    {
    }

    ~GMWorker()
    {
    }

    void load_data(char* file_path)
    {
        std::string fp = std::string(file_path);
        data_graph.loadGraphFromFile(fp);
    }
};
