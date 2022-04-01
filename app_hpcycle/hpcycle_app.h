#pragma once

#define MAX_STAGE_NUM 4
  
#define STAGE_1 1 // dfs in G -> LeftTable
#define STAGE_2 2 // dfs in G -> RightTable
#define STAGE_3 3 // left+right Join
#define STAGE_4 4 // dfs in G_idx, left+mid+right Join

#define CASE_1 1
#define CASE_2 2
#define CASE_3 3
#define CASE_4 4

// only used for STAGE_4
#define TYPE_1 1
#define TYPE_2 2

#define TIME_THRESHOLD 0.1
#define NPATH_THRESH 500000
#define RECORD_RESULT

#include "../system/workerOL.h"
#include "../system/task.h"
#include <fstream>
#include <assert.h>
#include <climits>
#include <algorithm>
#include <boost/functional/hash.hpp>

#include "../system/TaskProgMap.h"
#include "../system/conmap.h"
#include "graph.h"

typedef typename conmap<int, PathVec>::bucket bucket;
typedef typename conmap<int, PathVec>::bucket_map bucket_map;

typedef unsigned long long int ULL;

Graph g;

bool comp_plen(Path p, Path q)
{
    return p.size()<q.size();
}

struct Range
{
    int lstart, lend, rstart, rend;
};

struct JoinWrapper
{   
    vector<pair<int,int> > range_sd; // {lhp:rhp}: Range
    vector<Range> range_map;

    vector<HPEdgeSet> mid_paths;
    vector<int> total_weights;

    int capacity()
    {   
        int capacity = 0;
        for(auto& range: range_map) {
            capacity += (range.lend-range.lstart)*(range.rend-range.rstart);
        }
        return capacity;
    }

    void clear()
    {
        range_sd.clear();
        range_map.clear();
        mid_paths.clear();
        total_weights.clear();
    }

    bool empty()
    {
        return range_map.empty();
    }
};

// adjusted length, threshold length
struct ContextValue
{
    int start, cur, end, adlength, thlength, total_weight, direction;
    hash_set<int> visited;

    Path gpath;
    HPEdgeSet path;

    bool case3_spawn_flag;
    bool case4_spawn_flag;

    JoinWrapper joinwrapper;
    int type;

    ContextValue() 
    {
        case4_spawn_flag = false;
        case3_spawn_flag = false;
	}
    ~ContextValue() 
    {
	}
};

ofbinstream & operator>>(ofbinstream & m, Range & c)
{
    m >> c.lstart;
    m >> c.lend;
    m >> c.rstart;
    m >> c.rend;
    return m;
}
ifbinstream & operator<<(ifbinstream & m, const Range & c)
{   
    m << c.lstart;
    m << c.lend;
    m << c.rstart;
    m << c.rend;
    return m;
}

ofbinstream & operator>>(ofbinstream & m, pair<int,int> & c)
{
    m >> c.first;
    m >> c.second;
    return m;
}
ifbinstream & operator<<(ifbinstream & m, const pair<int,int> & c)
{
    m << c.first;
    m << c.second;
    return m;
}

ofbinstream & operator>>(ofbinstream & m, JoinWrapper & c) 
{   
    m >> c.range_sd;
    m >> c.range_map;
    m >> c.mid_paths;
    m >> c.total_weights;

    return m;
}
ifbinstream & operator<<(ifbinstream & m, const JoinWrapper & c)
{
    m << c.range_sd;
    m << c.range_map;
    m << c.mid_paths;
    m << c.total_weights;

    return m;
}

ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{   
	m >> c.start;
    m >> c.cur;
	m >> c.end;
	m >> c.adlength;
    m >> c.thlength;
    m >> c.total_weight;
    m >> c.direction;
	m >> c.visited;
    m >> c.gpath;
    m >> c.path;

    m >> c.case4_spawn_flag;
    m >> c.case3_spawn_flag;
    m >> c.joinwrapper;

    m >> c.type;

    return m;
}

ifbinstream & operator<<(ifbinstream & m, const ContextValue & c)
{
	m << c.start;
    m << c.cur;
	m << c.end;
	m << c.adlength;
    m << c.thlength;
    m << c.total_weight;
    m << c.direction;
	m << c.visited;
    m << c.gpath;
    m << c.path;

    m << c.case4_spawn_flag;
    m << c.case3_spawn_flag;
    m << c.joinwrapper;

    m << c.type;

    return m;
}

// check if input file is empty
bool empty(ifstream &pFile)
{
    return pFile.peek() == ifstream::traits_type::eof();
}

bool check_hppath_validity(HPPath &hp_paths, vector<bool> &visited)
{
    for(int i=0; i<hp_paths.hppath.size(); i++) {
        for(int j=0; j<hp_paths.hppath[i].path.size(); j++) {
            int ve = hp_paths.hppath[i].path[j];
            if(i!=0 && j==0) continue;
             if(visited[ve]) {
                return false;
            }
            visited[ve] = true;
        }
    }
    return true;
}

// hot-point index cycle-enum query
struct HCEQuery
{   
    // given by user
    int start, end, length;

    struct timeb start_t, end_t;

    // ====== resource held for each query ====
    typedef conmap<int, PathVec> TableType; // i.e. concurrent PathTable
    TableType ltable, rtable;
    
    int left_min, right_min;
    mutex left_mtx, right_mtx;

    int cur_stage_num, category;

    vector<ULL> counters;
    // =======================================

    HCEQuery() 
    {
        left_min = g.mnum_of_vertices;
        right_min = g.mnum_of_vertices;

        cur_stage_num = STAGE_1;

        category = -1;

        counters.assign(32, 0);
    }
};

typedef Task<ContextValue> HCETask;

class HCEComper: public Comper<HCETask, HCEQuery>
{
public:

    ULL counter;
    vector<bool> cvisited, cvisited2;

    HCEComper()
    {
        cvisited.assign(g.mnum_of_vertices, false);
        cvisited2.assign(g.mnum_of_vertices, false);
    }
    
    void dfs_gidx_inside_1(int s, int cur, hash_set<int> &visited, int adlength, int thlength, 
                            int total_weight, HPEdgeSet &path, int direction, JoinWrapper &jw,
                            HCEQuery &q)
    {
        double drun_time;
        timeb cur_time;

        visited.insert(cur);

        bucket & rbucket = q.rtable.get_bucket(cur);
    	bucket_map & rkvmap = rbucket.get_map();

        if(s!=cur && rkvmap.find(cur)!=rkvmap.end()) { // prevent self-loop
            // case 1
            bucket & lbucket = q.ltable.get_bucket(s);
            bucket_map & lkvmap = lbucket.get_map();

            int lsize = lkvmap[s].size(), rsize = rkvmap[cur].size();
            int batch_sz, nbatch, nlpath, nrpath, lastend, cap, leftover;
            Range rg;

            // Step 1: merge
            if(lsize < rsize) 
            {
                cap = jw.capacity();
                leftover = NPATH_THRESH-cap;
                nrpath = leftover%lsize==0? leftover/lsize: leftover/lsize+1;
                nrpath = nrpath>rsize?rsize:nrpath;
                rg.lstart = 0;
                rg.lend = lsize;
                rg.rstart = 0;
                rg.rend = nrpath;
                if(nrpath > 0) {
                    jw.range_sd.push_back({s,cur});
                    jw.range_map.push_back(rg);
                    jw.mid_paths.push_back(path);
                    jw.total_weights.push_back(total_weight);
                }
                if(rg.rend == rsize) {
                    if(jw.capacity() >= NPATH_THRESH) {
                        HCETask *task = new HCETask();
                        task->context.thlength = thlength;
                        task->context.joinwrapper = jw;
                        task->context.case4_spawn_flag = true;
                        task->context.type = TYPE_2;
                        add_task(task);
                        jw.clear();
                    }
                    // else continue adding
                }
                else {
                    HCETask *task = new HCETask();
                    task->context.thlength = thlength;
                    task->context.joinwrapper = jw;
                    task->context.case4_spawn_flag = true;
                    task->context.type = TYPE_2;
                    add_task(task);
                    jw.clear();
                
                    batch_sz = NPATH_THRESH%lsize==0? NPATH_THRESH/lsize: NPATH_THRESH/lsize+1;
                    nbatch = (rsize-nrpath)/batch_sz;
                    lastend = nrpath;
                    for(int j=0; j<nbatch; j++) {
                        rg.lstart = 0;
                        rg.lend = lsize;
                        rg.rstart = j*batch_sz+nrpath;
                        lastend = rg.rend = (j+1)*batch_sz+nrpath;
                        jw.range_sd.push_back({s,cur});
                        jw.range_map.push_back(rg);
                        jw.mid_paths.push_back(path);
                        jw.total_weights.push_back(total_weight);
                        
                        // Step 2: add a TYPE_2 task
                        HCETask *t = new HCETask();
                        t->context.thlength = thlength;
                        t->context.joinwrapper = jw;
                        t->context.case4_spawn_flag = true;
                        t->context.type = TYPE_2;
                        add_task(t);
                        jw.clear();
                    }
                    rg.lstart = 0;
                    rg.lend = lsize;
                    rg.rstart = lastend;
                    rg.rend = rsize;
                    if(lastend != rsize) {
                        jw.range_sd.push_back({s,cur});
                        jw.range_map.push_back(rg);
                        jw.mid_paths.push_back(path);
                        jw.total_weights.push_back(total_weight);
                    }
                }
            } 
            else {  // rsize<lsize
                cap = jw.capacity();
                leftover = NPATH_THRESH-cap;
                nlpath = leftover%rsize==0? leftover/rsize: leftover/rsize+1;
                nlpath = nlpath>lsize?lsize:nlpath;
                rg.lstart = 0;
                rg.lend = nlpath;
                rg.rstart = 0;
                rg.rend = rsize;
                if(nlpath > 0) {
                    jw.range_sd.push_back({s,cur});
                    jw.range_map.push_back(rg);
                    jw.mid_paths.push_back(path);
                    jw.total_weights.push_back(total_weight);
                }
                if(rg.lend == lsize) {
                    if(jw.capacity() >= NPATH_THRESH) {
                        HCETask *task = new HCETask();
                        task->context.thlength = thlength;
                        task->context.joinwrapper = jw;
                        task->context.case4_spawn_flag = true;
                        task->context.type = TYPE_2;
                        add_task(task);
                        jw.clear();
                    }
                    // else continue adding
                }
                else {
                    HCETask *task = new HCETask();
                    task->context.thlength = thlength;
                    task->context.joinwrapper = jw;
                    task->context.case4_spawn_flag = true;
                    task->context.type = TYPE_2;
                    add_task(task);
                    jw.clear();
                
                    batch_sz = NPATH_THRESH%rsize==0? NPATH_THRESH/rsize: NPATH_THRESH/rsize+1;
                    nbatch = (lsize-nlpath)/batch_sz;
                    lastend = nlpath;
                    for(int j=0; j<nbatch; j++) {
                        rg.lstart = j*batch_sz+nlpath;
                        lastend = rg.lend = (j+1)*batch_sz+nlpath;
                        rg.rstart = 0;
                        rg.rend = rsize;
                        jw.range_sd.push_back({s,cur});
                        jw.range_map.push_back(rg);
                        jw.mid_paths.push_back(path);
                        jw.total_weights.push_back(total_weight);

                        // Step 2: add a TYPE_2 task
                        HCETask *t = new HCETask(); 
                        t->context.thlength = thlength;
                        t->context.joinwrapper = jw;
                        t->context.case4_spawn_flag = true;
                        t->context.type = TYPE_2;
                        add_task(t);
                        jw.clear();
                    }
                    rg.lstart = lastend;
                    rg.lend = lsize;
                    rg.rstart = 0;
                    rg.rend = rsize;
                    if(lastend != lsize) {
                        jw.range_sd.push_back({s,cur});
                        jw.range_map.push_back(rg);
                        jw.mid_paths.push_back(path);
                        jw.total_weights.push_back(total_weight);
                    }
                }
            }
        }

        int new_weight;
        if(direction == OUT) {
            for(auto& ne: g.mhp_adjlist[cur]) {
                if(visited.find(ne.v) == visited.end()) {
                    path.push_back(ne);
                    if(total_weight == 0)
                        new_weight = ne.weight();
                    else
                        new_weight = total_weight+ne.weight()-1;

                    if(new_weight <= adlength) {
                        ftime(&cur_time);
					    drun_time = cur_time.time-g.gtime_start[thread_id].time+(double)(cur_time.millitm-g.gtime_start[thread_id].millitm)/1000;
                        if(drun_time < 10*TIME_THRESHOLD) {
                            dfs_gidx_inside_1(s, ne.v, visited, adlength, thlength, new_weight, path, OUT, jw, q); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
                        } else {
                            // to split
                            HCETask *t = new HCETask();
                            t->context.start = s;
                            t->context.cur = ne.v;
                            t->context.thlength = thlength;
                            t->context.adlength = adlength;
                            t->context.total_weight = new_weight;
                            t->context.path = path;
                            t->context.direction = OUT;
                            t->context.visited = visited;
                            t->context.case4_spawn_flag = true; // !!!
                            t->context.type = TYPE_1;
                            add_task(t);
                        }
                        path.pop_back();
                    } else {
                        path.pop_back();
                        break; 
                    }
                }
            }
        } else {
            for(auto& ne: g.mhp_adjlist_i[cur]) {
                if(visited.find(ne.v) == visited.end()) {
                    path.push_back(ne);
                    if(total_weight == 0)
                        new_weight = ne.weight();
                    else
                        new_weight = total_weight+ne.weight()-1;
                    if(new_weight <= adlength) {
                        dfs_gidx_inside_1(s, ne.v, visited, adlength, thlength, new_weight, path, IN, jw, q); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
                        path.pop_back();
                    } else {
                        path.pop_back();
                        break;
                    }
                }
            }
        }
        visited.erase(cur);
    }

    void dfs_gidx_inside_2(int s, int cur, hash_set<int> &visited, 
                            int adlength, int thlength, int total_weight, HPEdgeSet &path, 
                            int direction, FILE *gfpout, JoinWrapper &jw, HCEQuery &q) 
    {
        double drun_time;
        struct timeb cur_time;

        visited.insert(cur);

        bucket & bucket = q.ltable.get_bucket(cur);
        bucket_map & kvmap = bucket.get_map();

        if(s!=cur && kvmap.find(cur)!=kvmap.end()) {

            int size = kvmap[cur].size();
            int cap, npath, batch_sz, nbatch, lastend;
            Range rg;

            cap = jw.capacity();
            npath = NPATH_THRESH-cap;
            npath = npath>size?size:npath;
            rg.lstart = 0;
            rg.lend = npath;
            rg.rstart = 0;
            rg.rend = 1;

            if(npath > 0) {
                jw.range_sd.push_back({cur,cur});
                jw.range_map.push_back(rg);
                jw.mid_paths.push_back(path);
                jw.total_weights.push_back(total_weight);
            }
            if(rg.lend == size) {
                if(jw.capacity() >= NPATH_THRESH) {
                    HCETask *task = new HCETask();
                    task->context.thlength = thlength;
                    task->context.joinwrapper = jw;
                    task->context.type = TYPE_2;
                    task->context.direction = direction;
                    add_task(task);
                    jw.clear();
                }
            } else {
                HCETask *task = new HCETask();
                task->context.thlength = thlength;
                task->context.joinwrapper = jw;
                task->context.type = TYPE_2;
                task->context.direction = direction;
                add_task(task);
                jw.clear();
                
                batch_sz = NPATH_THRESH;
                nbatch = (size-npath)/batch_sz;
                lastend = npath;
                for(int j=0; j<nbatch; j++) {
                    rg.lstart = j*batch_sz+npath;
                    lastend = rg.lend = (j+1)*batch_sz+npath;
                    rg.rstart = 0;
                    rg.rend = 1;

                    jw.range_sd.push_back({cur,cur});
                    jw.range_map.push_back(rg);
                    jw.mid_paths.push_back(path);
                    jw.total_weights.push_back(total_weight);

                    // Step 2: add a TYPE_2 task
                    HCETask *t = new HCETask();
                    t->context.thlength = thlength;
                    t->context.joinwrapper = jw;
                    t->context.type = TYPE_2;
                    t->context.direction = direction;
                    add_task(t);
                    jw.clear();
                }

                rg.lstart = lastend;
                rg.lend = size;
                rg.rstart = 0;
                rg.rend = 1;
                if(lastend != size) {
                    jw.range_sd.push_back({cur,cur});
                    jw.range_map.push_back(rg);
                    jw.mid_paths.push_back(path);
                    jw.total_weights.push_back(total_weight);
                }
            }
        }
        
        int new_weight;
        for(auto &ne: g.mhp_adjlist_i[cur]) {
            if(visited.find(ne.v) == visited.end()) {
                path.push_back(ne);
                if(total_weight == 0)
                    new_weight = ne.weight();
                else
                    new_weight = total_weight+ne.weight()-1;

                if(new_weight <= adlength) {
                    ftime(&cur_time);
                    drun_time = cur_time.time-g.gtime_start[thread_id].time+(double)(cur_time.millitm-g.gtime_start[thread_id].millitm)/1000;
                    if(drun_time < 10*TIME_THRESHOLD) {
                        dfs_gidx_inside_2(s, ne.v, visited, adlength, thlength, new_weight, path, IN, gfpout, jw, q); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
                    } else {
                        // to split
                        HCETask *t = new HCETask();
                        t->context.start = s;
                        t->context.cur = ne.v;
                        t->context.thlength = thlength;
                        t->context.adlength = adlength;
                        t->context.total_weight = new_weight;
                        t->context.path = path;
                        t->context.direction = IN;
                        t->context.visited = visited;
                        t->context.type = TYPE_1;
                        add_task(t);
                    }
                    path.pop_back();
                } else {
                    path.pop_back();
                    break;
                }
            }
        }
        visited.erase(cur);
    }

    void dfs_gidx_inside_3(int s, int cur, hash_set<int> &visited, 
                            int adlength, int thlength, int total_weight, HPEdgeSet &path, 
                            int direction, FILE *gfpout, JoinWrapper &jw, HCEQuery &q)
    {   
        double drun_time;
        struct timeb cur_time;

        visited.insert(cur);

        bucket & bucket = q.rtable.get_bucket(cur);
        bucket_map & kvmap = bucket.get_map();
        
        if(s!=cur && kvmap.find(cur)!=kvmap.end()) { // prevent self-loop
            
            int size = kvmap[cur].size();
            int cap, npath, batch_sz, nbatch, lastend;
            Range rg;

            cap = jw.capacity();
            npath = NPATH_THRESH-cap;
            npath = npath>size?size:npath;
            
            rg.lstart = 0;
            rg.lend = 1;
            rg.rstart = 0;
            rg.rend = npath;
            
            if(npath > 0) {
                jw.range_sd.push_back({cur,cur});
                jw.range_map.push_back(rg);
                jw.mid_paths.push_back(path);
                jw.total_weights.push_back(total_weight);
            }
            
            if(rg.rend == size) {
                if(jw.capacity() >= NPATH_THRESH) {
                    HCETask *task = new HCETask();
                    task->context.thlength = thlength;
                    task->context.joinwrapper = jw;
                    task->context.type = TYPE_2;
                    task->context.direction = direction;
                    add_task(task);
                    jw.clear();
                }
            }
            else {
                HCETask *task = new HCETask();
                task->context.thlength = thlength;
                task->context.joinwrapper = jw;
                task->context.type = TYPE_2;
                task->context.direction = direction;
                add_task(task);
                jw.clear();

                batch_sz = NPATH_THRESH;
                nbatch = (size-npath)/batch_sz;
                lastend = npath;
                for(int j=0; j<nbatch; j++) {

                    rg.lstart = 0;
                    rg.lend = 1;
                    rg.rstart = j*batch_sz+npath;
                    lastend = rg.rend = (j+1)*batch_sz+npath;
                    jw.range_sd.push_back({cur,cur});
                    jw.range_map.push_back(rg);
                    jw.mid_paths.push_back(path);
                    jw.total_weights.push_back(total_weight);

                    // Step 2: add a TYPE_2 task
                    HCETask *t = new HCETask();
                    t->context.thlength = thlength;
                    t->context.joinwrapper = jw;
                    t->context.type = TYPE_2;
                    t->context.direction = direction;
                    add_task(t);
                    jw.clear();
                }

                rg.lstart = 0;
                rg.lend = 1;
                rg.rstart = lastend;
                rg.rend = size;
                if(lastend != size) {
                    jw.range_sd.push_back({cur,cur});
                    jw.range_map.push_back(rg);
                    jw.mid_paths.push_back(path);
                    jw.total_weights.push_back(total_weight);
                }
            }
        }
        int new_weight;
        if(direction == OUT) {
            for(auto& ne: g.mhp_adjlist[cur]) {
                if(visited.find(ne.v) == visited.end()) {
                    path.push_back(ne);
                    if(total_weight == 0)
                        new_weight = ne.weight();
                    else
                        new_weight = total_weight+ne.weight()-1;

                    if(new_weight <= adlength) {
                        ftime(&cur_time);
					    drun_time = cur_time.time-g.gtime_start[thread_id].time+(double)(cur_time.millitm-g.gtime_start[thread_id].millitm)/1000;
                        if(drun_time < 10*TIME_THRESHOLD) {
                            dfs_gidx_inside_3(s, ne.v, visited, adlength, thlength, new_weight, path, OUT, gfpout, jw, q); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
                        } else {
                            // to split
                            HCETask *t = new HCETask();
                            t->context.start = s;
                            t->context.cur = ne.v;
                            t->context.thlength = thlength;
                            t->context.adlength = adlength;
                            t->context.total_weight = new_weight;
                            t->context.path = path;
                            t->context.direction = OUT;
                            t->context.visited = visited;
                            t->context.type = TYPE_1;
                            add_task(t);
                        }
                        path.pop_back();
                    } else {
                        path.pop_back();
                        break; 
                    }
                }
            }
        } 
        visited.erase(cur);
    }

    void dfs_gidx_inside_4(int s, int d, hash_set<int> &visited, 
                            int adlength, int thlength, int total_weight, HPEdgeSet &path, 
                            FILE *gfpout, HCEQuery &q) 
    {      
        double drun_time;
        struct timeb cur_time;

        visited.insert(s);
        if(s == d) {
            HPPath hpp(path, total_weight);
            bool valid = check_hppath_validity(hpp, cvisited);    
            if(valid) {
                // satisfy length requirement
                if(total_weight<=thlength) {

#ifdef RECORD_RESULT
                    fprintf(gfpout, "@ "); // * means 2 Join
                    for(int i=0; i<hpp.hppath.size(); i++) {
                        for(int j=0; j<hpp.hppath[i].path.size(); j++) {
                            if(i!=0 && j==0) continue; 
                            fprintf(gfpout, "%d ", hpp.hppath[i].path[j]);
                        }
                    }
                    fprintf(gfpout, "\n");
#endif

                    counter++;
                }
            }
            for(auto& hpedge: hpp.hppath) {
                for(auto& ve: hpedge.path) {
                    cvisited[ve] = false;
                }
            }
        } else {
            int new_weight;
            for(auto& ne: g.mhp_adjlist[s]) {
                if(visited.find(ne.v) == visited.end()) {
                    path.push_back(ne);
                    if(total_weight == 0)
                        new_weight = ne.weight();
                    else 
                        new_weight = total_weight+ne.weight()-1;
                    if(new_weight <= adlength) {
                        ftime(&cur_time);
					    drun_time = cur_time.time-g.gtime_start[thread_id].time+(double)(cur_time.millitm-g.gtime_start[thread_id].millitm)/1000;
                        if(drun_time < 10*TIME_THRESHOLD) {
                            dfs_gidx_inside_4(ne.v, d, visited, adlength, thlength, new_weight, path, gfpout, q); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
                        } else {
                            // to split
                            HCETask *t = new HCETask();
                            t->context.cur = ne.v;
                            t->context.end = d;
                            t->context.thlength = thlength;
                            t->context.adlength = adlength;
                            t->context.total_weight = new_weight;
                            t->context.path = path;
                            t->context.visited = visited;
                            add_task(t);
                        }
                        path.pop_back();
                    } else {
                        path.pop_back();
                        break; 
                    }
                }
            }  
        }
        visited.erase(s);
    }


    void dfs_with_hp(int s, int d, int length, hash_set<int> &visited, Path &path, FILE *gfpout, 
                        int direction, bool bwrite, HCEQuery &q)
    {   
        double drun_time;
        timeb cur_time;

        visited.insert(s);
        path.push_back(s);
        
        if(s == d && bwrite) {
            // report simple path, directly path find, w/o HP-index

#ifdef RECORD_RESULT
            for(int i=0; i<path.size(); i++)
            {
                fprintf(gfpout, "%d ", path[i]);
            }
            fprintf(gfpout, "\n");
#endif

            counter++;

        } else {
            if(path.size() <= length) {
                if(direction == OUT) {
                    for(int i=1; i<=g.mppadj_lists_o[s][0]; i++) {
                        int ne = g.mppadj_lists_o[s][i];
                        // check whether ne is HP
                        if(g.mphp[ne]) {

                            bucket & bucket = q.ltable.get_bucket(ne);
                            bucket.lock();
                            bucket_map & kvmap = bucket.get_map(); //kvmap is a hash_map
                            auto it = kvmap.find(ne);
                            if(it == kvmap.end()) {
                                kvmap[ne] = PathVec(); // same as kvmap.insert
                            }

                            kvmap[ne].push_back(path);

                            bucket.unlock();
                            // tracks minimum length
                            q.left_mtx.lock();
                            if(path.size()<q.left_min) {
                                q.left_min = path.size();
                            }
                            q.left_mtx.unlock();

                            continue;
                        }
                        if(visited.find(ne) == visited.end()) {
                            ftime(&cur_time);
					        drun_time = cur_time.time-g.gtime_start[thread_id].time+(double)(cur_time.millitm-g.gtime_start[thread_id].millitm)/1000;
                            if(drun_time < TIME_THRESHOLD)  {
                                dfs_with_hp(ne, d, length, visited, path, gfpout, OUT, bwrite, q);
                            } else {
                                // to split
                                HCETask *t = new HCETask();
                                t->context.start = ne;
                                t->context.end = d;
                                t->context.thlength = length;
                                t->context.visited = visited;
                                t->context.gpath = path;
                                add_task(t);
                            }
                        }
                    }
                } else {
                    for(int i=1; i<=g.mppadj_lists_i[s][0]; i++) {
                        int ne = g.mppadj_lists_i[s][i];
                        // check whether ne is HP
                        if(g.mphp[ne]) {

                            bucket & bucket = q.rtable.get_bucket(ne);
                            bucket.lock();
                            bucket_map & kvmap = bucket.get_map(); //kvmap is a hash_map
                            auto it = kvmap.find(ne);
                            if(it == kvmap.end()) {
                                kvmap[ne] = PathVec(); // same as kvmap.insert
                            }
                            kvmap[ne].push_back(path);

                            bucket.unlock();
                            // tracks minimum length
                            q.right_mtx.lock();
                            if(path.size()<q.right_min) {
                                q.right_min = path.size();
                            }
                            q.right_mtx.unlock();

                            continue;
                        }
                        if(visited.find(ne) == visited.end()) {
                            ftime(&cur_time);
					        drun_time = cur_time.time-g.gtime_start[thread_id].time+(double)(cur_time.millitm-g.gtime_start[thread_id].millitm)/1000;
                            if(drun_time < TIME_THRESHOLD)  {
                                dfs_with_hp(ne, d, length, visited, path, gfpout, IN, bwrite, q);
                            } else {
                                // to split
                                HCETask *t = new HCETask();
                                t->context.start = ne;
                                t->context.end = d;
                                t->context.thlength = length;
                                t->context.visited = visited;
                                t->context.gpath = path;
                                add_task(t);
                            }
                        }
                    }
                } 
            }
        }
        path.pop_back();
        visited.erase(s);
    }

    // coarsely Join left and right table
    void TableJoinCoarse(int lhp, int thlength, vector<int> &visited, HCEQuery &q) 
    {
        bucket & lbucket = q.ltable.get_bucket(lhp);
        bucket_map & lkvmap = lbucket.get_map();
        bucket & rbucket = q.rtable.get_bucket(lhp);
        bucket_map & rkvmap = rbucket.get_map();
        if(rkvmap.find(lhp) != rkvmap.end()) {
            for(auto& lp: lkvmap[lhp]) {
                // set buffer
                for(auto& lv: lp) {
                    cvisited[lv] = true;
                }
                for(auto& rp: rkvmap[lhp]) {
                    bool valid = true;
                    for(auto& rv: rp) {
                        if(cvisited[rv]) { 
                            valid = false;
                            break;
                        }
                    }
                    // satisfy length requirement
                    if(valid && lp.size()+rp.size() <= thlength) {

#ifdef RECORD_RESULT
                        fprintf(gfpout, "* "); // * means 2 Join
                        for(int i=0; i<lp.size(); i++)
                            fprintf(gfpout, "%d ", lp[i]);
                        fprintf(gfpout, "%d ", lhp);
                        for(int i=rp.size()-1; i>=0; i--) 
                            fprintf(gfpout, "%d ", rp[i]);
                        fprintf(gfpout, "\n");
#endif

                        counter++;
                    }
                }
                // reset buffer
                for(auto& lv: lp) {
                    cvisited[lv] = false;
                }
            }
        }
    }

    // fine grained Join left and right table
    void TableJoinFineGrained(int thlength, JoinWrapper &jw, HCEQuery &q)
    {
        for(int cnt=0; cnt<jw.range_map.size(); cnt++) {
            int hp = jw.range_sd[cnt].first;
            auto& range = jw.range_map[cnt];

            bucket & lbucket = q.ltable.get_bucket(hp);
            bucket_map & lkvmap = lbucket.get_map();
            bucket & rbucket = q.rtable.get_bucket(hp);
            bucket_map & rkvmap = rbucket.get_map();

            // sort the vector with less elements, use deep copy
            if(range.lend-range.lstart>range.rend-range.rstart) {
                PathVec rpathvec;
                rpathvec.assign(rkvmap[hp].begin()+range.rstart, rkvmap[hp].begin()+range.rend);
                sort(rpathvec.begin(), rpathvec.end(), comp_plen);
                for(int j=range.lstart; j<range.lend; j++) {
                    // set buffer
                    Path &lp = lkvmap[hp][j];

                    for(auto& lv: lp) {
                        cvisited[lv] = true;
                    }
                    for(int k=0; k<rpathvec.size(); k++) { // sorted
                        Path &rp = rpathvec[k];
                        if(lp.size()+rp.size() > thlength) break;

                        bool valid = true;
                        for(auto& rv: rp) {
                            if(cvisited[rv]) {
                                valid = false;
                                break;
                            }
                        }
                        // satisfy length requirement
                        if(valid) {

#ifdef RECORD_RESULT
                            fprintf(gfpout, "* "); // * means 2 Join
                            for(int i=0; i<lp.size(); i++)
                                fprintf(gfpout, "%d ", lp[i]);
                            fprintf(gfpout, "%d ", hp);
                            for(int i=rp.size()-1; i>=0; i--) 
                                fprintf(gfpout, "%d ", rp[i]);
                            fprintf(gfpout, "\n");
#endif

                            counter++;
                        }
                    }
                    // reset buffer
                    for(auto& lv: lp) {
                        cvisited[lv] = false;
                    }
                }
            } else {
                PathVec lpathvec;
                lpathvec.assign(lkvmap[hp].begin()+range.lstart, lkvmap[hp].begin()+range.lend);
                sort(lpathvec.begin(), lpathvec.end(), comp_plen);

                for(int j=range.rstart; j<range.rend; j++) {
                    // set buffer
                    Path &rp = rkvmap[hp][j];
                    for(auto& rv: rp) {
                        cvisited[rv] = true;
                    }
                    for(int k=0; k<lpathvec.size(); k++) { // sorted
                        Path &lp = lpathvec[k];

                        if(lp.size()+rp.size() > thlength) break;

                        bool valid = true;
                        for(auto& lv: lp) {
                            if(cvisited[lv]) {
                                valid = false;
                                break;
                            }
                        }
                        // satisfy length requirement
                        if(valid) {

#ifdef RECORD_RESULT
                            fprintf(gfpout, "* "); // * means 2 Join
                            for(int i=0; i<lp.size(); i++)
                                fprintf(gfpout, "%d ", lp[i]);
                            fprintf(gfpout, "%d ", hp);
                            for(int i=rp.size()-1; i>=0; i--) 
                                fprintf(gfpout, "%d ", rp[i]);
                            fprintf(gfpout, "\n");
#endif

                            counter++;
                        }
                    }
                    // reset buffer
                    for(auto& rv: rp) {
                        cvisited[rv] = false;
                    }
                }
            }
        }
    }

    void TableJoinFineGrained2(JoinWrapper &jw, int thlength, FILE *gfpout, HCEQuery &q)
    {   
        for(int cnt=0; cnt<jw.range_map.size(); cnt++) {
            int s = jw.range_sd[cnt].first;
            int cur = jw.range_sd[cnt].second;
            auto& range = jw.range_map[cnt];
            auto& path = jw.mid_paths[cnt];

            int total_weight = jw.total_weights[cnt];

            HPPath hpp(path, total_weight);
            bool valid = check_hppath_validity(hpp, cvisited);
            if(valid) {
                bucket & lbucket = q.ltable.get_bucket(s);
                bucket_map & lkvmap = lbucket.get_map();
                bucket & rbucket = q.rtable.get_bucket(cur);
                bucket_map & rkvmap = rbucket.get_map();

                if(range.lend-range.lstart>range.rend-range.rstart) 
                {
                    PathVec rpathvec;
                    rpathvec.assign(rkvmap[cur].begin()+range.rstart, rkvmap[cur].begin()+range.rend);
                    sort(rpathvec.begin(), rpathvec.end(), comp_plen);

                    for(int it=range.lstart; it<range.lend; it++) {

                        // cout << "it: " << it << ", "<< range.lstart << ", " << range.lend << ", lkvmap[s].size() = " << lkvmap[s].size() << endl;
                    
                        Path &lp = lkvmap[s][it];
                        
                        if(total_weight+lp.size()>=thlength) continue;

                        bool valid2 = true;
                        for(auto& lv: lp) {
                            if(cvisited[lv]) {
                                valid2 = false;
                                break;
                            }
                            cvisited2[lv] = true;
                        }
                        if(valid2) {

                            for(int jt=0; jt<rpathvec.size(); jt++) {
                           
                                Path &rp = rpathvec[jt];

                                if(total_weight+lp.size()+rp.size()>thlength) break;

                                valid2 = true;
                                for(auto& rv: rp) {
                                    if(cvisited[rv] || cvisited2[rv]) {
                                        valid2 = false;
                                        break;
                                    }
                                }
                                if(valid2) {
                                    
#ifdef RECORD_RESULT
                                    fprintf(gfpout, "@ "); // @ means 3 Join

                                    for(int i=0; i<lp.size(); i++) {
                                        fprintf(gfpout, "%d ", lp[i]);    
                                    }
                                    
                                    for(int i=0; i<hpp.hppath.size(); i++) {
                                        for(int j=0; j<hpp.hppath[i].path.size(); j++) {
                                            if(i!=0 && j==0) continue;
                                            fprintf(gfpout, "%d ", hpp.hppath[i].path[j]);
                                        }
                                    }  

                                    for(int i=rp.size()-1; i>=0; i--) 
                                        fprintf(gfpout, "%d ", rp[i]);
                                    fprintf(gfpout, "\n");
#endif
                                    counter++;
                                } 
                            }
                        }
                        // reset buffer 
                        for(auto& lv: lp) {
                            cvisited2[lv] = false;
                        }
                    }
                }   
                else {

                    PathVec lpathvec;
                    lpathvec.assign(lkvmap[s].begin()+range.lstart, lkvmap[s].begin()+range.lend);
                    sort(lpathvec.begin(), lpathvec.end(), comp_plen);

                    for(int it=range.rstart; it<range.rend; it++) {
                    
                        Path &rp = rkvmap[cur][it];

                        if(total_weight+rp.size()>=thlength) continue;

                        bool valid2 = true;
                        for(auto& rv: rp) {
                            if(cvisited[rv]) {
                                valid2 = false;
                                break;
                            }
                            cvisited2[rv] = true;
                        }
                        if(valid2) {

                            for(int jt=0; jt<lpathvec.size(); jt++) { // sorted
                                Path &lp = lpathvec[jt];

                                if(total_weight+lp.size()+rp.size()>thlength) break;

                                valid2 = true;
                                for(auto& lv: lp) {
                                    if(cvisited[lv] || cvisited2[lv]) {
                                        valid2 = false;
                                        break;
                                    }
                                }
                                if(valid2) {

#ifdef RECORD_RESULT                    
                                    fprintf(gfpout, "@ "); // @ means 3 Join
                                    for(int i=0; i<lp.size(); i++)
                                        fprintf(gfpout, "%d ", lp[i]);

                                    for(int i=0; i<hpp.hppath.size(); i++) {
                                        for(int j=0; j<hpp.hppath[i].path.size(); j++) {
                                            if(i!=0 && j==0) continue;
                                            fprintf(gfpout, "%d ", hpp.hppath[i].path[j]);
                                        }
                                    }
                                    for(int i=rp.size()-1; i>=0; i--) 
                                        fprintf(gfpout, "%d ", rp[i]);
                                    fprintf(gfpout, "\n");
#endif

                                    counter++;
                                }
                            }
                        }
                        // reset buffer 
                        for(auto& rv: rp) {
                            cvisited2[rv] = false;
                        }
                    }
                }
            }
            for(auto& hpedge: hpp.hppath) {
                for(auto& ve: hpedge.path) {
                    cvisited[ve] = false;
                }
            }
        }
    }

    void TableJoinFineGrained3(int direction, JoinWrapper &jw, int thlength, FILE *gfpout, HCEQuery &q) 
    {   
        for(int cnt=0; cnt<jw.range_map.size(); cnt++) {
            int s = jw.range_sd[cnt].first;
            auto& range = jw.range_map[cnt];
            auto& path = jw.mid_paths[cnt];
            int total_weight = jw.total_weights[cnt];

            // case 2 & 3
            HPPath hpp(path, total_weight);
            if(direction == IN) {
                reverse(hpp.hppath.begin(), hpp.hppath.end());
            }
            bool valid = check_hppath_validity(hpp, cvisited);
            if(valid) {
                if(direction == OUT) {
                    bucket & bucket = q.rtable.get_bucket(s);
                    bucket_map & kvmap = bucket.get_map();

                    // cout<< range.rstart << " " << range.rend <<" "<< kvmap[s].size() << endl;
                    for(int it=range.rstart; it<range.rend; it++) {
                        Path &p = kvmap[s][it];

                        bool valid2 = true;
                        for(auto& ve: p) {
                            if(cvisited[ve])  {
                                valid2 = false;
                                break;
                            }
                        }
                        if(valid2) {
                            // check if satisfy length requirement
                            if(total_weight+p.size()<=thlength) {

#ifdef RECORD_RESULT
                                fprintf(gfpout, "@ ");
                                
                                for(int i=0; i<hpp.hppath.size(); i++) {
                                    for(int j=0; j<hpp.hppath[i].path.size(); j++) {
                                        if(i!=0 && j==0) continue; 
                                        fprintf(gfpout, "%d ", hpp.hppath[i].path[j]);
                                    }
                                }
                                for(int i=p.size()-1; i>=0; i--) 
                                    fprintf(gfpout, "%d ", p[i]);
                                fprintf(gfpout, "\n");
#endif

                                counter++;
                                
                            }
                        }
                    }
                } 
                else 
                {
                    bucket & bucket = q.ltable.get_bucket(s);
                    bucket_map & kvmap = bucket.get_map();
                    for(int it=range.lstart; it<range.lend; it++) {
                        Path &p = kvmap[s][it];

                        bool valid2 = true;
                        for(auto& ve: p) {
                            if(cvisited[ve])  {
                                valid2 = false;
                                break;
                            }
                        }
                        if(valid2) {
                            // check if satisfy length requirement
                            if(total_weight+p.size()<=thlength) {

#ifdef RECORD_RESULT                              
                                fprintf(gfpout, "@ ");

                                for(int i=0; i<p.size(); i++) 
                                    fprintf(gfpout, "%d ", p[i]);

                                for(int i=0; i<hpp.hppath.size(); i++) {
                                    for(int j=0; j<hpp.hppath[i].path.size(); j++) {
                                        if(i!=0 && j==0) continue; 
                                        fprintf(gfpout, "%d ", hpp.hppath[i].path[j]);
                                    }
                                } 
                                fprintf(gfpout, "\n");
#endif

                                counter++;
                            }
                        }
                    }
                }
            }
            for(auto& hpedge: hpp.hppath) {
                for(auto& ve: hpedge.path) {
                    cvisited[ve] = false;
                }
            }
        }
    }

    virtual bool toQuery(string& line, HCEQuery& q)
	{   
        ftime(&q.start_t);
        istringstream istrStream(line);
        istrStream>>q.start>>q.end>>q.length;
        if(!g.mphp[q.start] && !g.mphp[q.end]) {
            q.category = CASE_1;
        } else if(!g.mphp[q.start] && g.mphp[q.end]) {
            q.category = CASE_2;
        } else if(g.mphp[q.start] && !g.mphp[q.end]) {
            q.category = CASE_3;
        } else {
            q.category = CASE_4;
        }
        return true;
	}

    virtual bool task_spawn(HCEQuery &q) // specify which stage of task we're spawning
	{   
        HCETask *task = new HCETask();
        if(q.category == CASE_1) 
        {
            if(q.cur_stage_num == STAGE_1) {
                task->context.start = q.start;
                task->context.end = q.end;
                task->context.thlength = q.length; 
                add_task(task);

            } else if(q.cur_stage_num == STAGE_2) {
                task->context.start = q.end;
                task->context.end = q.start;
                task->context.thlength = q.length;
                add_task(task);

            } else if(q.cur_stage_num == STAGE_3) {
                task->context.thlength = q.length;
                add_task(task);

            } else {
                // q.cur_stage_num == STAGE_4
                int tmp = task->context.thlength = q.length+1;
                task->context.adlength = tmp-q.left_min-q.right_min;
                add_task(task); 
            }
        }
        else if(q.category == CASE_2)
        {
            if(q.cur_stage_num == STAGE_1) {
                task->context.start = q.start;
                task->context.end = q.end;
                task->context.thlength = q.length;
                add_task(task);
            }
            else if(q.cur_stage_num == STAGE_2) {
                add_task(task);
            } else if(q.cur_stage_num == STAGE_3) {
                task->context.end = q.end;
                task->context.thlength = q.length;
                add_task(task);
            } else if(q.cur_stage_num == STAGE_4){
                task->context.start = q.end;
                task->context.cur = q.end;
                task->context.thlength = q.length+1;
                task->context.adlength = q.length-q.left_min+1;
                task->context.total_weight = 0;
                task->context.direction = IN;
                task->context.type = TYPE_1;
                add_task(task);
            }
        }
        else if(q.category == CASE_3)  // CASE_3
        {   
            if(q.cur_stage_num == STAGE_1) {
                add_task(task);
            }
            else if(q.cur_stage_num == STAGE_2) {
                task->context.start = q.end;
                task->context.end = q.start;
                task->context.thlength = q.length;
                add_task(task);

            } else if(q.cur_stage_num == STAGE_3) {
                task->context.start = q.start;
                task->context.thlength = q.length;
                add_task(task);

            } else if(q.cur_stage_num == STAGE_4) {
                task->context.start = q.start;
                task->context.cur = q.start;
                task->context.thlength = q.length+1;
                task->context.adlength = q.length-q.right_min+1;
                task->context.total_weight = 0;
                task->context.direction = OUT;
                task->context.type = TYPE_1;
                add_task(task);
            }
        } 
        else {
            // q.category == CASE_4
            if(q.cur_stage_num == STAGE_1) {    
                add_task(task);
            }
            else if(q.cur_stage_num == STAGE_2) {
                add_task(task);
                
            } else if(q.cur_stage_num == STAGE_3) {
                add_task(task);

            } else if(q.cur_stage_num == STAGE_4) {
                task->context.cur = q.start;
                task->context.end = q.end;
                task->context.thlength = q.length+1;
                task->context.adlength = q.length+1;
                task->context.total_weight = 0;
                
                add_task(task); 
            }
        }
        return true;
	}

    virtual void compute(ContextT &context, HCEQuery &q)
    {   
        if(q.category == CASE_1) 
        {
            // Case 1
            if(q.cur_stage_num == STAGE_1) 
            {   
                counter = q.counters[thread_id];
                ftime(&g.gtime_start[thread_id]);
                dfs_with_hp(context.start, context.end, context.thlength, context.visited, context.gpath, gfpout, OUT, true, q);
                q.counters[thread_id] = counter;

            } else if(q.cur_stage_num == STAGE_2)
            {
                ftime(&g.gtime_start[thread_id]);
                dfs_with_hp(context.start, context.end, context.thlength, context.visited, context.gpath, gfpout, IN, false, q); 

            } else if(q.cur_stage_num == STAGE_3)
            { 
                
                if(context.case3_spawn_flag) {
                    //TableJoinCoarse(context.start, context.thlength, context.visited2, q);
                    counter = q.counters[thread_id];
                    TableJoinFineGrained(context.thlength, context.joinwrapper, q);
                    q.counters[thread_id] = counter;
                }
                else {
                    // std::cout<<"EXECUTE ONLY ONCE"<<std::endl;
                    /*
                    // =========== Coarse tasks ===========
                    for(int i=0; i<CONMAP_BUCKET_NUM; i++) {
                        bucket & lbucket = q.ltable.get_bucket(i);
                        bucket_map & lkvmap = lbucket.get_map();
                        for(auto& lpair: lkvmap) {
                            auto lhp = lpair.first;
                            std::cout<<"lhp: "<<lhp<<" has "<<lkvmap[lhp].size()<<" paths"<<std::endl;
                            HCETask *task = new HCETask();
                            task->context.case3_spawn_flag = true;
                            task->context.start = lhp;
                            task->context.thlength = context.thlength;
                            task->context.visited2.resize(g.mnum_of_vertices, 0);
                            add_task(task);
                        }
                    }
                    // ===================================
                    //*/

                    //*
                    // ========= Fine Grained tasks =============
                    //  1. start from min(ltable, rtable)
                    //  2. append and check >= NPATH_THRESH
                    int lsize, rsize, batch_sz, nbatch, nlpath, nrpath, lastend, cap, leftover;
                    JoinWrapper jw;
                    Range rg;

                    for(int i=0; i<CONMAP_BUCKET_NUM; i++) {
                        bucket & lbucket = q.ltable.get_bucket(i);
                        bucket_map & lkvmap = lbucket.get_map();
                        for(auto& lpair: lkvmap) {
                            auto lhp = lpair.first;
                            lsize = lpair.second.size();
                            bucket & rbucket = q.rtable.get_bucket(lhp);
                            bucket_map & rkvmap = rbucket.get_map();
                            if(rkvmap.find(lhp) != rkvmap.end()) {
                                rsize = rkvmap[lhp].size();
                                // std::cout<<"hp: "<<lhp<<" has "<<lsize<<" left paths and "<<rsize<<" right paths."<<std::endl;

                                if(lsize < rsize) {
                                    // Step 1: merge leftover
                                    cap = jw.capacity(); 
                                    leftover = NPATH_THRESH-cap;
                                    nrpath = leftover%lsize==0? leftover/lsize: leftover/lsize+1; 
                                    nrpath = nrpath>rsize?rsize:nrpath;
                                    rg.lstart = 0;
                                    rg.lend = lsize;
                                    rg.rstart = 0;
                                    rg.rend = nrpath;
                                    if(nrpath != 0) {
                                        jw.range_sd.push_back({lhp, lhp});
                                        jw.range_map.push_back(rg);
                                    }
                                    if(rg.rend == rsize) {
                                        if(jw.capacity() >= NPATH_THRESH) {
                                            HCETask *task = new HCETask();
                                            task->context.case3_spawn_flag = true;
                                            task->context.thlength = context.thlength;
                                            task->context.joinwrapper = jw;
                                            add_task(task);
                                            jw.clear();
                                        }
                                    } else {

                                        HCETask *task = new HCETask();
                                        task->context.case3_spawn_flag = true;
                                        task->context.thlength = context.thlength;
                                        task->context.joinwrapper = jw;
                                        add_task(task);
                                        jw.clear();

                                        // Step 2: handle paths of current HP
                                        batch_sz = NPATH_THRESH%lsize==0? NPATH_THRESH/lsize: NPATH_THRESH/lsize+1;
                                        nbatch = (rsize-nrpath)/batch_sz;
                                        lastend = nrpath;
                                        for(int j=0; j<nbatch; j++) {
                                            // wrap into task, NO LEFTOVER
                                            rg.lstart = 0;
                                            rg.lend = lsize;
                                            rg.rstart = j*batch_sz+nrpath;
                                            lastend = rg.rend = (j+1)*batch_sz+nrpath;
                                            jw.range_sd.push_back({lhp, lhp});
                                            jw.range_map.push_back(rg);

                                            HCETask *task = new HCETask();
                                            task->context.case3_spawn_flag = true;
                                            task->context.thlength = context.thlength;
                                            task->context.joinwrapper = jw;
                                            add_task(task);
                                            jw.clear();
                                        }
                                        // Step 3: handle leftover for current HP
                                        rg.lstart = 0;
                                        rg.lend = lsize;
                                        rg.rstart = lastend; 
                                        rg.rend = rsize;
                                        if(lastend != rsize) {
                                            jw.range_sd.push_back({lhp,lhp});
                                            jw.range_map.push_back(rg);
                                        }
                                    }
                                } else {
                                    // Step 1: merge leftover
                                    cap = jw.capacity();
                                    leftover = NPATH_THRESH-cap;
                                    nlpath = leftover%rsize==0? leftover/rsize: leftover/rsize+1;
                                    nlpath = nlpath>lsize?lsize:nlpath;
                                    rg.lstart = 0;
                                    rg.lend = nlpath;
                                    rg.rstart = 0;
                                    rg.rend = rsize;
                                    if(nlpath != 0) {
                                        jw.range_sd.push_back({lhp, lhp});
                                        jw.range_map.push_back(rg);
                                    }
                                    if(rg.lend == lsize) {
                                        if(jw.capacity() >= NPATH_THRESH) {
                                            HCETask *task = new HCETask();
                                            task->context.case3_spawn_flag = true;
                                            task->context.thlength = context.thlength;
                                            task->context.joinwrapper = jw;
                                            add_task(task);
                                            jw.clear();
                                        }
                                    }
                                    else {
                                        HCETask *task = new HCETask();
                                        task->context.case3_spawn_flag = true;
                                        task->context.thlength = context.thlength;
                                        task->context.joinwrapper = jw;
                                        add_task(task);
                                        jw.clear();

                                        // Step 2: handle paths of current HP
                                        batch_sz = NPATH_THRESH%rsize==0? NPATH_THRESH/rsize: NPATH_THRESH/rsize+1;
                                        nbatch = (lsize-nlpath)/batch_sz;
                                        lastend = nlpath;
                                        for(int j=0; j<nbatch; j++) {
                                            rg.lstart = j*batch_sz+nlpath;
                                            lastend = rg.lend = (j+1)*batch_sz+nlpath;
                                            rg.rstart = 0;
                                            rg.rend = rsize;
                                            jw.range_sd.push_back({lhp,lhp});
                                            jw.range_map.push_back(rg);

                                            HCETask *task = new HCETask();
                                            task->context.case3_spawn_flag = true;
                                            task->context.thlength = context.thlength;
                                            task->context.joinwrapper = jw;
                                            add_task(task);
                                            jw.clear();
                                        }
                                        // Step 3: handle leftover for current HP
                                        rg.lstart = lastend;
                                        rg.lend = lsize;
                                        rg.rstart = 0;
                                        rg.rend = rsize;
                                        if(lastend != lsize) {
                                            jw.range_sd.push_back({lhp,lhp});
                                            jw.range_map.push_back(rg);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // add remaining task
                    if(!jw.empty()) {
                        HCETask *task = new HCETask();
                        task->context.case3_spawn_flag = true;
                        task->context.thlength = context.thlength;
                        task->context.joinwrapper = jw;
                        add_task(task);
                    }

                    // ================================

                }  
            } else // q.cur_stage_num = 4
            {
                if(context.case4_spawn_flag) {
                    if(context.type == TYPE_1) {
                        ftime(&g.gtime_start[thread_id]);
                        dfs_gidx_inside_1(context.start, context.cur, context.visited, context.adlength, context.thlength, context.total_weight, context.path, context.direction, context.joinwrapper, q);
                        if(!context.joinwrapper.empty()) {
                            HCETask *task = new HCETask();
                            task->context.thlength = context.thlength;
                            task->context.joinwrapper = context.joinwrapper;
                            task->context.case4_spawn_flag = true;
                            task->context.type = TYPE_2;
                            add_task(task);
                        }
                    } else {
                        // context.type == TYPE_2
                        counter = q.counters[thread_id];
                        TableJoinFineGrained2(context.joinwrapper, context.thlength, gfpout, q);
                        q.counters[thread_id] = counter;
                    }
                } else {
                    for(int i=0; i<CONMAP_BUCKET_NUM; i++) {
                        bucket & lbucket = q.ltable.get_bucket(i);
                        bucket_map & lkvmap = lbucket.get_map();
                        for(auto& lpair: lkvmap) {
                            auto lhp = lpair.first; 
                            HCETask *task = new HCETask();
                            task->context.case4_spawn_flag = true;
                            task->context.start = lhp;
                            task->context.cur = lhp;
                            task->context.thlength = context.thlength;
                            task->context.adlength = context.adlength;
                            task->context.total_weight = 0;
                            task->context.direction = OUT;
                            task->context.type = TYPE_1;
                            add_task(task);
                        }
                    }
                }
            }
        }
        else if(q.category == CASE_2) {
            // Case 2, no STAGE_2
            if(q.cur_stage_num == STAGE_1) {
                counter = q.counters[thread_id];
                ftime(&g.gtime_start[thread_id]);
                dfs_with_hp(context.start, context.end, context.thlength, context.visited, context.gpath, gfpout, OUT, true, q);
                q.counters[thread_id] = counter;

            } else if(q.cur_stage_num == STAGE_3){
                // ===== let 1 thread do the following =====
                counter = q.counters[thread_id];
                bucket & bucket = q.ltable.get_bucket(context.end);
                bucket_map &kvmap = bucket.get_map();
                if(kvmap.find(context.end) != kvmap.end()) {
                    
                    for(auto& lp: kvmap[context.end]) {
                        if(lp.size() <= context.thlength) {

#ifdef RECORD_RESULT
                            fprintf(gfpout, "* "); // * means 2 Join
                            for(int i=0; i<lp.size(); i++)
                                fprintf(gfpout, "%d ", lp[i]);
                            fprintf(gfpout, "%d ", context.end);
                            fprintf(gfpout, "\n");
#endif
                            counter++;
                        }
                    }
                }
                q.counters[thread_id] = counter;
                //============================================
            } else if(q.cur_stage_num == STAGE_4) {
                if(context.type == TYPE_1) {
                    ftime(&g.gtime_start[thread_id]);
                    dfs_gidx_inside_2(context.start, context.cur, context.visited, context.adlength, context.thlength, context.total_weight, context.path, context.direction, gfpout, context.joinwrapper, q);
                    if(!context.joinwrapper.empty()) {
                        HCETask *task = new HCETask();
                        task->context.thlength = context.thlength;
                        task->context.joinwrapper = context.joinwrapper;
                        task->context.type = TYPE_2;
                        task->context.direction = context.direction;
                        add_task(task);
                    }
                } else {
                    counter = q.counters[thread_id];
                    TableJoinFineGrained3(context.direction, context.joinwrapper, context.thlength, gfpout, q);
                    q.counters[thread_id] = counter;
                }
            }
        }
        else if(q.category == CASE_3) {
            // Case 3, no STAGE_1
            if(q.cur_stage_num == STAGE_2) // flag for stage 1
            {
                counter = q.counters[thread_id];
                ftime(&g.gtime_start[thread_id]);
                dfs_with_hp(context.start, context.end, context.thlength, context.visited, context.gpath, gfpout, IN, true, q); 
                q.counters[thread_id] = counter;

            } else if(q.cur_stage_num == STAGE_3) {
                // ===== let 1 thread do the following =====
                counter = q.counters[thread_id];
                bucket & bucket = q.rtable.get_bucket(context.start);
                bucket_map & kvmap = bucket.get_map(); // don't have to lock since we just read.
                if(kvmap.find(context.start) != kvmap.end()) {
                    for(auto& rp: kvmap[context.start]) {
                        // satisfy length requirement
                        if(rp.size() <= context.thlength) {

#ifdef RECORD_RESULT
                            fprintf(gfpout, "* "); // * means 2 Join
                            fprintf(gfpout, "%d ", context.start);
                            for(int i=rp.size()-1; i>=0; i--)
                                fprintf(gfpout, "%d ", rp[i]);
                            fprintf(gfpout, "\n");
#endif

                            counter++;
                        }
                    }
                }
                q.counters[thread_id] = counter;
                // ==========================================
            } else if(q.cur_stage_num == STAGE_4) {
                if(context.type == TYPE_1) {
                    ftime(&g.gtime_start[thread_id]);
                    dfs_gidx_inside_3(context.start, context.cur, context.visited, context.adlength, context.thlength, context.total_weight, context.path, context.direction, gfpout, context.joinwrapper, q);
                    if(!context.joinwrapper.empty()) {
                        HCETask *task = new HCETask();
                        task->context.thlength = context.thlength;
                        task->context.joinwrapper = context.joinwrapper;
                        task->context.type = TYPE_2;
                        task->context.direction = context.direction;
                        add_task(task);
                    }
                } else {
                    counter = q.counters[thread_id];
                    TableJoinFineGrained3(context.direction, context.joinwrapper, context.thlength, gfpout, q);
                    q.counters[thread_id] = counter;
                }
            }
        } else {
            // q.category == CASE_4
            if(q.cur_stage_num == STAGE_4) {
                counter = q.counters[thread_id];
                ftime(&g.gtime_start[thread_id]);
                dfs_gidx_inside_4(context.cur, context.end, context.visited, context.adlength, context.thlength, context.total_weight, context.path, gfpout, q);
                q.counters[thread_id] = counter;
            }
        }
    }

    virtual bool postprocess(HCEQuery &q) 
    {
		cout<<"[INFO] Query "<<get_queryID()<<", Stage "<<q.cur_stage_num<<"/"<<MAX_STAGE_NUM<<" is done."<<endl;
        
        if(q.cur_stage_num < MAX_STAGE_NUM) {
            q.cur_stage_num++;
            task_spawn(q);
            return true;

        } else {
            ftime(&q.end_t);
            double totaltime = q.end_t.time-q.start_t.time+(double)(q.end_t.millitm-q.start_t.millitm)/1000;
            cout<<"Query "<<get_queryID()<<" total time: "<<totaltime<<endl;

            ULL total_results = 0;
            for(int i=0; i<32; ++i)
            {
                total_results += q.counters[i];
            }
            cout<<"Query "<<get_queryID()<<" total count: "<< total_results <<endl;
            return false;
        }
	}

    virtual bool is_bigTask(ContextValue &context)
	{
        return false;
	}
};

class HCEWorker : public Worker<HCEComper>
{
public:
    HCEWorker(int num_compers) : Worker(num_compers)
    {
    }

    ~HCEWorker()
    {
    }

    void load_data(const char* file_path, int num_compers)
    {
        g.loadGraphFromFile(file_path);
        g.build_graph_idx(4, 60, 60, num_compers); // User defined parameters
    }
};