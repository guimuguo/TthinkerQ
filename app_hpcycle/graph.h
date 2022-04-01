#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <assert.h>
#include <sys/timeb.h>
#include <omp.h>

#include "data.h"

#define OUT 0
#define IN 1

using namespace std;

typedef vector<int> Path;
typedef vector<Path> PathVec;

int comp_int(const void *e1, const void *e2);

struct HPEdge { // there could be multiple edges between u and v
    int u, v; // u and v are both hot points
    Path path; // represent path in G

    HPEdge() {}

    HPEdge(int u_, int v_, Path path_) {
        u = u_;
        v = v_;
        path = path_;
    }
    // # of vertices in path
    inline int weight() {
        return path.size();
    } 
};

bool vec_sort(HPEdge i, HPEdge j);
mutex HPlock_i;
mutex HPlock_o;

typedef vector<HPEdge> HPEdgeSet;

struct HPPath {
    HPEdgeSet hppath;
    int weight_sum;

    HPPath(HPEdgeSet hppath_, int weight_sum_) {
        hppath = hppath_;
        weight_sum = weight_sum_;
    }
};

/**
 |-- End --|--------------- Paths ---------------------|
 |   HP1   |   path 1, path 2, path 3, ...             |
 -------------------------------------------------------
 |   HP2   |   path 1', path 2', ...                   |
 -------------------------------------------------------
 |   HP3   |   ...                                     |
 -------------------------------------------------------
**/
typedef unordered_map<int, PathVec> PathTable; // Left or Right Table

typedef vector<HPPath> HPPathVec; // Middle table with single end point.
typedef unordered_map<int, HPPathVec> HPPathVecMap; // Middle table with mutiple end points.

typedef unordered_map<int, HPEdgeSet> HpAdj;


ofbinstream & operator>>(ofbinstream & m, HPEdgeSet & c)
{   
    size_t size;
    m >> size;
    c.resize(size);
    for (auto it = c.begin(); it != c.end(); ++it) {
        m >> it->u;
        m >> it->v;
        m >> it->path;
    }
    return m;
}
ifbinstream & operator<<(ifbinstream & m, const HPEdgeSet & c) 
{
    m << c.size();
    for (auto it = c.begin(); it != c.end(); ++it) {
        m << it->u;
        m << it->v;
        m << it->path;
    }
    return m;
}

ofbinstream & operator>>(ofbinstream & m, HpAdj & c) 
{
    size_t size;
    m >> size;
    for (size_t i = 0; i < size; i++) {
        int key;
        m >> key;
        m >> c[key];
    }
    return m;
}

ifbinstream & operator<<(ifbinstream & m, const HpAdj & c) 
{
    m << c.size();
    for (auto it = c.begin(); it != c.end(); ++it) {
        m << it->first;
        m << it->second;
    }
    return m;
}


class Graph {
public:
    int mhp_deg_threshold, mhp_deg_threshold_i, mmaxlen_threshold;
    int mnum_of_vertices, mnum_hp;
    int **mppadj_lists_o, *mpadj_list_buf_o, **mppadj_lists_i, *mpadj_list_buf_i; 
    vector<bool> mphp, visited, visited2, visited3; // mphp: flag array. If mphp[v]=true, then v is a hot point, vice versa.
    int *mpdeg_dist; // outdegree distribution in G
    vector<int> mvec_hp; // contains hot point vertex IDs
    unordered_set<int> mset_hp; // contains hot point vertex IDs
    Path path;

    HpAdj mhp_adjlist, mhp_adjlist_i; // adjacent table of G_idx

    // ======== Add from TthinkerQ =============

    timeb gtime_start[40];
    int num_compers;
    
    // =========================================

    // Graph(char *filepath) 
    // {
    //     loadGraphFromFile(filepath);
    // }
    Graph() 
    {
    }

    ~Graph() 
    {
        delete []mppadj_lists_i;
        delete []mpadj_list_buf_i;
        delete []mppadj_lists_o;
	    delete []mpadj_list_buf_o;
        delete[] mpdeg_dist;
    }

    void generateQueries(int num_of_queries, int len);

    int loadGraphFromFile(const char *file_path);
    void ReverseAdj(int num_of_vertices, int nbuf_size);
    // void findAllSimplePathsBase(int start, int end, int length);

    // void findAllSimplePathsHP(int start, int end, int length);
    // void findAllSimplePathsHPCase1(int start, int end, int length, ofstream &fout);
    // void findAllSimplePathsHPCase2(int start, int end, int length, ofstream &fout);
    // void findAllSimplePathsHPCase3(int start, int end, int length, ofstream &fout);
    // void findAllSimplePathsHPCase4(int start, int end, int length, ofstream &fout);

    void build_graph_idx(int hp_deg_threshold, int maxlen_threshold, int hp_deg_threshold_i, int num_compers);

    void dfs(int s, int d, int length, ofstream &fout, int direction);
    void dfs_with_hp(int s, int d, int length, PathTable &ptable, FILE *gfpout, int &min_len, int direction);
    void dfs_build_gidx(int s, int cur, unordered_set<int> &HP_set, unordered_set<int> &cvisited, Path &cpath);

    // dfs: single source single destination
    // void dfs_gidx_inside(int s, int d, int length, int ttllength, int total_weight, HPEdgeSet &path, ofstream &fout);
    // void dfs_gidx_inside(int s, int d, int npoints, int ttllength, ofstream &fout); // wrapper function
    // // dfs: single source mutiple destinations
    // void dfs_gidx_inside(int s, PathTable &table, int npoints, int ttllength, int direction, ofstream &fout); // wrapper function
    // void dfs_gidx_inside(int s, int cur, PathTable &table, int adlength, int ttllength, int total_weight, HPEdgeSet &path, int direction, ofstream &fout);

    // void dfs_gidx_inside(int s, int cur, PathTable &table, PathVec &pvec, int adlength, int ttllength, int total_weight, HPEdgeSet &path, int direction, ofstream &fout);

    // // utitlity
    // bool check_hppath_validity(HPPath &hp_paths, vector<bool> &visited);
    // void print_left_path(Path &lp, ofstream &fout);
    // void print_right_path(Path &rp, ofstream &fout);
    // void print_mid_path(HPPath &hp_paths, ofstream &fout);
};


void Graph::generateQueries(int num_of_queries, int len) 
{
    // generate queries with deg(s)>0 and deg(t)>0
    ofstream qfout("../client/batchFile_version/batch_in_grqc.txt");

    // do k-step random walk
    for(int i=0; i<num_of_queries; i++) {
        int s = rand() % mnum_of_vertices;
        int k = 0;
        int cur = s;
        bool succ = true;
        while(k++ < len) {
            int deg = mppadj_lists_o[cur][0];
            if(deg == 0) {
                succ = false;
                break;
            }
            cur = mppadj_lists_o[cur][rand()%deg+1];
        }
        if(!succ) {
            i--;
            continue;
        }
        qfout << s << " " << cur << " " << len << '\n';
    }
    qfout << "server_exit" << endl;
    qfout.close();
}


int Graph::loadGraphFromFile(const char *file_path)
{
    Data *pdata;
	Transaction *ptransaction;
	int num_of_vertices, nbuf_size, nvertex_no, nbuf_pos;
	int *padj_lens, nsize, *ptemp, nlist_len, i, nmax_deg;

	pdata = new Data(file_path);

	nsize = 5000; // capacity of padj_lens
	padj_lens = new int[nsize]; // buffer for keeping each transaction for parsing
	memset(padj_lens, 0, sizeof(int)*nsize);

	num_of_vertices = 0;
	nbuf_size = 0; // ---> all transactions are saved in one buffer.... (concat adj-lists)
	nmax_deg = 0; // track longest transaction, just for printing at Line 1855...
	ptransaction = pdata->getNextTransaction();
	while(ptransaction)
	{
		if(nmax_deg<ptransaction->length)
			nmax_deg = ptransaction->length;
		if(num_of_vertices>=nsize)
		{
			ptemp = new int[2*nsize]; // ptemp is used to double the capacity of padj_lens, (vector)
			memcpy(ptemp, padj_lens, sizeof(int)*nsize);
			memset(&ptemp[nsize], 0, sizeof(int)*nsize);
			delete []padj_lens;
			padj_lens = ptemp;
			nsize *= 2;
		}
		padj_lens[num_of_vertices] = ptransaction->length;
		num_of_vertices++;
		nbuf_size += ptransaction->length; // ---> all transactions are saved in one buffer....
		ptransaction = pdata->getNextTransaction(); // line by line transaction reading
	}

	mppadj_lists_o = new int*[num_of_vertices]; // mppadj_lists[] keeps the start position (pointer) of each vertex's adj-list in mpadj_list_buf
	mpadj_list_buf_o = new int[num_of_vertices+nbuf_size];
	nvertex_no = 0;
	nbuf_pos = 0;

    visited.resize(num_of_vertices, false);
    visited2.resize(num_of_vertices, false);
    visited3.resize(num_of_vertices, false);

    mpdeg_dist = new int[num_of_vertices];
    memset(mpdeg_dist, 0, sizeof(int)*num_of_vertices);


	ptransaction = pdata->getNextTransaction(); // note the rewind operation inside getNextTransaction()
	while(ptransaction)
	{
		nlist_len = 1;
		mppadj_lists_o[nvertex_no] = &mpadj_list_buf_o[nbuf_pos];
		for(i=0; i<ptransaction->length; i++)
		{
			if(ptransaction->t[i]!=nvertex_no) // remove self-loop
				mppadj_lists_o[nvertex_no][nlist_len++] = ptransaction->t[i];
		}
		nbuf_pos += nlist_len;
		mppadj_lists_o[nvertex_no][0] = nlist_len-1; // from position 1 onwards, neighbors are kept; position 0 keeps the number of neighbors (adj-list length)
        // temporarially don't need to sort.
		qsort(&mppadj_lists_o[nvertex_no][1], mppadj_lists_o[nvertex_no][0], sizeof(int), comp_int); // adj-lists are sorted !!!

        mpdeg_dist[nlist_len-1]++;
		nvertex_no++;
		ptransaction = pdata->getNextTransaction();
	}
	delete pdata;
	delete []padj_lens;

    printf("maximum vertex outdegree: %d\n", nmax_deg); // max_deg just for printing...

    // opposite direction
    ReverseAdj(num_of_vertices, nbuf_size);

	mnum_of_vertices = num_of_vertices;

    cout<<"# of vertices : "<<mnum_of_vertices<<endl;
    // for(int i=0; i<num_of_vertices; i++) {
    //     if(mpdeg_dist[i]>0) {
    //         cout<<"OutDegree: "<<i<<" has "<<mpdeg_dist[i]<<" vertices."<<endl;
    //     }
    // }
	return num_of_vertices;
}

void Graph::ReverseAdj(int num_of_vertices, int nbuf_size) 
{
	int nlist_len, indeg_max = 0;
	PathVec in_vec(num_of_vertices);
	for(int i=0;i<num_of_vertices;i++)
	{
		nlist_len = 0;
		if(mppadj_lists_o[i]!= NULL) {
			for(int j=1; j<=mppadj_lists_o[i][0];j++)
				in_vec[mppadj_lists_o[i][j]].push_back(i);
        }
	}

	mppadj_lists_i = new int*[num_of_vertices]; // mppadj_lists[] keeps the start position (pointer) of each vertex's adj-list in mpadj_list_buf
	mpadj_list_buf_i = new int[num_of_vertices+nbuf_size];

	int nbuf_pos = 0;
	for(int i=0;i<num_of_vertices;i++)
	{
		nlist_len = 1;
		mppadj_lists_i[i] = &mpadj_list_buf_i[nbuf_pos];
		for(int j=0;j<in_vec[i].size();j++)
		{
			mppadj_lists_i[i][nlist_len++] = in_vec[i][j];
		}
		nbuf_pos += nlist_len;
        // track max indegree
        if(indeg_max < nlist_len-1)
            indeg_max = nlist_len-1;

		mppadj_lists_i[i][0] = in_vec[i].size(); // from position 1 onwards, neighbors are kept; position 0 keeps the number of neighbors (adj-list length)
		qsort(&mppadj_lists_i[i][1], mppadj_lists_i[i][0], sizeof(int), comp_int); // adj-lists are sorted !!!
	}
    printf("maximum vertex indegree: %d\n", indeg_max);
}

/*
void Graph::findAllSimplePathsBase(int start, int end, int length)
{   
    ofstream fout;
    char file[100], no[100];
	strcpy(file, "path");
	sprintf(no, "_%d_%d", start, end);
	strcat(file, no);
	fout.open(file);
    dfs(start, end, length, fout, OUT); 
    fout.close();
}

void Graph::findAllSimplePathsHP(int start, int end, int length) 
{   
    assert(length <= mmaxlen_threshold);
    ofstream fout;
    char file[100], no[100];
	strcpy(file, "path_hp");
	sprintf(no, "_%d_%d", start, end);
	strcat(file, no);
	fout.open(file);

    if(!mphp[start] && !mphp[end])  // case 1: neither s nor d is HP
        findAllSimplePathsHPCase1(start, end, length, fout);
    else if(!mphp[start] && mphp[end]) // case 2: s is not HP, d is HP      
        findAllSimplePathsHPCase2(start, end, length, fout);
    else if(mphp[start] && !mphp[end]) // case 3: s is HP, d is not HP
        findAllSimplePathsHPCase3(start, end, length, fout);
    else  // case 4: s and d are both HP
        findAllSimplePathsHPCase4(start, end, length, fout);
    fout.close();
}


void Graph::findAllSimplePathsHPCase1(int start, int end, int length, ofstream &fout) 
{
    cout << "This is Case 1 ..."<<endl;
    bool valid;
    int left_min=mnum_of_vertices, right_min=mnum_of_vertices; // records the minimum length in left(right) table

    PathTable ltable, rtable;
    dfs_with_hp(start, end, length, ltable, fout, left_min, OUT);
    dfs_with_hp(end, start, length, rtable, fout, right_min, IN);

    // Step 2, left+right table join
    for(auto& lpair: ltable) { 
        auto lhp = lpair.first; 
        auto lpaths = lpair.second;
        if(rtable.find(lhp) != rtable.end()) {
            for(auto& lp: lpaths) {
                // set buffer
                for(auto& lv: lp) {
                    visited2[lv] = true;
                }
                for(auto& rp: rtable[lhp]) {
                    valid = true;
                    for(auto& rv: rp) {
                        if(visited2[rv]) {
                            valid = false;
                            break;
                        }
                    }
                    // satisfy length requirement
                    if(valid && lp.size()+rp.size()+1 <= length+1) {
                        print_left_path(lp, fout);
                        fout<<lhp<<" ";
                        print_right_path(rp, fout);
                        fout<<"\n";
                    }
                }
                // reset buffer
                for(auto& lv: lp) {
                    visited2[lv] = false;
                }
            }
        }
    }

    // left+mid+right table join with one for-loop
    for(auto& lpair: ltable) {
        auto lhp = lpair.first;
        auto lpaths = lpair.second;
        HPEdgeSet path; 
        dfs_gidx_inside(lhp, lhp, rtable, lpaths, length-left_min-right_min+1, length+1, 0, path, OUT, fout);
    }
}

// case 2: start is not HP, end is HP
void Graph::findAllSimplePathsHPCase2(int start, int end, int length, ofstream &fout)  
{   
    cout << "This is Case 2 ..."<<endl;
    int left_min=mnum_of_vertices;
    PathTable ltable;
    dfs_with_hp(start, end, length, ltable, fout, left_min, OUT);

    if(ltable.find(end) != ltable.end()) {
        for(auto& lp: ltable[end]) {
            // satisfy length requirement
            if(lp.size() <= length) { // lp.size()+1 <= length+1
                print_left_path(lp, fout);
                fout<<end<<" ";
                fout<<"\n";
            }
        }
    }
    // [hp2->e][hp1->hp2]
    // [hp1->hp2][hp2->e]
    dfs_gidx_inside(end, ltable, length-left_min+1, length+1, IN, fout);
}

// case 3: start is HP, end is not HP
void Graph::findAllSimplePathsHPCase3(int start, int end, int length, ofstream &fout)  
{
    // ========== serial ==========
    cout << "This is Case 3 ..."<<endl;
    int right_min=mnum_of_vertices;
    PathTable rtable;
    dfs_with_hp(end, start, length, rtable, fout, right_min, IN); 

    if(rtable.find(start) != rtable.end()) {
        for(auto& rp: rtable[start]) {
            // satisfy length requirement
            if(rp.size() <= length) {
                fout<<start<<" ";
                print_right_path(rp, fout);
                fout<<"\n";
            }
        }
    }
    // ======= parallel ============
    dfs_gidx_inside(start, rtable, length-right_min+1, length+1, OUT, fout);
}


void Graph::findAllSimplePathsHPCase4(int start, int end, int length, ofstream &fout)
{
    cout << "This is Case 4 ..."<<endl;
    dfs_gidx_inside(start, end, length+1, length+1, fout);
}


bool Graph::check_hppath_validity(HPPath &hp_paths, vector<bool>& visited)
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

void Graph::print_mid_path(HPPath &hp_paths, ofstream &fout)
{
    for(int i=0; i<hp_paths.hppath.size(); i++) {
        for(int j=0; j<hp_paths.hppath[i].path.size(); j++) {
            if(i!=0 && j==0) continue; 
            fout<<hp_paths.hppath[i].path[j]<<" ";
        }
    }
}

void Graph::print_left_path(Path &lp, ofstream &fout)
{
    for(int i=0; i<lp.size(); i++) fout<<lp[i]<<" ";
}

void Graph::print_right_path(Path &rp, ofstream &fout) 
{
    for(int i=rp.size()-1; i>=0; i--) fout<<rp[i]<<" ";
}

// total_num += # of newly found paths
void Graph::dfs(int s, int d, int length, ofstream &fout, int direction)
{
    visited[s] = true; 
    path.push_back(s);
    if(s == d) {
        // report simple path
        for(int i=0; i<path.size(); i++)
        {   
            fout<<path[i]<<" ";
        }
        fout<<"\n";
    } else {
        if(path.size() <= length) { 
            if(direction == OUT) {
                for(int i=1; i<=mppadj_lists_o[s][0]; i++) {
                    int ne = mppadj_lists_o[s][i];
                    if(!visited[ne]) {
                        dfs(ne, d, length, fout, OUT);
                    }
                }
            } else {
                for(int i=1; i<=mppadj_lists_i[s][0]; i++) {
                    int ne = mppadj_lists_i[s][i];
                    if(!visited[ne]) {
                        dfs(ne, d, length, fout, IN);
                    }
                }
            }
        }
    }
    path.pop_back(); 
    visited[s] = false;
}
*/

// ptable: L/R and p(h)
void Graph::dfs_with_hp(int s, int d, int length, PathTable &ptable, FILE *gfpout, int &min_len, int direction)
{
    visited[s] = true;
    path.push_back(s);
    if(s == d && direction == OUT) { 
        // report simple path, directly path find, w/o HP-index
        for(int i=0; i<path.size(); i++)
        {   
            fprintf(gfpout, "%d ", path[i]);
        }
        fprintf(gfpout, "\n");
        fflush(gfpout);
    } else {
        if(path.size() <= length) {
            if(direction == OUT) {
                for(int i=1; i<=mppadj_lists_o[s][0]; i++) {
                    int ne = mppadj_lists_o[s][i];
                    // check whether ne is HP
                    if(mphp[ne]) { 
                        if(ptable.find(ne) == ptable.end()) {
                            ptable[ne] = PathVec();
                        } 
                        ptable[ne].push_back(path);
                        // tracks minimum length
                        if(path.size()<min_len) min_len = path.size();
                        continue;
                    }
                    if(!visited[ne]) {
                        dfs_with_hp(ne, d, length, ptable, gfpout, min_len, OUT);
                    }
                }
            } else {
                for(int i=1; i<=mppadj_lists_i[s][0]; i++) {
                    int ne = mppadj_lists_i[s][i];
                    // check whether ne is HP
                    if(mphp[ne]) { 
                        if(ptable.find(ne) == ptable.end()) {
                            ptable[ne] = PathVec();
                        } 
                        ptable[ne].push_back(path);
                        // tracks minimum length
                        if(path.size()<min_len) min_len = path.size();
                        continue;
                    }
                    if(!visited[ne]) {
                        dfs_with_hp(ne, d, length, ptable, gfpout, min_len, IN);
                    }
                }
            } 
        }
    }
    path.pop_back(); 
    visited[s] = false;
}

void Graph::dfs_build_gidx(int s, int cur, unordered_set<int> &HP_set, unordered_set<int> &cvisited,
                            Path &cpath)
{
    cvisited.insert(cur);
    cpath.push_back(cur);
    // s, d, middle: 6-2=4=5-1, [s, x1,x2,x3,x4, d] is allowed
    // [s, hp, x1, x2, hp, d] 
    // [ltable=null, s=hp, x1, x2, x3, hp, d]
    // [s, hp, x1, x2, x3, d=hp, rtable=null] 
    // [ltable=null, s=hp, x1, x2, x3, x4, d=hp, rtable=null] 

    // if start and end aren't HP, total_weight <= mlen_threshold-1
    // if either start or end is HP(not both),  total_weight <= mlen_threshold
    // if start and end both are HP, total_weight <= mlen_threshold+1
    // For ease of computation, we can just set total_weight <= mlen_threshold-1

    int threshold = mmaxlen_threshold; // TODO:+1/0/-1
    if(s!=cur && HP_set.find(cur)!=HP_set.end() && cpath.size()<=threshold) { // prevent self-loop
        HPEdge e(s, cur, cpath);
        HPlock_o.lock();
        if(mhp_adjlist.find(s) == mhp_adjlist.end()) {
            mhp_adjlist[s] = HPEdgeSet();
        }
        mhp_adjlist[s].push_back(e);
        HPlock_o.unlock();
        // ===== add reverse HPEdge ======
        HPEdge er(cur, s, cpath);
        HPlock_i.lock();
        if(mhp_adjlist_i.find(cur) == mhp_adjlist_i.end()) {
            mhp_adjlist_i[cur] = HPEdgeSet();
        }
        mhp_adjlist_i[cur].push_back(er);
        HPlock_i.unlock();
    } else {
        if(cpath.size() <= threshold) {
            for(int i=1; i<=mppadj_lists_o[cur][0]; i++) {
                int ne = mppadj_lists_o[cur][i];
                if(cvisited.find(ne) == cvisited.end()) {
                    dfs_build_gidx(s, ne, HP_set, cvisited, cpath);
                }
            }
        }
    }
    cpath.pop_back();
    cvisited.erase(cur);
}

void Graph::build_graph_idx(int maxlen_threshold, int hp_deg_threshold, int hp_deg_threshold_i, int num_compers)
{   
    mmaxlen_threshold = maxlen_threshold;
    mhp_deg_threshold = hp_deg_threshold;
    mhp_deg_threshold_i = hp_deg_threshold_i;

    mphp.resize(mnum_of_vertices, false);
    mnum_hp = 0;

    // mark hot points by their outdegree
    for(int i=0; i<mnum_of_vertices; i++) {
        // consider both in-degree and out-degree
        if(mppadj_lists_o[i][0] >= mhp_deg_threshold && mppadj_lists_i[i][0] >= mhp_deg_threshold_i) {
            mphp[i] = true;
            mvec_hp.push_back(i);
            mset_hp.insert(i);
            mnum_hp++;
        }
    }
    sort(mvec_hp.begin(), mvec_hp.end());

    cout<<"# of hot points: "<<mnum_hp<<endl;
    cout<<"HPs: ";
    for(auto& el: mvec_hp)
        cout<<el<<" ";
    cout<<"\n";
    
    // build G_idx

    unordered_set<int> cvisited[32];
    Path cpath[32];

    // always use 32 threads here
    //**
#pragma omp parallel for schedule(dynamic, 1) num_threads(32)
    for(int i=0; i<mvec_hp.size(); i++) {
        int tid = omp_get_thread_num();
        dfs_build_gidx(mvec_hp[i], mvec_hp[i], mset_hp, cvisited[tid], cpath[tid]); 
    }

    // sort mhp_adjlist_i w.r.t. weight
#pragma omp parallel for schedule(dynamic, 1) num_threads(32)
    for(int i=0; i<mvec_hp.size(); i++) {
        sort(mhp_adjlist[mvec_hp[i]].begin(), mhp_adjlist[mvec_hp[i]].end(), vec_sort);
        sort(mhp_adjlist_i[mvec_hp[i]].begin(), mhp_adjlist_i[mvec_hp[i]].end(), vec_sort);
    }

    // do serialization on disk

    // ifbinstream hpadj("./baidu_hpadj");
    // hpadj << mhp_adjlist;
    // hpadj.close();

    // ifbinstream hpadj_i("./baidu_hpadj_i");
    // hpadj_i << mhp_adjlist_i;
    // hpadj_i.close();
    /**/

    // ==================================
    // directly loading from disk

    /**
    ofbinstream in("./baidu_hpadj");
    in >> mhp_adjlist;
    in.close();

    ofbinstream in_i("./baidu_hpadj_i");
    in_i >> mhp_adjlist_i;
    in_i.close();
    /**/

    // =================================

    // cout<<"*** Info of G_idx ***" << endl;
    int total_p = 0;
    for(auto& pair: mhp_adjlist_i) {
        auto hp = pair.first;
        
        // for(int j=0; j<mhp_adjlist_i[hp].size(); j++) {
        //     int prev = 0;
        //     if(prev > mhp_adjlist_i[hp][j].weight()) cout << "NOT SORT" << endl;
        //     prev = mhp_adjlist_i[hp][j];
        // }
        /* a switch
        for(int j=0; j<mhp_adjlist[hp].size(); j++) {
            cout<<"From "<<mhp_adjlist[hp][j].u<<" to "<<mhp_adjlist[hp][j].v<<", weight: "<<mhp_adjlist[hp][j].weight()<<", path: ";
            total_p++;
            for(int k=0; k<mhp_adjlist[hp][j].path.size(); k++) {
                cout<<mhp_adjlist[hp][j].path[k]<<" ";
            }
            cout<<endl;
        }
        //*/
        //*
        total_p += mhp_adjlist_i[hp].size();
        //*/
    }
    cout<<"# edges in G_idx (in): "<<total_p<<endl;

    total_p = 0;
    for(auto& pair: mhp_adjlist) {
        auto hp = pair.first;
        
        /* a switch
        for(int j=0; j<mhp_adjlist[hp].size(); j++) {
            cout<<"From "<<mhp_adjlist[hp][j].u<<" to "<<mhp_adjlist[hp][j].v<<", weight: "<<mhp_adjlist[hp][j].weight()<<", path: ";
            total_p++;
            for(int k=0; k<mhp_adjlist[hp][j].path.size(); k++) {
                cout<<mhp_adjlist[hp][j].path[k]<<" ";
            }
            cout<<endl;
        }
        //*/
        //*
        total_p += mhp_adjlist[hp].size();
        //*/ 
    }
    cout<<"# edges in G_idx (out): "<<total_p<<endl;

    // generateQueries(1000, 5);
}

/*
// wrapper function
void Graph::dfs_gidx_inside(int s, int d, int npoints, int ttllength, ofstream &fout) // npoints: # of points
{
    HPEdgeSet path;
    dfs_gidx_inside(s, d, npoints, ttllength, 0, path, fout); 
}

// output: plist
// adlength: adjusted length
void Graph::dfs_gidx_inside(int s, int d, int adlength, int ttllength, int total_weight, HPEdgeSet &path, ofstream &fout) 
{   
    visited[s] = true;
    if(s == d) {
        HPPath hpp(path, total_weight);
        bool valid = check_hppath_validity(hpp, visited2);    
        if(valid) {
            // satisfy length requirement
            if(total_weight<=ttllength) {
                print_mid_path(hpp, fout);
                fout<<"\n";
            }
        }
        for(auto& hpedge: hpp.hppath) {
            for(auto& ve: hpedge.path) {
                visited2[ve] = false;
            }
        }
    } else {
        int new_weight;
        for(auto& ne: mhp_adjlist[s]) {
            if(!visited[ne.v]) {
                path.push_back(ne);
                if(total_weight == 0) 
                    new_weight = ne.weight();
                else 
                    new_weight = total_weight+ne.weight()-1;
                if(new_weight <= adlength) { 
                    dfs_gidx_inside(ne.v, d, adlength, ttllength, new_weight, path, fout); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
                    path.pop_back();
                } else {
                    path.pop_back();
                    break; 
                }
            }
        }
        
    }
    visited[s] = false;
}

void Graph::dfs_gidx_inside(int s, PathTable &table, int npoints, int ttllength, int direction, ofstream &fout)
{
    HPEdgeSet path;
    dfs_gidx_inside(s, s, table, npoints, ttllength, 0, path, direction, fout); 
}

void Graph::dfs_gidx_inside(int s, int cur, PathTable &table, int adlength, int ttllength, int total_weight, HPEdgeSet &path, int direction, ofstream &fout)
{
    visited[cur] = true;
    if(s!=cur && table.find(cur)!=table.end()) { // prevent self-loop
        // case 2 & 3
        HPPath hpp(path, total_weight);
        if(direction == IN) {
            reverse(hpp.hppath.begin(), hpp.hppath.end());
        }
        bool valid = check_hppath_validity(hpp, visited2);
        if(valid) {
            for(auto& p: table[cur]) {
                bool valid2 = true;
                for(auto& ve: p) {
                    if(visited2[ve])  {
                        valid2 = false;
                        break;
                    }
                }
                if(valid2) {
                    // check if satisfy length requirement
                    if(total_weight+p.size()<=ttllength) {
                        if(direction == IN) {
                            print_left_path(p, fout);
                            print_mid_path(hpp, fout);
                        } else {
                            print_mid_path(hpp, fout);
                            print_right_path(p, fout);
                        }
                        fout<<"\n";
                    }
                }
            }
        }
        for(auto& hpedge: hpp.hppath) {
            for(auto& ve: hpedge.path) {
                visited2[ve] = false;
            }
        }
    }
    int new_weight;
    if(direction == OUT) {
        for(auto& ne: mhp_adjlist[cur]) {
            if(!visited[ne.v]) {
                path.push_back(ne);
                if(total_weight == 0)
                    new_weight = ne.weight();
                else
                    new_weight = total_weight+ne.weight()-1;
                if(new_weight <= adlength) {
                    dfs_gidx_inside(s, ne.v, table, adlength, ttllength, new_weight, path, OUT, fout); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
                    path.pop_back();
                } else {
                    path.pop_back();
                    break; 
                }
            }
        }
    } else {
        for(auto& ne: mhp_adjlist_i[cur]) {
            if(!visited[ne.v]) {
                path.push_back(ne);
                if(total_weight == 0)
                    new_weight = ne.weight();
                else 
                    new_weight = total_weight+ne.weight()-1;
                if(new_weight <= adlength) {
                    dfs_gidx_inside(s, ne.v, table, adlength, ttllength, new_weight, path, IN, fout); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
                    path.pop_back();
                } else {
                    path.pop_back();
                    break; 
                }
            }
        }
    }
    
    visited[cur] = false;
}

// ctable: counter table, table: rtable
void Graph::dfs_gidx_inside(int s, int cur, PathTable &table, PathVec &pvec, int adlength, int ttllength, int total_weight, HPEdgeSet &path, int direction, ofstream &fout)
{
    visited[cur] = true;
    if(s!=cur && table.find(cur)!=table.end()) { // prevent self-loop
        HPPath hpp(path, total_weight);
        bool valid = check_hppath_validity(hpp, visited2);
        if(valid) {
            for(auto& lp: pvec) {
                bool valid2 = true;
                for(auto& lv: lp) {
                    if(visited2[lv]) {
                        valid2 = false;
                        break;
                    }
                    visited3[lv] = true;
                }
                if(valid2) {
                    for(auto& rp: table[cur]) {
                        // first check lp and rp don't have repetitive vertices
                        valid2 = true;
                        for(auto& rv: rp) {
                            if(visited2[rv] || visited3[rv])  {
                                valid2 = false;
                                break;
                            }
                        }
                        if(valid2) {
                            // check if satisfy length requirement
                            if(total_weight+lp.size()+rp.size()<=ttllength) {
                                print_left_path(lp, fout);
                                print_mid_path(hpp, fout);
                                print_right_path(rp, fout);
                                fout<<"\n";
                            }
                        }
                    }
                }
                // reset buffer 
                for(auto& lv: lp) {
                    visited3[lv] = false;
                }
            }
        }
        for(auto& hpedge: hpp.hppath) {
            for(auto& ve: hpedge.path) {
                visited2[ve] = false;
            }
        }
    }
    int new_weight;
    if(direction == OUT) {
        for(auto& ne: mhp_adjlist[cur]) {
            if(!visited[ne.v]) {
                path.push_back(ne);
                if(total_weight == 0)
                    new_weight = ne.weight();
                else
                    new_weight = total_weight+ne.weight()-1;
                if(new_weight <= adlength) {
                    dfs_gidx_inside(s, ne.v, table, pvec, adlength, ttllength, new_weight, path, OUT, fout); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
                    path.pop_back();
                } else {
                    path.pop_back();
                    break; 
                }
            }
        }
    } else {
        for(auto& ne: mhp_adjlist_i[cur]) {
            if(!visited[ne.v]) {
                path.push_back(ne);
                if(total_weight == 0)
                    new_weight = ne.weight();
                else 
                    new_weight = total_weight+ne.weight()-1;
                if(new_weight <= adlength) {
                    dfs_gidx_inside(s, ne.v, table, pvec, adlength, ttllength, new_weight, path, IN, fout); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
                    path.pop_back();
                } else {
                    path.pop_back();
                    break; 
                }
            }
        }
    }
    visited[cur] = false;
}
*/

int comp_int(const void *e1, const void *e2)
{
    int n1, n2;
    n1 = *(int *) e1;
    n2 = *(int *) e2;
    if (n1>n2)
        return 1;
    else if (n1<n2)
        return -1;
    else
        return 0;
}

bool vec_sort(HPEdge i, HPEdge j)
{
    return i.weight()<j.weight();
}