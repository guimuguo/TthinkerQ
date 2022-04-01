#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <assert.h>
#include <sys/timeb.h>

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

class Graph {
public:
    int mhp_deg_threshold, mmaxlen_threshold;
    int mnum_of_vertices, mnum_hp;
    int **mppadj_lists_o, *mpadj_list_buf_o, **mppadj_lists_i, *mpadj_list_buf_i; 
    vector<bool> mphp, visited, visited2; // mphp: flag array. If mphp[v]=true, then v is a hot point, vice versa.
    int *mpdeg_dist; // outdegree distribution in G
    vector<int> mvec_hp; // contains hot point vertex IDs
    unordered_set<int> mset_hp; // contains hot point vertex IDs
    Path path; 

    HpAdj mhp_adjlist, mhp_adjlist_i; // adjacent table of G_idx

    // ======== Add from TthinkerQ =============

    timeb gtime_start[40];
    

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
        delete[] mpadj_list_buf_o;
        delete[] mpadj_list_buf_i;
        delete[] mppadj_lists_o;
        delete[] mppadj_lists_i;
        delete[] mpdeg_dist;
    }

    int loadGraphFromFile(const char *file_path);
    void ReverseAdj(int num_of_vertices, int nbuf_size);
    void findAllSimplePathsBase(int start, int end, int length);

    void findAllSimplePathsHP(int start, int end, int length);
    void findAllSimplePathsHPCase1(int start, int end, int length, ofstream &fout, int &total_num);
    void findAllSimplePathsHPCase2(int start, int end, int length, ofstream &fout, int &total_num);
    void findAllSimplePathsHPCase3(int start, int end, int length, ofstream &fout, int &total_num);
    void findAllSimplePathsHPCase4(int start, int end, int length, ofstream &fout, int &total_num);

    void build_graph_idx(int hp_deg_threshold, int maxlen_threshold);

    void dfs(int s, int d, int length, ofstream &fout, int &total, int direction);
    void dfs_with_hp(int s, int d, int length, PathTable &ptable, ofstream &fout, int &total, int &min_len, int direction);
    void dfs_build_gidx(int s, int cur, unordered_set<int> &HP_set);

    // dfs: single source single destination
    void dfs_gidx_inside(int s, int d, int length, int total_weight, HPEdgeSet &path, HPPathVec &plist);
    void dfs_gidx_inside(int s, int d, int npoints, HPPathVec &plist); // wrapper function
    // dfs: single source mutiple destinations
    void dfs_gidx_inside(int s, PathTable &table, int npoints, HPPathVecMap &plist, int direction); // wrapper function
    void dfs_gidx_inside(int s, int cur, PathTable &table, int adlength, int total_weight, HPEdgeSet &path, HPPathVecMap &plist, int direction);

    // utitlity
    bool check_hppath_validity(HPPath &hp_paths, vector<bool> &visited);
    void print_left_path(Path &lp, ofstream &fout);
    void print_right_path(Path &rp, ofstream &fout);
    void print_mid_path(HPPath &hp_paths, ofstream &fout);
};


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


void Graph::findAllSimplePathsBase(int start, int end, int length)
{   
    ofstream fout;
    char file[100], no[100];
	strcpy(file, "path");
	sprintf(no, "_%d_%d", start, end);
	strcat(file, no);
	fout.open(file);
    int total_num = 0; // total # results found
    dfs(start, end, length, fout, total_num, OUT);
    cout<<"Total # of results: "<<total_num<<endl; 
    fout.close();
}

void Graph::findAllSimplePathsHP(int start, int end, int length) 
{   
    assert(length <= mmaxlen_threshold);
    int total_num = 0;
    ofstream fout;
    char file[100], no[100];
	strcpy(file, "path_hp");
	sprintf(no, "_%d_%d", start, end);
	strcat(file, no);
	fout.open(file);

    if(!mphp[start] && !mphp[end])  // case 1: neither s nor d is HP
        findAllSimplePathsHPCase1(start, end, length, fout, total_num);
    else if(!mphp[start] && mphp[end]) // case 2: s is not HP, d is HP      
        findAllSimplePathsHPCase2(start, end, length, fout, total_num);
    else if(mphp[start] && !mphp[end]) // case 3: s is HP, d is not HP
        findAllSimplePathsHPCase3(start, end, length, fout, total_num);
    else  // case 4: s and d are both HP
        findAllSimplePathsHPCase4(start, end, length, fout, total_num);
    fout.close();
}

void Graph::findAllSimplePathsHPCase1(int start, int end, int length, ofstream &fout, int &total_num) 
{   
    int tmp, cur_length; 
    int left_min=mnum_of_vertices, right_min=mnum_of_vertices; // records the minimum length in left(right) table
    // valid: flag indicating path within mid table is valid
    // valid2: flag indicating path across left, mid, right is valid.
    bool valid, valid2;

    PathTable ltable, rtable;
    dfs_with_hp(start, end, length, ltable, fout, total_num, left_min, OUT);
    dfs_with_hp(end, start, length, rtable, fout, tmp, right_min, IN);

    cout << "# path with no HP: "<<total_num<<endl;

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
                        total_num += 1;
                    }
                }
                // reset buffer
                for(auto& lv: lp) {
                    visited2[lv] = false;
                }
            }
        }
    }
    cout << "# path with at most 1 HP: "<<total_num<<endl;

    // left+mid+right table join with one for-loop
    for(auto& lpair: ltable) {
        auto lhp = lpair.first; 
        auto lpaths = lpair.second;
        HPPathVecMap ptable_idx;
        dfs_gidx_inside(lhp, rtable, length-left_min-right_min+1, ptable_idx, OUT);
        for(auto& rpair: rtable) {
            auto rhp = rpair.first; 
            auto rpaths = rpair.second;
            for(auto& hp_paths: ptable_idx[rhp]) {
                cur_length = hp_paths.weight_sum;
                valid = check_hppath_validity(hp_paths, visited);
                if(valid) {
                    // continue to check repetition with lpaths and rpaths
                    for(auto& lp: lpaths) {
                        valid2 = true;
                        for(auto& lv: lp) {
                            if(visited[lv]) {
                                valid2 = false;
                                break;
                            }
                            visited2[lv] = true;
                        }
                        if(valid2) {
                            for(auto& rp: rpaths) {
                                // first check lp and rp don't have repetitive vertices
                                valid2 = true;
                                for(auto& rv: rp) {
                                    if(visited2[rv] || visited[rv])  {
                                        valid2 = false;
                                        break;
                                    }
                                }
                                if(valid2) {
                                    // check if satisfy length requirement
                                    if(cur_length+lp.size()+rp.size()<=length+1) {
                                        print_left_path(lp, fout);
                                        print_mid_path(hp_paths, fout);
                                        print_right_path(rp, fout);
                                        fout<<"\n";
                                        total_num += 1;
                                    }
                                }
                            }
                        }
                        // reset buffer 
                        for(auto& lv: lp) {
                            visited2[lv] = false;
                        }
                    }
                }
                for(auto& hpedge: hp_paths.hppath) {
                    for(auto& ve: hpedge.path) {
                        visited[ve] = false;
                    }
                }
            }
        }
    }
    cout<<"Total # of results: "<<total_num<<endl;
}

// case 2: start is not HP, end is HP
void Graph::findAllSimplePathsHPCase2(int start, int end, int length, ofstream &fout, int &total_num)  
{   
    int cur_length, left_min=mnum_of_vertices;
    bool valid, valid2;
    PathTable ltable;
    dfs_with_hp(start, end, length, ltable, fout, total_num, left_min, OUT);
    cout << "# path with no HP: "<<total_num<<endl; // total_num = 0;
     
    if(ltable.find(end) != ltable.end()) { 
        for(auto& lp: ltable[end]) {
            // satisfy length requirement
            if(lp.size() <= length) { // lp.size()+1 <= length+1
                print_left_path(lp, fout);
                fout<<end<<" ";
                fout<<"\n";
                total_num += 1;
            }
        }   
    }
    cout << "# path with at most 1 HP: "<<total_num<<endl;

    // [hp2->e][hp1->hp2] 
    // [hp1->hp2][hp2->e]
    HPPathVecMap ptable_idx;
    dfs_gidx_inside(end, ltable, length-left_min+1, ptable_idx, IN);
    for(auto& lpair: ltable) {
        auto lhp = lpair.first; 
        auto lpaths = lpair.second;
        if(lhp==end) continue; // already done by previous steps
        // check whether these subpaths have repetitive vertices
        for(auto& hp_paths: ptable_idx[lhp]) {
            reverse(hp_paths.hppath.begin(), hp_paths.hppath.end());
            cur_length = hp_paths.weight_sum;
            valid = check_hppath_validity(hp_paths, visited);
            if(valid) {
                // continue to check repetition with lpaths
                for(auto& lp: lpaths) {
                    valid2 = true;
                    for(auto& lv: lp) {
                        if(visited[lv]) {
                            valid2 = false;
                            break;
                        }
                    }
                    if(valid2) {
                        // check if satisfy length requirement
                        if(cur_length+lp.size()<=length+1) {
                            print_left_path(lp, fout);
                            print_mid_path(hp_paths, fout);
                            fout<<"\n";
                            total_num += 1;
                        }  
                    }
                }
            }
            for(auto& hpedge: hp_paths.hppath) {
                for(auto& ve: hpedge.path) {
                    visited[ve] = false;
                }
            }
        }
    }
    cout<<"Total # of results: "<<total_num<<endl;
}

// case 3: start is HP, end is not HP
void Graph::findAllSimplePathsHPCase3(int start, int end, int length, ofstream &fout, int &total_num)  
{   
    int cur_length, right_min=mnum_of_vertices;
    bool valid, valid2;
    PathTable rtable;
    dfs_with_hp(end, start, length, rtable, fout, total_num, right_min, IN); 
    cout << "# path with no HP: "<<total_num<<endl; // total_num = 0

    if(rtable.find(start) != rtable.end()) {
        for(auto& rp: rtable[start]) {
            // satisfy length requirement
            if(rp.size() <= length) {
                fout<<start<<" ";
                print_right_path(rp, fout);
                fout<<"\n"; 
                total_num += 1;
            }
        }
    }
    cout << "# path with at most 1 HP:  "<<total_num<<endl;

    HPPathVecMap ptable_idx;
    dfs_gidx_inside(start, rtable, length-right_min+1, ptable_idx, OUT);
    for(auto& rpair: rtable) {
        auto rhp = rpair.first; 
        auto rpaths = rpair.second;
        if(start==rhp) continue; // already done by previous steps
        // check whether these subpaths have repetitive vertices
        for(auto& hp_paths: ptable_idx[rhp]) {
            cur_length = hp_paths.weight_sum;
            valid = check_hppath_validity(hp_paths, visited);
            if(valid) {   
                for(auto& rp: rpaths) {
                    valid2 = true;
                    for(auto& rv: rp) {
                        if(visited[rv])  {
                            valid2 = false;
                            break;
                        }
                    }
                    if(valid2) {
                        // check if satisfy length requirement
                        if(cur_length+rp.size()<=length+1) {
                            print_mid_path(hp_paths, fout);
                            print_right_path(rp, fout);
                            fout<<"\n";
                            total_num += 1;
                        }
                    }
                }
            }
            for(auto& hpedge: hp_paths.hppath) {
                for(auto& ve: hpedge.path) {
                    visited[ve] = false;
                }
            }
        }
    }
    cout<<"Total # of results: "<<total_num<<endl;
}


void Graph::findAllSimplePathsHPCase4(int start, int end, int length, ofstream &fout, int &total_num)
{
    int cur_length;
    bool valid;

    cout << "This is Case 4."<<total_num<<endl; // total_num = 0;

    HPPathVec ptable_idx;
    dfs_gidx_inside(start, end, length+1, ptable_idx);
    // check whether these subpaths have repetitive vertices
    for(auto& hp_paths: ptable_idx) {
        cur_length = hp_paths.weight_sum;
        valid = check_hppath_validity(hp_paths, visited);    
        if(valid) {        
            // satisfy length requirement
            if(cur_length<=length+1) {
                print_mid_path(hp_paths, fout);
                fout<<"\n";
                total_num += 1;
            }  
        }
        for(auto& hpedge: hp_paths.hppath) {
            for(auto& ve: hpedge.path) {
                visited[ve] = false;
            }
        }
    }
    cout<<"Total # of results: "<<total_num<<endl;
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

/**
void Graph::findAllSimplePathsHP(int start, int end, int length) 
{   
    // if length > mmaxlen_threshold, we might miss some results. 
    assert(length <= mmaxlen_threshold);
    bool valid, valid2; 
    int cur_length = 0, total_num = 0, tmp; 
    ofstream fout;
    char file[100], no[100];
	strcpy(file, "path");
	sprintf(no, "_%d_%d", start, end);
	strcat(file, no);
	fout.open(file);
    
    PathTable ltable, rtable;
    // if start or end itself is HP
    if(mphp[start]) {
        PathVec vec(1);
        ltable[start] = vec;
    } else
        dfs_with_hp(start, end, length, ltable, fout, total_num, OUT);

    if(mphp[end]) {
        PathVec vec(1);
        rtable[end] = vec;
    } else
        dfs_with_hp(end, start, length, rtable, fout, tmp, IN); // path contains no HP
    
    cout << "# path with no HP: "<<total_num<<endl;

    // check intersection of ltable and rtable
    bool *check_dup = new bool[mnum_of_vertices];
    memset(check_dup, false, sizeof(bool)*mnum_of_vertices);
    for(auto& [lhp, lpaths]: ltable) {
        for(auto& [rhp, rpaths]: rtable) {
            if(lhp==rhp) {
                for(auto& lp: lpaths) {
                    // set buffer
                    for(auto& lv: lp) {
                        check_dup[lv] = true;
                    }
                    for(auto& rp: rpaths) {
                        valid = true;
                        for(auto& rv: rp) {
                            if(check_dup[rv]) {
                                valid = false;
                                break;
                            }
                        }
                        // satisfy length requirement
                        if(valid && lp.size()+rp.size()+1 <= length+1) {
                            for(int i=0; i<lp.size(); i++) fout<<lp[i]<<" ";
                            fout<<lhp<<" ";
                            for(int i=rp.size()-1; i>=0; i--) fout<<rp[i]<<" ";
                            fout<<"\n"; 
                            total_num += 1;
                        }
                    }
                    // reset buffer
                    for(auto& lv: lp) {
                        check_dup[lv] = false;
                    }
                }
            }
        } 
    }

    cout << "# path with at most 1 HP :  "<<total_num<<endl;


    // Step 3, left+mid+right table join
    memset(check_dup, false, sizeof(bool)*mnum_of_vertices); 
    for(auto& [lhp, lpaths]: ltable) {
        for(auto& [rhp, rpaths]: rtable) {
            if(lhp==rhp) continue; // already done by previous steps
            
            HPPathVec ptable_idx;
            dfs_gidx_inside(lhp, rhp, length, ptable_idx); 
            // check whether these subpaths have repetitive vertices
            for(auto& hp_paths: ptable_idx) {
                valid = true;
                cur_length = hp_paths.weight_sum;
                for(int i=0; i<hp_paths.path.size(); i++) {
                    for(int j=0; j<hp_paths.path[i].path.size(); j++) {
                        int ve = hp_paths.path[i].path[j];
                        if(i!=0 && j==0) continue;
                        if(visited[ve]) {
                            valid = false;
                            break;
                        }
                        visited[ve] = true;
                    }
                    if(!valid) break;
                }
                if(valid) {
                    // continue to check repetition with lpaths and rpaths
                    for(auto& lp: lpaths) {
                        valid2 = true;
                        for(auto& lv: lp) {
                            if(visited[lv]) {
                                valid2 = false;
                                break;
                            }
                            check_dup[lv] = true;
                        }
                        if(valid2) {
                            for(auto& rp: rpaths) {
                                // first check lp and rp don't have repetitive vertices
                                valid2 = true;
                                for(auto& rv: rp) {
                                    if(check_dup[rv] || visited[rv])  {
                                        valid2 = false;
                                        break;
                                    }
                                }
                                if(valid2) {
                                    // satisfy length requirement
                                    if(cur_length+lp.size()+rp.size()<=length+1) {
                                        //cout<<"Found One ..."<<endl;
                                        for(int i=0; i<lp.size(); i++) fout<<lp[i]<<" ";
                                        for(int i=0; i<hp_paths.path.size(); i++) {
                                            for(int j=0; j<hp_paths.path[i].path.size(); j++) {
                                                if(i!=0 && j==0) continue; 
                                                fout<<hp_paths.path[i].path[j]<<" ";
                                            }
                                        }
                                        for(int i=rp.size()-1; i>=0; i--) fout<<rp[i]<<" ";
                                        fout<<endl;
                                        total_num += 1;
                                    }
                                }
                            }
                        }
                        // reset buffer 
                        for(auto& lv: lp) {
                            check_dup[lv] = false;
                        }
                    }
                }
                for(auto& hpedge: hp_paths.path) {
                    for(auto& ve: hpedge.path) {
                        visited[ve] = false;
                    }
                }
            }
        }
    }
    cout<<"Total # of results: "<<total_num<<endl;
    fout.close();
}
**/

// total_num += # of newly found paths
void Graph::dfs(int s, int d, int length, ofstream &fout, int &total, int direction)
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
        total += 1;
    } else {
        if(path.size() <= length) { 
            if(direction == OUT) {
                for(int i=1; i<=mppadj_lists_o[s][0]; i++) {
                    int ne = mppadj_lists_o[s][i];
                    if(!visited[ne]) {
                        dfs(ne, d, length, fout, total, OUT);
                    }
                }
            } else {
                for(int i=1; i<=mppadj_lists_i[s][0]; i++) {
                    int ne = mppadj_lists_i[s][i];
                    if(!visited[ne]) {
                        dfs(ne, d, length, fout, total, IN);
                    }
                }
            }
        }
    }
    path.pop_back(); 
    visited[s] = false;
}


// ptable: L/R and p(h)
void Graph::dfs_with_hp(int s, int d, int length, PathTable &ptable, ofstream &fout, int &total, int &min_len, int direction)
{
    visited[s] = true; 
    path.push_back(s);
    if(s == d && direction == OUT) { 
        // report simple path, directly path find, w/o HP-index
        for(int i=0; i<path.size(); i++)
        {   
            fout<<path[i]<<" ";
        }
        fout<<"\n";
        total+=1;
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
                        dfs_with_hp(ne, d, length, ptable, fout, total, min_len, OUT);
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
                        dfs_with_hp(ne, d, length, ptable, fout, total, min_len, IN);
                    }
                }
            } 
        }
    }
    path.pop_back(); 
    visited[s] = false;
}

void Graph::dfs_build_gidx(int s, int cur, unordered_set<int> &HP_set) 
{
    visited[cur] = true; 
    path.push_back(cur);
    // s, d, middle: 6-2=4=5-1, [s, x1,x2,x3,x4, d] is allowed
    // [s, hp, x1, x2, hp, d] 
    // [ltable=null, s=hp, x1, x2, x3, hp, d]
    // [s, hp, x1, x2, x3, d=hp, rtable=null] 
    // [ltable=null, s=hp, x1, x2, x3, x4, d=hp, rtable=null] 

    // if start and end aren't HP, total_weight <= mlen_threshold-1
    // if either start or end is HP(not both),  total_weight <= mlen_threshold
    // if start and end both are HP, total_weight <= mlen_threshold+1
    // For ease of computation, we can just set total_weight <= mlen_threshold-1
    // for line 493 and 501
    int threshold = mmaxlen_threshold; // TODO:+1/-1/0
    if(s!=cur && HP_set.find(cur)!=HP_set.end() && path.size()<=threshold) { // prevent self-loop
        HPEdge e(s, cur, path); 
        if(mhp_adjlist.find(s) == mhp_adjlist.end()) {
            mhp_adjlist[s] = HPEdgeSet();
        }
        mhp_adjlist[s].push_back(e);
        // add reverse HPEdge
        HPEdge er(cur, s, path); 
        if(mhp_adjlist_i.find(cur) == mhp_adjlist_i.end()) {
            mhp_adjlist_i[cur] = HPEdgeSet();
        }
        mhp_adjlist_i[cur].push_back(er);
    } else {
        if(path.size() <= threshold) {
            for(int i=1; i<=mppadj_lists_o[cur][0]; i++) {
                int ne = mppadj_lists_o[cur][i];
                if(!visited[ne]) {
                    dfs_build_gidx(s, ne, HP_set);
                }
            }
        }
    }
    path.pop_back(); 
    visited[cur] = false;
}

void Graph::build_graph_idx(int maxlen_threshold, int hp_deg_threshold)
{   
    mmaxlen_threshold = maxlen_threshold;
    mhp_deg_threshold = hp_deg_threshold;
    mphp.resize(mnum_of_vertices, false);
    mnum_hp = 0;

    // mark hot points by their outdegree
    for(int i=0; i<mnum_of_vertices; i++) {
        if(mppadj_lists_o[i][0] >= mhp_deg_threshold) {
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
    for(int i=0; i<mvec_hp.size(); i++) {
        dfs_build_gidx(mvec_hp[i], mvec_hp[i], mset_hp); 
        sort(mhp_adjlist[mvec_hp[i]].begin(), mhp_adjlist[mvec_hp[i]].end(), vec_sort);
    }
    // sort mhp_adjlist_i w.r.t. weight
    for(int i=0; i<mvec_hp.size(); i++) {
        sort(mhp_adjlist_i[mvec_hp[i]].begin(), mhp_adjlist_i[mvec_hp[i]].end(), vec_sort);
    }

    // cout<<"*** Info of G_idx ***" << endl;
    int total_p = 0;
    for(auto& pair: mhp_adjlist) {
        auto hp = pair.first;
        auto adj = pair.second;
        /*
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
    cout<<"# edges in G_idx: "<<total_p<<endl;
}

// wrapper function
void Graph::dfs_gidx_inside(int s, int d, int npoints, HPPathVec &plist) // npoints: # of points
{
    HPEdgeSet path;
    dfs_gidx_inside(s, d, npoints, 0, path, plist); 
}

// output: plist
// adlength: adjusted length
void Graph::dfs_gidx_inside(int s, int d, int adlength, int total_weight, HPEdgeSet &path, HPPathVec &plist) 
{   
    visited[s] = true;
    if(s == d) {
        HPPath hpp(path, total_weight);
        plist.push_back(hpp);
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
                    dfs_gidx_inside(ne.v, d, adlength, new_weight, path, plist); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
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

void Graph::dfs_gidx_inside(int s, PathTable &table, int npoints, HPPathVecMap &plist, int direction)
{
    HPEdgeSet path;
    dfs_gidx_inside(s, s, table, npoints, 0, path, plist, direction); 
}

void Graph::dfs_gidx_inside(int s, int cur, PathTable &table, int adlength, int total_weight, HPEdgeSet &path, HPPathVecMap &plist, int direction)
{
    visited[cur] = true;
    if(s!=cur && table.find(cur)!=table.end()) { // prevent self-loop
        HPPath hpp(path, total_weight);
        if(plist.find(cur) == plist.end()) 
            plist[cur] = HPPathVec();
        plist[cur].push_back(hpp); // divided by end point
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
                    dfs_gidx_inside(s, ne.v, table, adlength, new_weight, path, plist, OUT); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
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
                    dfs_gidx_inside(s, ne.v, table, adlength, new_weight, path, plist, IN); // path1:[0,2,4] path2:[4,6,10]->[0,2,4,6,10]
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