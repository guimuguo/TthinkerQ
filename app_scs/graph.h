#include "Utility.h"
#include "ListLinearHeap.h"
#include <fstream>
#include <time.h>
#include <sys/timeb.h>
#include "../system/serialization.h"
#include "../system/rwlock.h"
#include "../system/conmap.h"
#include "data/format_graph.h"
#include<map>

// global variables
string graph_file;
map<string, string> config;

vector<ui> G; // vector for mapping index back to id
hash_map<ui, ui> outer_id2index; // id2index during loading the graph

bool EXE_heu2;
bool EXE_heu3;
bool EXE_heu4;

bool EXE_ub1;
bool EXE_ub2;
bool EXE_ub3;
bool EXE_ub3_optimization;
bool EXE_core_maintenance;
bool EXE_new2VI; // inclusion-based reduction
bool EXE_del_from_VR; // degree-based reduction
bool EXE_dom_ustar; // domination based branching
ui domS_Threshold;
double MaxTime;
ui srch_ord;
bool preprocess_graph;


ofbinstream & operator>>(ofbinstream & m, ui & i)
{
    i = *(ui*)m.raw_bytes(sizeof(ui));
    return m;
}

ifbinstream & operator<<(ifbinstream & m, ui i)
{
    m.raw_bytes(&i, sizeof(ui));
    return m;
}

// forward declaring ContextValue
struct ContextValue;


struct SCSQuery{
    ui QID;
	ui N1;
	ui N2;
	ui kl;
	vector<ui> H;
	ui ubD;
    rwlock hlock; // lock for updating query parameters
    vector<ui> counters;
    struct timeb start_t, end_t;

    SCSQuery(){
        counters.assign(32, 0);
    }

    ~SCSQuery(){
    
    }
};


// Size-Bounded Community Search
class SCS{ 
    public:

        ui num_vertices;
        ui num_edges;
        ui dMAX;
        ui ku;
        bool gUB;

        ui num_VIVR;

        ui * pstart;
        ui * peel_sequence;
        ui * degree;
        ui * core;
        ui * q_dist;

        vector<ui> G0;
        ui * G0_edges;
        ui * G0_x;

        hash_map<ui, ui> id2index; // gives new id based on key(old id)

        struct timeb gtime_start;

        void read_config(char * config_file);
        void load_graph();
        
        void core_decomposition_linear_list();

        void cal_query_dist(ui QID);

        void heu2(vector<ui> & H2, ui & kl2, SCSQuery & q);
        void heu3(vector<ui> & H3, ui & kl3, SCSQuery & q);
        void heu4(vector<ui> & H4, ui & kl4, SCSQuery & q);
        void CSSC_heu(SCSQuery & q);

        void reduction_g(SCSQuery & q, SCS &split_g);

        int find_ustar(ContextValue &t);
        int find_ustar_mindeg(ContextValue &t);
        int find_ustar_2phase(ContextValue &t, ui N2);
        int find_ustar_link(ContextValue &t);
        int find_ustar_random(ContextValue &t);

        ui get_ub1(ContextValue &t, ui N2);
        ui get_ub2(ContextValue &t, ui N2);
        ui get_ub3(ContextValue &t, ui kl, ui N2);
        bool compute_ub(ContextValue &t, ui kl, ui N2);

        void get_rVR(ContextValue &t, ui kl, bool & del_v_in_VI, vector<ui> & rVR, queue<ui> & Q, bool * inQ);
        void get_NEI_of_VI(ContextValue &t);
        void add_to_VI_from_VR(ContextValue &t, ui vertex);
        void pop_from_VI_to_VR(ContextValue &t, ui vertex);
        void pop_from_VI(ContextValue &t, ui vertex);
        void restore_in_VR(ContextValue &t, ui vertex);
        void remove_domS_from_VR(ContextValue &t, vector<pair<double, ui>> domS);
        void find_domS_of_ustar(ui ustar, ContextValue &t, vector<pair<double, ui>> & domS);

        void core_maintenance(ContextValue &t, ui kl, bool & del_v_in_VI, vector<ui> & rVR, ui ustar);
        void core_maintenance(ContextValue &t, ui kl, bool & del_v_in_VI, vector<ui> & rVR, vector<ui> & del_vec);
        void core_maintenance(ContextValue &t, ui kl, bool & del_v_in_VI, vector<ui> & rVR, vector<pair<double, ui>> & del_vec);

        ui find_ubD(ui kl, ui ubD, ui N2);
        void inclusion_based_reduction(ContextValue &t, unordered_set<ui> & new2VI, ui kl);
        void degree_based_reduction(ContextValue &t, vector<ui> & del_from_VR, ui kl, ui N2);
        void undo_degree_based_reduction(vector<ui> & del_from_VR, ContextValue &t);

        void ShrinkGraph(ContextValue &t, SCS & split_g);
        void DestroySplitGraph();
    
        // for SCS object
        friend ifbinstream & operator<<(ifbinstream & m, const SCS & g){
            ui num_vertices = g.num_vertices;
            m << num_vertices;

            ui num_edges = g.num_edges;
            m << num_edges;

            m << g.gUB;

            // m << g.id2index;

            m << g.G0;
            m << g.num_VIVR;

            for(ui i=0; i < num_vertices+1; i++){
                m << g.pstart[i];
            }

            for(ui i=0; i < 2*num_edges; i++){
                m << g.G0_edges[i];
            }

            for(ui i=0; i < num_vertices; i++){
                m << g.G0_x[i];
            }

            for(ui i=0; i < num_vertices; i++){
                m << g.q_dist[i];
            }

            return m;
        }


        friend ofbinstream & operator>>(ofbinstream & m, SCS & g){
            ui num_vertices;
            m >> num_vertices;
            g.num_vertices = num_vertices;

            ui num_edges;
            m >> num_edges;
            g.num_edges = num_edges;

            m >> g.gUB;

            // hash_map<ui, ui> id2index;
            // m >> id2index;
            // g.id2index = id2index;

            m >> g.G0;
            m >> g.num_VIVR;

            ui * m_pstart = new ui [num_vertices+1];
            for(ui i=0; i < num_vertices+1; i++){
                m >> m_pstart[i];
            }
            g.pstart = m_pstart;

            ui * m_G0_edges = new ui [2*num_edges];
            for(ui i=0; i < 2*num_edges; i++){
                m >> m_G0_edges[i];
            }
            g.G0_edges = m_G0_edges;

            ui * m_G0_x = new ui [num_vertices];
            for(ui i=0; i < num_vertices; i++){
                m >> m_G0_x[i];
            }
            g.G0_x = m_G0_x;

            ui * m_q_dist = new ui [num_vertices];
            for(ui i=0; i < num_vertices; i++){
                m >> m_q_dist[i];
            }
            g.q_dist = m_q_dist;

            return m;
        }

};


struct ContextValue{
    SCS split_g;

    vector<ui> VI;
    vector<ui> NEI;

    hash_set<ui> inVI;
    hash_set<ui> inVR;

    ui * degVI;
    ui * degVIVR;
    ui * inNEI;

    ContextValue(){
	}

	~ContextValue(){
        split_g.DestroySplitGraph();

        if(degVI != NULL) delete [] degVI;
        if(degVIVR != NULL) delete [] degVIVR;
        if(inNEI != NULL) delete [] inNEI;
	}
};


ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{
	m >> c.split_g;
	m >> c.VI;

    m >> c.inVI;
    m >> c.inVR;

    ui * m_degVI = new ui [c.split_g.num_vertices];
    for(ui i=0; i < c.split_g.num_vertices; i++){
        m >> m_degVI[i];
    }
    c.degVI = m_degVI;

    ui * m_degVIVR = new ui [c.split_g.num_vertices];
    for(ui i=0; i < c.split_g.num_vertices; i++){
        m >> m_degVIVR[i];
    }
    c.degVIVR = m_degVIVR;

    ui * m_inNEI = new ui [c.split_g.num_vertices];
    for(ui i=0; i < c.split_g.num_vertices; i++){
        m >> m_inNEI[i];
    }
    c.inNEI = m_inNEI;

    m >> c.NEI;

    return m;
}

ifbinstream & operator<<(ifbinstream & m, const ContextValue & c)
{
	m << c.split_g;
	m << c.VI;

    m << c.inVI;
    m << c.inVR;

    for(ui i=0; i < c.split_g.num_vertices; i++){
        m << c.degVI[i];
    }

    for(ui i=0; i < c.split_g.num_vertices; i++){
        m << c.degVIVR[i];
    }

    for(ui i=0; i < c.split_g.num_vertices; i++){
        m << c.inNEI[i];
    }

    m << c.NEI;

    return m;
}


void SCS::read_config(char * config_file){
    ifstream configFile(config_file, ios::in);
    
    string args;
    string val;

    if (!configFile.is_open()){
        cout << "Graph file Open Failed "<<endl;
        exit(1);
    }
    else{
        while (configFile >> args >> val){
            config[args] = val;
        }
        configFile.close();
    }

    graph_file = config["Graph_file"];
    
    EXE_heu2 = stoi(config["EXE_heu2"]); //Heuristic strategy 1
    EXE_heu3 = stoi(config["EXE_heu3"]); //Heuristic strategy 2
    EXE_heu4 = stoi(config["EXE_heu4"]); //Heuristic strategy 3

    EXE_ub1 = stoi(config["EXE_ub1"]); //UB1
    EXE_ub2 = stoi(config["EXE_ub2"]); //UB2
    EXE_ub3 = stoi(config["EXE_ub3"]); //UB3
    EXE_ub3_optimization = stoi(config["EXE_ub3_optimization"]); //UB3 optimization
    
    EXE_core_maintenance = stoi(config["EXE_core_maintenance"]); //reduction rule 1
    EXE_new2VI = stoi(config["EXE_new2VI"]); //reduction rule 2: Inclusion based reduction
    EXE_del_from_VR = stoi(config["EXE_del_from_VR"]); //reduction rule 3: Degree based reduction
    
    EXE_dom_ustar = stoi(config["EXE_dom_ustar"]); //Dominating based branching rule
    
    domS_Threshold = stoi(config["domS_Threshold"]); //Dom pair threshold, upper bound for vertices being dominated by a particular vertex
    MaxTime = stoi(config["MaxTime"]); // maximum running time for timeout
    
    srch_ord = stoi(config["srch_ord"]); //Branching order, selecting which vertex to branch on

    preprocess_graph = stoi(config["preprocess_graph"]); // whether we need to convert the graph format or not

    cout<<"    Heu: "<<EXE_heu2<<","<<EXE_heu3<<","<<EXE_heu4;
    cout<<"    UBs: "<<EXE_ub1<<","<<EXE_ub2<<","<<EXE_ub3<<","<<EXE_ub3_optimization;
    cout<<"    Rdt: "<<EXE_core_maintenance<<","<<EXE_new2VI<<","<<EXE_del_from_VR;
    cout<<"    Dom: "<<EXE_dom_ustar<<","<<domS_Threshold;
    cout<<"    MaxTime: "<<MaxTime;
    cout<<"    Ord: "<<srch_ord;
    cout << "   preprocess_graph: " << preprocess_graph;
    cout<<endl;
    
}


void SCS::load_graph(){
    // preprocess graph file and convert into required format
    if(preprocess_graph){
        FormatGraph fg;
        fg.loadGraphFromFile(const_cast<char*>(graph_file.c_str()));
        string outfile_path = graph_file + "_formatted";
        fg.graphDump(const_cast<char*>(outfile_path.c_str()));
        graph_file = outfile_path;
        cout << "Graph file format converted!" << endl;
    }

    map<ui, set<ui>> adj_map;
    set<ui> vertices;
    ifstream inputFile(graph_file, ios::in);

    if (!inputFile.is_open())
    {
        cout << "Graph file Open Failed "<<endl;
        exit(1);
    }
    else
    {
        ui tu,tv;
        while (inputFile >> tu >> tv)
        {
            vertices.insert(tu);
            vertices.insert(tv);

            if(tu==tv) continue;
            
            // count the number of G0_edges
            if(adj_map[tv].find(tu) == adj_map[tv].end()) num_edges++;

            adj_map[tu].insert(tv);
            adj_map[tv].insert(tu);
        }
        inputFile.close();
    }

    num_vertices = vertices.size();
    
    pstart = new ui[num_vertices+1];
    G0_edges = new ui[2*num_edges];

    pstart[0] = 0;

    ui * temp_edges = new ui[2*num_edges];
    ui new_index = 0;

    for(auto g: adj_map){
        ui id = g.first;
        G.push_back(id);

        ui j = 0;
        outer_id2index[id] = new_index;
        for(ui nei : g.second){
            temp_edges[pstart[new_index]+j] = nei;j++;
        }
        pstart[new_index+1] = pstart[new_index] + g.second.size();
        new_index++;
    }

    // get the new index for edges after recoding
    for(int i=0; i<2*num_edges; i++){
        G0_edges[i] = outer_id2index[temp_edges[i]];
    }

    // for(ui i =0;i<num_vertices;i++){
    //     if(adj_map.find(i)!=adj_map.end()){
    //         ui j = 0;
    //         for(ui nei : adj_map[i]){ //// nei: neighbors

    //             G0_edges[pstart[i]+j] = nei;j++;
    //         }
    //         pstart[i+1] = pstart[i] + adj_map[i].size();
    //     }
    //     else{
    //         pstart[i+1] = pstart[i];
    //     }
    // }

    peel_sequence = new ui[num_vertices];
    degree = new ui[num_vertices];
    core = new ui[num_vertices];
    dMAX = 0;

    for(auto g: adj_map){
        ui id = g.first;
        ui new_index = outer_id2index[id];
        peel_sequence[new_index] = new_index;

        degree[new_index] = g.second.size();
        if(degree[new_index] > dMAX) dMAX = degree[new_index];
    }

    // for(ui i = 0; i<num_vertices; i++){
    //     peel_sequence[i] = i;
    //     if(adj_map.find(i)!=adj_map.end()){
    //         degree[i] = adj_map.find(i)->second.size(); //// value of G[i]'s size
    //         if(degree[i]>dMAX) dMAX = degree[i];
    //     }
    //     else degree[i] = 0;
    // }

    cout<<"num_vertices="<<num_vertices<<",num_edges="<<num_edges<<",dMAX="<<dMAX<<endl;

    delete [] temp_edges;
}


void SCS::core_decomposition_linear_list(){
    ui max_core = 0;
    ListLinearHeap *linear_heap = new ListLinearHeap(num_vertices, num_vertices-1);
    linear_heap->init(num_vertices, num_vertices-1, peel_sequence, degree); //// initialize doubly linked list and insert initial elements
    memset(core, 0, sizeof(ui)*num_vertices);
    for(ui i = 0;i < num_vertices;i ++) {
        ui u, key;
        linear_heap->pop_min(u, key);
        if(key > max_core) max_core = key;
        peel_sequence[i] = u;
        core[u] = max_core; //// get core number of each vertex

        for(ui j = pstart[u];j < pstart[u+1];j ++) if(core[G0_edges[j]] == 0) {
            linear_heap->decrement(G0_edges[j]);
        }
    }
    delete linear_heap;
}


// Calculate hop distance from query node to each node
void SCS::cal_query_dist(ui QID)
{
    q_dist = new ui[num_vertices];
    for(ui i =0;i<num_vertices;i++)
        q_dist[i] = INF;

    queue<ui> Q;
    q_dist[QID] = 0;
    Q.push(QID);
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        for(ui i = pstart[v]; i < pstart[v+1]; i++){
            ui w = G0_edges[i];
            if(q_dist[w] == INF){
                q_dist[w] = q_dist[v] + 1;
                Q.push(w);
            }
        }
    }
}


void SCS::heu2(vector<ui> & H2, ui & kl2, SCSQuery & q)
{
    H2.clear();
    kl2 = 0;
    
    vector<ui> tH;
    ui tkl = 0;
    ui tH_size = 0;
    
    ui * sta = new ui[num_vertices];
    memset(sta, 0, sizeof(ui)*num_vertices);
    
    ui * deg = new ui[num_vertices];
    memset(deg, 0, sizeof(ui)*num_vertices);
    
    priority_queue<pair<ui, ui>> Q;
    Q.push(make_pair(0, q.QID));
    sta[q.QID] = 1;
    
    while(!Q.empty()){
        ui v = Q.top().second;

        Q.pop();
        if(sta[v] == 2) continue;
        tH.push_back(v);

        sta[v] = 2;
        for(ui j = pstart[v]; j < pstart[v+1]; j++){
            ui nei = G0_edges[j]; // neighbor of v(v is in tH)
            if(sta[nei] == 0){
                ui d = 0;
                for(ui w = pstart[nei]; w < pstart[nei+1]; w++){
                    if(sta[G0_edges[w]] == 2) ++ d; // if neighbor of neighbor is in tH, increase degree in tH
                }
                Q.push(make_pair(d, nei));
                sta[nei] = 1;
            }
            else{
                if(sta[nei] == 1){
                    ui new_d = 0;
                    for(ui w = pstart[nei]; w < pstart[nei+1]; w++){
                        if(sta[G0_edges[w]] == 2) ++ new_d;
                    }
                    Q.push(make_pair(new_d, nei));
                }
                else{
                    ++ deg[nei]; // if the neighbor is already in tH(sta[2] then just increase the degree of nei and v.
                    ++ deg[v];
                }
            }
        }
        if(tH.size()>=q.N1){
            ui mindeg = INF;
            for(ui i = 0; i < tH.size(); i++){
                if(deg[tH[i]] < mindeg)
                    mindeg = deg[tH[i]];
            }
            if(mindeg >= tkl){
                tkl = mindeg;
                tH_size = (ui)tH.size();
            }
        }
        if(tH.size()==q.N2)
            break;
    }
    
    for(ui i=0; i<tH_size; i++){
        H2.push_back(tH[i]);
    }
    kl2 = tkl;
    
    delete [] sta;
    delete [] deg;
}


//// uses connection score in VI and degree in original graph to select the vertices
void SCS::heu3(vector<ui> & H3, ui & kl3, SCSQuery & q)
{
    H3.clear();
    kl3 = 0;
    
    vector<ui> tH;
    ui tkl = 0;
    ui tH_size = 0;
    
    ui * sta = new ui[num_vertices];
    memset(sta, 0, sizeof(ui)*num_vertices);
    
    ui * deg = new ui[num_vertices];
    memset(deg, 0, sizeof(ui)*num_vertices);
    
    priority_queue<pair<double, ui>> Q;
    Q.push(make_pair(0, q.QID));
    sta[q.QID] = 1;
    
    while (!Q.empty()) {
        ui v = Q.top().second;
        Q.pop();
        if(sta[v] == 2) continue;
        tH.push_back(v);
        sta[v] = 2;
        for(ui nei = pstart[v]; nei<pstart[v+1]; nei++){
            if(sta[G0_edges[nei]] == 2){
                ++ deg[G0_edges[nei]]; // degree in tH
                ++ deg[v];
            }
        }
        for(ui nei = pstart[v]; nei<pstart[v+1]; nei++){
            if(!sta[G0_edges[nei]]){
                double score = 0;
                for(ui w = pstart[G0_edges[nei]]; w < pstart[G0_edges[nei]+1]; w++){
                    if(sta[G0_edges[w]] == 2 && deg[G0_edges[w]] != 0){
                        score += (double) 1/deg[G0_edges[w]]; // degree in tH
                    }
                }
                score += (double) degree[G0_edges[nei]]/dMAX; // degree in original graph
                Q.push(make_pair(score, G0_edges[nei]));
                sta[G0_edges[nei]] = 1;
            }
            else{
                if(sta[G0_edges[nei]] == 1){
                    double new_score = 0;
                    for(ui w = pstart[G0_edges[nei]]; w < pstart[G0_edges[nei]+1]; w++){
                        if(sta[G0_edges[w]] == 2 && deg[G0_edges[w]] != 0){
                            new_score += (double) 1/deg[G0_edges[w]];
                        }
                    }
                    new_score += (double) degree[G0_edges[nei]]/dMAX;
                    Q.push(make_pair(new_score, G0_edges[nei]));
                }
            }
        }
        if(tH.size()>=q.N1){
            ui mindeg = INF;
            for(ui i = 0; i < tH.size(); i++){
                if(deg[tH[i]] < mindeg)
                    mindeg = deg[tH[i]];
            }
            if(mindeg >= tkl){
                tkl = mindeg;
                tH_size = (ui)tH.size();
            }
        }
        if(tH.size()==q.N2)
            break;
    }
    
    for(ui i=0; i<tH_size; i++){
        H3.push_back(tH[i]);
    }
    kl3 = tkl;
    
    delete [] sta;
    delete [] deg;
}



//// core decomposition on the subgraph induced by QID and its 1-hop neighbors
void SCS::heu4(vector<ui> & H4, ui & kl4, SCSQuery & q)
{
    H4.clear();
    kl4 = 0;
    
    if( degree[q.QID] < q.N1-1){
        return;
    }
    
    vector<ui> S;
    bool * inS = new bool[num_vertices];
    memset(inS, 0, sizeof(bool)*num_vertices);
    
    ui * deg = new ui[num_vertices];
    memset(deg, 0, sizeof(ui)*num_vertices);
    
    S.push_back(q.QID);
    inS[q.QID] = 1;
    ++ deg[q.QID]; // degree in S
    
    for(ui j = pstart[q.QID]; j < pstart[q.QID+1]; j++){
        S.push_back(G0_edges[j]); // S is the subgraph induced by QID and its 1-hop neighbors
        inS[G0_edges[j]] = 1;
        for(ui w = pstart[G0_edges[j]]; w < pstart[G0_edges[j]+1]; w++){
            if(inS[G0_edges[w]]){
                ++ deg[G0_edges[w]];
                ++ deg[G0_edges[j]];
            }
        }
    }

    priority_queue<pair<ui,ui>, vector<pair<ui,ui>>, greater<pair<ui, ui>>> Q;
    for(auto e : S)
        Q.push(make_pair(deg[e], e));
    vector<ui> rmv;
    ui rmv_idx = 0;
    ui mindeg = 0;
    while (!Q.empty()) {
        ui d = Q.top().first;
        ui v = Q.top().second;
        Q.pop();
        if(inS[v] == 0) continue;
        ui remain = (ui)S.size() - (ui)rmv.size();
        if(remain >= q.N1 && remain <= q.N2){
            if(d > mindeg){
                mindeg = d;
                rmv_idx = (ui)rmv.size();
            }
        }
        if(remain == q.N1) break;

        inS[v] = 0;
        rmv.push_back(v);
        for(ui j = pstart[v]; j < pstart[v+1]; j++){
            if(inS[G0_edges[j]]){
                -- deg[G0_edges[j]];
                Q.push(make_pair(deg[G0_edges[j]], G0_edges[j]));
            }
        }
    }
    
    memset(inS, 0, sizeof(bool)*num_vertices);
    for(auto e : S)
        inS[e] = 1;
    for(ui i = 0; i < rmv_idx; i++)
        inS[rmv[i]] = 0; // remove the vertices from S
    
    vector<ui> tH;
    ui tkl = mindeg;
    for(auto e : S)
        if(inS[e]) tH.push_back(e);
    
    for(ui i=0; i<tH.size(); i++){
        H4.push_back(tH[i]);
    }
    kl4 = tkl;
    
    delete [] inS;
    delete [] deg;
    
}

void SCS::CSSC_heu(SCSQuery & q){
    q.H.clear();
    q.kl = 0;
    
    vector<ui> H2; ui kl2 = 0;
    
    if(EXE_heu2) heu2(H2, kl2, q);

    vector<ui> H3; ui kl3 = 0;
    if(EXE_heu3) heu3(H3, kl3, q);

    vector<ui> H4; ui kl4 = 0;
    if(EXE_heu4) heu4(H4, kl4, q);
    
    if(kl4 >= kl3 && kl4 >= kl2){
        // q.H = H4;
        for(auto v: H4){
            ui original_id = G[v];
            q.H.push_back(original_id);
        }
        q.kl = kl4;
    }
    else{
        if(kl3 >= kl2){
            // q.H = H3;
            for(auto v: H3){
                ui original_id = G[v];
                q.H.push_back(original_id);
            }
            q.kl = kl3;
        }
        else{
            // q.H = H2;
            for(auto v: H2){
                ui original_id = G[v];
                q.H.push_back(original_id);
            }
            q.kl = kl2;
        }
    }
}


void SCS::reduction_g(SCSQuery & q, SCS &split_g){
    cout << "Number of vertices before reduction: " << num_vertices << endl;
    cout << "Number of edges before reduction: " << num_edges << endl;

    vector<ui> G0_edges_temp; // we don't know what size these should be of yet
    vector<ui> G0_x_temp;
    vector<ui> pstart_temp; 
    
    bool * inQ = new bool[num_vertices];
    memset(inQ, 0, sizeof(bool)*num_vertices);
    
    queue<ui> Q;
    Q.push(q.QID);
    inQ[q.QID] = 1;

    ui new_index = 0;
    ui shrinked_total_edges = 0;

    pstart_temp.push_back(0);

    split_g.G0.clear();
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        split_g.G0.push_back(v);

        split_g.id2index[v] = new_index; // v is the original id

        ui new_deg = 0;

        for(ui i = pstart[v]; i < pstart[v+1]; i++){
            ui w = G0_edges[i]; // here w is original id
            if(core[w] > q.kl && q_dist[w] <= q.ubD){
                G0_edges_temp.push_back(w);
                ++new_deg;
                ++shrinked_total_edges;

                if(!inQ[w]){
                    Q.push(w);
                    inQ[w] = 1;
                }
            }
        }

        pstart_temp.push_back(pstart_temp[new_index] + new_deg);
        G0_x_temp.push_back(new_deg);
        ++ new_index;
    }

    // update num_vertices and num_edges after reduction
    split_g.num_vertices = new_index;
    split_g.num_edges = shrinked_total_edges / 2;
    split_g.num_VIVR = new_index;

    // update query id
    q.QID = split_g.id2index[q.QID];

    split_g.pstart = new ui[new_index+1]; // new_index is the total number of vertices after reduction
    split_g.G0_edges = new ui[shrinked_total_edges];
    split_g.G0_x = new ui[new_index];

    // copy from vector to array
    memcpy(split_g.pstart, pstart_temp.data(), sizeof(ui)*(split_g.num_vertices+1));
    memcpy(split_g.G0_x, G0_x_temp.data(), sizeof(ui)*(split_g.num_vertices));

    // get new index of id in G0_edges after recoding
    for(int i=0; i<shrinked_total_edges; i++){
        split_g.G0_edges[i] = split_g.id2index[G0_edges_temp[i]];
    }

    cout << "Number of vertices after reduction: " << split_g.num_vertices << endl;
    cout << "Number of edges after reduction: " << split_g.num_edges << endl;

    delete [] inQ;
}

// find the upper bound of Diameter based on optimal min-degree
ui SCS::find_ubD(ui kl, ui ubD, ui N2){
    for(ui d = 1; d <= N2; d++){
        if(d == 1 || d == 2){
            if(kl + d > N2){
                ubD = d-1;
                return ubD;
            }
        }
        else{
            ui min_n = kl + d + 1 + floor(d/3) * (kl - 2);
            if(N2 < min_n){
                ubD = d - 1;
                return ubD;
            }
        }
    }
    return ubD;
}


void SCS::inclusion_based_reduction(ContextValue &t, unordered_set<ui> & new2VI, ui kl){
    for(auto e : t.VI){
        if(t.degVIVR[e] == kl+1){
            vector<ui> its_nei;
            for(ui i = pstart[e]; i < pstart[e] + G0_x[e]; i++){
                ui w = G0_edges[i];
                if(t.inVR.find(w) != t.inVR.end()){
                    its_nei.push_back(w); //// add neighbors to its_nei vector
                }
            }
            if(its_nei.size() != 0){
                for(auto x : its_nei){
                    new2VI.insert(x);
                }
            }
        }
    }
    for(auto e : new2VI){
        if(t.inVR.find(e) != t.inVR.end()){
            t.inVI.insert(e);
            t.inVR.erase(e);
            t.VI.push_back(e); //// greedily move the neighbors to VI
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++ ){
                if(t.inVI.find(G0_edges[i]) != t.inVI.end()){
                    ++ t.degVI[G0_edges[i]]; //// increase the degree of neighbors of e in VI
                    ++ t.degVI[e]; //// increase the degree of e in VI
                }
            }
        }
    }

    return;
}


void SCS::degree_based_reduction(ContextValue &t, vector<ui> & del_from_VR, ui kl, ui N2){
    for(auto e : t.NEI){

        if(t.inNEI[e] < kl+1){
            int lack = kl + 1 - t.inNEI[e]; //// lacking degree to be a better degree vertex in VI
            int bugt = N2 - (int)t.VI.size() - 1; //// the vertices adjacent to e that can still be added to VI
            if( lack > bugt ){ //// even after adding all the possible vertices adjacent to e, its not any better
                del_from_VR.push_back(e);
                t.inVR.erase(e);
                t.inNEI[e] = 0;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if((t.inVR.find(w) != t.inVR.end()) || (t.inVI.find(w) != t.inVI.end())){
                        -- t.degVIVR[w];
                        -- t.degVIVR[e];
                    }
                }
            }
        }
    }

    return;
}


void SCS::undo_degree_based_reduction(vector<ui> & del_from_VR, ContextValue &t){
    for(auto e : del_from_VR){
        t.inVR.insert(e);
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            ui w = G0_edges[i];
            if((t.inVR.find(w) != t.inVR.end()) || (t.inVI.find(w) != t.inVI.end())){
                ++ t.degVIVR[w];
                ++ t.degVIVR[e];
            }
        }
    }

    return;
}


void SCS::get_NEI_of_VI(ContextValue &t){
    for(auto e : t.VI){
        for(ui i = pstart[e]; i < pstart[e] + G0_x[e]; i++){
            if(t.inVR.find(G0_edges[i]) != t.inVR.end()){
                if(t.inNEI[G0_edges[i]] == 0){
                    t.NEI.push_back(G0_edges[i]);
                    t.inNEI[G0_edges[i]] = 1;
                }
                else{
                    ++ t.inNEI[G0_edges[i]];
                }
            }
        }
    }

    return;
}

// function to add a vertex to VI from VR
void SCS::add_to_VI_from_VR(ContextValue &t, ui vertex){
    t.VI.push_back(vertex);

    t.inVI.insert(vertex);
    t.inVR.erase(vertex);
    for(ui i = pstart[vertex]; i < pstart[vertex]+G0_x[vertex]; i++){
        if(t.inVI.find(G0_edges[i]) != t.inVI.end()){
            ++ t.degVI[vertex]; //// increase degree of ustar in VI
            ++ t.degVI[G0_edges[i]]; //// increase degree of neighbors of ustar in VI
        }
    }

    return;
}

// function to pop vertex from VI and put into VR(update inVR)
void SCS::pop_from_VI_to_VR(ContextValue &t, ui vertex){
    t.VI.pop_back();

    t.inVI.erase(vertex);
    t.inVR.insert(vertex);

    for(ui i = pstart[vertex]; i < pstart[vertex]+G0_x[vertex]; i++){
        ui w = G0_edges[i];
        if(t.inVI.find(w) != t.inVI.end()){
            -- t.degVI[vertex];
            -- t.degVI[w];
        }
    }

    return;
}


void SCS::pop_from_VI(ContextValue &t, ui vertex){
    t.VI.pop_back();
    t.inVI.erase(vertex);

    for(ui i = pstart[vertex]; i < pstart[vertex]+G0_x[vertex]; i++){
        ui w = G0_edges[i];
        if(t.inVI.find(w) != t.inVI.end()){
            -- t.degVI[vertex];
            -- t.degVI[w];
        }

        if((t.inVI.find(w) != t.inVI.end()) || (t.inVR.find(w) != t.inVR.end())){
            -- t.degVIVR[w];
            -- t.degVIVR[vertex];
        }
    }

    return;
}

// function to restore vertices in VR
void SCS::restore_in_VR(ContextValue &t, ui vertex){
    t.inVR.insert(vertex);
    for(ui j = pstart[vertex]; j < pstart[vertex]+G0_x[vertex]; j++){
        ui w = G0_edges[j];
        if((t.inVI.find(w) != t.inVI.end()) || (t.inVR.find(w) != t.inVR.end())){
            ++ t.degVIVR[w];
            ++ t.degVIVR[vertex];
        }
    }

    return;
}


// function to remove vertices being dominated from VR
void SCS::remove_domS_from_VR(ContextValue &t, vector<pair<double, ui>> domS){
    for(ui i = 0; i < domS.size(); i++){
        ui dv = domS[i].second;
        t.inVR.erase(dv);
        for(ui j = pstart[dv]; j < pstart[dv]+G0_x[dv]; j++){
            ui w = G0_edges[j];
            if((t.inVI.find(w) != t.inVI.end()) || (t.inVR.find(w) != t.inVR.end())){
                -- t.degVIVR[w];
                -- t.degVIVR[dv];
            }
        }
    }
    
    return;
}


// find all the vertices being dominated by ustar    
void SCS::find_domS_of_ustar(ui ustar, ContextValue &t, vector<pair<double, ui>> & domS){
    if(domS_Threshold == 0) return;

    for(auto e : t.NEI){
        if( t.inNEI[e] != 0 && t.degVI[e] <= t.degVI[ustar] && t.degVIVR[e] <= t.degVIVR[ustar] && e != ustar){
            ui ustar_sidx = pstart[ustar];
            ui ustar_eidx = pstart[ustar] + G0_x[ustar];
            bool be_dom = true;
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                ui w = G0_edges[i];
                if( ((t.inVR.find(w) != t.inVR.end()) || (t.inVI.find(w) != t.inVI.end())) && w != ustar){
                    while(G0_edges[ustar_sidx] < w && ustar_sidx < ustar_eidx)
                        ++ ustar_sidx;
                    if(G0_edges[ustar_sidx] == w) continue;
                    else{
                        be_dom = false;
                        break;
                    }
                }
            }
            if(be_dom){
                double its_score = 0;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){

                    if((t.inVI.find(G0_edges[i]) != t.inVI.end()) && (t.degVI[G0_edges[i]] != 0))
                      its_score += (double) 1/t.degVI[G0_edges[i]];
                }
                its_score += (double)t.degVIVR[e]/dMAX;
                domS.push_back(make_pair(its_score, e));
            }
        }
        if(domS.size() >= domS_Threshold) break;
    }
    if(domS.size() > 1) sort(domS.begin(), domS.end(), greater<pair<double, ui>>());
}

//// find the vertex to branch on
int SCS::find_ustar(ContextValue &t)
{
    int uid = -1;
    double best_score = 0;
    for(auto e : t.NEI){
        if(t.inVR.find(e) != t.inVR.end()){
            double its_score = 0;
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                if((t.inVI.find(G0_edges[i]) != t.inVI.end()) && (t.degVI[G0_edges[i]] != 0))
                  its_score += (double) 1/t.degVI[G0_edges[i]]; //// connection score
            }
            its_score += (double)t.degVIVR[e]/dMAX;
            if(its_score > best_score){
                best_score = its_score;
                uid = (int) e;
            }
        }
    }
    return uid;
}


//// finding ustar giving priority to low degree vertices, using connection score
int SCS::find_ustar_mindeg(ContextValue &t)
{
    double best_score;
    int uid = -1;
    set<ui> dict_deg;
    for(auto e : t.VI){
        dict_deg.insert(t.degVI[e]);
    }

    for(auto deg : dict_deg){
        vector<ui> vt;
        for(auto e : t.VI){
            if(t.degVI[e]==deg){
                vt.push_back(e); // get all the vertices in VI with degree deg
            }
        }

        best_score = 0;
        for(auto e : vt){
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                ui w = G0_edges[i]; //// w is a vertex in VR
                if(t.inVR.find(w) != t.inVR.end()){
                    double its_score = 0;
                    for(ui j = pstart[w]; j < pstart[w]+G0_x[w]; j++){
                        if((t.inVI.find(G0_edges[j]) != t.inVI.end()) && (t.degVI[G0_edges[j]] != 0)){
                            its_score += (double) 1/t.degVI[G0_edges[j]]; //// w's connection score in VI
                        }
                    }
                    its_score += (double)t.degVIVR[w]/dMAX; //// add connection score in entire VIVR as well
                    if(its_score > best_score){  // get the best score among the vertices with same degree only
                        best_score = its_score; // prioritize low degree vertices even more
                        uid = (int) w;
                    }
                }
            }
        }
        if(uid != -1)
            break;
    }
    return uid;
}

int SCS::find_ustar_2phase(ContextValue &t, ui N2)
{
    if(t.VI.size() > (N2*2)/5 )
        return find_ustar_mindeg(t);
    else
        return find_ustar(t);
}

//// score is just a connection to vertices in VI
int SCS::find_ustar_link(ContextValue &t)
{
    int uid = -1;
    double best_score = 0;
    for(auto e : t.NEI){
        if(t.inVR.find(e) != t.inVR.end()){
            double its_score = 0;
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                if(t.inVI.find(G0_edges[i]) != t.inVI.end())
                    its_score ++;
            }
            if(its_score >= best_score){
                best_score = its_score;
                uid = e;
            }
        }
    }
    return uid;
}


//// find the vertex to branch on randomly
int SCS::find_ustar_random(ContextValue &t)
{
    int uid = -1;
    for(auto e : t.NEI){
        if(t.inVR.find(e) != t.inVR.end()){
            return e;
        }
    }
    return uid;
}

ui SCS::get_ub1(ContextValue &t, ui N2)
{
    ui min_deg = INF;
    ui r = N2 - (ui)t.VI.size(); //// vertices that can still be added to VI, h - |C|
    if(r<0){
        exit(1);
    }
    for(auto e : t.VI){
        ui its_deg_ub = 0;
        ui cands = 0;
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            if(t.inVR.find(G0_edges[i]) != t.inVR.end()) ++ cands;
        }
        its_deg_ub = t.degVI[e] + min(r, cands); //// degree based upper bound, 𝑑𝐶(𝑢) +min{ℎ − |𝐶|, 𝑑𝑅∪{𝑢 } (𝑢)},
        if(its_deg_ub < min_deg) min_deg = its_deg_ub; //// taking the minimum over all the upper bounds
    }
    return min_deg;
}

ui SCS::get_ub2(ContextValue &t, ui N2)
{
    vector<ui> nei;
    bool * innei = new bool[num_vertices];
    memset(innei, 0, sizeof(bool)*num_vertices);
    for(auto e : t.VI){
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            if(t.inVR.find(G0_edges[i]) != t.inVR.end() && !innei[G0_edges[i]]){
                nei.push_back(G0_edges[i]); //// add the neighbors of VI's vertices in VR
                innei[G0_edges[i]] = 1;
            }
        }
    }
    delete [] innei;
    
    ui r = N2 - (ui)t.VI.size(); //// h- |C| vertices
    vector<ui> cov_power;
    for(auto e : nei){
        ui its_power = 0;
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            if(t.inVI.find(G0_edges[i]) != t.inVI.end())
                ++ its_power; //// get the degree of VR's vertices in VI
        }
        cov_power.push_back(its_power);
    }
    ui r1 = (ui)cov_power.size();
    ui r2 = min(r, r1); //// r2 vertices will be in R' now
    sort(cov_power.begin(), cov_power.end(), greater<ui>());//decreasing order, to get vertices with greater degree
    vector<ui> interm_deg;
    for(auto e : t.VI){
        interm_deg.push_back(t.degVI[e]); //// add degree of each vertices in VI to interm_deg
    }
    sort(interm_deg.begin(), interm_deg.end(), less<ui>());//increasing order, to get vertices with less degree
    
    for(ui i = 0; i < r2; i++){
        ui budget = cov_power[i];
        for(ui j = 0; j < budget; j++){
            ++ interm_deg[j]; //// increase degree of low degree vertices in VI
        }
        sort(interm_deg.begin(), interm_deg.end(), less<ui>());//increasing order, we need to sort in every iteration
    }
    
    return interm_deg[0]; //// the first element is the one with min degree i.e, Unr
}

ui SCS::get_ub3(ContextValue &t, ui kl, ui N2)
{
    vector<ui> nei;
    bool * innei = new bool[num_vertices];
    memset(innei, 0, sizeof(bool)*num_vertices);
    
    for(auto e : t.VI){
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            if(t.inVR.find(G0_edges[i]) != t.inVR.end() && !innei[G0_edges[i]]){
                nei.push_back(G0_edges[i]);
                innei[G0_edges[i]] = 1;
            }
        }
    }
    delete [] innei;
    
    set<ui> ditc_deg;
    for(auto e : t.VI){
        ditc_deg.insert(t.degVI[e]); //// add degree of VI's vertices
    }
    vector<ui> deg_lev;
    for(auto e : ditc_deg)
        deg_lev.push_back(e); //// add the degree levels, 1,2,3,...
    
    ui ditc_deg_num = (ui)deg_lev.size();
    
    ui r = N2 - (ui)t.VI.size();
    
    ui min_deg = INF;
    ui rcd_i = 0;
    
    for(ui i = 0; i < ditc_deg_num; i++){ //// iterate from dmin to dmax
        ui t_deg = deg_lev[i];
        vector<ui> interm_deg;
        for(auto e : t.VI){
            if(t.degVI[e] <= t_deg){
                interm_deg.push_back(t.degVI[e]); //// add degree of each vertices in VI with degree <= deg_lev to interm_deg, C<=1 ...
            }
        }
        sort(interm_deg.begin(), interm_deg.end(), less<ui>()); //increasing order
        
        vector<ui> cov_power;
        for(auto e : nei){
            ui its_power = 0;
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                if(t.inVI.find(G0_edges[i]) != t.inVI.end() && t.degVI[G0_edges[i]] <= t_deg){
                    ++ its_power; //// get the degree of VR's vertices whose degrees are <= deg_lev in VI
                }
            }
            if(its_power > 0){
                cov_power.push_back(its_power);
            }
        }
        sort(cov_power.begin(), cov_power.end(), greater<ui>()); //decreasing order
        
        ui r1 = (ui)cov_power.size();
        ui r2 = min(r, r1);
        
        for(ui i = 0; i < r2; i++){
            ui budget = cov_power[i];
            for(ui j = 0; j < budget; j++){
                ++ interm_deg[j];
            }
            sort(interm_deg.begin(), interm_deg.end(), less<ui>());//increasing order
        }
        if(min_deg > interm_deg[0]){ //// the first element is the one with min degree i.e, Udc
            min_deg = interm_deg[0]; //// take minimum over each iteration
            rcd_i = i;
        }
        if(min_deg <= kl)
            return min_deg;
        
        if(EXE_ub3_optimization){
            if(i < ditc_deg_num-1 && min_deg <= deg_lev[i+1]){
                return min_deg;
            }
        }
    }
    return min_deg;
}


// bool SCS::compute_ub(bool * inVI, bool * inVR, vector<ui> VI, ui * degVI, ui kl, ui N2){
bool SCS::compute_ub(ContextValue &t, ui kl, ui N2){

    ui ub1 = INF;

    if(EXE_ub1) ub1 = get_ub1(t, N2);
    if(ub1 <= kl) return false;
    
    ui ub2 = INF;
    if(EXE_ub2) ub2 = get_ub2(t, N2);
    if(ub2 <= kl) return false;  

    ui ub3 = INF;
    if(EXE_ub3) ub3 = get_ub3(t, kl, N2);
    if(ub3 <= kl) return false;
    
    return true;
}

// get the vertices to be restored to VR
void SCS::get_rVR(ContextValue &t, ui kl, bool & del_v_in_VI, vector<ui> & rVR, queue<ui> & Q, bool * inQ){
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        if(t.inVI.find(v) != t.inVI.end()){
            del_v_in_VI = true; //// remove vertex from VI with low degree <= kl
            return;
        }
        else if(t.inVR.find(v) != t.inVR.end()){
            t.inVR.erase(v);
            rVR.push_back(v); //// add to rVR for restoration later on
        }
        for(ui i = pstart[v]; i < pstart[v]+G0_x[v]; i++){
            ui w = G0_edges[i]; // neighbors of v in VR
            if((t.inVI.find(w) != t.inVI.end()) || (t.inVR.find(w) != t.inVR.end())){
                -- t.degVIVR[w]; //// decrease the degree of neighbor in VIVR
                -- t.degVIVR[v]; //// decrease the degree of vertex in VIVR

                if(t.degVIVR[w] <= kl && !inQ[w]){
                    Q.push(w);
                    inQ[w] = 1;
                }
            }
        }
    }
    return;
}

void SCS::core_maintenance(ContextValue &t, ui kl, bool & del_v_in_VI, vector<ui> & rVR, ui ustar){
    queue<ui> Q;
    bool * inQ = new bool[num_vertices];
    memset(inQ, 0, sizeof(bool)*num_vertices);
    
    for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
        ui w = G0_edges[i];
        if(((t.inVI.find(w) != t.inVI.end()) || (t.inVR.find(w) != t.inVR.end())) && (t.degVIVR[w] <= kl)){ // unpromising vertex
            Q.push(w);
            inQ[w] = 1;
        }
    }

    get_rVR(t, kl, del_v_in_VI, rVR, Q, inQ);
    delete [] inQ;
}

void SCS::core_maintenance(ContextValue &t, ui kl, bool & del_v_in_VI, vector<ui> & rVR, vector<ui> & del_vec){
    queue<ui> Q;
    bool * inQ = new bool[num_vertices];
    memset(inQ, 0, sizeof(bool)*num_vertices);
    
    for(auto e : del_vec){
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            ui w = G0_edges[i];
            if( ((t.inVI.find(w) != t.inVI.end()) || (t.inVR.find(w) != t.inVR.end())) && t.degVIVR[w] <= kl && !inQ[w]){
                Q.push(w);
                inQ[w] = 1;
            }
        }
    }

    get_rVR(t, kl, del_v_in_VI, rVR, Q, inQ);
    delete [] inQ;
}

void SCS::core_maintenance(ContextValue &t, ui kl, bool & del_v_in_VI, vector<ui> & rVR, vector<pair<double, ui>> & del_vec){
    queue<ui> Q;
    bool * inQ = new bool[num_vertices];
    memset(inQ, 0, sizeof(bool)*num_vertices);
    
    for(auto ee : del_vec){
        auto e = ee.second;
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            ui w = G0_edges[i];
            if( ((t.inVI.find(w) != t.inVI.end()) || (t.inVR.find(w) != t.inVR.end())) && t.degVIVR[w] <= kl && !inQ[w]){
                Q.push(w);
                inQ[w] = 1;
            }
        }
    }

    get_rVR(t, kl, del_v_in_VI, rVR, Q, inQ);
    delete [] inQ;
}

// This function only shrinks G0_edges(adjacency list) for each vertices, not the entire graph
void SCS::ShrinkGraph(ContextValue &t, SCS & split_g){
    vector<ui> G0_edges_temp;
    ui shrinked_total_edges = 0;

    split_g.pstart = new ui [num_vertices+1];
    split_g.G0_x = new ui [num_vertices]; // degree of a vertex

    split_g.pstart[0] = 0;
    split_g.num_VIVR = 0;

    // check every vertices and their neighbor to see if they are either in VI or in VR
    for(ui v=0; v<num_vertices; v++){
        if((t.inVI.find(v) != t.inVI.end()) || (t.inVR.find(v) != t.inVR.end())){ // vertex should be in VI or in VR
            
            ++split_g.num_VIVR;
            ui new_deg = 0;
            for(ui i = pstart[v]; i < pstart[v] + G0_x[v]; i++){
                ui w = G0_edges[i];
                if((t.inVI.find(w) != t.inVI.end()) || (t.inVR.find(w) != t.inVR.end())){
                    // count new G0_edges after shrinking
                    shrinked_total_edges ++;

                    // update G0_edges_temp
                    G0_edges_temp.push_back(w);
                    ++new_deg;
                }
            }

            // update degree of current vertex and pstart of next vertex
            split_g.G0_x[v] = new_deg;
            split_g.pstart[v+1] = split_g.pstart[v] + new_deg; 
        }
        else{
            // update degree of current vertex and pstart of next vertex
            split_g.G0_x[v] = 0;
            split_g.pstart[v+1] = split_g.pstart[v];
        }
    } 

    // copy the vector into original array only
    split_g.G0_edges = new ui [shrinked_total_edges];
    memcpy(split_g.G0_edges, G0_edges_temp.data(), sizeof(ui)*(shrinked_total_edges));


    // update shrinked number of G0_edges
    split_g.num_edges = shrinked_total_edges/2;
}


// // TODO: 1. shrink pstart and G0_edges, copy to a new array 2. use vector
// // This function shrinks both pstart and G0_edges using vector
// void SCS::ShrinkGraph(ContextValue &t, SCS & split_g){

//     vector<ui> G0_edges_temp;
//     vector<ui> G0_x_temp;
//     vector<ui> pstart_temp;

//     // recoding map
//     // TODO: change to hash_map, one direction should be hash_map, another should be array
//     map<ui, ui> &index2id_temp = split_g.index2id; // get original id based on new index
//     map<ui, ui> &id2index_temp = split_g.id2index; // get new index based on original id

//     ui new_index = 0;
//     ui edges_counter = 0;
//     pstart_temp.push_back(0);

//     // check every vertices and their neighbor to see if they are either in VI or in VR
//     for(ui v=0; v<num_vertices; v++){
//         if(t.inVI[v] || t.inVR[v]){ // vertex should be in VI or in VR
//             // update recoding map

//             index2id_temp[new_index] = index2id[v]; // TODO: check
//             id2index_temp[index2id[v]] = new_index; // TODO: check
            
//             ui n_deg = 0; // degree
//             G0_x_temp.push_back(n_deg);
//             for(ui i = pstart[v]; i < pstart[v] + G0_x[v]; i++){
//                 ui w = G0_edges[i]; // nei --> here it is index
//                 // ui w = index2id[G0_edges[i]]; // here we get id
//                 if(t.inVI[w] || t.inVR[w]){
//                     // count new G0_edges after shrinking
//                     edges_counter ++;

//                     G0_edges_temp.push_back(w);
//                     G0_x_temp.insert(G0_x_temp.begin()+new_index, ++n_deg);
//                 }
//             }

//             // get new pstart and update counter
//             pstart_temp.push_back(pstart_temp[new_index] + G0_x_temp[new_index]);
//             new_index++;
//         }
//     }

//     // slice the pstart_new, G0_edges_new and G0_x_new array
//     split_g.pstart = new ui[new_index + 1];
//     split_g.G0_edges = new ui[2*edges_counter];
//     split_g.G0_x = new ui[new_index];

//     // copy the part of array only
//     memcpy(split_g.pstart, pstart_temp.data(), sizeof(ui)*(new_index+1));
//     // memcpy(split_g.G0_edges, G0_edges_temp, sizeof(ui)*(2*edges_counter));
//     memcpy(split_g.G0_x, G0_x_temp.data(), sizeof(ui)*(new_index));

//     // get new index of id in G0_edges
//     for(int i=0; i<2*edges_counter; i++){
//         split_g.G0_edges[i] = id2index[G0_edges_temp[i]];
//     }

//     // update shrinked number of vertices and G0_edges
//     split_g.num_vertices = new_index;
//     split_g.num_edges = edges_counter;
// }


void SCS::DestroySplitGraph(){
    if(G0_edges != NULL) delete [] G0_edges;
    if(G0_x != NULL) delete [] G0_x;
    if(pstart != NULL) delete [] pstart; 
    if(q_dist != NULL) delete [] q_dist;
}

