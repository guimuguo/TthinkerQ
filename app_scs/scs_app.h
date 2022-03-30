
#ifndef SCS_APP_H_
#define SCS_APP_H_

#define SPAWNED_TASK 1
#define DECOMPOSED_TASK 2

#include "../system/workerOL.h"
#include "../system/task.h"
#include <fstream>
#include <assert.h>
#include "graph.h"
#include <climits>
#include "../system/TaskProgMap.h"
#include <sstream>
#include "../system/rwlock.h"
#include "../system/conmap.h"


float TIME_THRESHOLD; // if running time >= TIME_THRESHOLD, split into subtasks

SCS global_g;

ui kl_temp;
ui ubD_temp;


//-------------------------------------------
// check if input file is empty
bool empty(ifstream &pFile)
{
    return pFile.peek() == ifstream::traits_type::eof();
}

typedef Task<ContextValue> SCSTask;


bool setup_task(SCSTask *task, SCSQuery &query, FILE * gfpout, int query_id){
	query.QID = outer_id2index[query.QID]; // get the index of QID
	global_g.CSSC_heu(query); // find initial feasible community H and min-degree kl 

    ui ku = miv(global_g.core[query.QID], query.N2-1); //UB ////upper bound, taking min val; min(cn(q), h-1), finding k-hat

    if(query.kl == ku){  // check if k-tilde == k-hat
		ftime(&query.end_t);
		double totaltime = query.end_t.time-query.start_t.time+(double)(query.end_t.millitm-query.start_t.millitm)/1000;

        cout<<"Heuristics found the Optimal Solution!"<<endl;
        cout<<"Optimal Minimum Degree: "<<query.kl<<endl;
		cout << "Optimal Community: ";
		for(auto h: query.H){
			cout << h << " ";
		}
		cout << "Size of Optimal Community: " << query.H.size() << endl;

		fprintf(gfpout, "[INFO] Query ");
		fprintf(gfpout, "%d", query_id);
		fprintf(gfpout, "\n");
		fprintf(gfpout, "Optimal Minimum Degree: ");
		fprintf(gfpout, "%d ", query.kl);
		fprintf(gfpout, "\n");
		fprintf(gfpout, "Optimal Community: ");
		for(auto h: query.H){
			fprintf(gfpout, "%d ", h);
			fprintf(gfpout, " ");
		}
		fprintf(gfpout, "\n");
		fprintf(gfpout, "Size of Optimal Community: ");
		fprintf(gfpout, "%d", query.H.size());
		fprintf(gfpout, "\n");
		fprintf(gfpout, "Total vertices explored: ");
		fprintf(gfpout, "%d", accumulate(query.counters.begin(), query.counters.end(), 0));
		fprintf(gfpout, "\n");
		fprintf(gfpout, "Total time: ");
		fprintf(gfpout, "%f", totaltime);
		fprintf(gfpout, "\n");

		fflush(gfpout);

        return false;		
    }

    // get the upper bound of diameter based on number of vertices and min-degree
    query.ubD = 0;
    if(query.kl<=1) query.ubD = query.N2-1;
    else query.ubD = global_g.find_ubD(query.kl, query.ubD, query.N2);

    global_g.cal_query_dist(query.QID); //// calculate hop distance from query vertex to each vertices // query dependent
    global_g.reduction_g(query, task->context.split_g); //// remove unpromising vertices, // query dependent

	// check if number of vertices after reduction is less than N1
	if(task->context.split_g.num_vertices < query.N1) return false;

	task->context.split_g.gUB = global_g.gUB;

	task->context.split_g.q_dist = new ui [task->context.split_g.num_vertices];
	memset(task->context.split_g.q_dist, 0, sizeof(ui)*task->context.split_g.num_vertices);

	// for(int i=0; i<global_g.num_vertices; i++){
	// 	if(find(task->context.split_g.G0.begin(), task->context.split_g.G0.end(), i) != task->context.split_g.G0.end()){
	// 		task->context.split_g.q_dist[task->context.split_g.id2index[i]] = global_g.q_dist[i];
	// 	}
    // }
    
    task->context.degVI = new ui[task->context.split_g.num_vertices];
    memset(task->context.degVI, 0, sizeof(ui)*task->context.split_g.num_vertices);
    task->context.degVIVR = new ui[task->context.split_g.num_vertices];
    memset(task->context.degVIVR, 0, sizeof(ui)*task->context.split_g.num_vertices);
    task->context.inNEI = new ui[task->context.split_g.num_vertices];
    memset(task->context.inNEI, 0, sizeof(ui)*task->context.split_g.num_vertices);

    for(auto e : task->context.split_g.G0){
		ui recoded_id = task->context.split_g.id2index[e];
		task->context.inVR.insert(recoded_id);
        task->context.degVIVR[recoded_id] = task->context.split_g.G0_x[recoded_id]; //// degree of vertices in VIVR after removing unpromising vertices
    }

	task->context.VI.clear();

    task->context.VI.push_back(query.QID);
	task->context.inVI.insert(query.QID);
	task->context.inVR.erase(query.QID);
    
    return true;
}


class SCSComper : public Comper<SCSTask, SCSQuery>
{
public:
	ui counter = 0;
	void BB_dom_ustar(ContextT &context, SCSQuery &query, SCS &gograph)
	{
		struct timeb cur_time;
		double drun_time;

		ftime(&query.end_t);

		query.hlock.rdlock();
		double DurTime = query.end_t.time-query.start_t.time+(double)(query.end_t.millitm-query.start_t.millitm)/1000;
		query.hlock.unlock();

		if(DurTime > MaxTime){
			cout << "Timeout occured at " << DurTime << " sec on " <<  "QueryId " << get_queryID() << endl;
			return;
		}

		if(context.VI.size() > query.N2){
			return;
		}

		// update the feasible community and upper bound of Diameter
		if(context.VI.size() >= query.N1 && context.VI.size() <= query.N2){
			ui cur_min_deg = INF;
			for(auto e : context.VI){
				if(context.degVI[e] < cur_min_deg)
					cur_min_deg = context.degVI[e];
			}

			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			if(cur_min_deg > kl_temp){
				query.hlock.rdlock();
				ubD_temp = query.ubD;
				query.hlock.unlock();

				ubD_temp = gograph.find_ubD(kl_temp, ubD_temp, query.N2);

				query.hlock.wrlock();
				if(cur_min_deg > query.kl){
					query.kl = cur_min_deg;

					// query.H = context.VI;

					// get reverse mapping of recoded id after pstart is shrinked
					query.H.clear();
					for(auto index: context.VI){
						ui intermediate_index = gograph.G0[index]; // get the intermediate index
						ui original_id = G[intermediate_index]; // get the original id
						query.H.push_back(original_id);
					}

					query.ubD = ubD_temp;
				}
				query.hlock.unlock();

			}
		}

		if(context.VI.size() == query.N2){
			return;
		}

		unordered_set<ui> new2VI;

		if(EXE_new2VI){ // Inclusion based reduction
			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			gograph.inclusion_based_reduction(context, new2VI, kl_temp);

			if(context.VI.size() > query.N2){
				// undo inclusion based reduction
				for(auto e : new2VI){
					gograph.pop_from_VI_to_VR(context, e);
				}
				return;
			}
		}

		//// update the feasible community and upper bound of Diameter
		if(context.VI.size() >= query.N1 && context.VI.size() <= query.N2){
			ui cur_min_deg = INF;
			for(auto e : context.VI){
				if(context.degVI[e] < cur_min_deg)
					cur_min_deg = context.degVI[e];
			}

			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			if(cur_min_deg > kl_temp){
				query.hlock.rdlock();
				ubD_temp = query.ubD;
				query.hlock.unlock();

				ubD_temp = gograph.find_ubD(kl_temp, ubD_temp, query.N2);

				query.hlock.wrlock();
				if(cur_min_deg > query.kl){
					query.kl = cur_min_deg;
					// query.H = context.VI;

					// get reverse mapping of recoded id after pstart is shrinked
					query.H.clear();
					for(auto index: context.VI){
						ui intermediate_index = gograph.G0[index]; // get the intermediate index
						ui original_id = G[intermediate_index]; // get the original id
						query.H.push_back(original_id);
					}

					query.ubD = ubD_temp;
				}
				query.hlock.unlock();

			}
		}

		if(context.VI.size() == query.N2){
			if(EXE_new2VI){
				// undo inclusion based reduction
				for(auto e : new2VI){
					gograph.pop_from_VI_to_VR(context, e);
				}
			}
			return;
		}
		
		context.NEI.clear(); //// NEI is used to find vertices being dominated later on
		memset(context.inNEI, 0, sizeof(ui)*gograph.num_vertices); // for array version
		
		// get all the neighbors of vertices in VI
		gograph.get_NEI_of_VI(context);

		vector<ui> del_from_VR;
		if(EXE_del_from_VR){ // degree-based reduction
			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			gograph.degree_based_reduction(context, del_from_VR, kl_temp, query.N2);
		}

		int ustar = -1;
		ustar = gograph.find_ustar_mindeg(context);
		
		if(ustar < 0){ //// return if no vertex to branch on is found and undo changes before returning(if any)
			if(EXE_del_from_VR){
				// undo degree based reduction
				for(auto e : del_from_VR){
					gograph.restore_in_VR(context, e);
				}
			}

			if(EXE_new2VI){
				// undo inclusion based reduction
				for(auto e : new2VI){
					gograph.pop_from_VI_to_VR(context, e);
				}
			}
			return;
		}

		counter++;

		vector<pair<double, ui>> domS;

		gograph.find_domS_of_ustar(ustar, context, domS); //// find all the vertices being dominated by ustar, returns a pair of <score, vertex>

		if(!domS.empty()){
			for(ui round = 0; round < domS.size(); round++){
				
				ui domv = domS[round].second;

				// add ustar to VI
				gograph.add_to_VI_from_VR(context, ustar);

				// add domv to VI
				gograph.add_to_VI_from_VR(context, domv);

				counter++;
				
				vector<ui> rmv_each_round;
				for(ui i = 0; i < round; i++){
					ui dv = domS[i].second;
					context.inVR.erase(dv);
					rmv_each_round.push_back(dv);
					for(ui j = gograph.pstart[dv]; j < gograph.pstart[dv]+gograph.G0_x[dv]; j++){
						if(context.inVI.find(gograph.G0_edges[j]) != context.inVI.end() || context.inVR.find(gograph.G0_edges[j]) != context.inVR.end()){
							-- context.degVIVR[gograph.G0_edges[j]];
							-- context.degVIVR[dv];
						}
					}
				}
				
				vector<ui> rVR;
				bool del_v_in_VI = false;


				// if the neighbors of vertices with higher connection score has degree <= kl, then remove them
				if(EXE_core_maintenance){
					query.hlock.rdlock();
					kl_temp = query.kl;
					query.hlock.unlock();

					gograph.core_maintenance(context, kl_temp, del_v_in_VI, rVR, rmv_each_round);
					if(del_v_in_VI){ //// if there are unpromising vertices in VI, then prune this branch
						// restore vertices in rVR
						for(ui i = 0; i < rVR.size(); i++){
							ui v = rVR[i];
							gograph.restore_in_VR(context, v);
						}

						// restore vertices being dominated
						for(ui i = 0; i < round; i++){
							ui dv = domS[i].second;
							gograph.restore_in_VR(context, dv);
						}
						
						// remove domv from VI
						gograph.pop_from_VI_to_VR(context, domv);
						
						// remove ustar from VI
						gograph.pop_from_VI_to_VR(context, ustar);

						continue;
					}
				}
				

				query.hlock.rdlock();
				kl_temp = query.kl;
				query.hlock.unlock();

				gograph.gUB = gograph.compute_ub(context, kl_temp, query.N2);
				
				if(gograph.gUB){
					ftime(&cur_time);
					drun_time = cur_time.time-gograph.gtime_start.time+(double)(cur_time.millitm-gograph.gtime_start.millitm)/1000;
					if(drun_time < TIME_THRESHOLD){
						BB_dom_ustar(context, query, gograph); //// C U {u*, v1}
					}
					else{
						SCSTask * t = new SCSTask();

						gograph.ShrinkGraph(context, t->context.split_g);

						// -----------------------
						// for map version
						t->context.VI = context.VI;
						t->context.NEI = context.NEI;
						t->context.inVI = context.inVI;
						t->context.inVR = context.inVR;

						t->context.split_g.q_dist = new ui [gograph.num_vertices];
						memcpy(t->context.split_g.q_dist, gograph.q_dist, sizeof(ui)*(gograph.num_vertices));

						t->context.degVI = new ui[gograph.num_vertices];
						memcpy(t->context.degVI, context.degVI, sizeof(ui)*(gograph.num_vertices));

						t->context.degVIVR = new ui[gograph.num_vertices];
						memcpy(t->context.degVIVR, context.degVIVR, sizeof(ui)*(gograph.num_vertices));

						t->context.inNEI = new ui[gograph.num_vertices];
						memcpy(t->context.inNEI, context.inNEI, sizeof(ui)*(gograph.num_vertices));

						t->context.split_g.num_vertices = gograph.num_vertices;
						t->context.split_g.gUB = gograph.gUB;

						t->context.split_g.G0 = gograph.G0;

						add_task(t);
					}
				}

				if(EXE_core_maintenance){
					// restore vertices in rVR
					for(ui i = 0; i < rVR.size(); i++){
						ui v = rVR[i];
						gograph.restore_in_VR(context, v);
					}
				}

				// restore vertices being dominated
				for(ui i = 0; i < round; i++){
					ui dv = domS[i].second;
					gograph.restore_in_VR(context, dv);
				}

				// remove domv from VI
				gograph.pop_from_VI_to_VR(context, domv);

				// remove ustar from VI
				gograph.pop_from_VI_to_VR(context, ustar);

			}
			
			// remove vertices being dominated from VR
			gograph.remove_domS_from_VR(context, domS);
			
			vector<ui> rVR;
			bool del_v_in_VI = false;

			if(EXE_core_maintenance){
				query.hlock.rdlock();
				kl_temp = query.kl;
				query.hlock.unlock();

				gograph.core_maintenance(context, kl_temp, del_v_in_VI, rVR, domS);
				if(del_v_in_VI){
					// restore vertices in rVR
					for(ui i = 0; i < rVR.size(); i++){
						ui v = rVR[i];
						gograph.restore_in_VR(context, v);
					}

					// restore vertices being dominated
					for(ui i = 0; i < domS.size(); i++){
						ui dv = domS[i].second;
						gograph.restore_in_VR(context, dv);
					}
					
					if(EXE_del_from_VR){
						// undo degree based reduction
						for(auto e : del_from_VR){
							gograph.restore_in_VR(context, e);
						}
					}

					if(EXE_new2VI){
						// undo inclusion based reduction
						for(auto e : new2VI){
							gograph.pop_from_VI_to_VR(context, e);
						}
					}

					return;
				}
			}

			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			gograph.gUB = gograph.compute_ub(context, kl_temp, query.N2);
			if(gograph.gUB){
				ftime(&cur_time);
				drun_time = cur_time.time-gograph.gtime_start.time+(double)(cur_time.millitm-gograph.gtime_start.millitm)/1000;
				if(drun_time < TIME_THRESHOLD){
					BB_dom_ustar(context, query, gograph); //// C U {u*}
				}
				else{
					SCSTask * t = new SCSTask();

					gograph.ShrinkGraph(context, t->context.split_g);

					// -----------------------
					// for map version
					t->context.VI = context.VI;
					t->context.NEI = context.NEI;
					t->context.inVI = context.inVI;
					t->context.inVR = context.inVR;

					t->context.degVI = new ui[gograph.num_vertices];
					memcpy(t->context.degVI, context.degVI, sizeof(ui)*(gograph.num_vertices));

					t->context.degVIVR = new ui[gograph.num_vertices];
					memcpy(t->context.degVIVR, context.degVIVR, sizeof(ui)*(gograph.num_vertices));

					t->context.inNEI = new ui[gograph.num_vertices];
					memcpy(t->context.inNEI, context.inNEI, sizeof(ui)*(gograph.num_vertices));

					t->context.split_g.q_dist = new ui [gograph.num_vertices];
					memcpy(t->context.split_g.q_dist, gograph.q_dist, sizeof(ui)*(gograph.num_vertices));

					t->context.split_g.num_vertices = gograph.num_vertices;
					t->context.split_g.gUB = gograph.gUB;

					t->context.split_g.G0 = gograph.G0;

					add_task(t);
				}
			}
			
			if(EXE_core_maintenance){
				// restore vertices in rVR
				for(ui i = 0; i < rVR.size(); i++){
					ui v = rVR[i];
					gograph.restore_in_VR(context, v);
				}
			}

			// restore vertices being dominated
			for(ui i = 0; i < domS.size(); i++){
				ui dv = domS[i].second;
				gograph.restore_in_VR(context, dv);
			}

			if(EXE_del_from_VR){
				// undo degree based reduction
				for(auto e : del_from_VR){
					gograph.restore_in_VR(context, e);
				}
			}
			if(EXE_new2VI){
				// undo inclusion based reduction
				for(auto e : new2VI){
					gograph.pop_from_VI_to_VR(context, e);
				}
			}
		}
		
		else{

			// add ustar to VI
			gograph.add_to_VI_from_VR(context, ustar);

			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			gograph.gUB = gograph.compute_ub(context, kl_temp, query.N2);
			if(gograph.gUB){
				ftime(&cur_time);
				drun_time = cur_time.time-gograph.gtime_start.time+(double)(cur_time.millitm-gograph.gtime_start.millitm)/1000;
				if(drun_time < TIME_THRESHOLD){
					BB_dom_ustar(context, query, gograph);
				}
				else{
					SCSTask * t = new SCSTask();

					gograph.ShrinkGraph(context, t->context.split_g);

					// for map version
					t->context.VI = context.VI;
					t->context.NEI = context.NEI;
					t->context.inVI = context.inVI;
					t->context.inVR = context.inVR;

					t->context.split_g.q_dist = new ui [gograph.num_vertices];
					memcpy(t->context.split_g.q_dist, gograph.q_dist, sizeof(ui)*(gograph.num_vertices));

					t->context.degVI = new ui[gograph.num_vertices];
					memcpy(t->context.degVI, context.degVI, sizeof(ui)*(gograph.num_vertices));

					t->context.degVIVR = new ui[gograph.num_vertices];
					memcpy(t->context.degVIVR, context.degVIVR, sizeof(ui)*(gograph.num_vertices));

					t->context.inNEI = new ui[gograph.num_vertices];
					memcpy(t->context.inNEI, context.inNEI, sizeof(ui)*(gograph.num_vertices));

					t->context.split_g.num_vertices = gograph.num_vertices;
					t->context.split_g.gUB = gograph.gUB;

					t->context.split_g.G0 = gograph.G0;

					add_task(t);
				}
			}
			
			// remove ustar from VI
			gograph.pop_from_VI(context, ustar);
			
			vector<ui> rVR;
			bool del_v_in_VI = false;
			if(EXE_core_maintenance){
				query.hlock.rdlock();
				kl_temp = query.kl;
				query.hlock.unlock();

				gograph.core_maintenance(context, kl_temp, del_v_in_VI, rVR, ustar);
				if(del_v_in_VI){
					// restore vertices in rVR
					for(ui i = 0; i < rVR.size(); i++){
						ui v = rVR[i];
						gograph.restore_in_VR(context, v);
					}
					
					// restore ustar in VR
					gograph.restore_in_VR(context, ustar);

					// undo degree based reduction
					for(auto e : del_from_VR){
						gograph.restore_in_VR(context, e);
					}
					
					// undo inclusion based reduction
					for(auto e : new2VI){
						gograph.pop_from_VI_to_VR(context, e);
					}

					return;
				}
			}

			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			gograph.gUB = gograph.compute_ub(context, kl_temp, query.N2);
			if(gograph.gUB){
				ftime(&cur_time);
				drun_time = cur_time.time-gograph.gtime_start.time+(double)(cur_time.millitm-gograph.gtime_start.millitm)/1000;
				if(drun_time < TIME_THRESHOLD){
					BB_dom_ustar(context, query, gograph);
				}
				else{
					//to split
					SCSTask * t = new SCSTask();

					gograph.ShrinkGraph(context,  t->context.split_g);

					// for map version
					t->context.VI = context.VI;
					t->context.NEI = context.NEI;
					t->context.inVI = context.inVI;
					t->context.inVR = context.inVR;

					t->context.split_g.q_dist = new ui [gograph.num_vertices];
					memcpy(t->context.split_g.q_dist, gograph.q_dist, sizeof(ui)*(gograph.num_vertices));

					t->context.degVI = new ui[gograph.num_vertices];
					memcpy(t->context.degVI, context.degVI, sizeof(ui)*(gograph.num_vertices));

					t->context.degVIVR = new ui[gograph.num_vertices];
					memcpy(t->context.degVIVR, context.degVIVR, sizeof(ui)*(gograph.num_vertices));

					t->context.inNEI = new ui[gograph.num_vertices];
					memcpy(t->context.inNEI, context.inNEI, sizeof(ui)*(gograph.num_vertices));

					t->context.split_g.num_vertices = gograph.num_vertices;
					t->context.split_g.gUB = gograph.gUB;

					t->context.split_g.G0 = gograph.G0;

					add_task(t);
				}
			}
			
			if(EXE_core_maintenance){
				// restore vertices in rVR
				for(ui i = 0; i < rVR.size(); i++){
					ui v = rVR[i];
					gograph.restore_in_VR(context, v);
				}
			}

			// restore ustar
			gograph.restore_in_VR(context, ustar);
			
			if(EXE_del_from_VR){
				// undo degree based reduction
				for(auto e : del_from_VR){
					gograph.restore_in_VR(context, e);
				}
			}
			if(EXE_new2VI){
				// undo inclusion based reduction
				for(auto e : new2VI){
					gograph.pop_from_VI_to_VR(context, e);
				}
			}
		}

		return;
	}


	void BB(ContextT &context, SCSQuery &query, SCS &gograph)
	{
		struct timeb cur_time;
		double drun_time;

		ftime(&query.end_t);
		double DurTime = query.end_t.time-query.start_t.time+(double)(query.end_t.millitm-query.start_t.millitm)/1000;
    
		if(DurTime > MaxTime){
			return;
		}

		if(context.VI.size() > query.N2){
			return;
		}

		// update the feasible community and upper bound of Diameter
		if(context.VI.size() >= query.N1 && context.VI.size() <= query.N2){
			ui cur_min_deg = INF;
			for(auto e : context.VI){
				if(context.degVI[e] < cur_min_deg)
					cur_min_deg = context.degVI[e];
			}

			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			if(cur_min_deg > kl_temp){
				query.hlock.rdlock();
				ubD_temp = query.ubD;
				query.hlock.unlock();

				ubD_temp = gograph.find_ubD(kl_temp, ubD_temp, query.N2);

				query.hlock.wrlock();
				if(cur_min_deg > query.kl){
					query.kl = cur_min_deg;

					// query.H = context.VI;

					// get reverse mapping of recoded id after pstart is shrinked
					query.H.clear();
					for(auto index: context.VI){
						ui intermediate_index = gograph.G0[index]; // get the intermediate index
						ui original_id = G[intermediate_index]; // get the original id
						query.H.push_back(original_id);
					}

					query.ubD = ubD_temp;
				}
				query.hlock.unlock();

			}
		}

		if(context.VI.size() == query.N2){
			return;
		}

		unordered_set<ui> new2VI;

		if(EXE_new2VI){ // Inclusion based reduction
			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			gograph.inclusion_based_reduction(context, new2VI, kl_temp);

			if(context.VI.size() > query.N2){
				// undo inclusion based reduction
				for(auto e : new2VI){
					gograph.pop_from_VI_to_VR(context, e);
				}
				return;
			}
		}

		//// update the feasible community and upper bound of Diameter
		if(context.VI.size() >= query.N1 && context.VI.size() <= query.N2){
			ui cur_min_deg = INF;
			for(auto e : context.VI){
				if(context.degVI[e] < cur_min_deg)
					cur_min_deg = context.degVI[e];
			}

			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			if(cur_min_deg > kl_temp){
				query.hlock.rdlock();
				ubD_temp = query.ubD;
				query.hlock.unlock();

				ubD_temp = gograph.find_ubD(kl_temp, ubD_temp, query.N2);

				query.hlock.wrlock();
				if(cur_min_deg > query.kl){
					query.kl = cur_min_deg;
					// query.H = context.VI;

					// get reverse mapping of recoded id after pstart is shrinked
					query.H.clear();
					for(auto index: context.VI){
						ui intermediate_index = gograph.G0[index]; // get the intermediate index
						ui original_id = G[intermediate_index]; // get the original id
						query.H.push_back(original_id);
					}

					query.ubD = ubD_temp;
				}
				query.hlock.unlock();

			}
		}

		if(context.VI.size() == query.N2){
			if(EXE_new2VI){
				// undo inclusion based reduction
				for(auto e : new2VI){
					gograph.pop_from_VI_to_VR(context, e);
				}
			}
			return;
		}
		
		context.NEI.clear(); //// NEI is used to find vertices being dominated later on
		memset(context.inNEI, 0, sizeof(ui)*gograph.num_vertices); // for array version
		
		// get all the neighbors of vertices in VI
		gograph.get_NEI_of_VI(context);

		vector<ui> del_from_VR;
		if(EXE_del_from_VR){ // degree-based reduction
			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			gograph.degree_based_reduction(context, del_from_VR, kl_temp, query.N2);
		}

		int ustar = -1;

		if(srch_ord == 1) ustar = gograph.find_ustar(context);
		else if(srch_ord == 2) ustar = gograph.find_ustar_2phase(context, query.N2);
		else if(srch_ord == 3) ustar = gograph.find_ustar_mindeg(context);
		else if(srch_ord == 4) ustar = gograph.find_ustar_link(context);
		else if(srch_ord == 5) ustar = gograph.find_ustar_random(context);
		
		if(ustar < 0){ //// return if no vertex to branch on is found and undo changes before returning(if any)
			if(EXE_del_from_VR){
				// undo degree based reduction
				for(auto e : del_from_VR){
					gograph.restore_in_VR(context, e);
				}
			}

			if(EXE_new2VI){
				// undo inclusion based reduction
				for(auto e : new2VI){
					gograph.pop_from_VI_to_VR(context, e);
				}
			}
			return;
		}

		counter++;
		
		// add ustar to VI
		gograph.add_to_VI_from_VR(context, ustar);

		query.hlock.rdlock();
		kl_temp = query.kl;
		ubD_temp = query.ubD;
		query.hlock.unlock();

		if(gograph.q_dist[ustar] <= ubD_temp){
			gograph.gUB = gograph.compute_ub(context, kl_temp, query.N2);
			if(gograph.gUB){
				ftime(&cur_time);
				drun_time = cur_time.time-gograph.gtime_start.time+(double)(cur_time.millitm-gograph.gtime_start.millitm)/1000;
				if(drun_time < TIME_THRESHOLD){
					BB(context, query, gograph);
				}
				else{
					SCSTask * t = new SCSTask();

					gograph.ShrinkGraph(context, t->context.split_g);

					// for map version
					t->context.VI = context.VI;
					t->context.NEI = context.NEI;
					t->context.inVI = context.inVI;
					t->context.inVR = context.inVR;

					t->context.split_g.q_dist = new ui [gograph.num_vertices];
					memcpy(t->context.split_g.q_dist, gograph.q_dist, sizeof(ui)*(gograph.num_vertices));

					t->context.degVI = new ui[gograph.num_vertices];
					memcpy(t->context.degVI, context.degVI, sizeof(ui)*(gograph.num_vertices));

					t->context.degVIVR = new ui[gograph.num_vertices];
					memcpy(t->context.degVIVR, context.degVIVR, sizeof(ui)*(gograph.num_vertices));

					t->context.inNEI = new ui[gograph.num_vertices];
					memcpy(t->context.inNEI, context.inNEI, sizeof(ui)*(gograph.num_vertices));

					t->context.split_g.num_vertices = gograph.num_vertices;
					t->context.split_g.gUB = gograph.gUB;

					t->context.split_g.G0 = gograph.G0;

					add_task(t);
				}
			}
		}
		
		// remove ustar from VI
		gograph.pop_from_VI(context, ustar);
		
		vector<ui> rVR;
		bool del_v_in_VI = false;
		if(EXE_core_maintenance){
			query.hlock.rdlock();
			kl_temp = query.kl;
			query.hlock.unlock();

			gograph.core_maintenance(context, kl_temp, del_v_in_VI, rVR, ustar);
			if(del_v_in_VI){
				// restore vertices in rVR
				for(ui i = 0; i < rVR.size(); i++){
					ui v = rVR[i];
					gograph.restore_in_VR(context, v);
				}
				
				// restore ustar in VR
				gograph.restore_in_VR(context, ustar);

				// undo degree based reduction
				for(auto e : del_from_VR){
					gograph.restore_in_VR(context, e);
				}
				
				// undo inclusion based reduction
				for(auto e : new2VI){
					gograph.pop_from_VI_to_VR(context, e);
				}

				return;
			}
		}

		query.hlock.rdlock();
		kl_temp = query.kl;
		query.hlock.unlock();

		gograph.gUB = gograph.compute_ub(context, kl_temp, query.N2);
		if(gograph.gUB){
			ftime(&cur_time);
			drun_time = cur_time.time-gograph.gtime_start.time+(double)(cur_time.millitm-gograph.gtime_start.millitm)/1000;
			if(drun_time < TIME_THRESHOLD){
				BB(context, query, gograph);
			}
			else{
				//to split
				SCSTask * t = new SCSTask();

				gograph.ShrinkGraph(context,  t->context.split_g);

				// for map version
				t->context.VI = context.VI;
				t->context.NEI = context.NEI;
				t->context.inVI = context.inVI;
				t->context.inVR = context.inVR;

				t->context.split_g.q_dist = new ui [gograph.num_vertices];
				memcpy(t->context.split_g.q_dist, gograph.q_dist, sizeof(ui)*(gograph.num_vertices));

				t->context.degVI = new ui[gograph.num_vertices];
				memcpy(t->context.degVI, context.degVI, sizeof(ui)*(gograph.num_vertices));

				t->context.degVIVR = new ui[gograph.num_vertices];
				memcpy(t->context.degVIVR, context.degVIVR, sizeof(ui)*(gograph.num_vertices));

				t->context.inNEI = new ui[gograph.num_vertices];
				memcpy(t->context.inNEI, context.inNEI, sizeof(ui)*(gograph.num_vertices));

				t->context.split_g.num_vertices = gograph.num_vertices;
				t->context.split_g.gUB = gograph.gUB;

				t->context.split_g.G0 = gograph.G0;

				add_task(t);
			}
		}
		
		if(EXE_core_maintenance){
			// restore vertices in rVR
			for(ui i = 0; i < rVR.size(); i++){
				ui v = rVR[i];
				gograph.restore_in_VR(context, v);
			}
		}

		// restore ustar
		gograph.restore_in_VR(context, ustar);
		
		if(EXE_del_from_VR){
			// undo degree based reduction
			for(auto e : del_from_VR){
				gograph.restore_in_VR(context, e);
			}
		}
		if(EXE_new2VI){
			// undo inclusion based reduction
			for(auto e : new2VI){
				gograph.pop_from_VI_to_VR(context, e);
			}
		}
		return;
	}

	virtual bool toQuery(string& line, SCSQuery& query){
		stringstream iss(line);
		iss >> query.QID;
		iss >> query.N1;
		iss >> query.N2;

		return true;
	}

	virtual bool postprocess(SCSQuery& q) {//called at the end of a stage, returns false if the query is done.
		ftime(&q.end_t);
		double totaltime = q.end_t.time-q.start_t.time+(double)(q.end_t.millitm-q.start_t.millitm)/1000;

		cout << "[INFO] Query " << get_queryID() << endl;
		cout <<"Optimal Minimum Degree: "<<q.kl<<endl;
		cout << "Optimal Community: ";
		for(auto h: q.H){
			cout << h << " ";
		}
		cout << "\n";
		cout << "Size of Optimal Community: " << q.H.size() << endl;
		cout << "Total vertices explored: " << accumulate(q.counters.begin(), q.counters.end(), 0) << endl;
		cout << "Total time: " << totaltime << endl;

		fprintf(gfpout, "[INFO] Query ");
		fprintf(gfpout, "%d", get_queryID());
		fprintf(gfpout, "\n");
		fprintf(gfpout, "Optimal Minimum Degree: ");
		fprintf(gfpout, "%d ", q.kl);
		fprintf(gfpout, "\n");
		fprintf(gfpout, "Optimal Community: ");
		for(auto h: q.H){
			fprintf(gfpout, "%d ", h);
			fprintf(gfpout, " ");
		}
		fprintf(gfpout, "\n");
		fprintf(gfpout, "Size of Optimal Community: ");
		fprintf(gfpout, "%d", q.H.size());
		fprintf(gfpout, "\n");
		fprintf(gfpout, "Total vertices explored: ");
		fprintf(gfpout, "%d", accumulate(q.counters.begin(), q.counters.end(), 0));
		fprintf(gfpout, "\n");
		fprintf(gfpout, "Total time: ");
		fprintf(gfpout, "%f", totaltime);
		fprintf(gfpout, "\n");

		fflush(gfpout);
		return false;
	}

	virtual bool task_spawn(SCSQuery& query){
		cout << "task spawn for QueryId " << get_queryID() << endl;
		ftime(&query.start_t);
		SCSTask *task = new SCSTask();

		int query_id = get_queryID(); // get the id of query sent by user
		if (setup_task(task, query, gfpout, query_id))
		{
			add_task(task);
			return true;
		}

		delete task;
		return false;
	}

    virtual void compute(ContextT &context, SCSQuery &query){	
    	SCS& split_g = context.split_g;
		ftime(&split_g.gtime_start);
		counter = query.counters[thread_id];
		if(EXE_dom_ustar) BB_dom_ustar(context, query, split_g);
		else BB(context, query, split_g);
		query.counters[thread_id] = counter;
    }

	virtual bool is_bigTask(ContextT &context){	
		if(context.split_g.num_VIVR > BIGTASK_THRESHOLD){
			return true;
		}
		return false;
	}
};

class SCSWorker : public Worker<SCSComper>{
public:
    SCSWorker(int num_compers) : Worker(num_compers)
    {}

    ~SCSWorker(){	
		delete [] global_g.peel_sequence;
		delete [] global_g.degree;
		delete [] global_g.core;
		delete [] global_g.q_dist;
		delete [] global_g.pstart;
		delete [] global_g.G0_edges;
		delete [] global_g.G0_x;
    }


    void load_data(char* file_path){
        global_g.read_config(file_path);
		global_g.load_graph();
		global_g.core_decomposition_linear_list();
		
    }

};

#endif /* SCS_APP_H_ */
