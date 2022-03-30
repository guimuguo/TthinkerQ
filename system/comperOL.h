#ifndef COMPER_H_
#define COMPER_H_

#include "global.h"
#include <iostream>
#include <string>
#include <deque>
#include <queue>
#include "ioser.h"
#include "TaskProgMap.h"

using namespace std;
using namespace std::chrono;

#define SEQNO_LIMIT 9223372036854775807

template <class TaskT, class QueryT>
class Comper
{
public:
    typedef Comper<TaskT, QueryT> ComperT;
    // Used in worker.h
    typedef TaskT TaskType;
    // Used in worker.h
    typedef typename TaskT::ContextType ContextT;

    bool add_root_task_tag = false; // flag will set to true when task_spawn or postprocess

    char outpath[1000];
	char qfile[50];
	char msg_temp[1000];

    int thread_id;
    FILE *gfpout;
    unsigned int seqno = 0;

    int big_tasks_count = 0;

    TaskID cur_tid;
    int cur_qid;
    string q_msg;

    typedef std::chrono::_V2::steady_clock::time_point timepoint;
    // Comper's start time
    timepoint start_time;
    // To save latest elapsed time in seconds, in order to find exact duration for each method.
    float latest_elapsed_time = 0;

    thread main_thread;

    typedef deque<TaskT *> TaskQ;

    struct task_container
    {
		int qid;
		mutex Qreg_mtx, Qbig_mtx;
		TaskQ Qreg, Qbig;
		conque<string> Lreg, Lbig;
		QueryT q;
		map<int, FILE *> fout_map;

		task_container(int id){
			qid = id;
		}
    };

    typedef list<task_container*> Qlist;
    Qlist & activeQ_list = *(Qlist *)global_activeQ_list;

    task_container* tc;

	TaskID get_next_taskID() ////take it
	{
		TaskID id = thread_id;
		id = (id << 48); //first 16-bit is thread_id
		id += seqno;
		seqno++;
		assert(seqno < SEQNO_LIMIT);
		return id;
	}
    
    int get_queryID() 
    {
        return cur_qid;
    }

	void backtrack(task_prog * prog){
		prog->prog_lock.lock();
		prog->completed = true;
		if(prog->t_counter == 0){
			prog->prog_lock.unlock();

			bool flag = true;
			while(flag)
			{
				TaskID parent_id = prog->parent_tid;
				int query_id = prog->query_id;

				global_prog_map->erase(prog->tid);
				delete prog;

				if(parent_id == -1)
				{

					flag = false;

					add_root_task_tag = true;
					if(!postprocess(tc->q)) {
						//remove the query from activeQ_list
						activeQ_lock.wrlock();
						auto it = activeQ_list.begin();
						for(; it != activeQ_list.end(); ++it){
							tc = *it;
							if(tc->qid == query_id){
								//--------
								for(int i=0; i<num_compers; i++)
									fclose(tc->fout_map[i]);
								delete tc;
								activeQ_list.erase(it);
								activeQ_num--;

								cout<<"[INFO] Query "<<query_id<<" is done."<<endl;
								sprintf(outpath, "%d", query_id);
								notifier->send_msg(type, outpath);

								break;
							}
						}
						assert(it != activeQ_list.end());
						activeQ_lock.unlock();
					}
					add_root_task_tag = false;

				} else {
					prog = global_prog_map->get(parent_id);
					prog->prog_lock.lock();
					prog->t_counter--;
					flag = (prog->completed && prog->t_counter == 0);
					prog->prog_lock.unlock();
				}
			}
		}
		else prog->prog_lock.unlock();
	}

    Comper()
    {
    }

    virtual ~Comper()
    {
//        fclose(fgfpout);
        main_thread.join();
    }

    // UDF1
    virtual bool task_spawn(QueryT &q) = 0;

    // UDF2
    virtual void compute(ContextT &context, QueryT &q) = 0;
    // UDF2 wrapper
    void compute(TaskT *task)
    {
    	gfpout = tc->fout_map[thread_id];
    	cur_tid = task->get_id();
    	cur_qid = task->get_qid();
        compute(task->context, tc->q);// compute(task->context, task);
        backtrack(task->prog);
    }

    // UDF3
    virtual bool is_bigTask(ContextT & context)
    {
        return false;
    }

    //UDF4: queryLoader ==============================
	virtual bool toQuery(string& line, QueryT& q) = 0; //this is what user specifies!!!!!!

	virtual bool postprocess(QueryT& q) {//called at the end of a stage, returns false if the query is done.
		return false;
	}

    void start(int thread_id)
    {
//        gfpout = fopen(("output_" + to_string(thread_id)).c_str(), "wt");
        this->thread_id = thread_id;

        start_time = steady_clock::now();

        main_thread = thread(&ComperT::run, this);
    }

    string fname;

    long long bigFileSeqNo = 1;
    void set_bigTask_fname()
    {
        fname = "";
        fname += TASK_DISK_BUFFER_DIR + "/" + to_string(thread_id) + "_" + to_string(bigFileSeqNo) + "_bt";
        bigFileSeqNo++;
    }

    long long regFileSeqNo = 1;
    void set_regTask_fname()
    {
        fname = "";
        fname += TASK_DISK_BUFFER_DIR + "/" + to_string(thread_id) + "_" + to_string(regFileSeqNo) + "_rt";
        regFileSeqNo++;
    }

	void create_output_path(const char* msg){
		//create folder outpath (qid_query)
		strcpy(outpath, out_path.c_str());
		sprintf(qfile, "/query%d_", cur_qid);
		strcat(outpath, qfile);
		//loop char* msg and convert '/' and ' ' to '_'
		strcpy(msg_temp, msg);
		for (int i=0; *(msg_temp+i) != '\0'; i++)
			if(*(msg_temp+i) == '/' || *(msg_temp+i) == ' ')
				msg_temp[i]='_';

		strcat(outpath, msg_temp);
		recursive_mkdir(outpath);

		//create file named comper_id
		for(int i=0; i<num_compers; i++)
		{
			char comper_file[200], comper_id[50];
			strcpy(comper_file, outpath);
			sprintf(comper_id, "/%d", i);
			strcat(comper_file, comper_id);

			tc->fout_map[i] = fopen(comper_file, "wt");
		}
	}

    bool refill_Qbig()
    {
        string file;
        bool succ = tc->Lbig.dequeue(file);
        if (!succ)
            // "global_Lbig" is empty
            return false; 
        else
        {
            ofbinstream in(file.c_str());
            while (!in.eof())
            {
                TaskT *task;
                in >> task;
                add_refilled_task(task);
            }
            in.close();

            if (remove(file.c_str()) != 0)
            {
                log("Error removing file: " + file);
                log("Error printed by perror");
            }
            return true;
        }
    }

    bool refill_Qreg()
    {
        string file;
        bool succ = tc->Lreg.dequeue(file);
        // Lreg is not empty, refill from disk.
        if (succ)
        {
            log("Comper:refill from Lreg !!");
            ofbinstream in(file.c_str());
            while (!in.eof())
            {
                TaskT *task;
                in >> task;
                add_refilled_task(task);
            }
            in.close();
            if (remove(file.c_str()) != 0)
            {
                log("Error removing file: " + file);
                log("Error printed by perror");
            }
            return true;
        }
        else return false;
    }

    void spill_Qbig()
    {
        int i = 0;
        // Fetch tasks into local queue first, to avoid locking while spilling to disk.
        queue<TaskT *> collector;
        while (i < BT_TASKS_PER_FILE && !tc->Qbig.empty())
        {
            // Get a task from the tail
            TaskT *t = tc->Qbig.back();
            tc->Qbig.pop_back();
            collector.push(t);
            i++;
        }
        tc->Qbig_mtx.unlock();

        if (!collector.empty())
        {
            set_bigTask_fname();
            ifbinstream bigTask_out(fname.c_str());

            while (!collector.empty())
            {

                TaskT *t = collector.front();
                collector.pop();
                // Stream to file
                bigTask_out << t;
                // Release from memory
                delete t;
            }

            bigTask_out.close();
            tc->Lbig.enqueue(fname);
        }
    }

    void spill_Qreg()
    {
        int i = 0;
        queue<TaskT *> collector;
        while (i < RT_TASKS_PER_FILE && !tc->Qreg.empty())
        {
            TaskT *t = tc->Qreg.back();
            tc->Qreg.pop_back();
            collector.push(t);
            i++;
        }
        tc->Qreg_mtx.unlock();

        if (!collector.empty())
        {
            set_regTask_fname();
            ifbinstream regTask_out(fname.c_str());

            while (!collector.empty())
            {
                TaskT *t = collector.front();
                collector.pop();
                // Stream to file
                regTask_out << t;
                // Release from memory
                delete t;
            }
            regTask_out.close();
            tc->Lreg.enqueue(fname);
        }
    }

//    //for task spawn
//    bool add_root_task(TaskT *task)
//    {
//    	task->prog = new task_prog(get_next_taskID(), -1, cur_qid);
//    	global_prog_map->insert(task->prog->tid, task->prog);
//
//        if (is_bigTask(task))
//        {
//            add_bigTask(task);
//            return true;
//        }
//
//        add_regTask(task);
//        return false;
//    }
//
//    //for task split
//    bool add_task(TaskT *task)
//    {
//    	task->prog = new task_prog(get_next_taskID(), cur_tid, cur_qid);
//    	global_prog_map->insert(task->prog->tid, task->prog);
//    	task_prog* parent_prog = global_prog_map->get(cur_tid);
//    	parent_prog->t_counter++; //no need lock since it happens before "completed -> ture"
//
//        if (is_bigTask(task))
//        {
//            add_bigTask(task);
//            return true;
//        }
//
//        add_regTask(task);
//        return false;
//    }


	void add_task(TaskT *task)
	{
		if (add_root_task_tag){
			//set parent prog id = -1 for root task.
			task->prog = new task_prog(get_next_taskID(), -1, cur_qid);
		} else {
			task->prog = new task_prog(get_next_taskID(), cur_tid, cur_qid);
			task_prog* parent_prog = global_prog_map->get(cur_tid);
			parent_prog->t_counter++; //no need lock since it happens before "completed -> ture"
		}
		global_prog_map->insert(task->prog->tid, task->prog);

		if (is_bigTask(task->context))
			add_bigTask(task);
		else
			add_regTask(task);
	}

    //for task refill; progress object is still in the global table
    bool add_refilled_task(TaskT *task)
	{
		if (is_bigTask(task->context))
		{
			add_bigTask(task);
			return true;
		}

		add_regTask(task);
		return false;
	}

    void add_bigTask(TaskT *task)
    {
        tc->Qbig_mtx.lock();
        // Check if spill is needed.
        while (tc->Qbig.size() >= Qbig_capacity)
        {
            spill_Qbig();
            tc->Qbig_mtx.lock();
        }
        tc->Qbig.push_back(task);
        tc->Qbig_mtx.unlock();
    }

    void add_regTask(TaskT *task)
    {
        tc->Qreg_mtx.lock();
        while (tc->Qreg.size() >= Qreg_capacity)
        {
            spill_Qreg();
            tc->Qreg_mtx.lock();
        }
        tc->Qreg.push_back(task);
        tc->Qreg_mtx.unlock();
    }

    bool get_and_process_tasks()
    {
    	//check activeQ_list:
		// 1- Check Qbig first
		// - if Qbig's size is less than Qbig_capacity, refill using Lbig
		// - if Qbig is not empty, pop a task and compute.
		// - if Qbig is empty, check Qreg
		// 2- Check Qreg
		// - if Qreg's size is less than Qreg_capacity, refill using Lreg
		// - if Qreg is not empty, pop a task
		// - if Qreg is empty, check active_Qlist
        TaskT *task = NULL;
        activeQ_lock.rdlock();
        if(!activeQ_list.empty())
        {
        	bool refilled = false;
        	for(auto it = activeQ_list.begin(); it != activeQ_list.end(); ++it)
			{
        		tc = *it;

				if (tc->Qbig_mtx.try_lock())
				{
					if (tc->Qbig.size() < BT_THRESHOLD_FOR_REFILL)
					{
						tc->Qbig_mtx.unlock();
						refilled = refill_Qbig();
					}
					else
						tc->Qbig_mtx.unlock();

					if (tc->Qbig_mtx.try_lock())
					{
						if (!tc->Qbig.empty())
						{
							task = tc->Qbig.front();
							tc->Qbig.pop_front();

							tc->Qbig_mtx.unlock();
							activeQ_lock.unlock();
							compute(task);

							delete task;
							return true;
						}
						else
							tc->Qbig_mtx.unlock();
					}
				}

				// Means Qbig is empty
				if (task == NULL)
				{
					//try lock until the last one in activeQ_list
					if (next(it) != activeQ_list.end())
					{
						if(tc->Qreg_mtx.try_lock())
						{
							// Refill Qreg using Lreg.
							if (tc->Qreg.size() < RT_THRESHOLD_FOR_REFILL)
							{
								tc->Qreg_mtx.unlock();
								refilled = refill_Qreg();
							}
							else
								tc->Qreg_mtx.unlock();

							queue<TaskT *> collector;
							size_t tasks_per_fetch = tasks_per_fetch_g; // _g: global variable
							if(tc->Qreg_mtx.try_lock())
							{
								while (!tc->Qreg.empty() && tasks_per_fetch > 0)
								{
									TaskT *task = tc->Qreg.front();
									tc->Qreg.pop_front();
									collector.push(task);
									tasks_per_fetch--;
								}
								tc->Qreg_mtx.unlock();

								if (collector.empty() && !refilled)
								{
									log("collector is empty...");
									continue;
								}

								activeQ_lock.unlock();
								// Process tasks in "collector"
								while (!collector.empty())
								{
									TaskT *task = collector.front();
									collector.pop();
									compute(task);
									delete task;
								}
								return true;
							}
						}
					}
					else
					{
						tc->Qreg_mtx.lock(); //try lock until the last one
						// Refill Qreg using Lreg.
						if (tc->Qreg.size() < RT_THRESHOLD_FOR_REFILL)
						{
							tc->Qreg_mtx.unlock();
							refilled = refill_Qreg();
						}
						else
							tc->Qreg_mtx.unlock();

						queue<TaskT *> collector;
						size_t tasks_per_fetch = tasks_per_fetch_g; // _g: global variable
						tc->Qreg_mtx.lock();
						while (!tc->Qreg.empty() && tasks_per_fetch > 0)
						{
							TaskT *task = tc->Qreg.front();
							tc->Qreg.pop_front();
							collector.push(task);
							tasks_per_fetch--;
						}
						tc->Qreg_mtx.unlock();

						if (collector.empty() && !refilled)
						{
							log("collector is empty...");
							continue;
						}

						activeQ_lock.unlock();
						// Process tasks in "collector"
						while (!collector.empty())
						{
							TaskT *task = collector.front();
							collector.pop();
							compute(task);
							delete task;
						}
						return true;
					}
				}
			}
        }
        activeQ_lock.unlock();
        //being here means did not get any task
        bool succ = false;
        if(!query_que.empty()){
        	activeQ_lock.wrlock();
			//case activeQ_list may be full but task num is less than comper num
			if(activeQ_num < activeQ_list_capacity && !query_que.empty()){
				activeQ_num++;
				activeQ_lock.unlock();
				qinfo spawn_qinfo;
				succ = query_que.dequeue(spawn_qinfo);
				// if query_que is not empty, spawn tasks from query_que.
				if (succ)
				{
					tc = new task_container(spawn_qinfo.qid);

					cur_qid = spawn_qinfo.qid;
					q_msg = spawn_qinfo.q;

					if(toQuery(spawn_qinfo.q, tc->q))
					{
						//create output file
						create_output_path(q_msg.c_str());
						gfpout = tc->fout_map[thread_id];

						add_root_task_tag = true;
						if(!task_spawn(tc->q))
						{
							//close output file
							for(int i=0; i<num_compers; i++)
								fclose(tc->fout_map[i]);

							delete tc;
							activeQ_lock.wrlock();
							activeQ_num--;
							activeQ_lock.unlock();

							cout<<"[INFO] Query "<<cur_qid<<" \""<<spawn_qinfo.q<<"\" spawns no task."<<endl;

							sprintf(outpath, "%d", cur_qid);
							notifier->send_msg(type, outpath);
						}
						else
						{
							activeQ_lock.wrlock();
							activeQ_list.push_back(tc);
							activeQ_lock.unlock();
						}
						add_root_task_tag = false;

					} else {
						delete tc;
						activeQ_lock.wrlock();
						activeQ_num--;
						activeQ_lock.unlock();

						cout<<"[INFO] Query "<<cur_qid<<" \""<<spawn_qinfo.q<<"\" is invalid."<<endl;

						sprintf(outpath, "%d", cur_qid);
						notifier->send_msg(type, outpath);
					}
				}
			} else activeQ_lock.unlock();
        }
        //being here means no task in activeQ_lock and can not spawn q in query_que
        return succ;
    }

    void run()
    {
        while (global_end_label == false) // Otherwise, thread terminates
        {
            // Process task or batch of tasks
            bool task_found = get_and_process_tasks();
            // Means that the Queues are empty
            if (!task_found)
            {
                unique_lock<mutex> lck(mtx_go);
                ready_go = false;
                global_num_idle++;
                while (!ready_go)
                {
                    log("wait!!!");
                    cv_go.wait(lck);
                }
                global_num_idle--;
            }
        }
    }
};

#endif
