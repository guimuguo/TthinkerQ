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

template <class TaskT, class DataT, class QueryT>
class Comper
{
public:
    typedef Comper<TaskT, DataT, QueryT> ComperT;
    // Used in worker.h
    typedef TaskT TaskType;
    // Used in worker.h
    typedef DataT DataType;
    typedef typename TaskT::ContextType ContextT;

    char outpath[200];
	char qfile[50];

    int thread_id;
    FILE *gfpout;
    unsigned int seqno = 0;

    int big_tasks_count = 0;

    TaskID cur_tid;
    int cur_qid;
    int spawn_qid;
    string q_msg;

    typedef stack<qinfo> DataStack;//todo delete

    typedef std::chrono::_V2::steady_clock::time_point timepoint;
    // Comper's start time
    timepoint start_time;
    // To save latest elapsed time in seconds, in order to find exact duration for each method.
    float latest_elapsed_time = 0;

    thread main_thread;

    DataStack &data_stack = *(DataStack *)global_data_stack;//todo delete
    typedef deque<TaskT *> TaskQ;
    TaskQ &Qreg = *(TaskQ *)global_Qreg;
    TaskQ &Qbig = *(TaskQ *)global_Qbig;

	TaskID get_next_taskID() ////take it
	{
		TaskID id = thread_id;
		id = (id << 48); //first 16-bit is thread_id
		id += seqno;
		seqno++;
		assert(seqno < SEQNO_LIMIT);
		return id;
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
					for(int i=0; i<num_compers; i++)
					{
						fclose(fout_map[query_id][i]);
					}
					cout<<"[INFO] Query "<<query_id<<" is done."<<endl;
					flag = false;
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
    virtual bool task_spawn(QueryT q) = 0;

    // UDF2
    virtual void compute(ContextT &context) = 0;
    // UDF2 wrapper
    void compute(TaskT *task)
    {
    	gfpout = fout_map[task->get_qid()][thread_id];
    	cur_tid = task->get_id();
    	cur_qid = task->get_qid();
        compute(task->context);// compute(task->context, task);
        backtrack(task->prog);
    }

    // UDF3
    virtual bool is_bigTask(TaskT *task)
    {
        return false;
    }

    //UDF4: qyeryLoader ==============================
	virtual QueryT toQuery(string& line) = 0; //this is what user specifies!!!!!!


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
		sprintf(qfile, "/query%d_", spawn_qid);
		strcat(outpath, qfile);
		strcat(outpath, msg);
		recursive_mkdir(outpath);

		//create file named comper_id
		for(int i=0; i<num_compers; i++)
		{
			char comper_file[200], comper_id[50];
			strcpy(comper_file, outpath);
			sprintf(comper_id, "/%d", i);
			strcat(comper_file, comper_id);

			fout_map[spawn_qid][i] = fopen(comper_file, "wt");
		}
	}

    bool refill_Qbig()
    {
        string file;
        bool succ = global_Lbig.dequeue(file);
        if (!succ)
            // "global_Lbig" is empty
            return false; 
        else
        {
            global_Lbig_num--;
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
        bool succ = global_Lreg.dequeue(file);
        // 1- Lreg is not empty, refill from disk.
        if (succ)
        {
            log("Comper:refill from Lreg !!");
            global_Lreg_num--;
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
        // 2- Check data_stack to refill
        else
        {
//            data_stack_mtx.lock();
//            if (!data_stack.empty())
//            {
//            	qinfo spawn_qinfo = data_stack.top();
//            	data_stack.pop();
//            	data_stack_mtx.unlock();
//            	QueryT q = toQuery(spawn_qinfo.q);
//            	spawn_qid = spawn_qinfo.qid;
//            	q_msg = spawn_qinfo.q;
//            	if(!task_spawn(q))
//            		cout<<"[INFO] Query "<<spawn_qinfo.qid<<" \""<<spawn_qinfo.q<<"\" spawns no task."<<endl;
//				return true;
//            }
//            data_stack_mtx.unlock();


        	qinfo spawn_qinfo;
        	bool succ = query_que.dequeue(spawn_qinfo);
			// if query_que is not empty, spawn tasks from query_que.
			if (succ)
			{
				QueryT q = toQuery(spawn_qinfo.q);
				spawn_qid = spawn_qinfo.qid;
				q_msg = spawn_qinfo.q;
				if(!task_spawn(q))
					cout<<"[INFO] Query "<<spawn_qinfo.qid<<" \""<<spawn_qinfo.q<<"\" spawns no task."<<endl;
				return true;
			}
            // Nothing to refill.
            return false;
        }
    }

    void spill_Qbig()
    {
        int i = 0;
        // Fetch tasks into local queue first, to avoid locking while spilling to disk.
        queue<TaskT *> collector;
        while (i < BT_TASKS_PER_FILE && !Qbig.empty())
        {
            // Get a task from the tail
            TaskT *t = Qbig.back();
            Qbig.pop_back();
            collector.push(t);
            i++;
        }
        Qbig_mtx.unlock();

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
            global_Lbig.enqueue(fname);
            global_Lbig_num++;
        }
    }

    void spill_Qreg()
    {
        int i = 0;
        queue<TaskT *> collector;
        while (i < RT_TASKS_PER_FILE && !Qreg.empty())
        {
            TaskT *t = Qreg.back();
            Qreg.pop_back();
            collector.push(t);
            i++;
        }
        Qreg_mtx.unlock();

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
            global_Lreg.enqueue(fname);
            global_Lreg_num++;
        }
    }

    //for task spawn
    bool add_root_task(TaskT *task)
    {
    	create_output_path(q_msg.c_str());
    	task->prog = new task_prog(get_next_taskID(), -1, spawn_qid);
    	global_prog_map->insert(task->prog->tid, task->prog);

        if (is_bigTask(task))
        {
            add_bigTask(task);
            return true;
        }

        add_regTask(task);
        return false;
    }

    //for task split
    bool add_task(TaskT *task)
    {
    	task->prog = new task_prog(get_next_taskID(), cur_tid, cur_qid);
    	global_prog_map->insert(task->prog->tid, task->prog);
    	task_prog* parent_prog = global_prog_map->get(cur_tid);
    	parent_prog->t_counter++; //no need lock since it happens before "completed -> ture"

        if (is_bigTask(task))
        {
            add_bigTask(task);
            return true;
        }

        add_regTask(task);
        return false;
    }

    //for task refill; progress object is still in the global table
    bool add_refilled_task(TaskT *task)
	{
		if (is_bigTask(task))
		{
			add_bigTask(task);
			return true;
		}

		add_regTask(task);
		return false;
	}

    void add_bigTask(TaskT *task)
    {
        Qbig_mtx.lock();
        // Check if spill is needed.
        while (Qbig.size() >= Qbig_capacity)
        {
            spill_Qbig();
            Qbig_mtx.lock();
        }
        Qbig.push_back(task);
        Qbig_mtx.unlock();
    }

    void add_regTask(TaskT *task)
    {
        Qreg_mtx.lock();
        while (Qreg.size() >= Qreg_capacity)
        {
            spill_Qreg();
            Qreg_mtx.lock();
        }
        Qreg.push_back(task);
        Qreg_mtx.unlock();
    }

    bool get_and_process_tasks()
    {
        // 1- Check Qbig first
        // - if Qbig's size is less than Qbig_capacity, refill using Lbig
        // - if Qbig is not empty, pop a task and compute.
        // - if Qbig is empty, check Qreg
        // 2- Check Qreg
        // - if Qreg's size is less than Qreg_capacity, refill using Lreg
        // - if Qreg is not empty, pop a task
        // - if Qreg is empty, check data_array
        TaskT *task = NULL;

        if (Qbig_mtx.try_lock())
        {
            if (Qbig.size() < BT_THRESHOLD_FOR_REFILL)
            {
                Qbig_mtx.unlock();
                refill_Qbig();
            }
            else
                Qbig_mtx.unlock();

            if (Qbig_mtx.try_lock())
            {
                if (!Qbig.empty())
                {
                    task = Qbig.front();
                    Qbig.pop_front();

                    Qbig_mtx.unlock();
                    compute(task);

                    delete task;
                    return true;
                }
                else
                    Qbig_mtx.unlock();
            }
        }

        // Means Qbig is empty
        if (task == NULL) 
        {
            Qreg_mtx.lock();
            // Refill Qreg using Lreg.
            // If Lreg is empty check data_array.
            bool refilled = false;
            if (Qreg.size() < RT_THRESHOLD_FOR_REFILL)
            {
                Qreg_mtx.unlock();
                refilled = refill_Qreg();
            }
            else
                Qreg_mtx.unlock();

            queue<TaskT *> collector;
            size_t tasks_per_fetch = tasks_per_fetch_g; // _g: global variable
            Qreg_mtx.lock();
            while (!Qreg.empty() && tasks_per_fetch > 0)
            {
                TaskT *task = Qreg.front();
                Qreg.pop_front();
                collector.push(task);
                tasks_per_fetch--;
            }
            Qreg_mtx.unlock();

            if (collector.empty() && !refilled)
            {
                log("collector is empty...");
                return false;
            }

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
