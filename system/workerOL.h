#ifndef WORKER_H_
#define WORKER_H_

#include "comperOL.h"
#include <unistd.h>
#include <omp.h>
#include "TaskProgMap.h"

using namespace std;

template <class ComperT>
class Worker
{
public:
    typedef typename ComperT::TaskType TaskT;
    typedef deque<TaskT *> TaskQ;
    typedef typename ComperT::Qlist Qlist;
//    typedef hash_map<int, TaskT> QMap;

    unsigned int seqno = 0;

    // Dynamic array of compers
    ComperT *compers = nullptr;
    // Contains all data loaded from file
    // Contains pointers to data-array without initialized tasks pointers
    // Regular tasks queue
    TaskQ *Qreg;
    // Big tasks queue
    TaskQ *Qbig;

    Qlist *activeQ_list;

    typedef std::chrono::_V2::steady_clock::time_point timepoint;
    // Worker's start time
    timepoint start_time = steady_clock::now();
    // To save latest elapsed time in seconds, in order to find exact duration for each method
    float latest_elapsed_time = 0;

    float init_time;
    float load_data_time;

    // Save files sequence number when spill big tasks to disk
    vector<size_t> big_files_seq;
    // Save files sequence number when spill reg tasks to disk
    vector<size_t> reg_files_seq;

    vector<string> big_files_names;
    vector<string> reg_files_names;

    //query Variables
//    QMap queries;
    msg_queue_server* server;

	size_t nxt_qid;

	const char* output_folder;
	char outpath[200];
	char qfile[50];

    // Init files seq with 1 for each thread.
    Worker(int comper_num) : big_files_seq(comper_num, 1), reg_files_seq(comper_num, 1), big_files_names(comper_num), reg_files_names(comper_num)
    {
        // Create disk buffer dir
        TASK_DISK_BUFFER_DIR = "buffered_tasks";
        recursive_mkdir(TASK_DISK_BUFFER_DIR.c_str());

//        global_Qreg = Qreg = new TaskQ;
//        global_Qbig = Qbig = new TaskQ;
        global_activeQ_list = activeQ_list = new Qlist;
        num_compers = comper_num;
        global_end_label = false;
        global_prog_map = new TaskProgMap();

        server=new msg_queue_server;
        notifier=new msg_queue_notifier;
//		notifier=new msg_queue_notifier;
		nxt_qid = 1;
		type = 1;
    }

    virtual ~Worker()
    {
        if (compers)
            delete[] compers;
//        delete Qreg;
//        delete Qbig;

        delete global_prog_map;
        delete activeQ_list;
        delete server;
        delete notifier;
    }

    // UDF1: read data from file_path
//    virtual void load_data(const string &file_path) {}

    // UDF2
//    virtual bool task_spawn(DataT &data) = 0;
//    virtual bool task_spawn() = 0;

    // UDF3
//    virtual bool is_bigTask(TaskT *task)
//    {
//        return false;
//    }

    TaskID get_next_taskID() ////take it
	{
    	TaskID id = num_compers;
		id = (id << 48); //first 16-bit is thread_id
		id += seqno;
		seqno++;
		assert(seqno < SEQNO_LIMIT);
		return id;
	}


    void set_bigTask_fname(const int thread_id)
    {
        // Reset filename
        big_files_names[thread_id] = "";
        big_files_names[thread_id] += TASK_DISK_BUFFER_DIR + "/w_" + to_string(thread_id) + "_" + to_string(big_files_seq[thread_id]) + "_bt";
        big_files_seq[thread_id]++;
    }

    void set_regTask_fname(const int thread_id)
    {
        // Reset filename
        reg_files_names[thread_id] = "";
        reg_files_names[thread_id] += TASK_DISK_BUFFER_DIR + "/w_" + to_string(thread_id) + "_" + to_string(reg_files_seq[thread_id]) + "_rt";
        reg_files_seq[thread_id]++;
    }

    void create_compers()
    {
        compers = new ComperT[num_compers];
        for (int i = 0; i < num_compers; i++)
        {
            compers[i].start(i);
        }
    }


    bool update_tasks()//return true when there is queries arrived
	{
		vector<qinfo> new_queries;

		while(server->recv_msg(type))
		{
			char* msg=server->get_msg();
			cout<<"Q"<<nxt_qid<<": "<<msg<<endl;
			if(strcmp(msg, "server_exit")==0)
			{
				server_exit = true;
				break;
			}
			else
			{
				qinfo q(nxt_qid, msg);
				new_queries.push_back(q);
				nxt_qid++;
			}
		}
		//add from vector to deque
		if(!new_queries.empty())
		{
			query_que.enqueue(new_queries);
			return true;
		}
		return false;
	}

    // Program entry point
    void run()
    {
        assert(RT_TASKS_PER_FILE <= Qreg_capacity);
        assert(BT_TASKS_PER_FILE <= Qbig_capacity);
        // Kernel_app's load_data method add tasks directly to two queues not to data_array.
//        if (data_array.size() > 0)
//        {
//            // Initialize some tasks
//            auto start = steady_clock::now();
//            initialize_tasks();
//            auto end = steady_clock::now();
//            init_time = (float)duration_cast<milliseconds>(end - start).count() / 1000;
//            cout << "initialize_tasks() execution time:" << init_time << endl;
//        }

//        output_folder = params.output_path.c_str();
        output_folder = out_path.c_str();

        // Setup computing threads
        create_compers();

        // Call status_sync() periodically
        while (global_end_label == false)
        {
            // Avoid busy-checking
            // usleep(WAIT_TIME_WHEN_IDLE);

            if(!server_exit)//will not recieve new queries after "server_exit"
            {
            	if(update_tasks())
            	{
            		mtx_go.lock();
					ready_go = true;
					// Release threads to spawn task from query
					cv_go.notify_all();
					mtx_go.unlock();
            	}
            }

            activeQ_lock.rdlock();
            if(activeQ_num > 0)
            {
            	activeQ_lock.unlock();
            	mtx_go.lock();
				ready_go = true;
				// Release threads to compute tasks
				cv_go.notify_all();
				mtx_go.unlock();
            }
            else {
            	activeQ_lock.unlock();
            	if(!query_que.empty()){
					mtx_go.lock();
					ready_go = true;
					// Release threads to compute tasks
					cv_go.notify_all();
					mtx_go.unlock();
				}
				else{
					mtx_go.lock();
					if (global_num_idle == num_compers && server_exit)
					{
						// every thread is waiting, guaranteed since mtx_go is locked
						// Since we are in else-branch, Qreg must be empty
						cout << "global_num_idle: " << global_num_idle << endl;
						global_end_label = true;
						ready_go = true;
						// Release threads
						cv_go.notify_all();
					}
					// else, some threads are still processing tasks, check in next round
					mtx_go.unlock();
				}
            }
        }
    }
};

#endif
