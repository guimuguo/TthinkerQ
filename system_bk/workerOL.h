#ifndef WORKER_H_
#define WORKER_H_

#include <iostream>
#include <string>
#include "comper.h"
#include <unistd.h>
#include <omp.h>
#include "TaskProgMap.h"
#include "msgtool.h"

using namespace std;

template <class ComperT, class QueryT>
class Worker
{
public:
    typedef typename ComperT::TaskType TaskT;
    typedef typename ComperT::DataType DataT;
//    typedef typename ComperT::QueryT QueryT;
    typedef deque<TaskT *> TaskQ;
    typedef stack<DataT *> DataStack;
//    typedef hash_map<int, TaskT> QMap;

    unsigned int seqno = 0;

    // Dynamic array of compers
    ComperT *compers = nullptr;
    // Contains all data loaded from file
    vector<DataT *> data_array;
    // Contains pointers to data-array without initialized tasks pointers
    DataStack *data_stack;
    // Regular tasks queue
    TaskQ *Qreg;
    // Big tasks queue
    TaskQ *Qbig;
    QueryT cur_q;
    int cur_qid;

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
//	msg_queue_notifier* notifier;//ADDED FOR AUTO
	long type;
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

        global_data_stack = data_stack = new stack<DataT *>;

        global_Qreg = Qreg = new TaskQ;
        global_Qbig = Qbig = new TaskQ;
        num_compers = comper_num;
        global_end_label = false;
        global_prog_map = new TaskProgMap();

        server=new msg_queue_server;
//		notifier=new msg_queue_notifier;
		nxt_qid = 1;
		type = 1;
    }

    virtual ~Worker()
    {
        for (int i = 0; i < data_array.size(); i++)
        {
            delete data_array[i];
        }

        if (compers)
            delete[] compers;
        delete Qreg;
        delete Qbig;
        delete global_prog_map;
    }

    // UDF1: read data from file_path, and insert into data_array
    virtual void load_data(const string &file_path) {}

    // UDF2
//    virtual bool task_spawn(DataT &data) = 0;
    virtual bool task_spawn() = 0;

    // UDF3
    virtual bool is_bigTask(TaskT *task)
    {
        return false;
    }

    TaskID get_next_taskID() ////take it
	{
    	TaskID id = num_compers;
		id = (id << 48); //first 16-bit is thread_id
		id += seqno;
		seqno++;
		assert(seqno < SEQNO_LIMIT);
		return id;
	}

    // Insert some tasks into Qreg before spawn compers, so compers have some tasks to work on
    void initialize_tasks()
    {

        // 1- Spawn some tasks from data_array
        size_t _size = min(Qreg_capacity, data_array.size());

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_compers)
        for (int i = 0; i < _size; i++)
        {
            task_spawn(*(data_array[i]));
        }
        // 2- Add the rest of data_array to a data stack, to be used by comper when spawn
        for (int i = data_array.size() - 1; i >= _size; i--)
        {
            data_stack->push(data_array[i]);
        }
    }

    bool add_task(TaskT *task, int query_id)
    {
    	task->prog = new task_prog(get_next_taskID(), -1, query_id);
    	global_prog_map->insert(task->prog->tid, task->prog);
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
        while (Qbig->size() >= Qbig_capacity)
        {
            spill_Qbig();
            Qbig_mtx.lock();
        }
        Qbig->push_back(task);
        Qbig_mtx.unlock();
    }

    void add_regTask(TaskT *task)
    {
        Qreg_mtx.lock();
        while (Qreg->size() >= Qreg_capacity)
        {
            spill_Qreg();
            Qreg_mtx.lock();
        }
        Qreg->push_back(task);
        Qreg_mtx.unlock();
    }

    void spill_Qbig()
    {
        int i = 0;
        queue<TaskT *> collector;
        while (i < BT_TASKS_PER_FILE && !Qbig->empty())
        {
            // Get task at the tail
            TaskT *t = Qbig->back();
            Qbig->pop_back();
            collector.push(t);
            i++;
        }
        Qbig_mtx.unlock();

        if (!collector.empty())
        {
            int thread_id = omp_get_thread_num();
            set_bigTask_fname(thread_id);
            ifbinstream bigTask_out(big_files_names[thread_id].c_str());

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
            global_Lbig.enqueue(big_files_names[thread_id]);
            global_Lbig_num++;
        }
    }

    void spill_Qreg()
    {
        int i = 0;
        queue<TaskT *> collector;
        while (i < RT_TASKS_PER_FILE && !Qreg->empty())
        {
            // Get task at the tail
            TaskT *t = Qreg->back();
            Qreg->pop_back();
            collector.push(t);
            i++;
        }
        Qreg_mtx.unlock();

        if (!collector.empty())
        {
            int thread_id = omp_get_thread_num();
            set_regTask_fname(thread_id);
            ifbinstream regTask_out(reg_files_names[thread_id].c_str());

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
            global_Lreg.enqueue(reg_files_names[thread_id]);
            global_Lreg_num++;
        }
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


    //user-defined qyeryLoader ==============================
	virtual QueryT toQuery(char* line) = 0; //this is what user specifies!!!!!!

	struct qinfo
	{
		int qid;
		char* q;

		qinfo(){}

		qinfo(int qid, char* q)
		{
			this->qid=qid;
			this->q=q;
		}

//		friend ibinstream& operator<<(ibinstream& m, const qinfo& v)
//		{
//			m << v.qid;
//			m << v.q;
//			return m;
//		}
//
//		friend obinstream& operator>>(obinstream& m, qinfo& v)
//		{
//			m >> v.qid;
//			m >> v.q;
//			return m;
//		}
	};

	void create_output_path(){
		strcpy(outpath, output_folder);
		sprintf(qfile, "/query%d", nxt_qid);
		strcat(outpath, qfile);
		//todo create folder outpath
		for(int i=0; i<num_compers; i++)
		{
			char comper_file[200], comper_id[50];
//			int tid = omp_get_thread_num();
			strcpy(comper_file, outpath);
			sprintf(comper_id, "/%d", i);
			strcat(comper_file, comper_id);
			//create file named comper_file
		}
	}

	//functions to be called in UDF init():
	//todo remove the 2 func
	QueryT get_query()
	{
		return cur_q;
	}
	int get_query_id()
	{
		return cur_qid;
	}

	//todo 1, 2, 3, 4, |server_exit|, 5, 6 => shut down the server when all queries are finished () //only the query after
    bool update_tasks()//return false to shut_down server
	{
//		vector<qinfo> new_queries;

		while(server->recv_msg(type))
		{
			char* msg=server->get_msg();
			cout<<"Q"<<nxt_qid<<": "<<msg<<endl;
			if(strcmp(msg, "server_exit")==0)
			{
//				new_queries.clear();
				return false;
			}
			else
			{
				//todo move toQuery to comper
				QueryT q=toQuery(msg);
				qinfo qentry(nxt_qid, q);//todo check ??? qinfo <id, string>
//				new_queries.push_back(qentry);
//				TaskT& task=queries[nxt_qid]=TaskT(q);
				//create empty output_folder {
				create_output_path(); //todo nxt_qid_str(QueryT q) as the filename
				//} create empty output_folder
				cur_q = q;// todo remove cur_q
				cur_qid = nxt_qid;
				//init active vertices
				//todo worker do not have task_spawn
				//todo add to vector
				task_spawn(); //todo add QueryT & q and nxt_qid as input //!!!!using "&"
				//----
				nxt_qid++; //todo change to size_t
			}
		}
		//add from vector to deque
		return true;
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
            usleep(WAIT_TIME_WHEN_IDLE);
            //recieve new queries
            if(!server_exit)
            	if(!update_tasks())
            		server_exit = true;

			Qbig_mtx.lock();
			if (!Qbig->empty())
			{
				// Case 1: there are big tasks to process, wake up threads
				Qbig_mtx.unlock();
				mtx_go.lock();
				ready_go = true;
				// Release threads to compute tasks
				cv_go.notify_all();
				mtx_go.unlock();
			}
			else
			{
				Qbig_mtx.unlock();

				Qreg_mtx.lock();
				if (!Qreg->empty())
				{
					// Case 1: there are reg tasks to process, wake up threads
					Qreg_mtx.unlock();
					mtx_go.lock();
					ready_go = true;
					// Release threads to compute tasks
					cv_go.notify_all();
					mtx_go.unlock();
				}
				else
				{
					Qreg_mtx.unlock();

					mtx_go.lock();
					if (global_num_idle == num_compers && server_exit)
					{
						// Case 2: every thread is waiting, guaranteed since mtx_go is locked
						// Since we are in else-branch, Qreg must be empty
						cout << "global_num_idle: " << global_num_idle << endl;
						global_end_label = true;
						ready_go = true;
						// Release threads
						cv_go.notify_all();
					}
					// Case 3: else, some threads are still processing tasks, check in next round
					mtx_go.unlock();
				}
			}


        }
    }
};

#endif
