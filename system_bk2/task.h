

#ifndef TASK_H_
#define TASK_H_
#include "TaskProgMap.h"

using namespace std;

template <class ContextT>
class Task
{
public:
	typedef ContextT ContextType; // used in comper.h
	ContextT context;
	task_prog* prog;

	TaskID get_id(){
		return prog->tid;
	}

	TaskID get_qid(){
		return prog->query_id;
	}

	//todo (delete) prog initialization is moved to add_task()
//	Task(TaskID tid, TaskID parent_tid, int query_id){
//		prog = new task_prog(tid, parent_tid, query_id);
//	}
//
//	Task(TaskID tid, int query_id)
//	{
//		prog = new task_prog(tid, -1, query_id);
//	}
//
//    virtual ~Task()
//    {
//        if(prog != NULL)
//        	delete prog;
//    }

	friend ibinstream &operator<<(ibinstream &m, const Task &t)
	{
		m << t.context;
		m << t.prog->tid;
		return m;
	}

	friend obinstream &operator>>(obinstream &m, Task &t)
	{
		m >> t.context;
		TaskID task_id;
		m >> task_id;
		t.prog = global_prog_map->get(task_id);
		return m;
	}

	friend ifbinstream &operator<<(ifbinstream &m, const Task &t)
	{
		m << t.context;
		m << t.prog->tid;
		return m;
	}

	friend ofbinstream &operator>>(ofbinstream &m, Task &t)
	{
		m >> t.context;
		TaskID task_id;
		m >> task_id;
		t.prog = global_prog_map->get(task_id);
		return m;
	}
};

#endif /* TASK_H_ */
