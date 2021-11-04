//########################################################################
//## Copyright 2018 Da Yan http://www.cs.uab.edu/yanda
//##
//## Licensed under the Apache License, Version 2.0 (the "License");
//## you may not use this file except in compliance with the License.
//## You may obtain a copy of the License at
//##
//## //http://www.apache.org/licenses/LICENSE-2.0
//##
//## Unless required by applicable law or agreed to in writing, software
//## distributed under the License is distributed on an "AS IS" BASIS,
//## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//## See the License for the specific language governing permissions and
//## limitations under the License.
//########################################################################
#ifndef TASKPROGMAP_H_
#define TASKPROGMAP_H_

#include "conmap.h"
#include <assert.h>

using namespace std;

//====== task progress object ======
typedef long long TaskID;

struct task_prog
{
	mutex prog_lock;
	atomic<int> t_counter; // could be updated by 2 different children
	TaskID tid; //task id
	TaskID parent_tid;
	int query_id;
	bool completed;

	task_prog(TaskID tid_, TaskID parent_tid_ ,int query_id_){
		t_counter = 0;
		completed = false;
		tid = tid_;
		parent_tid = parent_tid_;
		query_id = query_id_;
	}
};


//====== progress progress table ======

class TaskProgMap
{
public:

    //internal conmap for cached objects
	typedef conmap<TaskID, task_prog*> ProgMap;
	typedef ProgMap::bucket bucket;
	typedef ProgMap::bucket_map bucket_map;
	ProgMap prog_map;
    
    ~TaskProgMap() //The prog should be deleted when corresponding task is finished
    {
    	for(int i=0; i<CONMAP_BUCKET_NUM; i++)
    	{
    		bucket & bucket = prog_map.pos(i);//todo auto
    		bucket_map & kvmap = bucket.get_map();
    		bucket.lock();
    		assert(kvmap.empty()); //todo remove assert when release the code???
//    		for(auto it = kvmap.begin(); it != kvmap.end(); it++)
//    		{
//    			delete it->second; //release task_prog
//    		}
    		bucket.unlock();
    	}
    }
    

    task_prog * get(TaskID & key)
	{
    	bucket & bucket = prog_map.get_bucket(key);
    	bucket_map & kvmap = bucket.get_map();
    	bucket.lock();
		auto it = kvmap.find(key);
		assert(it != kvmap.end()); //todo remove assert when release the code
		task_prog* prog = it->second;
		bucket.unlock();
		return prog;
	}


    void insert(TaskID key, task_prog* value)
	{
    	bucket & bucket = prog_map.get_bucket(key);
		bucket.lock();
		bool inserted = bucket.insert(key, value);
		assert(inserted); //todo remove assert when release the code
		bucket.unlock();
	}
    

    void erase(TaskID key){
    	bucket & bucket = prog_map.get_bucket(key);
		bucket.lock();
		bool erased = bucket.erase(key);
		assert(erased); //todo remove assert when release the code
    	bucket.unlock();
    }


};

#endif
