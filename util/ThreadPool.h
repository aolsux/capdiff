/*
 * ThreadPool.h
 *
 *  Created on: 01.12.2010
 *      Author: Martin Rueckl
 */

#ifndef THREADPOOL_H_
#define THREADPOOL_H_

struct Task;

class ThreadPool
{
	public:
		ThreadPool();
		virtual ~ThreadPool();
		void schedule(const Task &task);
};

#endif /* THREADPOOL_H_ */
