/*
 * Progress.cpp
 *
 *  Created on: 16.11.2010
 *      Author: Martin Rueckl
 */
#include "Progress.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <exception>

#include <boost/thread.hpp>
#include "vt100.h"

ProgressMonitor::ProgressMonitor():
	workTime	(1),
	active		(true),
	eraseLines	(0),
	maxLength	(0)
{
	thread = new boost::thread(boost::bind(&ProgressMonitor::run, this));
}

ProgressMonitor::~ProgressMonitor()
{
	active = false;
	thread->join();
	{
		boost::mutex::scoped_lock l(mutex);
		for (unsigned u = 0; u < ProgressBars.size(); u++)
			delete ProgressBars[u];
		ProgressBars.clear();
	}
	delete thread;
}

void ProgressMonitor::addProgressBar(ProgressBar *pbar)
{
	boost::mutex::scoped_lock l(mutex);
	ProgressBars.push_back(pbar);
	maxLength = maxLength < pbar->title.size() ? pbar->title.size() : maxLength;
}

void ProgressMonitor::removeProgressBar(ProgressBar *pbar)
{
//	boost::mutex::scoped_lock l(mutex);
//	for (unsigned u = 0; u < ProgressBars.size(); u++)
//		if (ProgressBars[u] == &pbar) {
//			ProgressBars.erase(ProgressBars.begin() + u);
//			break;
//		}
}

void ProgressMonitor::run()
{
	cout << show_cursor( FALSE );
	while (active) {
		{
			boost::mutex::scoped_lock l(mutex);
			if(ProgressBars.size()>0)
			{
				cout << clear_to_eos;
				for (unsigned u = 0; u < ProgressBars.size(); u++) {
					cout << setw(maxLength + 1) << setfill(' ') << left << ProgressBars[u]->title << ": " << progressString(ProgressBars[u]) << endl;
				}
				removeFinished();
				cout << cursor_up(ProgressBars.size());
			}
		}
		boost::this_thread::sleep(workTime);
	}
	cout << finalize;
}

unsigned ProgressMonitor::removeFinished()
{
	unsigned elems	 = ProgressBars.size();
	unsigned removes = 0;
	if(elems < 15) return 0;
	for(vector<ProgressBar *>::iterator it = ProgressBars.begin(); it < ProgressBars.end(); it++)
	{
		ProgressBar *pbar = *it;
		if(pbar->finished() && removes + 15 < elems){
			delete pbar;
			ProgressBars.erase(it);
			removes++;
		}
	}
	return removes;
}

string ProgressMonitor::progressString(const ProgressBar * pbar)
{
	stringstream ss;
	ss << setfill(' ') << setw(10) << pbar->done << "/" << setfill(' ') << setw(10) << pbar->work << "(" << setfill(' ') << setw(6) << setprecision(4)
			<< (double) pbar->done * 100 / (double) pbar->work << "%) ";
	for (double u = 0; u <= (double) pbar->done / (pbar->percent * 2.) - 1.; u++)
		ss << '#';
	return ss.str();
}

