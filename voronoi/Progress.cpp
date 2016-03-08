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

ProgressMonitor* ProgressMonitor::instance = NULL;
boost::mutex ProgressMonitor::mutex;

ProgressMonitor::ProgressMonitor() {
	thread = new boost::thread(boost::bind(&ProgressMonitor::run, this));
	eraseLines = 0;
	if (instance != NULL)
		cerr << "Singleton invalidity!!!" << endl;
}

void ProgressMonitor::addProgressBar(const ProgressBar &pbar) {
	boost::mutex::scoped_lock l(mutex);
	if (instance == NULL)
		instance = new ProgressMonitor();
	instance->ProgressBars.push_back(&pbar);
}

void ProgressMonitor::removeProgressBar(const ProgressBar &pbar) {
	boost::mutex::scoped_lock l(mutex);
	for (unsigned u = 0; u < instance->ProgressBars.size(); u++)
		if (instance->ProgressBars[u] == &pbar) {
			instance->ProgressBars.erase(instance->ProgressBars.begin() + u);
			return;
		}
}

void ProgressMonitor::run() {
	cout << clear_to_eos;
	while (active()) {
		cout << show_cursor( FALSE );
		boost::mutex::scoped_lock l(mutex);
		unsigned maxLength = 0;
		for (unsigned u = 0; u < ProgressBars.size(); u++)
			maxLength = maxLength < ProgressBars[u]->title.size() ? ProgressBars[u]->title.size() : maxLength;
		cout << cursor_up(eraseLines);
		for (unsigned u = 0; u < ProgressBars.size(); u++) {
			cout << setw(maxLength + 1) << setfill(' ') << left << ProgressBars[u]->title << ": " << progressString(
					ProgressBars[u]) << endl;
		}
		eraseLines = ProgressBars.size();
		boost::posix_time::seconds workTime(1);
		boost::this_thread::sleep(workTime);
	}
	cout << finalize;
}

bool ProgressMonitor::active() {
	boost::mutex::scoped_lock l(mutex);
	for (unsigned u = 0; u < ProgressBars.size(); u++)
		if (!ProgressBars[u]->finished())
			return true;
	return false;
}

string ProgressMonitor::progressString(const ProgressBar * pbar) {
	stringstream ss;
	ss << setfill(' ') << setw(10) << pbar->done << "/" << setfill(' ') << setw(10) << pbar->work << "("
			<< setfill(' ') << setw(6) << setprecision(4) << (double) pbar->done * 100 / (double) pbar->work << "%) ";
	for (double u = 0; u <= (double) pbar->done / (pbar->percent * 2) - 1; u++)
		ss << '#';
	return ss.str();
}

