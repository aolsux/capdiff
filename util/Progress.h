/*
 * ProgressBar.h
 *
 *  Created on: 14.04.2010
 *      Author: Enigma
 */

#ifndef PROGRESSBAR_H_
#define PROGRESSBAR_H_

/*
 *
 */
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

#include <boost/thread.hpp>

using namespace std;
class ProgressBar;

class ProgressMonitor {
private:
	ProgressMonitor(const ProgressMonitor &ref);//noncopyable
	ProgressMonitor& operator=(const ProgressMonitor &ref);//noncopyable

public:
	ProgressMonitor();
	virtual ~ProgressMonitor();
	void addProgressBar(ProgressBar *pbar);
	void removeProgressBar(ProgressBar *pbar);

private:


	void run();
	string progressString(const ProgressBar * pbar);
	unsigned removeFinished();

	//	bool active();
	boost::mutex mutex;
	vector<ProgressBar *> ProgressBars;
	boost::thread *thread;
	boost::posix_time::seconds workTime;
	bool active;
	unsigned eraseLines;
	unsigned maxLength;
};

class ProgressBar {
public:

	ProgressBar(string text, unsigned _work) :
		title(text), work(_work), done(0), percent(max(1., (double) work / 100)) {
	}

	~ProgressBar() {
	}

	void progress(unsigned s = 1) {
		done += s;
	}

	void set(unsigned s) {
		done = s;
	}

	ProgressBar &operator+(unsigned s) {
		done += s;
		return *this;
	}

	ProgressBar &operator++() {
		done++;
		return *this;
	}

	string getTitle() {
		return title;
	}

	bool finished() const {
		return done == work;
	}

	void finish() {
		done = work;
	}

	string title;
	unsigned work, done;
	double percent;

};
#endif /* PROGRESSBAR_H_ */
