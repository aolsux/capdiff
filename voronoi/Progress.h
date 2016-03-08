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
	ProgressMonitor();//singleton
	ProgressMonitor(const ProgressMonitor &ref);//singleton->noncopyable
	ProgressMonitor& operator=(const ProgressMonitor &ref);//singleton->noncopyable
public:
	static void addProgressBar(const ProgressBar &pbar);
	static void removeProgressBar(const ProgressBar &pbar);

private:
	static ProgressMonitor* instance;
	static boost::mutex mutex;

	void run();

	bool active();

	string progressString(const ProgressBar * pbar);

	vector<const ProgressBar *> ProgressBars;
	boost::thread *thread;
	unsigned eraseLines;
};

class ProgressBar {
public:

	ProgressBar(string text, unsigned _work) :
		title(text), work(_work), done(0), percent(max(1., (double) work / 100)) {
		ProgressMonitor::addProgressBar(*this);
	}

	~ProgressBar(){
		ProgressMonitor::removeProgressBar(*this);
	}

	void progress(unsigned s = 1) {
		done += s;
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

	string title;
	unsigned work, done, pdone;
	double percent;

};
#endif /* PROGRESSBAR_H_ */
