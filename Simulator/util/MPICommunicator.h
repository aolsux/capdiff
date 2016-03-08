/*
 * MPICommunicator.h
 *
 *  Created on: Aug 24, 2009
 *      Author: gstu0908
 */

#ifndef MPICOMMUNICATOR_H_
#define MPICOMMUNICATOR_H_

//#include "boost/mpi.hpp"

using namespace std;
class MPICommunicator {
public:
	static int rank;
	static int size;
	static bool isMainRank;

	MPICommunicator(int argc, char* argv[]) {
		//MPI_Init(&argc, &argv);
		//MPI_Comm_size(MPI_COMM_WORLD, &size);
		//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		isMainRank = true;//(rank == 0);
		rank = 0;
		size = 1;
	}

	virtual ~MPICommunicator() {
		//MPI_Finalize();
	}
};

#endif /* MPICOMMUNICATOR_H_ */
