#include <list>
#include <vector>
#include <ilcplex/cplex.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <string.h> 
#include <algorithm>
#include <set>


using namespace std;

#define MAX_INTEGER INT_MAX
#define MAX_DOUBLE CPX_INFBOUND
#define S_VALUE_HOLZMANN 1


typedef  vector<double>  point;

struct solution { //Save the solutions (x,z) where z = f(x) and the time when obtaining the solution
	point *x;
	point *z;
	double time_point;
};

struct Box {
	int box_id = -1;					//Box id 
	point *lb = NULL;					//Lower bound of the box
	point *ub = NULL;					//Upper bound of the box
	double value = -1;					//Value of the box
	int pos_sibling = -1;				//Position of the box refering to its parent
	Box* parent;
};

struct Counters {
	double TSolver = 0.0;				//Total time used by CPLEX
	double TSearching_reachable = 0.0;	//Total time used searching reachable boxes by a point
	double TSplit = 0.0;				//Total time used splitting boxes
	double TAddList = 0.0;				//Total time used adding new points to the Pareto front
	double TSearching_next_box = 0.0;	//Total time used searching next box to analyze
	int Remainingboxestoexplore = 0;	//Number of boxes to explore when algorithm finishes
	int nosolution = 0;					//Number of empty boxes
	double total_time = 0;				//Total time used by the algorithm
};

struct Input {
	//main archive1 % archive2 % tmax % pointlimit % partition % parameterization_model % box_value % set_of_boxes % filtering % cpx_param_parallel % cpx_param_threads (11 input parameters)
	const char *path_obj_file;
	string obj_file;
	const char *path_lp_model;
	string lp_model;
	double tmax;
	int pointlimit;
	string partition;
	string parameterization_model;
	string box_value;
	string set_of_boxes;
	string filtering;
	int cpx_param_parallel;
	int cpx_param_threads;
};

struct MOILP { //Multi-objective integer linear problem
	CPXENVptr *env;				//Pointer to CPLEX environtment
	CPXLPptr *lp;				//Pointer to CPLEX problem
	list<solution> *PF;			//Pareto front
	int dimension;				//Number of objectives
	double **F;					//Matrix of objective coefficients
	point Ideal;				//Ideal point
	point BoundforNadir;		//Bound for Nadir point
	int n_var;					//Number of variables in the original model
	int n_const;				//Number of constraints in the original model
	int n_iterations;			//CPLEX iterations
	double hypervolume;			//Total hypervolume
	double spread;				//Total spread
	double max_time;			//Max defined-time for the user
	int pointlimit;				//Limit size of PF(Pareto Front)
	Counters Count;				//Counting the times
	
	int original_problem_type;	//-1 (MAX)  or 1 (MIN)
	std::string name;			//Problem name
};

#ifdef _WIN32

	#ifdef _WIN64
		#include <direct.h>
		//WINDOWS 64 bits
	#elif
		//Exlusive windows 32 bits
	#endif

class TIEMPO {

public:
	clock_t  t, t2;
	int tipo;

public:
	TIEMPO() {	
		tipo = 0;
	} //Constructor
	
	void TIEMPO::init() { 
		t = clock();
		t2 = 0;
	}

	void TIEMPO::acum() { //Acummulated time 
		clock_t t3 = clock();
		t2 += clock() - t;
		t = t3;
	}

	double TIEMPO::value() { //Transform into seconds
		return ((double)(t2 / double(CLOCKS_PER_SEC)));
	}
	int gettype() {
		return tipo;
	}
	void chgtype(int newt) {
		tipo = newt;
	}

	~TIEMPO() { } //Destructor
};

#elif __linux__
#include <sys/stat.h> //Para mkdir en Linux

class TIEMPO {
	timespec t1, t2;
	int tipo;
	double time;

private:
	timespec diff(timespec start, timespec end) {
		timespec temp;
		if ((end.tv_nsec - start.tv_nsec)<0) {
			temp.tv_sec = end.tv_sec - start.tv_sec - 1;
			temp.tv_nsec = 1000000000.0 + end.tv_nsec - start.tv_nsec;
		}
		else {
			temp.tv_sec = end.tv_sec - start.tv_sec;
			temp.tv_nsec = end.tv_nsec - start.tv_nsec;
		}
		return temp;
	}

	timespec add(timespec accum, timespec now) {
		timespec temp;
		if (accum.tv_nsec + now.tv_nsec >= 1E9) {
			temp.tv_sec = accum.tv_sec + now.tv_sec + 1;
			temp.tv_nsec = accum.tv_nsec + now.tv_nsec - 1E9;
		}
		else {
			temp.tv_sec = accum.tv_sec + now.tv_sec;
			temp.tv_nsec = accum.tv_nsec + now.tv_nsec;
		}
		return temp;
	}


public:

	TIEMPO() {	//Create object TIEMPO and initializes to one of the four types
		tipo = CLOCK_REALTIME;
		//tipo = CLOCK_MONOTONIC;
		//tipo = CLOCK_PROCESS_CPUTIME_ID;
		//tipo = CLOCK_THREAD_CPUTIME_ID;
	} //Constructor


	void init() { //Restart
		time = 0.0;

		clock_gettime(tipo, &t1);

		if (tipo == CLOCK_THREAD_CPUTIME_ID) {
			t1.tv_sec = 0.0;
			t1.tv_nsec = 0.0;
		}
	}

	void acum() { //Acummulate 
		timespec aux;

		clock_gettime(tipo, &t2);

		if (tipo == CLOCK_THREAD_CPUTIME_ID) {
			aux = diff(t1, t2);
			t1 = add(t1, aux);

			time += t1.tv_sec + t1.tv_nsec / 1E9;
		}
		else if (tipo == CLOCK_REALTIME) {
			time += (double)(1.0*(1.0*t2.tv_nsec - t1.tv_nsec*1.0)*1e-9 + 1.0*t2.tv_sec - 1.0*t1.tv_sec);
		}
		t1 = t2;
	}


	double value() { //Transform into seconds
		return (time);
	}

	int gettype() {
		return tipo;
	}

	void chgtype(int newt) {
		//0 = CLOCK_REALTIME
		//1 = CLOCK_MONOTONIC
		//2 = CLOCK_PROCESS_CPUTIME_ID
		//3 = CLOCK_THREAD_CPUTIME_ID
		tipo = newt;
	}

	~TIEMPO() { } //Destructor
};

#elif __APPLE__
//APPLE
#elif __unix__
//Unix
#elif defined(_POSIX_VERSION)
// POSIX
#else
#   error "Unknown compiler"
#endif

