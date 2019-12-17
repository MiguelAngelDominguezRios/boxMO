//Algorithm which calculates the Pareto front of an MOCO (multiobjective combinatorial optimization) problem
//
//Author: Miguel Angel Dominguez Rios
//Date of last modification : 09/12/2019


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
#include "boxes2.h"


//GLOBAL VARIABLES
int box_id = 0;
int boxes_destroyed = 0;
int p = 0;	//Dimension of the objective space


struct Container {
	Box *VBox;
	int position;
};

class Compare_value {
public:
	bool operator()(Box *x, Box *y) const {
		return x->value > y->value;
	}
};

void get_name_string(string *name) {
	//Take the name of an archive, ommiting the path
	//Example: Input   ./KP/data/hello.lp		Output	hello.lp

	std::size_t pos;

	pos = name->find("/");
	while (pos != -1) {
		*name = name->substr(pos + 1);
		pos = name->find("/");
	}
}

void error_input() {
	printf("\nTo execute BOXES write ./BOXES (arg1) (arg2) (arg3) (arg4) (arg5) (arg6) (arg7) (arg8) (arg9) (arg10) (arg11)");
	printf("\n(arg1) is the objective costs file");
	printf("\n(arg2) is the .lp file");
	printf("\n(arg3) is the maximum total execution time in seconds (0 for unlimited time)");
	printf("\n(arg4) is the maximum size of the Pareto front (0 for unlimited size)");
	printf("\n(arg5) is the type of partition. Type 1(full), 2(p-partition)");
	printf("\n(arg6) is the parameterization model. Type 1(chalmet), 2(tchebycheff) or 3(benson)");
	printf("\n(arg7) is the box value. Type 11(volume), 12(scaled), 13(reduced), 14(reduced_scaled)");
	printf("\n(arg8) is the number of set of boxes used. Type an integer positive number, or (inf) or (alternate).");
	printf("\n(arg9) is the filtering proccess. Type 1(RE).");
	printf("\n(arg10) is the value of CPX_PARAM_PARALLEL. Type -1, 0 or 1");
	printf("\n(arg11) is the value of CPX_PARAM_THREADS. Type 0 or a(positive integer)\n");
	exit(EXIT_FAILURE);
}

void error_open_file(const char *name) {
	printf("\nImpossible to open file %s\n", name);
	exit(EXIT_FAILURE);
}

void Input_control(int argc, char **argv, Input *input) {
	//Controlling the input parameters. In case of error, exit.
	FILE *fp;

	if (argc != 12) { error_input(); }

	//Data archive obj_file
	const char *arch1 = argv[1];
	fp = fopen(arch1, "r");
	if (fp == NULL)
		error_open_file(arch1);
	else {
		fclose(fp);
		string a1 = argv[1];
		input->path_obj_file = argv[1];
		get_name_string(&a1);
		input->obj_file = a1;
	}

	//Data archive model_file
	const char *arch2 = argv[2];
	fp = fopen(arch2, "r");
	if (fp == NULL)
		error_open_file(arch2);
	else {
		fclose(fp);
		string a2 = argv[2];
		input->path_lp_model = argv[2];
		get_name_string(&a2);
		input->lp_model = a2;
	}

	//Maximum execution time
	std::string time = argv[3];
	double tmax = atof(argv[3]);
	if ((time != "0") && (time != "0.0")) {
		if (tmax <= 0.0)
			error_input();
	}
	input->tmax = tmax;

	//Maximum size of Pareto front
	std::string pointlimit = argv[4];
	int maxsize = atoi(argv[4]);
	if (pointlimit != "0") {
		if (maxsize <= 0)
			error_input();
	}
	input->pointlimit = maxsize;

	//Partition
	std::string type = argv[5];
	if ((type == "1") || (type == "full"))				input->partition = "full";
	else if ((type == "2") || (type == "p-partition")) input->partition = "p-partition";
	else error_input();

	//Problem model
	std::string model = argv[6];
	if ((model == "1") || (model == "chalmet"))		input->parameterization_model = "chalmet";
	else if ((model == "2") || (model == "tchebycheff"))	input->parameterization_model = "tchebycheff";
	else if ((model == "3") || (model == "benson"))		input->parameterization_model = "benson";
	else error_input();

	//Box value
	std::string box = argv[7];
	if ((box == "11") || (box == "volume"))			input->box_value = "volume";
	else if ((box == "12") || (box == "scaled"))			input->box_value = "scaled";
	else if ((box == "13") || (box == "reduced"))			input->box_value = "reduced";
	else if ((box == "14") || (box == "reduced_scaled"))	input->box_value = "reduced_scaled";
	else error_input();

	//Set of boxes
	std::string set_b = argv[8];
	if (set_b == "inf")		input->set_of_boxes = "inf";
	else if (set_b == "alternate")	input->set_of_boxes = "alternate";
	else {
		input->set_of_boxes = argv[8];
		int set_b = atoi(argv[8]);
		if (set_b < 1) error_input();
	}

	//Filtering process
	std::string filtering = argv[9];
	if ((filtering == "1") || (filtering == "RE"))			input->filtering = "RE";
	else { error_input(); }

	//CPX_PARAM_PARALLEL
	std::string sparallel = argv[10];
	if ((sparallel != "-1") && (sparallel != "0") && (sparallel != "1")) error_input();
	input->cpx_param_parallel = atoi(argv[10]);

	//CPX_PARAM_THREADS
	std::string sthreads = argv[11];
	if (sthreads != "0") {
		int x = atoi(argv[11]);
		if ((x <= 0) || (x > 32)) error_input();
	}
	input->cpx_param_threads = atoi(argv[11]);

	//Incompatibilities
	if ((input->partition == "p-partition") && (input->filtering == "RA")) error_input();
}

void error_CPLEX_status(int code) {
	printf("\nCPLEX error status = %d", code);
	exit(EXIT_FAILURE);
}

void Create_CPLEX_object_and_set_CPLEX_parameters(MOILP *P, string name, int parallel, int threads) {
	CPXENVptr *env = new(CPXENVptr); //Create CPLEX environment
	CPXLPptr *lp = new(CPXLPptr);	//Create associated problem to CPLEX environment
	int status = 0;
	const char *namec = name.c_str();

	*env = CPXopenCPLEX(&status);
	*lp = CPXcreateprob(*env, &status, namec);
	status = CPXsetintparam(*env, CPX_PARAM_PARALLELMODE, parallel); //-1 CPX_PARALLEL_OPPORTUNISTIC;   0 (default) CPX_PARALLEL_AUTO  1 CPX_PARALLEL_DETERMINISTIC
	status = CPXsetintparam(*env, CPX_PARAM_THREADS, threads);	//0 (default, CPLEX decides) ;	1 (sequential, single Thread) ;		N (uses up to N threads)
	status = CPXsetdblparam(*env, CPX_PARAM_EPGAP, 0.0); //Default 1e-04. Value between 0 and 1. Sets a relative tolerance on the gap between the best integer objective and the objective of the best node remaining
	status = CPXsetdblparam(*env, CPX_PARAM_EPAGAP, 0.0); //Default 1e-06
	status = CPXsetdblparam(*env, CPX_PARAM_EPINT, 0.0); //Default 1e-05. CPLEX tolerance for integer variables
	P->env = env;
	P->lp = lp;
}

int compare_vectors(vector<double> *v1, vector<double> *v2, int p) {
	//Compare two p-dimensional vectors v1 and v2
	//Return 1 if v1 <= v2 
	//Return 0 otherwise

	for (int index = 0; index < p; index++) {
		if (v1->at(index) > v2->at(index)) {
			return 0;
		}
	}
	return 1;
}

void Set_MOILP_structure(MOILP *P, Input *input) {
	//Read .lp CPLEX model and objectives archive. If the problem is MAX, convert it into MIN

	FILE *fp;
	int i, j, n, m, status, original_sense, pointlimit;
	double max_time, bound;
	string name;
	list<solution> *PF = new(list<solution>);

	//Read CPLEX problem
	status = CPXreadcopyprob(*P->env, *P->lp, input->path_lp_model, NULL);
	if (status != 0) error_CPLEX_status(status);

	fp = fopen(input->path_obj_file, "r");
	if (fp == NULL) error_open_file(input->path_obj_file);

	fscanf(fp, "%d", &p); //Read number of objectives
	fscanf(fp, "%d", &n); //Read number of variables
	fscanf(fp, "%d", &m); //Read number of constraints

	//Creating the objective costs matrix
	double **F = (double **)malloc(p * sizeof(double *));
	for (i = 0; i < p; i++) {
		F[i] = (double *)malloc(n * sizeof(double));
		for (j = 0; j < n; j++) {
			fscanf(fp, "%lf", &F[i][j]);	//Read costs for every objective
		}
	}
	fclose(fp);

	//MAX problems are transformed into MIN problems
	original_sense = CPXgetobjsen(*P->env, *P->lp);
	if (original_sense == CPX_MAX) {
		for (i = 0; i < p; i++) {
			for (j = 0; j < n; j++) {
				F[i][j] = -F[i][j];
			}
		}
		CPXchgobjsen(*P->env, *P->lp, CPX_MIN);
	}

	//Max time and point limit
	max_time = input->tmax;
	if (input->tmax == 0) max_time = MAX_DOUBLE;
	pointlimit = input->pointlimit;
	if (input->pointlimit == 0) pointlimit = MAX_INTEGER;

	//Bound for Nadir and Ideal point	
	bound = (original_sense == CPX_MIN) ? -CPX_INFBOUND : CPX_INFBOUND;
	P->Ideal.resize(p);
	P->BoundforNadir.resize(p);
	for (i = 0; i < p; i++) {
		P->Ideal.at(i) = bound;
		P->BoundforNadir.at(i) = -bound;
	}
	name = input->lp_model;
	name.resize(name.length() - 3); //Eliminate  ".lp" in the name 

	//Initialize MOILP structure
	//P->env, P->lp, P->Ideal, P->BoundforNadir were previously assigned
	P->PF = PF;
	P->F = F;
	P->dimension = p;
	P->n_var = n;
	P->n_const = m;
	P->n_iterations = 0;
	P->hypervolume = 0.0;
	P->spread = -1.0;
	P->max_time = max_time;
	P->pointlimit = pointlimit;
	P->original_problem_type = original_sense; // CPXgetobjsen(*P->env, *P->lp);
	P->name = name;
}

void Solve_z(MOILP *P, int *stat, double *z) {
	//Solve a MIP problem, returning the stat and objective value
	CPXmipopt(*P->env, *P->lp);
	CPXgetobjval(*P->env, *P->lp, z);
	*stat = CPXgetstat(*P->env, *P->lp);
}

bool Calculate_problem_bounds(MOILP *P, int *indices, double *scaling) {
	//Calculate Ideal point and a bound for Nadir point of the problem. We always consider MIN problems.
	//In vector scaling we save the range (bound_of_nadir[i] - ideal[i]) for every component

	TIEMPO t;
	t.init();
	int i, stat;
	double obj, t_ref;
	double tmax = P->max_time;
	int cur_objsense = CPXgetobjsen(*P->env, *P->lp); //CPX_MIN 1	//CPX_MAX -1

	printf("\nCalculating bounds...");

	//Ideal point
	for (i = 0; i < p; i++) {
		t.acum();	t_ref = tmax - t.value();
		if (t_ref < 0) t_ref = 0.0;
		CPXsetdblparam(*P->env, CPX_PARAM_TILIM, t_ref);
		CPXchgobj(*P->env, *P->lp, P->n_var, indices, P->F[i]);
		Solve_z(P, &stat, &obj);
		P->Ideal.at(i) = obj;
		if (stat != CPXMIP_OPTIMAL) {
			printf("\nNON OPTIMAL SOLUTION IS FOUND");
			P->Count.total_time = P->max_time;
			return false;
		}
	}

	//Upper bound for Nadir point 
	CPXchgobjsen(*P->env, *P->lp, -cur_objsense);
	for (i = 0; i < p; i++) {
		t.acum();	t_ref = tmax - t.value();
		if (t_ref < 0) t_ref = 0.0;
		CPXsetdblparam(*P->env, CPX_PARAM_TILIM, t_ref);
		CPXchgobj(*P->env, *P->lp, P->n_var, indices, P->F[i]);
		Solve_z(P, &stat, &obj);
		P->BoundforNadir.at(i) = obj + 1;
		if (stat != CPXMIP_OPTIMAL) {
			printf("\nNON OPTIMAL SOLUTION IS FOUND");
			P->Count.total_time = P->max_time;
			return false;
		}
	}

	//Restore
	CPXchgobjsen(*P->env, *P->lp, cur_objsense);

	//Calculating ranges for every component
	for (i = 0; i < P->dimension; i++) {
		scaling[i] = P->BoundforNadir.at(i) - P->Ideal.at(i);
	}
	t.acum();
	printf("DONE after %lf seconds", t.value());
	return true;
}

Box *Create_box(point *LB, point *UB) {
	//Given the opposite extreme points for the box, we create the corresponding structure (except its value)

	Box *B = new(Box);

	//Create a pointer to a point with dimension dim
	point *lb = new(point);	lb->resize(p);
	point *ub = new(point);	ub->resize(p);

	for (int i = 0; i < p; i++) {
		lb->at(i) = LB->at(i);
		ub->at(i) = UB->at(i);
	}

	//Creating the Box
	B->box_id = box_id++;
	B->lb = lb;
	B->ub = ub;
	
	return B;
}

void Destroy_box(Box *B) {
	//Destroy box and free allocated memory
	boxes_destroyed++;
	free(B->lb);
	free(B->ub);
	free(B);
}

double get_volume(point *lb, point *ub) {
	//Get the volume of the box
	double volume = 1.0;
	for (int i = 0; i < p; i++) {
		volume = volume * (ub->at(i) - lb->at(i));
	}
	return volume;
}

double get_scaled_volume(point *lb, point *ub, double *scaling) {
	//Get scaled volume of the box
	double volume = 1.0;
	for (int i = 0; i < p; i++) {
		volume = volume * ((ub->at(i) - lb->at(i)) / scaling[i]);
	}
	return volume;
}

void get_value(Box *B, string value, double *scaling) {
	//Get the value of the box
	if ((value == "volume") || (value == "reduced"))				B->value = get_volume(B->lb, B->ub);
	else if ((value == "scaled") || (value == "reduced_scaled"))	B->value = get_scaled_volume(B->lb, B->ub, scaling);
}

Box *CREATE_BOX_AND_VALUE(point *LB, point *UB, point *z, Input *input, double *scaling, int coordinate) {
	string X = input->box_value;
	Box *B = Create_box(LB, UB);

	if (B != NULL) {
		get_value(B, X, scaling);
		if (input->partition == "full") {
			if (X == "reduced")					B->value -= get_volume(LB, z);
			else if (X == "reduced_scaled")	B->value -= get_scaled_volume(LB, z, scaling);
		}
	}
	return B;
}

void CREATE_INITIAL_BOX(Box **B0, point *LB, point *UB, Input *input, double *scaling) {
	//The initial box is never empty, because it is suppossed that there is at least one non-dominated point
	
	//Create the box
	*B0 = Create_box(LB, UB);

	//Assign the value of the box
	get_value(*B0, input->box_value, scaling);		//Assign the value of the box

	//Initial box is considered as the first sibling. In fact, it has no parent, but we need this
	(*B0)->pos_sibling = 0;
}

void INSERT_INTO_LIST_2(MOILP *P, vector<multiset<Box *, Compare_value>> *L, Box *B, int position) {
	//We insert in L[position] the new box. If that element of the vector does not exist, we create it
	//Every element of vector L is a tree, ordered by its value

	TIEMPO t;
	t.init();

	if (position == int(L->size())) { //Create new tree 
		multiset<Box *, Compare_value> new_tree;
		L->push_back(new_tree);
	}
	L->at(position).insert(B);

	t.acum();
	P->Count.TAddList += t.value();
}

void INSERT_INTO_LIST_1(MOILP *P, vector<multiset<Box *, Compare_value>> *L, Box *B, int *RemainingBoxesToExplore) {
	//We insert the new box in the corresponding cell of L, according to its sibling's position

	TIEMPO t;
	t.init();

	L->at(B->pos_sibling).insert(B);
	(*RemainingBoxesToExplore)++;

	t.acum();
	P->Count.TAddList += t.value();
}

void Create_chalmet_model(MOILP *P, int *indice_col) {
	//Given the original model, create Chalmet et al model
	// MIN (f1 + .... + fp)
	// s.t.		x in X
	//			fi(x) <= ki		for i = 1,...,p
	//
	int i, j, matbeg = 0;
	double *new_c = (double *)malloc(P->n_var * sizeof(double));

	//The sense of the new constraints are <= because we have a MIN problem
	const char *sense = "L";

	for (i = 0; i < P->n_var; i++) {
		new_c[i] = 0;
		for (j = 0; j < P->dimension; j++) {
			new_c[i] += P->F[j][i];
		}
	}
	CPXchgobj(*P->env, *P->lp, P->n_var, indice_col, new_c); //min sum(fi)
	for (i = 0; i < P->dimension; i++) {
		CPXaddrows(*P->env, *P->lp, 0, 1, P->n_var, 0, sense, &matbeg, indice_col, P->F[i], NULL, NULL); //new constraint
	}
	free(new_c);
}

void Create_tchebycheff_model(MOILP *P, int *indice_col) {
	//Given the original model, create the Tchebycheff model
	// MIN ((f1 + .... + fp) - alfa)
	// s.t.		x in X
	//			fi(x) <= ki - alfa		for i = 1,...,p
	//			alfa >= 0
	//

	int i, j, matbeg = 0;
	double *new_c = (double *)malloc((P->n_var + 1) * sizeof(double));


	//The sense of the new constraints are <=, since we have a MIN problem
	const char *sense = "L";

	for (i = 0; i < P->n_var; i++) {
		new_c[i] = 0;
		for (j = 0; j < p; j++) {
			new_c[i] += P->F[j][i];
		}
	}

	int index_alfa = P->n_var;
	double coef_alfa = -1.0;
	CPXaddcols(*P->env, *P->lp, 1, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
	CPXchgobj(*P->env, *P->lp, P->n_var, indice_col, new_c); //min (eps*sum(fi))
	CPXchgobj(*P->env, *P->lp, 1, &index_alfa, &coef_alfa); //min (eps*sum(fi) - alfa)

	int numrows = CPXgetnumrows(*P->env, *P->lp);


	for (i = 0; i < p; i++) {
		string rowname = "c" + to_string(P->n_const + i + 1);
		char *rownamechar = new char[rowname.length() + 1];
		strcpy(rownamechar, rowname.c_str());
		CPXaddrows(*P->env, *P->lp, 0, 1, P->n_var, 0, sense, &matbeg, indice_col, P->F[i], NULL, &rownamechar); //new constraint   
		numrows++;
		CPXchgcoef(*P->env, *P->lp, numrows - 1, P->n_var, 1.0);	//fi(x)+alfa <= 0
		free(rownamechar);
	}

	free(new_c);
}

void Create_benson_model(MOILP *P) {
	//Given the original MIN model, create the Benson model
	// MAX (l1 + l2 + ... + lp)
	// s.t.		x in X
	//			fi(x) + li = ui		for i = 1,...,p
	//			li >= 0				for i = 1,...,p
	//
	int i, status, matbeg = 0;
	int *indice_col = (int *)malloc((P->n_var + p) * sizeof(int));
	double *new_c = (double *)malloc((P->n_var + p) * sizeof(double));
	for (i = 0; i < P->n_var; i++) {
		indice_col[i] = i;
		new_c[i] = 0.0;
	}
	for (i = P->n_var; i < P->n_var + p; i++) {
		indice_col[i] = i;
		new_c[i] = 1.0;
	}

	//The sense of the new constraints are "="
	const char *sense = "E";

	status = CPXaddcols(*P->env, *P->lp, p, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);	//Add new variables li
	status = CPXchgobjsen(*P->env, *P->lp, CPX_MAX);
	status = CPXchgobj(*P->env, *P->lp, P->n_var + p, indice_col, new_c); //MAX (l1 + l2 +... +lp)

	int numrows = CPXgetnumrows(*P->env, *P->lp);
	for (i = 0; i < p; i++) {
		status = CPXaddrows(*P->env, *P->lp, 0, 1, P->n_var, 0, sense, &matbeg, indice_col, P->F[i], NULL, NULL); //new constraint   
		numrows++;
		status = CPXchgcoef(*P->env, *P->lp, numrows - 1, P->n_var + i, 1.0);	//fi(x)+ li = 0
	}
	free(new_c);
}

void CREATE_MODEL(MOILP *P, int *indices, Input *input) {
	if (input->parameterization_model == "chalmet") {
		Create_chalmet_model(P, indices);
	}
	else if (input->parameterization_model == "tchebycheff") {
		Create_tchebycheff_model(P, indices);
	}
	else if (input->parameterization_model == "benson") {
		Create_benson_model(P);
	}
}

Box *SELECT_NEXT_BOX_2(MOILP *P, vector<multiset<Box *, Compare_value>> *L) {
	//The next box to select is in the top of L[0] (maximum value in the current list of boxes)
	TIEMPO tt;
	tt.init();

	multiset<Box *>::iterator it = L->at(0).begin();
	tt.acum(); P->Count.TSearching_next_box += tt.value();

	return *it;
}

Box *SELECT_NEXT_BOX_1(MOILP *P, vector<multiset<Box *, Compare_value>> *L, int *counter) {
	//The next box to select varies in every component of L. We use a counter
	TIEMPO tt;
	tt.init();
	int index;

	index = *counter % p;

	while (L->at(index).size() == 0) {	//In case of L(index) been empty
		(*counter)++;
		index = *counter % p;
	}
	multiset<Box *>::iterator it = L->at(index).begin();
	(*counter)++;
	tt.acum(); P->Count.TSearching_next_box += tt.value();

	return *it;
}

bool check_time_and_pointlimit(TIEMPO *t, double maxtime, size_t size_PF, int pointlimit) {
	//Return false if time limit or the desired number of solutions are reached. Return true, otherwise
	t->acum();
	return ((t->value() < maxtime) && (size_PF < pointlimit));
}

void set_time_for_solver(MOILP *P, TIEMPO *t) {
	//Set the maximum time available for CPLEX
	t->acum();
	double tt = P->max_time - t->value();
	if (tt < 0) tt = 0;
	CPXsetdblparam(*P->env, CPX_PARAM_TILIM, tt);
}

void Solve(MOILP *P, int *stat, double *x, double *obj) {
	//Solve a MIP problem, returning the stat, solution, objective value and execution time
	TIEMPO t;
	t.init();
	CPXmipopt(*P->env, *P->lp);
	CPXgetobjval(*P->env, *P->lp, obj);
	*stat = CPXgetx(*P->env, *P->lp, x, 0, P->n_var - 1);
	*stat = CPXgetstat(*P->env, *P->lp);
	P->n_iterations++;
	t.acum();
	P->Count.TSolver += t.value();
}

void Set_constraint_extreme_box(MOILP *P, point *LB, point *UB) {
	//Set the rhs values for the new constraints (fi(x) <= ki) in the current iteration
	//ki = Ni - 1, being Ni the bound of the local nadir point (MIN problems)
	point *N = UB;
	int *indices = (int *)malloc(p * sizeof(int));
	double *values = (double *)malloc(p * sizeof(double));

	for (int i = 0; i < p; i++) {
		indices[i] = P->n_const + i;
		values[i] = N->at(i) - 0.5;	//Because of rounding errors when using CPXPARAMEPINT = 0 , we substract a small value 0<alfa<1
		
	}
	CPXchgrhs(*P->env, *P->lp, p, indices, values); //s.t.  fi <= ki,  i = 1,...,p 
	free(indices); free(values);
}

void Calculate_Image_of_x(double **F, double *x, int n, point *z, int p) {
	//Given a vector-solution x, return the point z = (f1(x), ... , fp(x))
	int i, j;
	for (i = 0; i < p; i++) {
		z->at(i) = 0.0;
		for (j = 0; j < n; j++) {
			z->at(i) += F[i][j] * x[j];
		}
	}
}

void Add_new_PF_point(double *x, int dim_x, point *z, std::list<solution> *PF, double time) {
	//Insert a new efficient solution, its Pareto Front point (x,z), and the time when the solution was obtained
	point *xx = new(point);		xx->resize(dim_x);
	point *zz = new(point);		zz->resize(p);
	solution *A = new(solution);

	for (int i = 0; i < dim_x; i++)		xx->at(i) = x[i];
	for (int i = 0; i < z->size(); i++) zz->at(i) = z->at(i);
	A->x = xx;
	A->z = zz;
	A->time_point = time;
	PF->push_back(*A);
}

void Find_upper_bounds_that_contains_z_2(MOILP *P, vector<list<Container>> *Bj, list<Container> *R, vector<multiset<Box *, Compare_value>> *L, point *z) {
	//In case of option "alternate RE", R will be a subset of L, which contains all boxes with upper bound dominating the point z.
	//Bj , j=1,..,p   are the set of points that have 1 component equals z (weakly dominates)
	TIEMPO t;
	t.init();
	int i, j;
	multiset<Box *>::iterator k;
	int coincidence_axis;
	int number_of_coincidences;
	bool dominates;

	for (i = 0; i < L->size(); i++) { //For every tree in vector L
		k = L->at(i).begin();
		while (k != L->at(i).end()) {
			coincidence_axis = -1;
			number_of_coincidences = 0;
			dominates = true;
			for (j = 0; j < p; j++) {
				if (z->at(j) >(*k)->ub->at(j)) {
					dominates = false;
					j = p;
				}
				else if (z->at(j) == (*k)->ub->at(j)) {
					coincidence_axis = j;
					number_of_coincidences++;
				}
			}
			if (dominates) {
				if (coincidence_axis == -1) { //Strictly dominates
					Container c;
					c.VBox = *k;	c.position = i;
					R->push_back(c);
					k = L->at(i).erase(k);

				}
				else if (number_of_coincidences == 1) {
					Container c;
					c.VBox = *k;	c.position = i;
					Bj->at(coincidence_axis).push_back(c);
					k = L->at(i).erase(k);
				}
				else {
					++k;
				}
			}
			else {
				++k;
			}
		}
	}
	t.acum();
	P->Count.TSearching_reachable += t.value();
}

void Find_upper_bounds_that_strictly_contains_z_2(MOILP *P, list<Container> *R, vector<multiset<Box *, Compare_value>> *L, point *z) {
	//In case of option "alternate RA", R will be a subset of L, which contains all boxes with upper bound dominating the point z.
	TIEMPO t;
	t.init();
	int i, j;
	multiset<Box *>::iterator k;
	int coincidence_axis;
	bool dominates;

	for (i = 0; i < L->size(); i++) { //For every tree in vector L
		k = L->at(i).begin();
		while (k != L->at(i).end()) {
			coincidence_axis = -1;
			dominates = true;
			for (j = 0; j < p; j++) {
				if (z->at(j) >(*k)->ub->at(j)) {
					dominates = false;
					j = p;
				}
				else if (z->at(j) == (*k)->ub->at(j)) {
					coincidence_axis = j;
				}
			}
			if (dominates) {
				if (coincidence_axis == -1) { //Strictly dominates
					Container c;
					c.VBox = *k;	c.position = i;
					R->push_back(c);
					k = L->at(i).erase(k);

				}
				else {
					++k;
				}
			}
			else {
				++k;
			}
		}
	}
	t.acum();
	P->Count.TSearching_reachable += t.value();
}

void Find_upper_bounds_that_contains_z_1(MOILP *P, vector<list<Box *>> *Bj, list<Box *> *R, vector<multiset<Box *, Compare_value>> *L, point *z, int *RemainingBoxesToExplore) {
	//In case of option of several lists and "RE", R will be a subset of L, which contains all boxes with upper bound dominating the point z.
	//Bj , j=1,..,p   are the set of points that have 1 component equals z (weakly dominates)
	TIEMPO t;
	t.init();
	int i, j;
	multiset<Box *>::iterator k;
	int coincidence_axis;
	int number_of_coincidences;
	bool dominates;

	for (i = 0; i < L->size(); i++) { //For every tree in vector L
		k = L->at(i).begin();
		while (k != L->at(i).end()) {
			coincidence_axis = -1;
			number_of_coincidences = 0;
			dominates = true;
			for (j = 0; j < p; j++) {
				if (z->at(j) >(*k)->ub->at(j)) {
					dominates = false;
					j = p;
				}
				else if (z->at(j) == (*k)->ub->at(j)) {
					coincidence_axis = j;
					number_of_coincidences++;
				}
			}
			if (dominates) {
				if (coincidence_axis == -1) { //Strictly dominates
					R->push_back(*k);
					k = L->at(i).erase(k);
					(*RemainingBoxesToExplore)--;

				}
				else if (number_of_coincidences == 1) {
					Bj->at(coincidence_axis).push_back(*k);
					k = L->at(i).erase(k);
					(*RemainingBoxesToExplore)--;
				}
				else {
					++k;
				}
			}
			else {
				++k;
			}
		}
	}
	t.acum();
	P->Count.TSearching_reachable += t.value();
}

void Find_upper_bounds_that_strictly_contains_z_1(MOILP *P, list<Box *> *R, vector<multiset<Box *, Compare_value>> *L, point *z, int *RemainingBoxesToExplore) {
	//In case of option of several lists and "RA", R will be a subset of L, which contains all boxes with upper bound dominating the point z.
	TIEMPO t;
	t.init();
	int i, j;
	multiset<Box *>::iterator k;
	int coincidence_axis;
	bool dominates;

	for (i = 0; i < L->size(); i++) { //For every tree in vector L
		k = L->at(i).begin();
		while (k != L->at(i).end()) {
			coincidence_axis = -1;
			dominates = true;
			for (j = 0; j < p; j++) {
				if (z->at(j) >(*k)->ub->at(j)) {
					dominates = false;
					j = p;
				}
				else if (z->at(j) == (*k)->ub->at(j)) {
					coincidence_axis = j;
				}
			}
			if (dominates) {
				if (coincidence_axis == -1) { //Strictly dominates
					R->push_back(*k);
					k = L->at(i).erase(k);
					(*RemainingBoxesToExplore)--;

				}
				else {
					++k;
				}
			}
			else {
				++k;
			}
		}
	}
	t.acum();
	P->Count.TSearching_reachable += t.value();
}

bool Calculate_point_t(point *t, Box *B, point *zz) {
	//When the point zz is outside box B, we project all the coordenates outside the box to a point to the box-frontier. Return that "new" point, t
	//t is the projection of zz to the box B
	bool is_inside = true;

	int i, p = int(zz->size());
	for (i = 0; i < p; i++) {
		if (zz->at(i) < B->lb->at(i)) {
			t->at(i) = B->lb->at(i);
			is_inside = false;
		}
		else {
			t->at(i) = zz->at(i);
		}
	}
	return (is_inside);
}

void create_u_j_bounds_p_partition(point **ib, point **ub, Box *B, point *z, point *tz, int coordinate) {
	//Creating the new bounds for the coordinate-sibling using p-partition	

	for (int i = 0; i < p; i++) {
		if (i < coordinate) {
			(*ib)->at(i) = B->lb->at(i);		(*ub)->at(i) = B->ub->at(i);
		}
		else if (i == coordinate) {
			(*ib)->at(i) = B->lb->at(i);		(*ub)->at(i) = tz->at(i);
		}
		else {
			(*ib)->at(i) = tz->at(i);				(*ub)->at(i) = B->ub->at(i);
		}
	}
	if (z->at(coordinate) < (*ib)->at(coordinate)) {
		*ub = NULL;
		return;
	}
}

void create_u_j_bounds_full(MOILP *P, point **ib, point **ub, Box *B, point *z, int coordinate) {
	//Creating the new bounds for the coordinate-sibling using full partition	
	for (int i = 0; i < p; i++) {
		(*ib)->at(i) = P->Ideal.at(i);
		if (i == coordinate) {
			(*ub)->at(i) = z->at(i);
		}
		else {
			(*ub)->at(i) = B->ub->at(i);
		}
		if ((*ib)->at(i) >= (*ub)->at(i)) {
			*ub = NULL;
			return;
		}
	}
}

Box *PARTITION_BOX_full(MOILP *P, Box *it, point *z, int coordinate, int p, Input *input, double *scaling) {
	//Create new box partitioning (it) respect to z, according the coordinate "coordinate"
	//We use full method
	TIEMPO t;
	t.init();

	Box *u_j = new(Box);

	point *ib = new(point), *ub = new(point);
	ib->resize(p); ub->resize(p);

	create_u_j_bounds_full(P, &ib, &ub, it, z, coordinate);
	if (ub != NULL) {
		u_j = CREATE_BOX_AND_VALUE(ib, ub, z, input, scaling, coordinate);
	}
	else {
		t.acum();
		P->Count.TSplit += t.value();
		return (NULL);
	}

	u_j->pos_sibling = coordinate;
	
	t.acum();
	P->Count.TSplit += t.value();
	return u_j;
}

Box *PARTITION_BOX_partition(MOILP *P, Box *it, point *z, point *tz, bool *is_inside, int coordinate, Input *input, double *scaling) {
	//Create new box partitioning (it) respect to tz (projection of z to it), according the coordinate "coordinate".
	//We use p-partition method
	TIEMPO t;
	t.init();

	Box *u_j = new(Box);

	point *ib = new(point), *ub = new(point);
	ib->resize(p); ub->resize(p);

	create_u_j_bounds_p_partition(&ib, &ub, it, z, tz, coordinate);
	if (ub != NULL) {
		u_j = CREATE_BOX_AND_VALUE(ib, ub, z, input, scaling, coordinate);

		if ((coordinate == p - 1) && (is_inside)) {
			if (input->box_value == "reduced") {
				double vol_dominated_part = get_volume(ib, z);
				u_j->value -= vol_dominated_part;
			}
			else if (input->box_value == "reduced_scaled") {
				double vol_dominated_part = get_scaled_volume(ib, tz, scaling);
				u_j->value -= vol_dominated_part;
			}
		}
	}
	else {
		t.acum();
		P->Count.TSplit += t.value();
		return (NULL);
	}


	u_j->pos_sibling = coordinate;

	t.acum();
	P->Count.TSplit += t.value();
	return u_j;
}

void update_to_min_lb(point *LB, point *lb2) {
	//Point LB = min (LB, lb2)
	int i;
	int p = int(LB->size());
	for (i = 0; i < p; i++) {
		if (lb2->at(i) < LB->at(i)) {
			LB->at(i) = lb2->at(i);
		}
	}
}

void UPDATE_BOXES_RE_2(MOILP *P, point *z, vector<multiset<Box *, Compare_value>> *L, Input *input, double *scaling, int *number_of_lists) {
	list<Container> R;
	list<Container>::iterator it, it2;
	int j;

	vector<list<Container>> Bj(p);
	vector<list<Container>> Pj(p);

	//First part. Calculate R an Bj
	Find_upper_bounds_that_contains_z_2(P, &Bj, &R, L, z);

	//Second part. Partition boxes and save in P_j
	if (input->partition == "full") {
		while (R.size() > 0) {
			it = R.begin();
			for (j = 0; j < p; j++) {
				Box *u_j = NULL;
				u_j = PARTITION_BOX_full(P, &(*it->VBox), z, j, p, input, scaling);

				if (u_j != NULL) {
					Container c;
					c.VBox = u_j;	c.position = it->position;
					Pj.at(j).push_back(c);

				}
			}
			Destroy_box(it->VBox);
			R.pop_front();
		}
	}
	else if (input->partition == "p-partition") {
		while (R.size() > 0) {
			it = R.begin();
			point tz(p);
			bool is_inside = Calculate_point_t(&tz, &(*it->VBox), z);

			for (j = 0; j < p; j++) {
				Box *u_j = NULL;
				u_j = PARTITION_BOX_partition(P, &(*it->VBox), z, &tz, &is_inside, j, input, scaling);

				if (u_j != NULL) {
					Container c;
					c.VBox = u_j;	c.position = it->position;
					Pj.at(j).push_back(c);

				}
			}
			Destroy_box(it->VBox);
			R.pop_front();
		}

	}

	//Third part. Filtering solutions
	for (j = 0; j < p; j++) { //Comparing Pj with elements of Pj 
		it = Pj.at(j).begin();
		bool flag = false;;
		while (it != Pj.at(j).end()) {
			it2 = it;
			++it2;
			if (it2 == Pj.at(j).end()) ++it;
			while (it2 != Pj.at(j).end()) {
				flag = false;
				if (compare_vectors(&(*it->VBox->ub), &(*it2->VBox->ub), p)) { //if  (it <= it2)
					update_to_min_lb(&(*it2->VBox->lb), &(*it->VBox->lb)); //Updating it2->lb
					get_value(it2->VBox, input->box_value, scaling);		//Recalculate the value when joining the two boxes

					Destroy_box(it->VBox);

					it = Pj.at(j).erase(it);
					it2 = Pj.at(j).end();
					flag = true;
				}
				else if (compare_vectors(&(*it2->VBox->ub), &(*it->VBox->ub), p)) { //it2 <= it
					update_to_min_lb(&(*it->VBox->lb), &(*it2->VBox->lb)); //Updating it->lb
					get_value(it->VBox, input->box_value, scaling);			//Recalculate the value when joining the two boxes

					Destroy_box(it2->VBox);

					it2 = Pj.at(j).erase(it2);
				}
				else {
					++it2;
				}
			}
			if ((!flag) && (it != Pj.at(j).end()))
				++it;
		}
	}

	for (j = 0; j < p; j++) { //Comparing Pj with elements of Bj 
		it = Pj.at(j).begin();
		bool flag = false;
		while (it != Pj.at(j).end()) {
			it2 = Bj.at(j).begin();
			if (it2 == Bj.at(j).end()) ++it;
			while (it2 != Bj.at(j).end()) {
				flag = false;
				if (compare_vectors(&(*it->VBox->ub), &(*it2->VBox->ub), p)) { //it <= it2
					update_to_min_lb(&(*it2->VBox->lb), &(*it->VBox->lb)); //Updating it->lb			
					get_value(it2->VBox, input->box_value, scaling);		//Recalculate the value when joining the two boxes

					Destroy_box(it->VBox);

					it = Pj.at(j).erase(it);
					it2 = Bj.at(j).end();
					flag = true;
				}
				else {
					++it2;
				}
			}
			if ((!flag) && (it != Pj.at(j).end()))
				++it;
		}
	}

	//Fourth part. Updating TREES
	for (j = 0; j < p; j++) {
		it = Pj.at(j).begin();
		while (it != Pj.at(j).end()) {
			if (it->position < *number_of_lists - 1) {
				INSERT_INTO_LIST_2(P, L, it->VBox, it->position + 1);
			}
			else {
				INSERT_INTO_LIST_2(P, L, it->VBox, it->position);
			}
			it = Pj.at(j).erase(it);
		}
		it = Bj.at(j).begin();
		while (it != Bj.at(j).end()) {
			INSERT_INTO_LIST_2(P, L, it->VBox, it->position);
			it = Bj.at(j).erase(it);
		}
	}

	//Removing all empty trees from L
	list<int> indexes_of_trees_to_remove;
	for (j = 0; j < L->size(); j++) {
		if (L->at(j).size() == 0) { //Empty tree
			indexes_of_trees_to_remove.push_front(j);
		}
	}
	while (indexes_of_trees_to_remove.size() > 0) {
		L->erase(L->begin() + indexes_of_trees_to_remove.front());
		indexes_of_trees_to_remove.pop_front();
	}
}

void UPDATE_BOXES_RE_1(MOILP *P, point *z, vector<multiset<Box *, Compare_value>> *L, Input *input, double *scaling, int *RemainingBoxesToExplore) {
	list<Box *> R;
	list<Box *>::iterator it, it2;
	int j;

	vector<list<Box *>> Bj(p);
	vector<list<Box *>> Pj(p);

	//First part. Calculate R an Bj
	Find_upper_bounds_that_contains_z_1(P, &Bj, &R, L, z, RemainingBoxesToExplore);

	//Second part. Partition boxes and save them in Pj
	if (input->partition == "full") {
		while (R.size() > 0) {
			it = R.begin();
			for (j = 0; j < p; j++) {
				Box *u_j = NULL;
				u_j = PARTITION_BOX_full(P, *it, z, j, p, input, scaling);

				if (u_j != NULL) {
					Pj.at(j).push_back(u_j);
				}
			}
			Destroy_box(*it);
			R.pop_front();
		}
	}
	else if (input->partition == "p-partition") {
		while (R.size() > 0) {
			it = R.begin();
			point tz(p);
			bool is_inside = Calculate_point_t(&tz, *it, z);

			for (j = 0; j < p; j++) {
				Box *u_j = NULL;
				u_j = PARTITION_BOX_partition(P, *it, z, &tz, &is_inside, j, input, scaling);

				if (u_j != NULL) {
					Pj.at(j).push_back(u_j);

				}
			}
			Destroy_box(*it);
			R.pop_front();
		}
	}

	//Third part. Filtering solutions
	for (j = 0; j < p; j++) { //Comparing Pj with elements of Pj 
		it = Pj.at(j).begin();
		bool flag = false;;
		while (it != Pj.at(j).end()) {
			it2 = it;
			++it2;
			if (it2 == Pj.at(j).end()) ++it;
			while (it2 != Pj.at(j).end()) {
				flag = false;
				if (compare_vectors((*it)->ub, (*it2)->ub, p)) { //if  (it <= it2)
					update_to_min_lb((*it2)->lb, (*it)->lb); //Updating it2->lb
					get_value((*it2), input->box_value, scaling);		//Recalculate the value when joining the two boxes

					Destroy_box(*it);
					it = Pj.at(j).erase(it);
					it2 = Pj.at(j).end();
					flag = true;
				}
				else if (compare_vectors((*it2)->ub, (*it)->ub, p)) { //it2 <= it
					update_to_min_lb((*it)->lb, (*it2)->lb); //Updating it2->lb
					get_value((*it), input->box_value, scaling);		//Recalculate the value when joining the two boxes

					Destroy_box(*it2);
					it2 = Pj.at(j).erase(it2);
				}
				else {
					++it2;
				}
			}
			if ((!flag) && (it != Pj.at(j).end()))
				++it;
		}
	}

	for (j = 0; j < p; j++) { //Comparing Pj with elements of Bj 
		it = Pj.at(j).begin();
		bool flag = false;
		while (it != Pj.at(j).end()) {
			it2 = Bj.at(j).begin();
			if (it2 == Bj.at(j).end()) ++it;
			while (it2 != Bj.at(j).end()) {
				flag = false;
				if (compare_vectors((*it)->ub, (*it2)->ub, p)) { //it <= it2
					update_to_min_lb((*it2)->lb, (*it)->lb); //Updating it2->lb
					get_value((*it2), input->box_value, scaling);		//Recalculate the value when joining the two boxes

					Destroy_box(*it);

					it = Pj.at(j).erase(it);
					it2 = Bj.at(j).end();
					flag = true;
				}
				else {
					++it2;
				}
			}
			if ((!flag) && (it != Pj.at(j).end()))
				++it;
		}
	}

	//Fourth part. Updating TREES
	for (j = 0; j < p; j++) {
		it = Pj.at(j).begin();
		while (it != Pj.at(j).end()) {
			INSERT_INTO_LIST_1(P, L, *it, RemainingBoxesToExplore);
			it = Pj.at(j).erase(it);
		}
		it = Bj.at(j).begin();
		while (it != Bj.at(j).end()) {
			INSERT_INTO_LIST_1(P, L, *it, RemainingBoxesToExplore);
			it = Bj.at(j).erase(it);
		}
	}
}

void RUN_ANYTIME_ALGORITHM_2(MOILP *P, Input *input) {
	TIEMPO t_ref;
	t_ref.init();

	Box *B0 = new (Box), *B = new(Box);
	vector<multiset<Box *, Compare_value>> L;

	int stat;
	double obj, tmax = P->max_time;
	bool intime = true;
	int *indices = (int *)malloc(P->n_var * sizeof(int));		for (int i = 0; i < P->n_var; i++) indices[i] = i;
	double *scaling = (double *)malloc(P->dimension * sizeof(double));
	double *x = (double *)malloc(P->n_var * sizeof(double));
	point z;	z.resize(P->dimension);
	int number_of_lists;
	if (input->set_of_boxes == "inf") number_of_lists = MAX_INTEGER;
	else number_of_lists = atoi(input->set_of_boxes.c_str());

	printf("\nExecuting problem %s with <%s %s %s %s %s>", P->name.c_str(), input->partition.c_str(), input->parameterization_model.c_str(), input->box_value.c_str(), input->set_of_boxes.c_str(), input->filtering.c_str());

	if (!Calculate_problem_bounds(P, indices, scaling)) return;			//Ideal, bound for Nadir point and range of objective functions

	printf("\nCalculating Pareto front...");

	CREATE_INITIAL_BOX(&B0, &P->Ideal, &P->BoundforNadir, input, scaling);
	INSERT_INTO_LIST_2(P, &L, B0, 0);
	CREATE_MODEL(P, indices, input);

	while ((L.size() > 0) && (intime)) {

		B = SELECT_NEXT_BOX_2(P, &L);

		Set_constraint_extreme_box(P, B->lb, B->ub);
		set_time_for_solver(P, &t_ref);
		
		Solve(P, &stat, x, &obj);

		if ((stat == CPXMIP_OPTIMAL) || (stat == CPXMIP_OPTIMAL_TOL)) { //stat 101 or 102
			Calculate_Image_of_x(P->F, x, P->n_var, &z, P->dimension);

			intime = check_time_and_pointlimit(&t_ref, tmax, P->PF->size(), P->pointlimit);

			if (intime) {
				Add_new_PF_point(x, P->n_var, &z, P->PF, t_ref.value());


				if (input->filtering == "RE") {
					UPDATE_BOXES_RE_2(P, &z, &L, input, scaling, &number_of_lists);
				}
			}
		}
		else { //Empty box or out of time 
			L.at(0).erase(L.at(0).begin());
			Destroy_box(B);
			if (L.at(0).size() == 0) {
				L.erase(L.begin());
			}
			P->Count.nosolution++;
		}
		intime = check_time_and_pointlimit(&t_ref, tmax, P->PF->size(), P->pointlimit);
	}

	//Finish execution
	t_ref.acum();
	P->Count.total_time = t_ref.value();
	printf("\nEnd of execution after %lf seconds\n", P->Count.total_time);

	//Free allocated memory
	multiset<Box *>::iterator it;
	P->Count.Remainingboxestoexplore = 0;
	while (L.size() > 0) {
		P->Count.Remainingboxestoexplore += (int)L.back().size();
		while (L.back().size() > 0) {
			it = L.back().begin();
			Destroy_box(*it);
			L.back().erase(L.back().begin());
		}
		L.pop_back();
	}
	free(indices);	free(scaling);	free(x);
}

void RUN_ANYTIME_ALGORITHM_1(MOILP *P, Input *input) {
	TIEMPO t_ref;
	t_ref.init();

	Box *B0 = new (Box), *B = new(Box);
	vector<multiset<Box *, Compare_value>> L(p);			//Vector of dimension p. 

	int stat;
	double obj, tmax = P->max_time;
	bool intime = true;
	int *indices = (int *)malloc(P->n_var * sizeof(int));		for (int i = 0; i < P->n_var; i++) indices[i] = i;
	double *scaling = (double *)malloc(p * sizeof(double));
	double *x = (double *)malloc(P->n_var * sizeof(double));
	point z;	z.resize(p);
	int RemainingBoxesToExplore = 0;
	int counter = 0;								//A counter for the current coordinate list to analyze

	printf("\nExecuting problem %s with <%s %s %s %s %s>", P->name.c_str(), input->partition.c_str(), input->parameterization_model.c_str(), input->box_value.c_str(), input->set_of_boxes.c_str(), input->filtering.c_str());

	if (!Calculate_problem_bounds(P, indices, scaling)) return;			//Ideal, bound for Nadir point and range of objective functions

	printf("\nCalculating Pareto front...");

	CREATE_INITIAL_BOX(&B0, &P->Ideal, &P->BoundforNadir, input, scaling);
	INSERT_INTO_LIST_1(P, &L, B0, &RemainingBoxesToExplore);
	CREATE_MODEL(P, indices, input);

	while ((RemainingBoxesToExplore > 0) && (intime)) {

		B = SELECT_NEXT_BOX_1(P, &L, &counter);

		Set_constraint_extreme_box(P, B->lb, B->ub);
		set_time_for_solver(P, &t_ref);
		Solve(P, &stat, x, &obj);

		if ((stat == CPXMIP_OPTIMAL) || (stat == CPXMIP_OPTIMAL_TOL)) { //stat 101 or 102
			Calculate_Image_of_x(P->F, x, P->n_var, &z, p);

			intime = check_time_and_pointlimit(&t_ref, tmax, P->PF->size(), P->pointlimit);

			if (intime) {
				Add_new_PF_point(x, P->n_var, &z, P->PF, t_ref.value());

				if (input->filtering == "RE") {
					UPDATE_BOXES_RE_1(P, &z, &L, input, scaling, &RemainingBoxesToExplore);
				}
			}
		}
		else {	//The box is empty
			L.at(B->pos_sibling).erase(L.at(B->pos_sibling).begin());
			Destroy_box(B);
			RemainingBoxesToExplore--;
			P->Count.nosolution++;
		}
		intime = check_time_and_pointlimit(&t_ref, tmax, P->PF->size(), P->pointlimit);
	}

	//Finish execution
	t_ref.acum();
	P->Count.total_time = t_ref.value();
	printf("\nEnd of execution after %lf seconds\n", P->Count.total_time);

	//Free allocate memory
	multiset<Box *>::iterator it;
	P->Count.Remainingboxestoexplore = RemainingBoxesToExplore;
	for (int i = 0; i < p; i++) {
		it = L.at(i).begin();
		while (it != L.at(i).end()) {
			Destroy_box(*it);
			it = L.at(i).erase(it);
		}
	}
	free(indices);	free(scaling);	free(x);
}

void free_problem_memory(MOILP *P) {
	free(P->lp);
	free(P->env);
	list<solution>::iterator it;
	it = P->PF->begin();
	while (it != P->PF->end()) {
		free(it->x);
		free(it->z);
		it = P->PF->erase(it);
	}
	free(P->PF);
	for (int i = 0; i < P->dimension; i++)
		free(P->F[i]);
	free(P->F);
	P->Ideal.clear();
	P->BoundforNadir.clear();
}

void Calculate_date_and_hour(std::string *fecha, std::string *hora) {
	struct tm *newtime;
	time_t long_time;
	time(&long_time);
	newtime = localtime(&long_time);

	if (newtime->tm_mday < 10) 		*fecha += ("0");
	*fecha += std::to_string(newtime->tm_mday);
	if (newtime->tm_mon < 9)		*fecha += ("0");
	*fecha += std::to_string(newtime->tm_mon + 1);		//Month
	*fecha += std::to_string(newtime->tm_year + 1900);	//Year

	if (newtime->tm_hour < 10) *hora += ("0");
	*hora += std::to_string(newtime->tm_hour);			//Hours
	if (newtime->tm_min < 10) *hora += ("0");
	*hora += std::to_string(newtime->tm_min);			//Minutes
	if (newtime->tm_sec < 10) *hora += ("0");
	*hora += std::to_string(newtime->tm_sec);			//Seconds	
}

void Transform_pareto_front_if_neccesary(MOILP *P) {
	//We have used a MIN problem in CPLEX. If the original problem was MAX, we transform the solution
	int i, j;

	if (P->original_problem_type == CPX_MAX) {
		for (list<solution>::iterator it = P->PF->begin(); it != P->PF->end(); ++it) {
			for (j = 0; j < p; j++) {
				it->z->at(j) = -(it->z->at(j));
			}
		}
		//Return to the original cost values
		for (i = 0; i < p; i++) {
			for (j = 0; j < P->n_var; j++) {
				P->F[i][j] = -(P->F[i][j]);
			}
			P->BoundforNadir.at(i) = -P->BoundforNadir.at(i);
			P->Ideal.at(i) = -P->Ideal.at(i);
		}
	}
}

int compare_lex_4(solution *v1, solution *v2, int index, int p) {
	//Compare two p-dimensional vectors v1 and v2 in lexicographic order.
	//Return -2 if v1 <_lex v2 and not v1 < v2
	//Return -1 if v1 < v2 
	//Return 0 if v1 = v2
	//Return 1  if v1 > v2 
	//Return 2 if v1 >_lex v2 and not v1 > v2

	if (v1->z->at(index) < v2->z->at(index)) {
		for (int i = index + 1; i < p; i++) {
			if (v1->z->at(i) > v2->z->at(i))
				return -2;
		}
		return -1;
	}
	else if (v1->z->at(index) > v2->z->at(index)) {
		for (int i = index + 1; i < p; i++) {
			if (v1->z->at(i) < v2->z->at(i))
				return 2;
		}
		return 1;
	}
	else {
		if (index < p - 1)
			return (compare_lex_4(v1, v2, index + 1, p));
		else
			return 0; //Equal vectors
	}
}

int insert_lexi_solution_without_domination(int type, solution *v, list<solution> *L, list<solution> *Eliminated_vectors, int p) {
	//Insert a new nondominated point z in lexicographic order into the list L
	//Return 0 if success
	//Return -1 if the point is dominated (thus, not included in the list)
	//Return a >= 1 , where a is the number of points which are deleted from the list, because they are dominated by the new point

	//Outputs of compare_lex_4.
	//Compare two p-dimensional vectors v1 and v2 in lexicographic order.
	//Return -2 if v1 < v2
	//Return -1 if v1 <_lex v2 and not v1 < v2
	//Return 0 if v1 = v2
	//Return 1  if v1 > v2 
	//Return 2 if v1 >_lex v2 and not v1 > v2

	std::list<solution>::iterator it1, it2, aux;
	int n_eliminated = 0;
	bool finish = false;

	//Introduce the new point in the front of the list. Then compare it with every element on the list, one by one until its lexicographical position is found.
	L->push_front(*v);
	it1 = it2 = L->begin();
	
	if (L->size() == 2) { //There is only one element to compare with
		std::advance(it2, 1);
		int aa = compare_lex_4(&(*it1), &(*it2), 0, p);
		if (type == CPX_MIN) {
			if (aa == 1) { //it1 = v >  it2
				L->pop_front();
				return -1;
			}
			else if (aa == 2) { //it1 = v >_lex it2
				L->push_back(*it1);				L->pop_front();
			}
			else if (aa == 0) { //The two vectors are equal   it1 = it2
				L->pop_front();
				return -1;
			}
			else if (aa == -2) {//it1 = v <_lex it2. They are ordered

			}
			else if (aa == -1) { //it1 <  it2. We eliminate it2
				Eliminated_vectors->push_back(L->back());
				L->pop_back();				n_eliminated++;
			}
		}
		else if (type == CPX_MAX) {
			if (aa == 1) { //it1 = v >  it2
				Eliminated_vectors->push_back(L->back());
				L->pop_back();				n_eliminated++;
			}
			else if (aa == 2) { //it1 = v >_lex it2. They are ordered

			}
			else if (aa == 0) { //The two vectors are equal   it1 = it2
				L->pop_front();
				return -1;
			}
			else if (aa == -2) {//it1 = v <_lex it2. 
				L->push_back(*it1);				L->pop_front();
			}
			else if (aa == -1) { //it1 <  it2. We eliminate it1
				L->pop_front();
				return -1;
			}
		}
	}
	else if (L->size() > 2) { //The list has at least two elements to compare with
		std::advance(it2, 1);
		int aa;

		if (type == CPX_MIN) {
			while (!finish) {
				aa = compare_lex_4(&(*it1), &(*it2), 0, p);
				if (aa == 1) { //v1 >  v2
					it1 = L->erase(it1);
					return -1;
				}
				else if (aa == 2) { //If v1 >_lex v2, go to the next element of the list
					++it2;
				}
				else if (aa == 0) {
					it1 = L->erase(it1);
					return -1;
				}
				else if (aa == -1) { //v1 < v2
					Eliminated_vectors->push_back(*it2);
					it2 = L->erase(it2);
					n_eliminated++;
					if (it2 == L->end()) {
						L->insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
						L->pop_front();
						return (n_eliminated);
					}
				}
				else { //Position found
					finish = true;
					bool finish2 = false;

					list<solution>::iterator it3 = it2;
					while (!finish2) {			//We look for more dominated points
						aa = compare_lex_4(&(*it1), &(*it3), 0, p);
						if (aa == -1) { 
							Eliminated_vectors->push_back(*it3);
							it3 = L->erase(it3);
							n_eliminated++;
						}
						else {
							++it3;
						}

						if (it3 == L->end()) {
							finish2 = true;
						}
					}
				}
				if (it2 == L->end())
					finish = true;
			}

			L->insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
			L->pop_front();
		}
		else if (type == CPX_MAX) {
			while (!finish) {
				aa = compare_lex_4(&(*it1), &(*it2), 0, p);

				if (aa == 1) { //v1 >  v2
					Eliminated_vectors->push_back(*it2);
					it2 = L->erase(it2);
					n_eliminated++;
					if (it2 == L->end()) {
						L->insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
						L->pop_front();
						return (n_eliminated);
					}
				}
				else if (aa == 2) { //Position found 
					finish = true;
					bool finish2 = false;

					list<solution>::iterator it3 = it2;
					while (!finish2) {			//We look for more dominated points
						aa = compare_lex_4(&(*it1), &(*it3), 0, p);
						if (aa == 1) {
							Eliminated_vectors->push_back(*it3);

							it3 = L->erase(it3);
							n_eliminated++;
						}
						else {
							++it3;
						}

						if (it3 == L->end()) {
							finish2 = true;
						}
					}
				}
				else if (aa == 0) {
					it1 = L->erase(it1);
					return -1;
				}
				else if (aa == -1) { //v1 < v2
					it1 = L->erase(it1);
					return -1;
				}
				else { //If v1 >_lex v2, go to the next element of the list
					++it2;
				}
				if (it2 == L->end())
					finish = true;
			}
			L->insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
			L->pop_front();
		}
	}
	else { //PF->size == 1
	}
	
	return (n_eliminated);
}

void convert_solution_to_a_normal_vector(list<solution> *PF, list<vector<double>> *newPF, int p) {
	//Creates a copy of the Pareto set in a list of vectors of size p.
	list<solution>::iterator it;
	newPF->clear();
	for (it = PF->begin(); it != PF->end(); ++it) {
		vector<double> new_point(p);
		for (int j = 0; j < p; j++) {
			new_point.at(j) = it->z->at(j);
		}
		newPF->push_back(new_point);
	}
}

bool Check_Pareto_front(list<solution> *PF, int original_sense, int p, int n_var) {
	//If the Pareto front have dominated points or a point dominates other points, return false. Otherwise return true
	//When the output is false, PF is modified and ordered in a lexicographic way
	//If Pareto front is correct, the output is not modified.
	list<solution> PF2, Eliminated;
	bool X = true;
	int out;

	for (list<solution>::iterator it = PF->begin(); it != PF->end(); ++it) {
		Eliminated.clear();

		out = insert_lexi_solution_without_domination(original_sense, &(*it), &PF2, &Eliminated, p);

		if (out != 0) {
			X = false;
			if (out == -1) {
				printf("\nRemoving vector "); for (int v = 0; v < p; v++) printf("%0.0f ", it->z->at(v)); printf(" because it is dominated.");
			}
			else { //out > 0
				for (list<solution>::iterator it2 = Eliminated.begin(); it2 != Eliminated.end(); ++it2) {
					printf("\nVector "); for (int v = 0; v < p; v++) printf("%0.0f ", it->z->at(v)); printf(" dominates the previous vector ");
					for (int v = 0; v < p; v++) printf("%0.0f ", it2->z->at(v));
				}
			}
		}

	}
	if (X == false) {
		PF->clear();
		for (list<solution>::iterator it = PF2.begin(); it != PF2.end(); ++it) {
			point *xx = new(point);		xx->resize(n_var);
			point *zz = new(point);		zz->resize(p);
			solution *A = new(solution);

			for (int i = 0; i < n_var; i++)		xx->at(i) = (*it).x->at(i);
			for (int i = 0; i < p; i++) zz->at(i) = (*it).z->at(i);
			A->x = xx;
			A->z = zz;
			A->time_point = (*it).time_point;
			PF->push_back(*A);
		}
		PF2.clear();
	}
	return (X);
}

void copy_PF_to_file(list<solution> *PF, string *file_name) {
	FILE *fp;
	int i, j, size;

#ifdef _WIN64
	#include<direct.h>
	std::string a = ".\\PARETO_FRONTS\\";
	_mkdir(a.c_str()); //Create a new folder. Return 0 if it does not exist, -1 otherwise
#elif __linux__
	std::string a = "./PARETO_FRONTS/";
	mkdir(a.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

	*file_name = a + *file_name;
	fp = fopen((*file_name).c_str(), "w");

	size = int(PF->size());

	if (size > 0) {
		list<solution>::iterator it = PF->begin();
		for (j = 0; j < size - 1; j++) {
			for (i = 0; i < p - 1; i++) fprintf(fp, "%1.0f ", it->z->at(i));
			fprintf(fp, "%1.0f\n", it->z->at(i));
			++it;
		}
		for (i = 0; i < p - 1; i++) fprintf(fp, "%1.0f ", it->z->at(i));
		fprintf(fp, "%1.0f", it->z->at(i));
	}
	fclose(fp);
}

void Print_solution_in_file(MOILP *P, string date, string hour, bool correctPF) {
	string file_name = P->name + '_' + date + hour + to_string(P->Count.total_time);

	copy_PF_to_file(P->PF, &file_name);

	FILE *fp;
	char *file = new char[file_name.length() + 1];
	strcpy(file, file_name.c_str());

	file_name = file_name + "_time";
	fp = fopen(file_name.c_str(), "w");

	if (P->PF->size() > 0) {
		list<solution>::iterator it = P->PF->begin();
		int i, j, size = int(P->PF->size());
		for (j = 0; j < size - 1; j++) {
			for (i = 0; i < p; i++) fprintf(fp, "%.0f ", it->z->at(i));
			fprintf(fp, "%0.5f\n", it->time_point);
			++it;
		}
		for (i = 0; i < p; i++) fprintf(fp, "%.0f ", it->z->at(i));
		fprintf(fp, "%0.5f", it->time_point);
	}
	fclose(fp);	
}

void Write_Results_file(MOILP *P, Input *input, string date, string hour, bool correctPF) {
	FILE *fp;
	std::string archive = P->name;
	std::string algorithm;
	const char *name = "Boxes_results.txt";

	fp = fopen(name, "a+");
	char c;
	fread(&c, sizeof(c), 1, fp);	//Read the first character to check if the file is empty

	if (feof(fp)) { //If the file is empty, create the header
		fprintf(fp, "Date\tHour\tInstance_name\tMaxt\tAlgorithm\t|N|\tIter\tTotal_Time\tPARALLELMODE\tTHREADS\tCorrect\n");
	}
	rewind(fp); //Rewind file
	
	//Write results in file
	int i_param[2];
	CPXgetintparam(*P->env, CPX_PARAM_PARALLELMODE, &i_param[0]);
	CPXgetintparam(*P->env, CPX_PARAM_THREADS, &i_param[1]);

	fprintf(fp, "%s\t%s\t", date.c_str(), hour.c_str());  //Date and hour
	fprintf(fp, "%s\t", archive.c_str());
	if (P->max_time == MAX_DOUBLE)
		fprintf(fp, "----\t");
	else	fprintf(fp, "%0.2f\t", P->max_time);

	algorithm = input->partition + "_" + input->parameterization_model + "_" + input->box_value + "_" + input->set_of_boxes + "_" + input->filtering;

	fprintf(fp, "%s\t%d\t%d\t", algorithm.c_str(), (int)P->PF->size(), P->n_iterations);
	fprintf(fp, "%0.6f\t%d\t%d\t", P->Count.total_time, i_param[0], i_param[1]);
	if (correctPF)
		fprintf(fp, "Yes");
	else
		fprintf(fp, "No");
	fprintf(fp, "\n");
	fclose(fp);
}

int main(int argc, char **argv) {
	//main archive1 % archive2 % tmax % pointlimit % partition % model % value % set_of_boxes % filtering % cpx_param_parallel % cpx_param_threads (11 input parameters)
	
	Input input;								//Input parameters
	MOILP P;									//New Model 
	std::string date, hour;						//Date and system time					
	bool correctPF = true;						//Flag

	Input_control(argc, argv, &input);			//Controlling that the input parameters are correct
	Create_CPLEX_object_and_set_CPLEX_parameters(&P, input.lp_model, input.cpx_param_parallel, input.cpx_param_threads);
	Set_MOILP_structure(&P, &input);

	if (input.set_of_boxes == "alternate") {
		RUN_ANYTIME_ALGORITHM_1(&P, &input);		//Consider p lists of boxes and alternating the lists for every iteration
	}
	else {
		RUN_ANYTIME_ALGORITHM_2(&P, &input);		//Consider a fixed list of boxes (or infinite). 
	}
	
	Calculate_date_and_hour(&date, &hour);
	Transform_pareto_front_if_neccesary(&P);		//If original model was MAX, transform the solution
	correctPF = Check_Pareto_front((&P)->PF, (&P)->original_problem_type, p, (&P)->n_var);
	Print_solution_in_file(&P, date, hour, correctPF);
	Write_Results_file(&P, &input, date, hour, correctPF);
		
	//Free MOILP 
	free_problem_memory(&P);

	return 0;
}
