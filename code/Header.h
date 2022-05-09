#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map> 
#include "math.h"
#include <queue>
#include <algorithm>
#include <set>
#include <iterator>
#include <chrono>
#include <random>
#include <thread>

#include <iomanip>
using namespace std::chrono;
using namespace std;
using std::string;

const double pi = 3.1415926;
const double RADIUS = 6381372;


/* 
the border of CD 
*/

const double lon1_hz = 103.2;
const double lon2_hz = 104.7;
const double lat1_hz = 30.29;
const double lat2_hz = 31.04;

/*
the border of AIS
*/
/*
const double lon1_hz = -180.0;
const double lon2_hz = 155.0;
const double lat1_hz = 3.0;
const double lat2_hz = 90;
*/



/*
the step size for adjusting varepsilon at each time step, i.e.,  \Delta\varepsilon,  (cf. Section 5.3)
*/
const double delta_epsilon = 50;  


const int vehicle_size = 30000; 

const int rho = 5;  /* \rho */
const int xi = 3; /* minPts */

const double alpha = 0.1; /* \alpha */
const double mu = 11.1; /* 
						speed constraint 
						mu=8.23 (16 knots) for AIS dataset
						*/

static double delta = 1000;  /* \delta*/


struct point {  /* a location */
	double x=INT_MIN;
	double y=INT_MIN;
};
struct grid_id { /* a grid cell */
	int u, v;
	vector<int> x_coordinate;
	vector<int> y_coordinate;
};

struct vehicle { /* a moving object */
	int id = -1;

	point loc;
	point loc_pre;
	point loc_ppre;
	point loc_copy;
	point loc_original;

	int no_C_unit = -1;

	grid_id g_id_int;
	string g_id_str;

	bool core = false;
	bool core_pre = false;
	bool seed = false;
	double d_seed = -1;

	bool stop = true;
	bool arrive = false;

	int original_id;
	string t1 = "-1";
	string t2 = "-1";
	string t3 = "-1";
	int delta_t;

	vector<int> neighbor;
	
};
struct cluster { /* a cluster */
	vector<vehicle> o;
	int cluster_id;
};

struct C_unit { /* a set of minimal group */
	vector<int> o;
	int seed;
	double d_max = INT_MAX;
};

/*
  f is the optimal value of f_k(\tilde(s))
  b is the corresponding set of b_opt
*/
struct solution { 
	double f;
	vector<double> b;
};

struct circle_t { /* a circle */
	point center;
	double r;
};

struct unit { /*a structure for cluster mapping */
	int id1;
	int id2;
	double mi;
};


