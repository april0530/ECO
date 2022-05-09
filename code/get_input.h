#pragma once
#include "Header.h"
point millerXY(point& l);
grid_id obtain_grid_id(point& p, double d);
int get_new_data_initial(vector<string>& str, vector<vehicle>& obj, 
	int* id_index, double epsilon);
void get_new_data(vector<string>& str, vector<vehicle>& obj, int* id_index,
	multimap<string, int>& veh, double epsilon);
void reset_stop_flag(vector<vehicle>& obj);

double get_fitness_qs(vector<vector<int>>& cluster, vector<vehicle>& o);
double get_low_qs(vector<vector<int>>& cluster, vector<vehicle>& o, double e);
double get_high_qs(vector<vector<int>>& cluster, vector<vehicle>& o, double e);

double metric_distance(point& loc1, point& loc2);
vector<vector<int>> optimal(vector<vehicle>& o, double& e, multimap<string, int>& grid_index_seed,
	multimap<string, int>& grid_index, vector<C_unit>& C);
void insert_object(std::multimap<string, int>& grid_index, vehicle& o, double epsilon);
string transfer_grid_id(grid_id g);

vector<grid_id> get_affected_grid_id(grid_id id, double d, double epsilon);
void get_C_unit(vector<vehicle>& o, multimap<string, int>& grid_index,
	multimap<string, int>& grid_index_seed, multimap<string, int>& veh,
	vector<C_unit>& C, double epsilon, double beta);
void clean_grid_index(multimap<string, int>& grid_index);
vector<int> count_neighbors(vehicle& o, double d, multimap<string, int>& grid_index_seed,
	multimap<string, int>& grid_index, vector<vehicle>& v, double epsilon);
vector<vector<int>> clustering(multimap<string, int>& grid_index_seed,
	multimap<string, int>& grid_index, vector<C_unit>& C, vector<vehicle>& o, double epsilon);
void smooth(vector<vehicle>& o, vector<C_unit>& C, double epsilon, 
	multimap<string, int>& veh);
long int StringToTime(char* time, int size);
bool CheckTime(char* time1, char* time2, int size1, int size2);
long int TimeInterval(char* time1, char* time2);
int double_equals(double const a, double const b);
double distance_sqr(point const a, point const b);
int intersect(struct circle_t circles[], struct point points[]);
vector<point> get_circles_interscts(struct circle_t circles[],
	struct point points[]);
vector<point> get_border_coordinate();
bool constrain_by_border(vector<point>& intersect);
void speed_cleaning(vehicle& o, vehicle& seed);
double compute_f_value(double b, double d, double eta);

solution return_f_value(vector<vehicle>& o, vehicle seed, double F, double eta);
vector<point>  getPoint(double cx, double cy, double r,
	double stx, double sty, double edx, double edy, int flag);

void reset_C_unit(vector<vehicle>& obj, vector<C_unit>& C);

vector<double> cluster_mapping(vector<vector<int>>& c1, vector<vector<int>>& c2, vector<vehicle> &o, map<int, int>& r);

