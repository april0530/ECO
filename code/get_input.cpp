#include "get_input.h"

/*
transform the geodetic coordinates (longitude and latitude) to rectangular coordinates (x, y)
*/
point millerXY(point& l) {
	double L = (RADIUS * pi) * 2;
	double W = L;
	double H = L / 2;
	double mill = 2.3;
	double x = 1.0 * l.x * pi / 180;
	double y = 1.0 * l.y * pi / 180;
	y = 1.25 * log(tan(0.25 * pi + 0.4 * y));

	x = (W / 2) + (W / (2 * pi)) * x;
	y = (H / 2) - (H / (2 * mill)) * y;
	point ll;
	ll.x = 1.0 * x;
	ll.y = 1.0 * y;
	return ll;
}

/*
process an arriving location that was not active at the last time step
*/
int get_new_data_initial(vector<string>& str, vector<vehicle>& obj,
	int* id_index, double epsilon) {

	point loc_temp;
	loc_temp.x = stod(str[2]);
	loc_temp.y = stod(str[3]);
	vehicle obj_temp;
	obj_temp.loc_original = loc_temp;
	loc_temp = millerXY(loc_temp);

	obj_temp.g_id_int = obtain_grid_id(loc_temp, 1.0 * epsilon / sqrt(2));

	obj_temp.id = obj.size();
	id_index[stoi(str[1])] = obj_temp.id;
	obj_temp.original_id = stoi(str[1]);
	obj_temp.loc.x = loc_temp.x;
	obj_temp.loc.y = loc_temp.y;
	obj_temp.arrive = true;
	obj_temp.stop = false;

	obj_temp.t1 = str[14];

	obj.push_back(obj_temp);

	return obj.size() - 1;

}
/*
process an arriving location
*/
void get_new_data(vector<string>& str, vector<vehicle>& obj, int* id_index,
	multimap<string, int>& veh, double epsilon) {

	int id = id_index[stoi(str[1])];
	int i;
	if (id >= 0)
	{
		obj[id].loc_pre = obj[id].loc;
		obj[id].loc_pre = obj[id].loc;

		point p;
		p.x = stod(str[2]);
		p.y = stod(str[3]);
		obj[id].loc_original = p;
		p = millerXY(p);

		point o = p;

		obj[id].loc.x = o.x;
		obj[id].loc.y = o.y;
		obj[id].loc_copy = obj[id].loc;
		obj[id].g_id_int = obtain_grid_id(obj[id].loc, 1.0 * epsilon / sqrt(2));
		obj[id].original_id = stoi(str[1]);
		obj[id].stop = false;
		obj[id].arrive = false;
		obj[id].t3 = obj[id].t2;
		obj[id].t2 = obj[id].t1;
		obj[id].t1 = str[14];
		i = id;
		int delta_t = TimeInterval((char*)obj[id].t2.data(), (char*)obj[id].t1.data());
		obj[id].delta_t = delta_t;
	}
	else {
		i = get_new_data_initial(str, obj, id_index, epsilon);
	}
}

/*
clear minimal groups at the curretn time step
*/
void reset_C_unit(vector<vehicle>& obj, vector<C_unit>& C) {
	vector<C_unit> c_clear;
	C.swap(c_clear);
	for (int i = 0; i < obj.size(); i++) {
		obj[i].no_C_unit = -1;
		if (obj[i].stop) obj[i].seed = false;
		else {
			if (obj[i].seed) {
				int c_id = C.size();
				obj[i].no_C_unit = c_id;
				C_unit C_temp;
				C_temp.seed = obj[i].id;
				C.push_back(C_temp);
			}
		}
	}
}
/*
reset flags of moving objects for processings at the next time step
*/
void reset_stop_flag(vector<vehicle>& obj) {
	for (int i = 0; i < obj.size(); i++) {
		obj[i].stop = true;
		obj[i].arrive = false;
		obj[i].core_pre = obj[i].core;
		obj[i].core = false;
	}
}
/*
clear grid indexes at the current time step
*/
void clean_grid_index(multimap<string, int>& grid_index)
{
	multimap<string, int> tmp1;
	grid_index.swap(tmp1);
}
/*
calculate Euclidean distance
*/
double metric_distance(point& loc1, point& loc2) {
	double x = 1.0 * (pow(loc1.x - loc2.x, 2) + pow(loc1.y - loc2.y, 2));
	return 1.0 * sqrt(x);
}

/*
calculate QS_m, i.e., QS
*/
double get_fitness_qs(vector<vector<int>>& cluster, vector<vehicle>& o) {
	double ts = 0.0;
	double is = 0.0;
	double ds = 0.0;
	double TS_DIS = 0.0;
	double IS_DIS = 0.0;
	double DS_DIS = 0.0;
	int TS_COUNT = 0;
	int IS_COUNT = 0;
	int DS_COUNT = 0;
	for (int i = 0; i < o.size(); i++) {
		if (!o[i].stop) {
			for (int j = i + 1; j < o.size(); j++) {
				if (!o[j].stop) {
					double d = metric_distance(o[i].loc, o[j].loc);
					TS_DIS += d;
					TS_COUNT++;
					if (d == 0) ts += 0.0;
					else ts += 1.0 / d;
				}
			}
		}
	}
	int n1 = 0;
	for (int i = 0; i < cluster.size(); i++) {
		for (int j = 0; j < cluster[i].size(); j++) {
			int a = cluster[i][j];
			if (!o[a].stop) {
				for (int k = j + 1; k < cluster[i].size(); k++)
				{
					int b = cluster[i][k];
					if (!o[b].stop) {
						double d = metric_distance(o[a].loc, o[b].loc);
						IS_DIS += d;
						IS_COUNT++;
						if (d == 0) is += 0.0;
						else is += 1.0 / d;
					}
					n1++;
				}
			}
		}
	}
	int n2 = 0;
	for (int i = 0; i < cluster.size(); i++) {
		double tmp = 0.0;
		for (int m = 0; m < cluster[i].size(); m++) {
			int a = cluster[i][m];
			if (!o[a].stop) {
				for (int k = i + 1; k < cluster.size(); k++) {
					for (int j = 0; j < cluster[k].size(); j++) {
						int b = cluster[k][j];
						{
							if (!o[b].stop) {
								double d = metric_distance(o[a].loc, o[b].loc);
								DS_DIS += d;
								DS_COUNT++;
								if (d == 0) ds += 0.0;
								else ds += (1.0 / d)*(1.0 / d);
							}
						}

					}
					n2++;
				}
			}
		}
	}
	double qs = 1.0 * is / ts - 1.0 * ds / (ts * ts);
	return qs;
}

/*
calculate QS_h
*/
double get_high_qs(vector<vector<int>>& cluster, vector<vehicle>& o, double e) {
	double ts = 0.0;
	double is = 0.0;
	double ds = 0.0;
	for (int i = 0; i < o.size(); i++) {
		if (!o[i].stop) {
			for (int j = i + 1; j < o.size(); j++) {
				if (!o[j].stop) {
					double d = metric_distance(o[i].loc, o[j].loc);
					if (d >= e && d <= e + delta_epsilon) {
						if (d == 0) ts += 0.0;
						else ts += 1.0 / d;
					}

				}
			}
		}
	}

	for (int i = 0; i < cluster.size(); i++) {
		for (int j = 0; j < cluster[i].size(); j++) {
			int a = cluster[i][j];
			if (!o[a].stop) {
				for (int k = j + 1; k < cluster[i].size(); k++)
				{
					int b = cluster[i][k];
					if (!o[b].stop) {

						double d = metric_distance(o[a].loc, o[b].loc);
						if (d >= e && d <= e + delta_epsilon) {
							if (d == 0) is += 0.0;
							else is += 1.0 / d;
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < cluster.size(); i++) {
		double tmp = 0.0;
		for (int m = 0; m < cluster[i].size(); m++) {
			int a = cluster[i][m];
			if (!o[a].stop) {
				for (int k = i + 1; k < cluster.size(); k++) {
					for (int j = 0; j < cluster[k].size(); j++) {
						int b = cluster[k][j];
						{
							if (!o[b].stop) {
								double d = metric_distance(o[a].loc, o[b].loc);
								if (d >= e && d <= e + delta_epsilon) {
									if (d == 0) ds += 0.0;
									else ds += (1.0 / d) * (1.0 / d);
								}
							}
						}
					}
				}
			}
		}
	}
	double qs = 1.0 * is / ts - 1.0 * ds / (ts * ts);

	return qs;
}
/*
calculate QS_l
*/
double get_low_qs(vector<vector<int>>& cluster, vector<vehicle>& o, double e) {
	double ts = 0.0;
	double is = 0.0;
	double ds = 0.0;
	for (int i = 0; i < o.size(); i++) {
		if (!o[i].stop) {
			for (int j = i + 1; j < o.size(); j++) {
				if (!o[j].stop) {
					double d = metric_distance(o[i].loc, o[j].loc);
					if (d >= e - delta_epsilon && d < e ) {
						if (d == 0) ts += 0.0;
						else ts += 1.0 / d;
					}

				}
			}
		}

	}

	for (int i = 0; i < cluster.size(); i++) {
		for (int j = 0; j < cluster[i].size(); j++) {
			int a = cluster[i][j];
			if (!o[a].stop) {
				for (int k = j + 1; k < cluster[i].size(); k++)
				{
					int b = cluster[i][k];
					if (!o[b].stop) {

						double d = metric_distance(o[a].loc, o[b].loc);
						if (d >= e - delta_epsilon && d < e) {
							if (d == 0) is += 0.0;
							else is += 1.0 / d;
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < cluster.size(); i++) {
		double tmp = 0.0;
		for (int m = 0; m < cluster[i].size(); m++) {
			int a = cluster[i][m];
			if (!o[a].stop) {
				for (int k = i + 1; k < cluster.size(); k++) {
					for (int j = 0; j < cluster[k].size(); j++) {
						int b = cluster[k][j];
						{
							if (!o[b].stop) {
								double d = metric_distance(o[a].loc, o[b].loc);
								if (d >= e - delta_epsilon && d < e) {
									if (d == 0) ds += 0.0;
									else ds += (1.0 / d) * (1.0 / d);
								}
							}
						}
					}
				}
			}
		}
	}
	double qs = 1.0 * is / ts - 1.0 * ds / (ts * ts);

	return qs;
}
/*
iteratively optimize QS
*/
vector<vector<int>> optimal(vector<vehicle>& o, double& e, multimap<string, int>& grid_index_seed,
	multimap<string, int>& grid_index, vector<C_unit>& C) {
	vector<vector<int>> c_mid;
	vector<vector<int>> pre;
	double q_mid, q_high, q_low;
	double q = 1.0 * INT_MIN;
	double e_pre = e;
	int xx = 0;
	while (true)
	{
		auto start = chrono::high_resolution_clock::now();
		c_mid = clustering(grid_index_seed,
			grid_index, C, o, e);

		auto end = chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end - start;

		q_mid = get_fitness_qs(c_mid, o);

		q_high = get_high_qs(c_mid, o, e);
		q_low = get_low_qs(c_mid, o, e);

		if (q_mid >= q_high && q_mid >= q_low) {
			break;
		}
		else if (q >= q_mid)
		{
			break;
		}
		else {
			if (q_high < q_low) {
				e += delta_epsilon;
			}
			else e -= delta_epsilon;
		}
		q = q_mid;
		pre = c_mid;
		xx++;
	
	}
	return pre;
}
/*
map a location to a grid cell
*/
grid_id obtain_grid_id(point& p, double d) {
	point l1;
	l1.x = lon1_hz;
	l1.y = lat1_hz;
	l1 = millerXY(l1);

	grid_id id;

	double a = 1.0 * (p.x - l1.x) / d;
	double b = 1.0 * (p.y - l1.y) / d;

	int aa = floor(a);
	int bb = floor(b);
	/*
	map the points located at the border of a grid cell to its top-right grid cell
	*/
	
	id.u = aa;
	id.v = bb;

	return id;
}
/*
transfer the two-dimensional ID of a grid cell to a one-dimensional ID
*/
string transfer_grid_id(grid_id g) {
	string a, b, c;
	if (g.u >= 0) a = to_string(g.u);
	else {
		a = "-";
		a += to_string(-g.u);
	}
	if (g.v >= 0) b = to_string(g.v);
	else {
		b = "-";
		b += to_string(-g.v);
	}
	c = a + "&" + b;
	return c;
}

/*
insert moving objects into grid indexes
*/
void insert_object(std::multimap<string, int>& grid_index, vehicle& o, double epsilon) {
	o.g_id_int = obtain_grid_id(o.loc, 1.0 * epsilon / sqrt(2));
	string id = transfer_grid_id(o.g_id_int);
	grid_index.insert(std::pair<string, int>(id, o.id));
	o.g_id_str = id;
}

/*
get d-close grid cells of a location
*/
vector<grid_id> get_affected_grid_id(grid_id id, double d, double epsilon) {
	vector<grid_id> affected_grid_id;
	if (d <= epsilon && d > epsilon / sqrt(2)) {
		for (int i = id.u - 2; i <= id.u + 2; i++) {
			for (int j = id.v - 2; j <= id.v + 2; j++) {
				if ((i == id.u - 2 && j == id.v - 2) || (i == id.u - 2 && j == id.v + 2) ||
					(i == id.u + 2 && j == id.v - 2) || (i == id.u + 2 && j == id.v + 2)) continue;
				grid_id id_temp;
				id_temp.u = i;
				id_temp.v = j;
				affected_grid_id.push_back(id_temp);
			}
		}
		
	}
	else {
		for (int i = id.u - 1; i <= id.u + 1; i++) {
			for (int j = id.v - 1; j <= id.v + 1; j++) {

				grid_id id_temp;
				id_temp.u = i;
				id_temp.v = j;
				affected_grid_id.push_back(id_temp);

			}
		}
	}

	return affected_grid_id;
}
/*
generate minimal groups
*/
void get_C_unit(vector<vehicle>& o, multimap<string, int>& grid_index,
	multimap<string, int>& grid_index_seed, multimap<string, int>& veh,
	vector<C_unit>& C, double epsilon, double beta) {

	srand((unsigned)time(NULL));
	int x = 0;
	/*
	get all seed points
	*/

	for (int u = 0; u < o.size(); u++) {
		if (!o[u].stop) {
			bool find = false;
			
			o[u].g_id_int = obtain_grid_id(o[u].loc, 1.0*epsilon/sqrt(2));

			vector<grid_id> g(get_affected_grid_id(o[u].g_id_int, beta, epsilon));

			for (int j = 0; j < g.size(); j++) {
				string id = transfer_grid_id(g[j]);
				auto res1 = grid_index_seed.equal_range(id);

				for (multimap<string, int>::iterator it = res1.first; it != res1.second; ++it)
				{
					int seed = it->second;
					if (o[seed].id != o[u].id ) {
						double d = metric_distance(o[seed].loc, o[u].loc);
						if (d <= beta) {
							find = true;
							break;
						}
					}
				}
			}
			
			if (!find)
			{
				x++;
				o[u].seed = true;
				int c_id = C.size();
				o[u].no_C_unit = c_id;
				C_unit C_temp;
				C_temp.seed = o[u].id;
				C.push_back(C_temp);
				/*
				grid_index_seed stores seed points
				*/
				insert_object(grid_index_seed, o[u], epsilon);
			}
		}
	}
	/*
	find the seed point for each non-seed point
	*/
	for (int u = 0; u < o.size(); u++) {

		if (!o[u].seed && !o[u].stop) {
			/*
			grid_index stores non-seed points
			*/
			insert_object(grid_index, o[u], epsilon);
			double min = 1.0 * INT_MAX;
			vehicle closest_seed;

			o[u].g_id_int = obtain_grid_id(o[u].loc, 1.0 * epsilon / sqrt(2));
			vector<grid_id> g(get_affected_grid_id(o[u].g_id_int, beta, epsilon));

			for (int j = 0; j < g.size(); j++) {
				string id = transfer_grid_id(g[j]);
				auto res1 = grid_index_seed.equal_range(id);

				for (multimap<string, int>::iterator it = res1.first; it != res1.second; ++it)
				{
					int seed = it->second;

					if (o[seed].id != o[u].id) {
						double d = metric_distance(o[seed].loc, o[u].loc);
						if (d <= beta) {
							if (min > d) {
								closest_seed = o[seed];
								min = d;
							}
						}

					}

				}
			}


			if (closest_seed.id >= 0) {

				x++;
				int c_id = closest_seed.no_C_unit;
				o[u].no_C_unit = c_id;

				vehicle tmp = o[C[c_id].seed];

				double d = metric_distance(o[u].loc, tmp.loc);

				o[u].d_seed = d;

				C[c_id].o.push_back(o[u].id);

				if (min > C[c_id].d_max) C[c_id].d_max = min;

			}
		}
	}
}
/*
get d-neighbors of a moving object
*/
vector<int> count_neighbors(vehicle& o, double d, multimap<string, int>& grid_index_seed,
	multimap<string, int>& grid_index, vector<vehicle>& v, double epsilon) {
	int count = 0;
	grid_id id = o.g_id_int;
	vector<grid_id> g(get_affected_grid_id(id, d, epsilon));
	vector<int> n;
	for (int i = 0; i < g.size(); i++) {
		string id_str = transfer_grid_id(g[i]);
		auto res = grid_index_seed.equal_range(id_str);
		for (auto iter2 = res.first; iter2 != res.second; ++iter2)
		{

			int id = iter2->second;
			double dist = metric_distance(v[id].loc, o.loc);
			if (dist <= d) n.push_back(id);

		}
		res = grid_index.equal_range(id_str);
		for (auto iter2 = res.first; iter2 != res.second; ++iter2)
		{
			int id = iter2->second;
			double dist = metric_distance(v[id].loc, o.loc);
			if (dist <= d) n.push_back(id);

		}
	}
	return n;
}
/*
DBSCAN clustering
*/
vector<vector<int>> clustering(multimap<string, int>& grid_index_seed,
	multimap<string, int>& grid_index, vector<C_unit>& C, vector<vehicle>& o, double epsilon) {

	for (int j = 0; j < o.size(); j++)
	{
		bool if_core = false;

		if (!o[j].stop && !o[j].core) {

			int a = grid_index_seed.count(o[j].g_id_str);
			int b = grid_index.count(o[j].g_id_str);
			
			if (a + b >= xi) { /* xi is minPts in the paper */
				o[j].core = true; /* grid index based accelration */

				auto it_tmp = grid_index_seed.equal_range(o[j].g_id_str);

				while (it_tmp.first != it_tmp.second) {
					if(!o[it_tmp.first->second].stop)
					o[it_tmp.first->second].core = true;

					it_tmp.first++;
				}
				it_tmp = grid_index.equal_range(o[j].g_id_str);

				while (it_tmp.first != it_tmp.second) {
					if (!o[it_tmp.first->second].stop)
					o[it_tmp.first->second].core = true;

					it_tmp.first++;
				}

				if_core = true;
			}
			else if (o[j].seed && (C[o[j].no_C_unit].o.size() + 1) >= xi) {
				/*
				minimal group based accelration
				*/
				
				o[j].core = true;
				if_core = true;
			}
			if (!if_core) {

				int neighbor = count_neighbors(o[j], epsilon, grid_index_seed, grid_index, o, epsilon).size();
				if (neighbor >= xi) {
					if_core = true;
					o[j].core = true;
				}
			}

		}
	}

	vector<int> c, ID;
	for (int i = 0; i < o.size(); i++) {
		if (o[i].core) {
			c.push_back(o[i].id);
		}
		ID.push_back(o[i].id);
	}


	vector<vector<int>> cl;

	while (c.size() > 0) {

		vector<int> ID_old(ID);
		queue<int> q;
		q.push(c[0]);

		std::pair<multimap<int, int>::iterator, multimap<int, int>::iterator> res1;
		vector<int>::iterator  it = find(ID.begin(), ID.end(), c[0]);
		if (it != ID.end()) {
			int d = distance(ID.begin(), it);
			ID.erase(ID.begin() + d);
		}

		while (q.size() > 0) {
			int tf = q.front();
			q.pop();
			if (o[tf].core) {
				
				vector<vehicle> Delta;

				vector<int> neighbor(count_neighbors(o[tf], epsilon, grid_index_seed, grid_index, o, epsilon));

				for (int i = 0; i < neighbor.size(); i++) {
					if (!o[neighbor[i]].stop){
					
						vector<int>::iterator  it = find(ID.begin(), ID.end(), o[neighbor[i]].id);
						if (it != ID.end()) {
							Delta.push_back(o[neighbor[i]]);
							q.push(o[neighbor[i]].id);

							int d = distance(ID.begin(), it);
							ID.erase(ID.begin() + d);
						}
					}
				}

			}
		}
		vector<int> cl_tmp;

		for (int i = 0; i < ID_old.size(); i++) {
			int id = ID_old[i];
			vector<int>::iterator it = find(ID.begin(), ID.end(), id);
			if (it == ID.end()) {
				cl_tmp.push_back(o[id].id);
				for (int j = 0; j < c.size(); j++) {
					if (c[j] == id)
					{
						c.erase(c.begin() + j);
					}
				}
			}
		}
		cl.push_back(cl_tmp);
	}

	return cl;
}

/*
smoothing
*/
void smooth(vector<vehicle>& o, vector<C_unit>& C, double epsilon, multimap<string, int>& veh) {


	for (int i = 0; i < C.size(); i++) {

		if ((C[i].o.size()+1) < rho) {
			continue;
		}
		int seed_index;
		solution s_min;
		s_min.f = 1.0 * INT_MAX;
		C[i].o.push_back(C[i].seed);

		vector<vehicle> o_temp;

		bool find = false;
		o[C[i].seed].seed = false;
		for (int j = 0; j < C[i].o.size(); j++) {

			int id = C[i].o[j];
			if (!o[id].stop && !o[id].arrive)
			{
				o_temp.push_back(o[id]);
			}
		}

		for (int j = o_temp.size() - 1; j >= 0; j--) {
			double d = metric_distance(o_temp[j].loc_copy, o_temp[j].loc_pre);
			if (d > o_temp[j].delta_t* mu) {
				vector<point> loc(getPoint(o_temp[j].loc_pre.x, o_temp[j].loc_pre.y, o_temp[j].delta_t * mu, o_temp[j].loc_pre.x, o_temp[j].loc_pre.y, o_temp[j].loc_copy.x, o_temp[j].loc_copy.y, 1));
				o_temp[j].loc = loc[0];
			}
			for (int w = 0; w < o_temp.size(); w++) {
				if (w != j) speed_cleaning(o_temp[w], o_temp[j]);
				/*
				speed-based preprocessing
				*/
			}
			solution s = return_f_value(o_temp, o_temp[j], s_min.f, delta);
			/*
			computing f(o_temp[j].l) (cf. Formula 7 in the paper)
			*/
			if (s.f < s_min.f) {
				s_min = s;
				/*
				store the new seed point, s_new
				*/
				seed_index = j;
			}

		}

		if (s_min.f > 0 && s_min.f < 1.0 * INT_MAX) {
			o_temp[seed_index].seed = true;
			vehicle seed = o_temp[seed_index];

			for (int j = 0; j < o_temp.size(); j++) {

				if (s_min.b[j] > -1) {
					int id_temp = o_temp[j].id;

					point tmp = o[id_temp].loc;

					vector<point> loc(getPoint(seed.loc.x, seed.loc.y, s_min.b[j] * delta, seed.loc.x, seed.loc.y, o_temp[j].loc.x, o_temp[j].loc.y, 1));

					o[id_temp].loc = loc[0];
				}
			}
		}
		C[i].o.pop_back();
	}
	for (int i = 0; i < o.size(); i++) {
		if (!o[i].stop) insert_object(veh, o[i], epsilon);
	}
	return;
}

/*
functions "StringToTime", "CheckTime", and "TimeInterval" transform the type of timestamp from string to int
*/
long int StringToTime(char* time, int size)
{
	string timestr;
	timestr = time;
	int hour, minute, second;
	long int timevalue;
	if (size == 7)
	{
		hour = atoi(timestr.substr(0, 1).c_str());
		minute = atoi(timestr.substr(2, 2).c_str());
		second = atoi(timestr.substr(5, 2).c_str());
		timevalue = second + minute * 60 + hour * 3600;

	}
	else if (size == 8)
	{
		hour = atoi(timestr.substr(0, 2).c_str());
		minute = atoi(timestr.substr(3, 2).c_str());
		second = atoi(timestr.substr(6, 2).c_str());
		timevalue = second + minute * 60 + hour * 3600;
	}
	else timevalue = 99999;
	return timevalue;
}
bool CheckTime(char* time1, char* time2, int size1, int size2)
{
	string timestr1, timestr2;
	timestr1 = time1;
	timestr2 = time2;
	int hour1, hour2;
	if (size1 == size2 && size1 == 7) {
		hour1 = atoi(timestr1.substr(0, 1).c_str());
		hour2 = atoi(timestr2.substr(0, 1).c_str());
	}
	else if (size1 == size2 && size1 == 8) {
		hour1 = atoi(timestr1.substr(0, 2).c_str());
		hour2 = atoi(timestr2.substr(0, 2).c_str());
	}
	else cout << "error formats of timestamps" << endl;
	if (hour1 > hour2) return 0; 
	else return 1;
}

long int TimeInterval(char* time1, char* time2)
{

	int size1 = strlen(time1);
	int size2 = strlen(time2);
	long int value1, value2;
	if (CheckTime(time1, time2, size1, size1))
	{
		value1 = StringToTime(time1, size1);
		value2 = StringToTime(time2, size2);

	}
	else
	{
		value1 = StringToTime(time1, size1);
		value2 = StringToTime(time2, size2) + 24 * 3600;

	}
	if (value2 - value1 < 0) cout << "error in transferring timestamps" << endl;
	return value2 - value1;

}
/*
calculate the intersection of a circle and a segment
*/
vector<point>  getPoint(double cx, double cy, double r,
	double stx, double sty, double edx, double edy, int flag) {
	int flagg = -1;
	point p1, p2;
	point p;
	p.x = edx;
	p.y = edy;
	if (abs(edx - stx) <= pow(10, -6)) {
		p1.x = r + cx;
		p1.y = cy;
		p2.x = r - cx;
		p2.y = cy;
		flagg = 0;
	}
	else if (abs(edy - sty) <= pow(10, -6)) {
		p1.x = cx;
		p1.y = r + cy;
		p2.x = cx;
		p2.y = r - cy;
		flagg = 2;

	}
	else {
		double k = (edy - sty) / (edx - stx);
		double b = edy - k * edx;

		double x1, y1, x2, y2;
		double c = cx * cx + (b - cy) * (b - cy) - r * r;
		double a = (1 + k * k);
		double b1 = -(2 * cx - 2 * k * (b - cy));
		double u = b1 * b1 - 4 * a * c;
		if (u < 0) {
			u = 0;
		}
		double tmp = sqrt(u);

		p1.x = (-b1 - tmp) / (2 * a);
		p1.y = k * p1.x + b;


		p2.x = (-b1 + tmp) / (2 * a);
		p2.y = k * p2.x + b;
		flagg = 2;

	}

	double d1 = metric_distance(p1, p);
	double d2 = metric_distance(p2, p);
	vector<point> re;
	if (flag == 1) {
		if (d1 > d2) re.push_back(p2);
		else re.push_back(p1);
	}
	else if (flag == 0) {
		if (d1 <= d2) re.push_back(p2);
		else re.push_back(p1);
	}
	else {
		re.push_back(p1);
		re.push_back(p2);
	}
	return re;
}
/*
functions "double_equals", "distance_sqr", "intersect" and "get_circles_interscts"
calculate the intersection between two circles
*/
int double_equals(double const a, double const b)
{
	static const double ZERO = 1e-9;
	return fabs(a - b) < ZERO;
}
double distance_sqr(point const a, point const b)
{
	return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}
int intersect(struct circle_t circles[], struct point points[])
{
	double d, a, b, c, p, q, r;
	double cos_value[2], sin_value[2];
	if (double_equals(circles[0].center.x, circles[1].center.x)
		&& double_equals(circles[0].center.y, circles[1].center.y)
		&& double_equals(circles[0].r, circles[1].r)) {
		return -1;
	}

	d = metric_distance(circles[0].center, circles[1].center);
	if (d > circles[0].r + circles[1].r
		|| d < fabs(circles[0].r - circles[1].r)) {
		return 0;
	}

	a = 2.0 * circles[0].r * (circles[0].center.x - circles[1].center.x);
	b = 2.0 * circles[0].r * (circles[0].center.y - circles[1].center.y);
	c = circles[1].r * circles[1].r - circles[0].r * circles[0].r
		- distance_sqr(circles[0].center, circles[1].center);
	p = a * a + b * b;
	q = -2.0 * a * c;

	if (double_equals(d, circles[0].r + circles[1].r)
		|| double_equals(d, fabs(circles[0].r - circles[1].r))) {
		cos_value[0] = -q / p / 2.0;
		sin_value[0] = sqrt(1 - cos_value[0] * cos_value[0]);

		points[0].x = circles[0].r * cos_value[0] + circles[0].center.x;
		points[0].y = circles[0].r * sin_value[0] + circles[0].center.y;

		if (!double_equals(distance_sqr(points[0], circles[1].center),
			circles[1].r * circles[1].r)) {
			points[0].y = circles[0].center.y - circles[0].r * sin_value[0];
		}
		return 1;
	}

	r = c * c - b * b;
	if (q * q - 4.0 * p * r < 0) return 0;
	cos_value[0] = (sqrt(q * q - 4.0 * p * r) - q) / p / 2.0;
	cos_value[1] = (-sqrt(q * q - 4.0 * p * r) - q) / p / 2.0;
	sin_value[0] = sqrt(1 - cos_value[0] * cos_value[0]);
	sin_value[1] = sqrt(1 - cos_value[1] * cos_value[1]);

	points[0].x = circles[0].r * cos_value[0] + circles[0].center.x;
	points[1].x = circles[0].r * cos_value[1] + circles[0].center.x;
	points[0].y = circles[0].r * sin_value[0] + circles[0].center.y;
	points[1].y = circles[0].r * sin_value[1] + circles[0].center.y;

	if (!double_equals(distance_sqr(points[0], circles[1].center),
		circles[1].r * circles[1].r)) {
		points[0].y = circles[0].center.y - circles[0].r * sin_value[0];
	}
	if (!double_equals(distance_sqr(points[1], circles[1].center),
		circles[1].r * circles[1].r)) {
		points[1].y = circles[0].center.y - circles[0].r * sin_value[1];
	}
	if (double_equals(points[0].y, points[1].y)
		&& double_equals(points[0].x, points[1].x)) {
		if (points[0].y > 0) {
			points[1].y = -points[1].y;
		}
		else {
			points[0].y = -points[0].y;
		}
	}
	return 2;
}
vector<point> get_circles_interscts(struct circle_t circles[],
	struct point points[])
{
	vector<point> p;
	point p1;

	switch (intersect(circles, points)) {
	case -1:
		p1.x = INT_MIN;
		p1.y = INT_MIN;
		p.push_back(p1);
		return p;
	case 0:
		p1.x = INT_MAX;
		p1.y = INT_MAX;
		p.push_back(p1);
		return p;
	case 1:
		p.push_back(points[0]);

		break;
	case 2:
		p.push_back(points[0]);
		p.push_back(points[1]);
	}
	return p;
}
/*
"get_border_coordinate" and "constrain_by_border" constrain that
the smoothed location does not exceed the border of the city
*/
vector<point> get_border_coordinate() {
	point l1;
	l1.x = lon1_hz;
	l1.y = lat1_hz;
	l1 = millerXY(l1);


	point l2;
	l2.x = lon2_hz;
	l2.y = lat2_hz;
	l2 = millerXY(l2);


	vector<point> p;
	p.push_back(l1);
	p.push_back(l2);
	return p;
}

bool constrain_by_border(vector<point>& intersect) {

	vector<point> border = get_border_coordinate();
	point p1 = border[0];
	point p2 = border[1];
	bool change = false;
	for (int i = 0; i < intersect.size(); i++) {
		if (intersect[i].x > std::max(p1.x, p2.x)) {
			intersect[i].x = std::max(p1.x, p2.x);
			change = true;
		}
		if (intersect[i].x < min(p1.x, p2.x))
		{
			intersect[i].x = min(p1.x, p2.x);
			change = true;
		}
		if (intersect[i].y > std::max(p1.y, p2.y)) {
			intersect[i].y = std::max(p1.y, p2.y);
			change = true;
		}
		if (intersect[i].y < min(p1.y, p2.y)) {
			intersect[i].y = min(p1.y, p2.y);
			change = true;
		}
	}
	return change;
}
/*
speed-based preprocessing
*/
void speed_cleaning(vehicle& o, vehicle& seed) {

	if (o.t2 != "-1") {
		int delta_t = o.delta_t;
		double max_dis = 1.0 * delta_t * mu;
		max_dis = 1.0 * max_dis;
		double d = metric_distance(o.loc_pre, o.loc_copy);
		struct circle_t circles[2];
		point initial = o.loc_copy;

		vector<point> intersect;

		if (d > max_dis)
		{/*
		 the original location does not follow the spped constraint
		 */
			int flag;
			double d1 = metric_distance(seed.loc, o.loc_pre);
			double d2 = metric_distance(o.loc_copy, seed.loc);
			if (d1 > d2 + max_dis || d1 < fabs(d2 - max_dis)) {
				if (abs(o.loc_pre.x - seed.loc.x) >= pow(10, -6) && abs(o.loc_pre.y - seed.loc.y) >= pow(10, -6)) {
					intersect = getPoint(o.loc_pre.x, o.loc_pre.y, max_dis, o.loc_pre.x, o.loc_pre.y, seed.loc.x, seed.loc.y, 0);
					flag = 0;
				}
				else {
					intersect = getPoint(o.loc_pre.x, o.loc_pre.y, max_dis, o.loc_pre.x, o.loc_pre.y, o.loc_copy.x, o.loc_copy.y, 1);
					flag = 2;
				}

			}
			else {

				circles[0].center = seed.loc;
				circles[0].r = d2;

				circles[1].center = o.loc_pre;
				circles[1].r = max_dis;

				struct point points[2];
				intersect = get_circles_interscts(circles, points);
				if (intersect[0].x == INT_MAX || intersect[0].x == INT_MIN) {
					intersect = getPoint(o.loc_pre.x, o.loc_pre.y, max_dis, o.loc_pre.x, o.loc_pre.y, seed.loc.x, seed.loc.y, 0);
				}
				flag = 1;
			}

			vector<point> tmp(intersect);
			bool change = constrain_by_border(intersect);

			if (intersect.size() > 1) {
				double a = metric_distance(o.loc_copy, intersect[0]);
				double b = metric_distance(o.loc_copy, intersect[1]);
				if (flag == 1) {
					if (a < b) o.loc = intersect[0];
					else o.loc = intersect[1];
				}
				else {
					if (d1 > d2 + max_dis) {
						if (a > b) o.loc = intersect[0];
						else o.loc = intersect[1];
					}
					else {
						if (a < b) o.loc = intersect[0];
						else o.loc = intersect[1];
					}
				}
			}
			change = false;

		}

	}
}
/*
functions "compute_f_value",  "return_f_value", and "return_ff_value" compute the minimal value of f_k(\tilde(s)) 
and the correpsonding set of b_opt (cf. Formulas 7-10)
*/
double compute_f_value(double b, double d, double eta) {
	double f = 1.0 * pow(eta * (b - 1), 2);

	f += 1.0 * alpha * pow(d - b * eta, 2);
	return f;
}

double compute_ff_value(double b, double d, double eta) {
	double f = 1.0 * pow(eta * (ceil(b) - 1), 2);

	f += 1.0 * alpha * pow(d - b * eta, 2);
	return f;
}

solution return_f_value(vector<vehicle>& o, vehicle seed, double F, double eta) {
	double f = 0;
	double b;

	solution s;
	for (int i = 0; i < o.size(); i++) {
		if (o[i].id == seed.id) {
			s.b.push_back(-1);
			continue;
		}
		else {

			double d = 1.0 * metric_distance(o[i].loc, seed.loc);

			if (o[i].loc.x < 0 || o[i].loc.y < 0 || seed.loc.x < 0 || seed.loc.y < 0) {
				cout << "error in point locations" << endl;
				getchar();
			}

			if (d <= eta) {
				s.b.push_back(-1);
				continue;
			}
			else {
				b = 1.0 * (1.0 * d / eta * alpha + 1) / (1 + alpha);

				int b_max = ceil(b);
				int b_min = floor(b);
				if (abs(b - b_max) > abs(b - b_min)) b = b_min;
				else b = b_max;

				double d2 = 1.0 * b * eta;
				double ratio = 1.0 * d2 / d;
				double x = ratio * (o[i].loc.x - seed.loc.x) + seed.loc.x;
				double y = ratio * (o[i].loc.y - seed.loc.y) + seed.loc.y;
				point r;
				r.x = x;
				r.y = y;

				if (metric_distance(r, o[i].loc_pre) <= mu * o[i].delta_t) {}
				else {

					int b_opt1 = -1, b_opt2 = -1;
					vector<point> p(getPoint(o[i].loc_pre.x, o[i].loc_pre.y, mu * o[i].delta_t, o[i].loc.x, o[i].loc.y, seed.loc.x, seed.loc.y, 2));
					/*
					approximately compute r(o) on the line segment se(\tilde{s},o)
					*/

					if (p[0].x >= min(o[i].loc.x, seed.loc.x) && p[0].x <= max(o[i].loc.x, seed.loc.x) &&
						p[0].y >= min(o[i].loc.y, seed.loc.y) && p[0].y <= max(o[i].loc.y, seed.loc.y)) {
						int bb = floor(1.0 * metric_distance(p[0], seed.loc) / eta);
						if (abs(bb - b) > abs(bb + 1 - b)) {
							b_opt1 = bb + 1;
						}
						else b_opt1 = bb;
					}

					if (p[1].x >= min(o[i].loc.x, seed.loc.x) && p[1].x <= max(o[i].loc.x, seed.loc.x) &&
						p[1].y >= min(o[i].loc.y, seed.loc.y) && p[1].y <= max(o[i].loc.y, seed.loc.y)) {
						int bb = floor(1.0 * metric_distance(p[0], seed.loc) / eta);
						if (abs(bb - b) > abs(bb + 1 - b)) {
							b_opt2 = bb + 1;
						}
						else b_opt2 = bb;
					}
					if (b_opt1 != b_opt2 && b_opt1 > 0 && b_opt2 > 0) {
						if (abs(b - b_opt1) < abs(b - b_opt2)) b = b_opt1;
						else b = b_opt2;
					}
					else if (b_opt1 > 0) {
						b = b_opt1;
					}
					else if (b_opt2 > 0) {
						b = b_opt2;
					}
					else {
						b = 1.0 * d / eta;
					}
				}
				s.b.push_back(b);
				double f_min;
				if (b == 1.0 * d / eta) f_min = compute_ff_value(b, d, eta);
				else f_min = compute_f_value(b, d, eta);
				double ff =
					f += f_min;
				if (f >= F) {
					s.f = 1.0 * INT_MAX;
					return s;
				}
			}
		}
	}
	s.f = f;

	return s;
}


/*
map the clusters during consecutive time steps and calculate NMI
*/
vector<double> cluster_mapping(vector<vector<int>>& c1, vector<vector<int>>& c2, vector<vehicle>& o, map<int, int>& r) {
	double n = 0.0;
	/*
	calculate NMI
	*/
	vector<vector<int>> matrix(c1.size());
	for (int i = 0; i < matrix.size(); i++) {
		matrix[i].resize(c2.size());
	}

	for (int i = 0; i < c1.size(); i++) {
		set<int> a;
		for (int k = 0; k < c1[i].size(); k++) {
			a.insert(c1[i][k]);

		}
		vector<int> cij;
		for (int j = 0; j < c2.size(); j++) {

			set<int> b;
			for (int k = 0; k < c2[j].size(); k++) {
				b.insert(c2[j][k]);
			}
			set<int> r3;
			set_intersection(a.begin(), a.end(), b.begin(), b.end(), inserter(r3, r3.begin()));
			matrix[i][j] = r3.size();

			n += r3.size();
		}
	}

	vector<int> row;
	double a = 0.0;
	for (int i = 0; i < matrix.size(); i++) {
		double tmp = 0.0;
		for (int j = 0; j < matrix[i].size(); j++) {
			tmp += matrix[i][j];
		}
		row.push_back(tmp);
		if (tmp == 0.0) a += 0.0;
		else {
			if (tmp == n) {
				a += 1.0 * tmp * log(1.0 * tmp / (n + 1));
			}
			else a += 1.0 * tmp * log(1.0 * tmp / n);
		}

	}

	vector<int> column;
	double b = 0.0;
	for (int i = 0; i < c2.size(); i++) {
		double tmp = 0.0;
		for (int j = 0; j < c1.size(); j++) {
			tmp += matrix[j][i];
		}
		column.push_back(tmp);
		if (tmp == 0.0) b += 0.0;
		else {
			if (tmp == n) {
				b += 1.0 * tmp * log(1.0 * tmp / (n + 1));
			}
			else b += 1.0 * tmp * log(1.0 * tmp / n);
		}
	}

	double c = 0.0;

	vector<double> u;
	vector<vector<double>> unit;
	vector<vector<bool>> flag;
	for (int i = 0; i < c1.size(); i++) {
		vector<double> tmp;
		vector<bool> tmp_f;
		for (int j = 0; j < c2.size(); j++) {
			double nmi;
			if (matrix[i][j] == 0 || row[i] == 0 || column[j] == 0) {
				nmi = 0.0;
			}
			else {
				nmi = 1.0 * matrix[i][j] * log(matrix[i][j] * n / (row[i] * column[j]));
			}
			c += nmi;

			tmp.push_back(nmi);
			u.push_back(nmi);
			tmp_f.push_back(false);
		}
		unit.push_back(tmp);
		flag.push_back(tmp_f);
	}
	/*
	map clusters
	*/
	vector<double> u_copy(u);
	std::sort(u.begin(), u.end());

	vector<bool> mm(c2.size());
	for (int i = 0; i < c2.size(); i++) mm[i] = false;

	for (int i = u.size() - 1; i >= 0; i--) {

		if (u[i] <= 0) break;
		auto it = find(u_copy.begin(), u_copy.end(), u[i]);
		int dist = distance(u_copy.begin(), it);

		int id1 = floor(1.0 * dist / c2.size());
		int id2 = dist - c2.size() * id1;

		if (!flag[id1][0] && !flag[0][id2]) {

			double tmp = 0.0;
			for (int p = 0; p < unit.size(); p++) {

				tmp += unit[p][id2];
			}

			if (tmp > 2 * unit[id1][id2]) continue;

			for (int p = 0; p < c2.size(); p++) {

				tmp += unit[id1][p];
			}
			if (tmp > 2 * unit[id1][id2]) continue;
			else
			{
				r.insert(make_pair(id1, id2));
				mm[id2] = true;
				for (int p = 0; p < flag.size(); p++) {
					flag[p][id2] = true;
					unit[p][id2] = 0;
				}
				for (int p = 0; p < c2.size(); p++) {
					flag[id1][p] = true;
					unit[id1][p] = 0;
				}
			}

		}
	}

	double NMI;
	if (c == 0 || (a == 0 && b == 0)) {
		NMI = 0;
	}
	else  NMI = c / sqrt(a * b);

	vector<double> result;
	result.push_back(NMI);
	return result;
}

