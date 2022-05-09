

#include "get_input.h"

const int a = 16;
static int id_index[vehicle_size];


int main()
{
	for (int i = 0; i < vehicle_size; i++) id_index[i] = -1;

	vector<string> str(a);
	vector<vehicle> obj;

	int xx = 0;
	
	vector<vector<int>> cl_pre, cl, cl_non_optimal, cl_non_optimal_pre;

	int s = 0;

	double qs = 0.0, qs_non_optimal = 0.0;
	double nmi = 0.0, nmi_non_optimal = 0.0;
	multimap<string, int> grid_index;
	multimap<string, int> grid_index_seed;
	multimap<string, int> veh;
	multimap<string, int> veh_eta;
	map<int, int> gid;
	vector<C_unit> C;
	auto start = chrono::high_resolution_clock::now();

	double total_qs = 0.0, total_qs_non_optimal=0.0;
	double total_nmi = 0.0, total_nmi_non_optimal = 0.0;
	double total_epsilon = 0.0;
	long long int total_time = 0;

	map<string, double> grid_eta;
	static double epsilon = 1200; /* initial value of \varepsilon */ 
	map<int, int> m1;

	int FILE_NUM = 1; 
	int file = 1;
	vector<int> detect;
	int temp = 0;
	int test = 0;
	ofstream out;
	while (file <= FILE_NUM) {
		ifstream in;
		in.open(".\\data\\CD_example.txt", ios::in); /* file route */
		file++;
		while (!in.eof()) {

			for (int i = 0; i < a; i++) {
				in >> str[i];			
			}
			test++;

			if (str[0] == "=") {
				if (temp == xx) {		
					break;
				}
				temp=xx;
			
				if (s > 0) {
					/* 
					smoothing 
					*/
					smooth(obj, C, epsilon, veh); 
					/*
					clear the minimal groups generated at the last time step
					*/
					reset_C_unit(obj, C);
				}
				/*
				generate the minimal groups at the current time step
				*/
				get_C_unit(obj, grid_index, grid_index_seed, veh, C, epsilon, delta);

				/*
				constrain \delta <= \varepsilon 
				*/
				if (epsilon < delta) epsilon=delta;
				if (s == 0)
				{
					/*
					iteratively optimizing \varepsilon only at the first time step
					*/
					cl = optimal(obj, epsilon, grid_index_seed, grid_index, C);
					
				}
				else
				{
					/*
					clustering using DBSCAN
					*/
					
					cl = clustering(grid_index_seed, grid_index, C, obj, epsilon);
				

					/*
					calculate QS_m, i.e., QS
					*/
					qs = get_fitness_qs(cl, obj);
					
					total_qs += qs;

					/*
				    calculate QS_h 
				    */
					double q_low = get_high_qs(cl, obj, epsilon);

					/*
				    calculate QS_l
				    */
				    double q_high = get_low_qs(cl, obj, epsilon);
				    
					if (qs > q_high&& qs > q_low) {}

					else {
						if (q_high > q_low) {
							epsilon -= delta_epsilon;
						}
						else epsilon += delta_epsilon;
					}
			
				}
				
				vector<double> result;
				if (s > 0) {
					/*
					map clusters and compute NMI
					*/
					result = cluster_mapping(cl, cl_pre, obj, m1);

					nmi = result[0];

					total_nmi += nmi;

				}
				
				auto end = chrono::high_resolution_clock::now();
				std::chrono::duration<double> diff = end - start;
				total_time = diff.count();

				m1.erase(m1.begin(), m1.end());		
				cl_pre = cl;

				s++;

				reset_stop_flag(obj);
				clean_grid_index(grid_index);
				clean_grid_index(grid_index_seed);
				clean_grid_index(veh);
			
			}
		
			else {
				
				/*
				process each arriving location
				*/			
				get_new_data(str, obj, id_index, veh, epsilon);
				xx++;

			}

		}	
	}
	std::cout << "================" << endl;
	std::cout << "time step: " << s-1 << endl;
	std::cout << "average processing time: " << 1.0 * total_time / xx << endl;
	std::cout << "QS: " << 1.0 * total_qs / (s-1) << endl;
	std::cout << "NMI: " << 1.0 * total_nmi / (s-1) << endl;
	std::cout << "================" << endl;

}



