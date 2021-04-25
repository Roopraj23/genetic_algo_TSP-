#include <bits/stdc++.h>
#include <limits.h>
using namespace std;

// Number of cities in TSP
#define V 5

// Initial population size for the algorithm
#define POP_SIZE 20

// Structure of a GNOME
// string defines the path traversed
// by the salesman while the fitness value
// of the path is stored in an integer

struct individual {
	string gnome;
	double fitness;
};

// Function to return a random number
// from start and end
int rand_num(int start, int end)
{
	int r = end - start;
	int rnum = start + rand() % r;
	return rnum;
}

bool repeat(string s, char ch) // To check weather city has already visited
{
	for (int i = 0; i < s.size(); i++) {
		if (s[i] == ch)
			return true;
	}
	return false;
}

// Function to return a mutated GNOME
string mutatedGene(string gnome)
{
	while (true) {
		int r = rand_num(1, V);
		int r1 = rand_num(1, V);
		if (r1 != r) {
			char temp = gnome[r];
			gnome[r] = gnome[r1];
			gnome[r1] = temp;
			break;
		}
	}
	return gnome;
}

// required to create the population
string create_gnome() //Random path
{
	string gnome = "0"; //Starting City
	while (true) {
		if (gnome.size() == V) {
			gnome += gnome[0];
			break;
		}
		int temp = rand_num(1, V);
		if (!repeat(gnome, (char)(temp + 48)))
			gnome += (char)(temp + 48);
	}
	return gnome;
}

// Function to return the fitness value of a gnome.
double cal_fitness(string gnome)
{
       double map[V][V] = 		{{0, 5.882329, 5.421476, 3.348194, 10.966859 },
					{ 5.882329, 0, 1.291898, 4.487427, 16.7559 },
					{ 5.421476, 1.291898, 0, 4.830114, 16.388252 },
					{ 3.348194, 4.487427, 4.830114, 0, 12.883509 },
					{ 10.966859, 16.7559, 16.388252, 12.883509, 0 } };
	double f = 0;
	for (int i = 0; i < gnome.size() - 1; i++) {
		if (map[gnome[i] - 48][gnome[i + 1] - 48] == INT_MAX)
			return INT_MAX;
		f += map[gnome[i] - 48][gnome[i + 1] - 48];
	}
	return f;
}

// Function to return the updated value
// of the cooling element.
int cooldown(int temp)
{
	return (90 * temp) / 100;
}

// Comparator for GNOME struct.
bool lessthan(struct individual t1,
			struct individual t2)
{
	return t1.fitness < t2.fitness;
}

// Utility function for TSP problem.
void TSPUtil(double map[V][V])
{
	// Generation Number
	int gen = 1;
	// Number of Gene Iterations
	int gen_thres = 20;

	vector<struct individual> population;
	struct individual temp;

	// Populating the GNOME pool.
	for (int i = 0; i < POP_SIZE; i++) {
		temp.gnome = create_gnome();
		temp.fitness = cal_fitness(temp.gnome);
		population.push_back(temp);
	}

	cout << "\nInitial population: " << endl << "GNOME	 FITNESS VALUE\n";
	for (int i = 0; i < POP_SIZE; i++)
		cout << population[i].gnome << " " << population[i].fitness << endl;
	cout << "\n";

	//bool found = false;
	int temperature = 10000;


	// Total time calculation start for all iterations taken together
	clock_t begin_all = clock();
	// Iteration to perform
	// population crossing and gene mutation.
	while (temperature > 100 && gen <= gen_thres)
    {
		sort(population.begin(), population.end(), lessthan);
		
		cout << "\nCurrent temp: " << temperature << "\n";

		vector<struct individual> new_population;
		// start time calculation
		clock_t begin = clock();
		for (int i = 0; i < POP_SIZE; i++) {
			struct individual p1 = population[i];

			while (true) {
				string new_g = mutatedGene(p1.gnome);
				struct individual new_gnome;
				new_gnome.gnome = new_g;
				new_gnome.fitness = cal_fitness(new_gnome.gnome);

				if (new_gnome.fitness <= population[i].fitness) {
					new_population.push_back(new_gnome);
					break;
				}
				else {
					// Accepting the rejected children at
					// a possible probablity above threshold.
					float prob = pow(2.7, -1 * ((float)(new_gnome.fitness - population[i].fitness) / temperature));
					if (prob > 0.5) {
						new_population.push_back(new_gnome);
						break;
					}
				}
			}
		}
		
		clock_t end = clock();
		double time_spent = (double)(end-begin); 	// clocks per sec
		cout<<"Total time spent: "<<time_spent<<" \n";

		temperature = cooldown(temperature);
		population = new_population;
		cout << "Generation " << gen << " \n";
		cout << "GNOME	 FITNESS VALUE\n";
        
		double higest = 100000;
		for (int i = 0; i < POP_SIZE; i++){
			cout << population[i].gnome << " "<< population[i].fitness << endl;
		    if(population[i].fitness < higest){
				higest = population[i].fitness;
			}
		}	
        cout<<"\nBest path legnth is :" << higest;

		gen++;
	}
	// Total time calculation end for all iterations taken together
	clock_t end_all = clock();
	double time_spent_all = (double)(end_all-begin_all); 	// clocks per sec
  	cout<<endl;
	cout<<"Total time spent all: "<<time_spent_all/3600<<" \n";
	
}

int main()
{

	double map[V][V] = 		{{0, 5.882329, 5.421476, 3.348194, 10.966859 },
					{ 5.882329, 0, 1.291898, 4.487427, 16.7559 },
					{ 5.421476, 1.291898, 0, 4.830114, 16.388252 },
					{ 3.348194, 4.487427, 4.830114, 0, 12.883509 },
					{ 10.966859, 16.7559, 16.388252, 12.883509, 0 } };
	TSPUtil(map);
}

