#include <bits/stdc++.h>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <climits>
#include <bitset>
#include <set>
#include <sys/time.h>
#include <omp.h>

using namespace std;

	int str2int(string str)
	{
		int value;
		istringstream (str)>>value;	
		return value;
	}

	
	string input_file="";
	string testcase_name="";
	string output_file="";
	int dimensions;
	int path_in_one_generation;
	int threads_number;

	double crossover_probability=0.8;
	double mutation_probability=0.1;

	vector< vector 	<double> > distance_bw;

	struct dna      // datastructure to define a city
	{
		string id;
		double x_coordinate;
		double y_coordinate;
	};


	struct trip		//datastructure to define a path between cities
	{
		vector<int> list;
		long double fitness;
	};

	vector<dna> gnome;//vector of struct dna
	
	vector<trip> new_population;//current generation

	std::vector<trip> shortest_fitness;

	trip overall_shortest_fitness;  //overall shortest path


	//function  to read the input file
	void readfile()
	{
		ifstream infile;
		char input_file_array[input_file.length()+1];
		
		for(int i=0;i<input_file.length();i++)
		{
			input_file_array[i] = input_file.at(i);
		}
		
		input_file_array[input_file.length()]='\0';

		infile.open(input_file_array);


		//TestCase Name Line
		string name_line;
		getline(infile,name_line);
		// print_to_err(name_line);
		int name_line_length=int(name_line.length());
		int start_name_line;
		for(int i=5;i<name_line_length;i++)
		{
			if(name_line.at(i)!=' ')
			{
				start_name_line=i;
				break;
			}
		}
		testcase_name=name_line.substr(start_name_line,(name_line_length- start_name_line));

	//Dimension Line
	string dimension_line;
	getline(infile,dimension_line);
	// print_to_err(dimension_line);
	int dimension_line_length=int(dimension_line.length());
	int start_dimension_line;
	for(int i=10;i<dimension_line_length;i++)
	{
		if(dimension_line.at(i)!=' ')
		{
			start_dimension_line=i;
			break;
		}
	}
	string dimension_string=dimension_line.substr(start_dimension_line,	(dimension_line_length- start_dimension_line));
	
	// print_to_err(dimension_string);
	dimensions=str2int(dimension_string);

	//NODE_COORD_SECTION garbage line
	string garbage_line;
	getline(infile,garbage_line);
	// print_to_err(garbage_line);

	//Below lines
	for(int i=0;i<dimensions;i++)
	{
		dna new_gnome;
		infile>>new_gnome.id>>new_gnome.x_coordinate>>new_gnome.y_coordinate;
		gnome.push_back(new_gnome);
	}
	infile.close();
}


	//function to direct the output to an external file
void outfile()
	{
		ofstream outfile;

		char output_file_array[output_file.length()+1];
		for(int i=0;i<output_file.length();i++)
		{
			output_file_array[i]=output_file.at(i);
		}
		output_file_array[output_file.length()]='\0';


		outfile.open (output_file_array);


		outfile<<"DIMENSION : "<<dimensions<<endl;
		outfile<<"TOUR_LENGTH : "<<overall_shortest_fitness.fitness<<endl;
		outfile<<"TOUR_SECTION"<<endl;

		for(int i=0;i<dimensions;i++)
		{
			outfile<< gnome[ overall_shortest_fitness.list[i] ].id<<" ";
		}
		outfile<<endl;
		outfile<<"-1"<<endl;
		outfile<<"End of file";
		outfile.close();
	}
	//function to make the distance vector
	void make_distance_vector()
	{
		distance_bw.resize(dimensions, vector<double>(dimensions,0));
		for(int i=0;i<dimensions;i++)
		{
			for(int j=i;j<dimensions;j++)
			{
				distance_bw[i][j]= pow((pow((gnome[i].x_coordinate - gnome[j].x_coordinate),2)+pow((gnome[i].y_coordinate - gnome[j].y_coordinate),2)),0.5);
				distance_bw[j][i]=distance_bw[i][j]; 
			}
		}
	}


	vector<int> crossover_genes(vector<int>& p1, vector<int>& p2)
	{

		int start_swath=(dimensions-1);
		int end_swath=(0);

		while(start_swath>=end_swath)
		{
			double r1 = ((double) rand() / (RAND_MAX));
			double r2 = ((double) rand() / (RAND_MAX));
			start_swath = int(floor(r1*dimensions))%dimensions;
			end_swath = int(floor(r2*dimensions))%dimensions;
		}



		set<int> swath;
		vector<int> child(int(p1.size()),-1);
		for(int i=start_swath;i<=end_swath;i++)
		{
			child[i]=p1[i];
			swath.insert(child[i]);
		}

		for(int i=start_swath;i<=end_swath;i++)
		{
			if(swath.find(p2[i])==swath.end())
			{
				int temp_pmx=p2[i];
				int same_in_p1=p1[i];
				int index_in_p2;
				LABEL_CASEY:
				for(int j=0;j<dimensions;j++)
				{
					if(p2[j]==same_in_p1)
					{
						index_in_p2=j;
						break;
					}
				}

				if(child[index_in_p2]==(-1))
				{
					
					child[index_in_p2]=p2[i];
					swath.insert(p2[i]);
				}
				else
				{
					temp_pmx=p2[index_in_p2];
					same_in_p1=p1[index_in_p2];
					goto LABEL_CASEY;
				}
			}
		}

		for(int i=0;i<dimensions;i++)
		{
			if(child[i]==(-1))
			{
				child[i]=p2[i];
			}
		}


		return child;
	}

  //greedy algorithm to make crossover to form new generations
	vector<int> greedy_crossover(vector<int>& p1, vector<int>& p2)
	{
		set<int> already_taken;
		set<int> not_taken;
		for(int i=0;i<dimensions;i++)
		{
			not_taken.insert(i);
		}
		vector<int> child(dimensions,(-1));
		child[0]=p1[0];
		not_taken.erase(child[0]);
		already_taken.insert(child[0]);

		int pos_in_p1;
		int pos_in_p2;
		bool already_p1;
		bool already_p2;
		bool random_help;
		int random_num;
		for(int i=1;i<dimensions;i++)
		{
			for(int j=0;j<dimensions;j++)
			{
				if(p1[j]==child[i-1])
				{
					pos_in_p1=j;
					break;
				}
			}
			for(int j=0;j<dimensions;j++)
			{
				if(p2[j]==child[i-1])
				{
					pos_in_p2=j;
					break;
				}
			}
		already_p1=(already_taken.find(p1[(pos_in_p1+1)%dimensions])!=already_taken.end());
		already_p2=(already_taken.find(p2[(pos_in_p2+1)%dimensions])!=already_taken.end());
		if(already_p1 && already_p2)
		{
			random_help=false;
			//select a random number from not_taken
			for(int j=0;j<dimensions;j++)
			{
				if(not_taken.find(j)!=not_taken.end())
				{
					random_num=j;
					double r3 = ((double) rand() / (RAND_MAX));
					if(r3<0.3)
					{
						random_help=true;
						child[i]=j;
						not_taken.erase(j);
						already_taken.insert(j);
						break;
					}
				}
			}

			if(!random_help)
			{
				child[i]=random_num;
				not_taken.erase(random_num);
				already_taken.insert(random_num);
			}

		}
		else if(already_p1 && !already_p2)
		{
			child[i]=p2[(pos_in_p2+1)%dimensions];
			not_taken.erase(child[i]);
			already_taken.insert(child[i]);
		}
		else if(!already_p1 && already_p2)
		{
			child[i]=p1[(pos_in_p1+1)%dimensions];
			not_taken.erase(child[i]);
			already_taken.insert(child[i]);
		}
		else if(!already_p1 && !already_p2)
		{
			if(distance_bw[child[i-1]][p1[(pos_in_p1+1)%dimensions]]<distance_bw[child[i-1]][p2[(pos_in_p2+1)%dimensions]])
			{
				child[i]=p1[(pos_in_p1+1)%dimensions];
				not_taken.erase(child[i]);
				already_taken.insert(child[i]);
			}
			else
			{
				child[i]=p2[(pos_in_p2+1)%dimensions];
				not_taken.erase(child[i]);
				already_taken.insert(child[i]);	
			}
		}		
	}

	return child;
}

//To generate initial Population
void generate_initial_population_and_calculate_distance()
{
	new_population.clear();
	new_population.resize(path_in_one_generation);

	//make permutations
	std::vector<int> initial_trip(dimensions);
	for(int i=0;i<dimensions;i++)
	{
		initial_trip[i]=i;
	}

	for(int i=0;i<path_in_one_generation;i++)
	{
		random_shuffle ( initial_trip.begin(), initial_trip.end() );

		new_population[i].list=initial_trip;
	}

	double temp_distance=0;
	for(int i=0;i<path_in_one_generation;i++)
	{
		temp_distance=0;
		for(int j=0;j<(dimensions);j++)
		{
			temp_distance+=distance_bw[new_population[i].list[j]][new_population[i].list[(j+1)%dimensions]];
		}
		new_population[i].fitness = temp_distance;
	}

}

//Function to sort according to fitness
bool sortbytriplength(const trip &a,const trip &b)
{
    return (a.fitness < b.fitness);
}

//function to sort the population
void sort_population()
{
	sort(new_population.begin(), new_population.end(),sortbytriplength);
}



int main(int argc, char const *argv[])
{
	srand (time (NULL));
	
	string temp_str(argv[1]);
	string temp_str_threads(argv[2]);
	input_file=temp_str;
	threads_number=str2int(temp_str_threads);
	readfile();
	
	output_file="output_"+testcase_name+".txt";
	make_distance_vector();

    //Population
	path_in_one_generation=20;
	
	//total generations to be made
	int total_generations=20;
	
	//temporary trips
	int temp_trips=threads_number*4;
	path_in_one_generation=(path_in_one_generation/temp_trips);
	path_in_one_generation*=temp_trips;

	int last_update_threshold=30;

	int half_trips=path_in_one_generation/2;
	int one_side=(path_in_one_generation/(2*threads_number));
	vector<int> shuffle_parent_vector(half_trips);

	for(int i=0;i<half_trips;i++)
	{
		shuffle_parent_vector[i]=i;
	}
	
	vector<trip> temp_store_parent(half_trips);

	int number_restarts=4;
	double start_time=omp_get_wtime();

	shortest_fitness.resize(number_restarts);

	for(int k=0; k<number_restarts; k++)
	{
		//Generate the initial population
		generate_initial_population_and_calculate_distance();
		sort_population();

		shortest_fitness[k].fitness=LDBL_MAX;

		if(new_population[0].fitness < shortest_fitness[k].fitness)
		{
			shortest_fitness[k]=new_population[0];
		}

		int last_time_updated=0;

		for(int i=0;i<total_generations;i++)
		{
			
			cout<<"generation number :"<<i<<endl;
			
			//parallelization started here
			#pragma omp parallel num_threads(threads_number)
			{
				int mutation_first_index;
				int mutation_second_index;
				int mutation_temp;
				int tid=omp_get_thread_num();
				//CrossOver happening here
				for(int j=(tid*one_side);j<((tid+1)*one_side);)
				{
					double r4 = ((double) rand() / (RAND_MAX));
					if(r4<crossover_probability)
					{
						new_population[half_trips+j].list = crossover_genes( new_population[j].list, new_population[j+1].list );				
					}

					double r5 = ((double) rand() / (RAND_MAX));
					if(r5<crossover_probability)
					{
						new_population[half_trips+j+1].list=greedy_crossover( new_population[j+1].list, new_population[j].list );
					}
					j+=2;
				}

				//Mutation in the gene pool
				for(int j=(half_trips+(tid*one_side));j<(half_trips+(tid+1)*one_side);j++)
				{
					double r6 = ((double) rand() / (RAND_MAX));

					if(r6 < mutation_probability)
					{
						do
						{
							double r7 = ((double) rand() / (RAND_MAX));
							double r8 = ((double) rand() / (RAND_MAX));
							mutation_first_index=int(floor(r7*dimensions))%dimensions;
							mutation_second_index=int(floor(r8*dimensions))%dimensions;
						}
						while(mutation_first_index==mutation_second_index);

						mutation_temp = new_population[j].list[mutation_first_index];
						new_population[j].list[mutation_first_index] = new_population[j].list[mutation_second_index];
						new_population[j].list[mutation_second_index] = mutation_temp;
					}
				}

				double temp_distance=0;
				for(int j=(half_trips+(tid*one_side));j<(half_trips+(tid+1)*one_side);j++)
				{
					temp_distance=0;
					for(int k=0;k<(dimensions);k++)
					{
						temp_distance+=distance_bw[new_population[j].list[k]][new_population[j].list[(k+1)%dimensions]];
					}
					new_population[j].fitness=temp_distance;
				}

				#pragma omp barrier  //so that all threads come and synchronize here

				#pragma omp single
				{
					sort_population();
					//Modify the current best trip if a better is found
					if(new_population[0].fitness < shortest_fitness[k].fitness) 
					{
						shortest_fitness[k]=new_population[0];
						last_time_updated=0;
					}
					else
					{
						last_time_updated+=1;
					}

					cout<<"Best path length is "<<(shortest_fitness[k].fitness)<<endl;
				}
			}
			if(last_time_updated >= last_update_threshold)
			{
				break;
			}
		}
	}

	overall_shortest_fitness=shortest_fitness[0]; //Final trip length 
	
	for(int k=0;k<number_restarts;k++)
	{
		if(overall_shortest_fitness.fitness > shortest_fitness[k].fitness)
		{
			overall_shortest_fitness=shortest_fitness[k];
		}
	}

	double end_time=omp_get_wtime();
    double time_taken = end_time - start_time;
    cout<<"Best path length is "<<(overall_shortest_fitness.fitness)<<endl;
    cout<<"\nTotal time spent all:"<< time_taken <<endl;
	outfile();
	return 0;
}