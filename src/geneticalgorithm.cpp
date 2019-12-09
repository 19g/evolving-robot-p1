#include <iostream>
#include <cmath>
#include <thread>
#include <random>
#include <vector>
#include <fstream>
#include "simulate.hpp"
#include "geneticalgorithm.hpp"

using namespace std;

int main() {
    srand(time(nullptr));

    // use threads for performance
    /*
    thread t[NUM_OF_TRIALS];
    for (int i = 0; i < NUM_OF_TRIALS; i++) {
        t[i] = thread(loop, i);
    }
    for (int i = 0; i < NUM_OF_TRIALS; i++) {
        t[i].join();
    }*/

    loop(0);

    return 0;
}

void loop(int thread_num) {
    // begin timer
    clock_t begin = clock();

    // initialize files
//    ofstream diversity_file;
//    diversity_file.open(DIVERSITY_TXT);
    ofstream learning_file;
    learning_file.open(to_string(thread_num) + LEARNING_TXT);
    ofstream robot_parameters_file;
    robot_parameters_file.open(to_string(thread_num) + "robot_parameters.txt");

    // setup gpu:
    Gpu_info* gi = (Gpu_info*)malloc(sizeof(Gpu_info));
    setup_gpu(gi);

    // initialize parent population randomly
    vector<Cube> parent(POP_SIZE);
    for (int i = 0; i < POP_SIZE; i++) {
        parent[i] = initialize_cube();
        simulation_loop(parent[i], thread_num, false, gi);
    }

//    for (int i = 0; i < POP_SIZE; i++) {
//        print_mass(parent[i].mass[0]);
//    }
    // evolutionary loop
    for (int eval = 0; eval <= NUM_OF_EVALS; eval+=POP_SIZE) {
        // get random order of individuals for crossover
        vector<int> order = randomize_array_of_springs();
        // initialize offspring
        vector<Cube> child(parent);

        cout << "\n";
        // crossover
        for (int i = 0; i < POP_SIZE; i += 2) {
            crossover(child[order[i]], child[order[i+1]]);
        }

        // mutation
        for (int i = 0; i < POP_SIZE; i++) {
            mutation(child[i]);
        }

        // get fitness of population
        for (int i = 0; i < POP_SIZE; i++) {
            simulation_loop(child[i], thread_num, false, gi);
        }

        // selection
        vector<Cube> all(POP_SIZE * 2);
        for (int i = 0; i < POP_SIZE; i++) {
            all[i] = parent[i];
        }
        for (int i = POP_SIZE; i < POP_SIZE * 2; i++) {
            all[i] = child[i - POP_SIZE];

        }
        tournament_selection(parent, child, all);

        // write diversity to a file
//        calculate_diversity(parent, diversity_file);

        // find most fit individual
        int max_fit_index = 0;
        for (int i = 0; i < POP_SIZE; i++) {
            if (parent[i].fitness > parent[max_fit_index].fitness) {
                max_fit_index = i;
            }
        }

        // write learning curve to file
        for (int i = 0; i < POP_SIZE; i++) {
            learning_file << parent[max_fit_index].distance << ",";
        }

        // print fitnesses of population
        if (eval % POP_SIZE == 0) {
            for (int i = 0; i < POP_SIZE; i++) {
                cout << eval << ": " << parent[i].distance << "\n";
            }
        }
    }

    // end timer
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    double iters_per_sec = NUM_OF_ITERATIONS / elapsed_secs;
    cout << "iter/sec: " << iters_per_sec << "\n";

    int max_fit_index = 0;
    for (int i = 0; i < POP_SIZE; i++) {
        if (parent[i].fitness > parent[max_fit_index].fitness) {
            max_fit_index = i;
        }
    }

    // output robot parameters to file
    for (int i = 0; i < NUM_OF_SPRINGS; i++) {
        robot_parameters_file << i << ":\n";
        robot_parameters_file << "k: " << parent[max_fit_index].spring[i].k << "\n";
        robot_parameters_file << "b: " << parent[max_fit_index].spring[i].b << "\n";
        robot_parameters_file << "c: " << parent[max_fit_index].spring[i].c << "\n\n";
    }

    cout << "MAX FITNESS: " << parent[max_fit_index].fitness << "\n";

    // output to file for opengl
    simulation_loop(parent[max_fit_index], thread_num, true, gi);

    // close files
    learning_file.close();
//    diversity_file.close();
    robot_parameters_file.close();
}

void setup_gpu(Gpu_info* info) {
    //info = malloc(sizeof(Gpu_info));
    clGetPlatformIDs(1, &(info->platform), NULL);

    info->err = clGetDeviceIDs(info->platform,
                               CL_DEVICE_TYPE_GPU, 
                               1, 
                               &(info->device_id), 
                               NULL
                              );
    cout << "err after get device ids: " << info->err << endl;

    info->context = clCreateContext(0,
                                    1, 
                                    &(info->device_id),
                                    NULL, 
                                    NULL, 
                                    &(info->err)
                                   );
    cout << "err after create context: " << info->err << endl;

    info->commands = clCreateCommandQueue(info->context,
                                          info->device_id, 
                                          0, 
                                          &(info->err)
                                         );

    cout << "err after create command: " << info->err << endl;

    // read program from file and convert to string:
    char* kernel_source = (char*)malloc(sizeof(char)*2048);
    parse_file_into_str("breathing_cube.cl", 
                        kernel_source, 
                        sizeof(char)*2048
                       );

    info->program = clCreateProgramWithSource(info->context, 
                                              1, 
                                              (const char**)&kernel_source, 
                                              NULL, 
                                              &(info->err)
                                             );
    cout << "err after create program: " << info->err << endl;
    free(kernel_source);
    info->err = clBuildProgram(info->program, 
                               0, 
                               NULL, 
                               "-cl-finite-math-only", 
                               NULL, 
                               NULL
                              );

    cout << "err after build program: " << info->err << endl;

    if (info->err != 0) {
		// Determine the size of the log
    	size_t log_size;
    	clGetProgramBuildInfo(info->program, info->device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

    	// Allocate memory for the log
    	char *log = (char *) malloc(log_size);

    	// Get the log
    	clGetProgramBuildInfo(info->program, info->device_id, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);

    	// Print the log
    	printf("%s\n", log);
		exit(1);
    }


    info->kernel_bc = clCreateKernel(info->program, 
                                     "breathing_cube", 
                                     &(info->err)
                                    );
    cout << "err after create kernel: " << info->err << endl;


    // cl_mem and friends are typedef pointers to structs
    info->input_bc = clCreateBuffer(info->context,  
                                    CL_MEM_READ_WRITE,  
                                    sizeof(Spring) * NUM_OF_SPRINGS, 
                                    NULL, 
                                    NULL
                                   );
    cout << "err after create buffer: " << info->err << endl;

    /*
    info->output_bc = clCreateBuffer(info->context, 
                                     CL_MEM_WRITE_ONLY, 
                                     sizeof(float) * NUM_OF_SPRINGS, 
                                     NULL, 
                                     NULL
                                    );*/
}

Cube initialize_cube() {
    // random variable generation
    random_device rd;
    mt19937 mt(rd());
//    uniform_real_distribution<double> a_val(MIN_A, MAX_A);
    uniform_real_distribution<float> b_val(MIN_B, MAX_B);
    uniform_real_distribution<float> c_val(MIN_C, MAX_C);
    uniform_real_distribution<double> k_spring_val(MIN_K_SPRING, MAX_K_SPRING);

    // temp objects
    vector<Mass> mass;
    vector<Spring> spring;

    // create masses
    Mass temp_mass = {WEIGHT_PER_MASS,
                      {0.0, 0.0, 0.0},
                      {0.0, 0.0, 0.0},
                      {0.0, 0.0, 0.0}};
    temp_mass.p = {0.0, 0.0, INITIAL_HEIGHT}; // 0
    mass.emplace_back(temp_mass);
    temp_mass.p = {L0_SIDE, 0.0, INITIAL_HEIGHT}; // 1
    mass.emplace_back(temp_mass);
    temp_mass.v = {0.0, 0.0, 0.0}; // 0
    temp_mass.p = {L0_SIDE, L0_SIDE, INITIAL_HEIGHT}; // 2
    mass.emplace_back(temp_mass);
    temp_mass.p = {0.0, L0_SIDE, INITIAL_HEIGHT}; // 3
    mass.emplace_back(temp_mass);
    temp_mass.p = {0.0, 0.0, INITIAL_HEIGHT + L0_SIDE}; // 4
    mass.emplace_back(temp_mass);
    temp_mass.p = {L0_SIDE, 0.0, INITIAL_HEIGHT + L0_SIDE}; // 5
    mass.emplace_back(temp_mass);
    temp_mass.p = {L0_SIDE, L0_SIDE, INITIAL_HEIGHT + L0_SIDE}; // 6
    mass.emplace_back(temp_mass);
    temp_mass.p = {0.0, L0_SIDE, INITIAL_HEIGHT + L0_SIDE}; // 7
    mass.emplace_back(temp_mass);

    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[1].p), 0, 1, dist(mass[0].p, mass[1].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[2].p), 0, 2, dist(mass[0].p, mass[2].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[3].p), 0, 3, dist(mass[0].p, mass[3].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[4].p), 0, 4, dist(mass[0].p, mass[4].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[5].p), 0, 5, dist(mass[0].p, mass[5].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[6].p), 0, 6, dist(mass[0].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[7].p), 0, 7, dist(mass[0].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[2].p), 1, 2, dist(mass[1].p, mass[2].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[3].p), 1, 3, dist(mass[1].p, mass[3].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[4].p), 1, 4, dist(mass[1].p, mass[4].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[5].p), 1, 5, dist(mass[1].p, mass[5].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[6].p), 1, 6, dist(mass[1].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[7].p), 1, 7, dist(mass[1].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[3].p), 2, 3, dist(mass[2].p, mass[3].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[4].p), 2, 4, dist(mass[2].p, mass[4].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[5].p), 2, 5, dist(mass[2].p, mass[5].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[6].p), 2, 6, dist(mass[2].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[7].p), 2, 7, dist(mass[2].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[4].p), 3, 4, dist(mass[3].p, mass[4].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[5].p), 3, 5, dist(mass[3].p, mass[5].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[6].p), 3, 6, dist(mass[3].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[7].p), 3, 7, dist(mass[3].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[4].p, mass[5].p), 4, 5, dist(mass[4].p, mass[5].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[4].p, mass[6].p), 4, 6, dist(mass[4].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[4].p, mass[7].p), 4, 7, dist(mass[4].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[5].p, mass[6].p), 5, 6, dist(mass[5].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[5].p, mass[7].p), 5, 7, dist(mass[5].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[6].p, mass[7].p), 6, 7, dist(mass[6].p, mass[7].p), b_val(mt), c_val(mt)});

    Cube cube = {mass, spring, 0};

    return cube;
}

vector<int> randomize_array_of_springs() {
    // create vector with order of parent to be crossed over
    vector<int> order(POP_SIZE);
    for (int i = 0; i < POP_SIZE; i++) {
        order[i] = i;
    }
    for (int i = order.size() - 1; i > 0; i--) {
        int j = (int)(rand() % (i+1));
        int temp = order[i];
        order[i] = order[j];
        order[j] = temp;
    }

    return order;
}

void crossover(Cube &A, Cube &B) {
    // choose random points for crossover
    int start = rand() % NUM_OF_SPRINGS;
    int end = rand() % NUM_OF_SPRINGS;
    while (end == start) {
        end = rand() % NUM_OF_SPRINGS;
    }

    // swap the spring values at that point
    for (int i = start; i < end; i++) {
        Spring temp = A.spring[i];
        A.spring[i] = B.spring[i];
        B.spring[i] = temp;
    }
}

void mutation(Cube &individual) {
    // random variable generation
    // TODO: move this outside of the mutation function becasue it's probably pretty slow
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> mut_chance(0, 1);

    for (int i = 0; i < NUM_OF_SPRINGS; i++) {
        if (mut_chance(mt) < PROB_OF_MUT) {
            // pick a spring and multiply its values by numbers between MIN_SWING and MAX_SWING
            uniform_int_distribution<> spring(0, NUM_OF_SPRINGS);
            uniform_real_distribution<double> swing(MIN_SWING, MAX_SWING);
//            individual.spring[spring(mt)].a =  individual.spring[spring(mt)].a * swing(mt);
            if (mut_chance(mt) < PROB_PER_PARAM) {
                uniform_real_distribution<double> b_val(0, individual.spring[i].a/2);
                individual.spring[spring(mt)].b = b_val(mt);
            }
            if (mut_chance(mt) < PROB_PER_PARAM) {
                uniform_real_distribution<double> c_val(MIN_C, MAX_C);
                individual.spring[spring(mt)].c = c_val(mt);
            }
            if (mut_chance(mt) < PROB_PER_PARAM) {
                uniform_real_distribution<double> k_val(MIN_K_SPRING, MAX_K_SPRING);
                individual.spring[spring(mt)].k = k_val(mt);
            }
        }
    }
}

void tournament_selection(vector<Cube> &parent, vector<Cube> &child, vector<Cube> &all) {
    // random variable generation
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<> compare(0, POP_SIZE * 2 - 1);

    // take elite child;
    int best_parent_index = 0;
    int best_child_index = 0;
    for (int i = 0; i < POP_SIZE; i++) {
        if (parent[i].fitness > parent[best_parent_index].fitness) {
            best_parent_index = i;
        }
        if (child[i].fitness > child[best_child_index].fitness) {
            best_child_index = i;
        }
    }
    parent[0] = parent[best_parent_index];
    parent[1] = child[best_child_index];

    for (int i = 2; i < POP_SIZE; i++) {
        int m = compare(mt);
        int n = compare(mt);
        while (m == n) {
            n = compare(mt);
        }

        if (all[m].fitness > all[n].fitness) {
            parent[i] = all[m];
        } else {
            parent[i] = all[n];
        }
    }
}

void calculate_diversity(vector<Cube> &population, ofstream &diversity_file) {
    double avg_k;
    double avg_b;
    double avg_c;
    double mse_k;
    double mse_b;
    double mse_c;
    double sum_k;
    double sum_b;
    double sum_c;
    // get avg value of parameters
    for (int i = 0; i < POP_SIZE; i++) {
        for (int j = 0; j < NUM_OF_SPRINGS; j++) {
            avg_k += population[i].spring[j].k;
            avg_b += population[i].spring[j].b;
            avg_c += population[i].spring[j].c;
        }
    }
    avg_k /= POP_SIZE;
    avg_b /= POP_SIZE;
    avg_c /= POP_SIZE;

    // calculate diversity form these averages
    for (int i = 0; i < POP_SIZE; i++) {
        for (int j = 0; j < NUM_OF_SPRINGS; j++) {
            sum_k += population[i].spring[j].k;
            sum_b += population[i].spring[j].b;
            sum_c += population[i].spring[j].c;
        }
        mse_k += pow(sum_k - avg_k, 2);
        mse_b += pow(sum_b - avg_b, 2);
        mse_c += pow(sum_c - avg_c, 2);
    }

    mse_k /= POP_SIZE;
    // scale up b and c so they have a meaningful contribution to diversity
    mse_b = (mse_b / POP_SIZE) * (MAX_K_SPRING - MIN_K_SPRING) / (MAX_B - MIN_B);
    mse_c =(mse_c / POP_SIZE) * (MAX_K_SPRING - MIN_K_SPRING) / (MAX_C- MIN_C);

    // write diversity to a file
    diversity_file << mse_k + mse_b + mse_c << ",";
}

bool parse_file_into_str(const char* filename, char* shader_str, int max_len) {
    FILE* fp;
    size_t cnt;
    fp = fopen(filename, "r");
    if (!fp) {
        cout << "COULD NOT OPEN FILE" << endl;
        return false;
    }

    cnt = fread(shader_str, sizeof(char), max_len-1, fp);
    if ((int)cnt >= max_len-1) {
        cout << "FILE WAS TRUNCATED" << endl;
    }

    if (ferror(fp)) {
        cout << "FILE ERROR" << endl;
        fclose(fp);
        return false;
    }

    // append \0 to end of file converted to string
    shader_str[cnt] = '\0';
    fclose(fp);
    return true;
}
