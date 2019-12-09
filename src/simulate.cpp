#include <iostream>
#include <fstream>
#include "simulate.hpp"
#include "geneticalgorithm.hpp"

using namespace std;


void simulation_loop(Cube &individual, 
                     int thread_num, 
                     bool opengl, 
                     Gpu_info* info) {
    // initialize files
    ofstream energy_file;
    energy_file.open(to_string(thread_num) + ENERGY_TXT);
    ofstream opengl_file;
    opengl_file.open(to_string(thread_num) + OPENGL_TXT);

    // declare variables
    float T = 0.0;
    vector<Mass> mass = individual.mass;
    vector<Spring> spring = individual.spring;
    vector<double> kinetic_energy;
    vector<double> potential_energy;
    vector<double> starting_position(DIMENSIONS);

//    for (int i = 0; i < NUM_OF_MASSES; i++) {
//        cout << mass[i].p[0] << "\t";
//        cout << mass[i].p[1] << "\t";
//        cout << mass[i].p[2] << "\t\n";
//    }

    // Calculate starting center of mass:
    vector<double> starting_com = calculate_center_of_mass(mass);
    individual.fitness = 0;

    // write to file for opengl
    if (opengl) {
        write_to_opengl_file(mass, opengl_file);
    }

    int num_iterations = NUM_OF_ITERATIONS;
    if (opengl) num_iterations *= 3;

    // simulation loop
    for (int iteration = 0; iteration < num_iterations; iteration++) {
        // initialize force vector
        vector<vector<double>> force(NUM_OF_MASSES, vector<double>(DIMENSIONS));

        // breathing cube
        //cout << spring[0].l0 << endl;
        //breathing_cube(spring, T);
        breathing_cube_gpu(&spring[0], T, info);
        //cout << "after: " << spring[0].l0 << endl;
        // calculate force on each spring
        calculate_force(mass, spring, force);
//        add_external_force(mass, spring, force);
        add_ground_force(mass, spring, force);

        // update position of cube
        update_position(mass, spring, force);

        // subtract height from fitness:
        for (int i=0; i<NUM_OF_MASSES; i++) {
            if (mass[i].p[2] > L0_SIDE)
                individual.fitness -= (mass[i].p[2])/NUM_OF_MASSES;
        }

        // calculate energy
        kinetic_energy.emplace_back(calculate_kinetic_energy(mass, spring));
        potential_energy.emplace_back(calculate_potential_energy(mass, spring));

//        cout << "T: " << T << "\n";
//        for (int i = 0; i < NUM_OF_MASSES; i++) {
//            print_mass(mass, 0);
//        }
//        cout << "e: " << kinetic_energy[iteration] + potential_energy[iteration] << "\n\n";

        // write to file for opengl
        if (opengl) {
            write_to_opengl_file(mass, opengl_file);
        }

        // update time
        T += DT;
    }

    // Calculate center of mass after NUM_OF_ITERATIONS iterations:
    vector<double> ending_com = calculate_center_of_mass(mass);

    // Calculate total distance travelled and its x componenet:
    double dist_travelled = dist(starting_com, ending_com);
    double dist_travelled_x = ending_com[0] - starting_com[0];
    double dist_travelled_y = ending_com[1] - starting_com[1];
    double dist_travelled_z = ending_com[2] - starting_com[2];

    // assign fitness equal to distance travelled in the positive x direction
    // TODO: maybe change later to be a function of dist_travelled as well?
    //individual.fitness = dist_travelled - dist_travelled_z;
    individual.fitness += dist_travelled_x * num_iterations - (abs(dist_travelled_z) * Z_PENALTY);
    individual.fitness -= abs(dist_travelled_y) * (num_iterations / 2);
//    cout << "dist travelled (x-dir): " << dist_travelled_x << endl;
    individual.distance = dist_travelled_x;

    // write energy to file
    if (opengl) {
        for (int i = 0; i < kinetic_energy.size(); i++) {
            energy_file << kinetic_energy[i];
            if (i + 1 < kinetic_energy.size()) {
                energy_file << ",";
            }
        }
        energy_file << "\n";
        for (int i = 0; i < potential_energy.size(); i++) {
            energy_file << potential_energy[i];
            if (i + 1 < potential_energy.size()) {
                energy_file << ",";
            }
        }
        energy_file << "\n";
        for (int i = 0; i < potential_energy.size(); i++) {
            energy_file << kinetic_energy[i] + potential_energy[i];
            if (i + 1 < potential_energy.size()) {
                energy_file << ",";
            }
        }
        
        cout << "Final distance (x direction): " << dist_travelled_x << endl;
    }
    opengl_file.close();
    energy_file.close();
}

float dist(vector<double> a, vector<double> b) {
    return sqrt(pow(b[0]-a[0],2) + pow(b[1]-a[1],2) + pow(b[2]-a[2],2));
}

void calculate_force(vector<Mass> &mass, vector<Spring> &spring, vector<vector<double>> &force) {
    // force due to spring
    for (int i = 0; i < NUM_OF_SPRINGS; i++) {
        // calculate force vector for spring
        double length = dist(mass[spring[i].m1].p, mass[spring[i].m2].p);

        double forceNormalized = spring[i].k * (length - spring[i].l0);
//        cout << "fN: " << forceNormalized << "\n";
        // not sure why, but even when there is no horizontal motion,
        // length becomes != spring[i].l0 at some seemingly random point
//        if (forceNormalized < 0.00000001) {
//            forceNormalized = 0;
//        }

        vector<double> forceVector = {forceNormalized * (mass[spring[i].m2].p[0] - mass[spring[i].m1].p[0]) / length,
                                      forceNormalized * (mass[spring[i].m2].p[1] - mass[spring[i].m1].p[1]) / length,
                                      forceNormalized * (mass[spring[i].m2].p[2] - mass[spring[i].m1].p[2]) / length};

        // now update force vector for masses that spring touches
        force[spring[i].m1][0] += forceVector[0];
        force[spring[i].m1][1] += forceVector[1];
        force[spring[i].m1][2] += forceVector[2];
        force[spring[i].m2][0] -= forceVector[0];
        force[spring[i].m2][1] -= forceVector[1];
        force[spring[i].m2][2] -= forceVector[2];
    }

    // force due to gravity
    for (int i = 0; i < NUM_OF_MASSES; i++) {
        force[i][2] -= mass[i].m * G; // note: acceleration due to gravity is defined as negative
    }

}

void add_external_force(vector<Mass> &mass, vector<Spring> &spring, vector<vector<double>> &force) {
    // this is where i would calculate the external forces ... if there were any
}

void add_ground_force(vector<Mass> &mass, vector<Spring> &spring, vector<vector<double>> &force) {
    for (int i = 0; i < NUM_OF_MASSES; i++) {
        // if mass is "under" ground
        if (mass[i].p[2] <= 0) {

            // calculate force due to friction:
            /**/

            double force_horizontal = sqrt(pow(force[i][0], 2) + pow(force[i][1], 2));
            // get angle of horizontal force:
            double cos_theta = force[i][0]/force_horizontal;  // F_x/F
            double sin_theta = force[i][1]/force_horizontal;  // F_y/F

            if (force[i][2] < 0) {
                if (force_horizontal < (-1) * force[i][2] * U_S) {
                    force[i][0] = 0;
                    force[i][1] = 0;
                    //cout << "STATIC\n";
                }
                else {
                    // Add force of kineteic friction in the opposite direction
                    // of the force. We know that force[i][2] is always negative
                    // here, so we add it to force_horizontal (because we want
                    // to actually subtract it):
                    double force_kinetic_friction = U_K * force[i][2];
                    //cout << "Fk = " << force_kinetic_friction << endl;
                    //cout << "Fh = " << force_horizontal << endl;
                    force_horizontal += force_kinetic_friction;
                    //cout << "Fh' = " << force_horizontal << endl;

                    // get back the x and y components of the horizontal force:
                    //cout << "Before: Fx = " << force[i][0] << ", Fy = " << force[i][1] << endl;
                    force[i][0] = force_horizontal*cos_theta;
                    force[i][1] = force_horizontal*sin_theta;
                    //cout << "After: Fx = " << force[i][0] << ", Fy = " << force[i][1] << endl;

                    
                    //force[i][0] += force[i][2] * U_K * (force[i][0] / force_horizontal);
                    //force[i][1] += force[i][2] * U_K * (force[i][1] / force_horizontal);
                    
                    
                }
            }/**/  // END FRICTION

            // Apply restorative force
            force[i][2] += K_GROUND * abs(mass[i].p[2]);
        }
    }
}

void update_position(vector<Mass> &mass, vector<Spring> &spring, vector<vector<double>> &force) {
    for (int i = 0; i < NUM_OF_MASSES; i++) {
        // acceleration, velocity, position calculation
        for (int j = 0; j < DIMENSIONS; j++) {
            mass[i].a[j] = force[i][j] / mass[i].m;
            mass[i].v[j] += mass[i].a[j] * DT;
            mass[i].v[j] = mass[i].v[j] * V_DAMP_CONST; // velocity dampening
            mass[i].p[j] += mass[i].v[j] * DT;
        }
    }
}

double calculate_potential_energy(vector<Mass> &mass, vector<Spring> &spring) {
    double ground_potential_energy = 0.0;
    double gravity_potential_energy = 0.0;
    double spring_potential_energy = 0.0;

    for (int i = 0; i < NUM_OF_MASSES; i++) {
        // potential energy due to gravity
        gravity_potential_energy += mass[i].m * abs(G) * mass[i].p[2];
        // energy due to ground
        if (mass[i].p[2] < 0) {
            ground_potential_energy += 0.5 * K_GROUND * pow(mass[i].p[2], 2); // maybe change to just a reverse kinetic energy
        }
    }

    // potential energy due to springs
    for (int i = 0; i < NUM_OF_SPRINGS; i++) {
        spring_potential_energy += 0.5 * spring[i].k * pow(spring[i].l0 - dist(mass[spring[i].m2].p, mass[spring[i].m1].p), 2);
    }

    return gravity_potential_energy + ground_potential_energy + spring_potential_energy;
}

double calculate_kinetic_energy(vector<Mass> &mass, vector<Spring> &spring) {
    double kinetic_energy = 0.0;

    for (int i = 0; i < NUM_OF_MASSES; i++) {
        // kinetic energy
        kinetic_energy += 0.5 * mass[i].m * pow(dist(mass[i].v, {0, 0, 0}), 2);
    }

    return kinetic_energy;
}

void print_mass (Mass &mass) {
    cout << "p: " << mass.p[0] << " " << mass.p[1] << " " << mass.p[2] << "\n";
    cout << "v: " << mass.v[0] << " " << mass.v[1] << " " << mass.v[2] << "\n";
    cout << "a: " << mass.a[0] << " " << mass.a[1] << " " << mass.a[2] << "\n";
}

void write_to_opengl_file(vector<Mass> &mass, ofstream &opengl_file) {
    vector<int> order = {4, 5, 7, 7, 5, 6, 7, 6, 3, 3, 6, 2, 3, 2, 0, 0, 2, 1, 0,
                         1, 4, 4, 1, 5, 5, 1, 6, 6, 1, 2, 0, 4, 3, 3, 4, 7};
    for (int i = 0; i < NUM_OPEN_GL_CUBE_PTS; i++) {
        opengl_file << mass[order[i]].p[0] << "," << mass[order[i]].p[1] << "," << mass[order[i]].p[2];
        if (i != NUM_OPEN_GL_CUBE_PTS - 1) {
            opengl_file << ",";
        }
    }
    opengl_file << "\n";
}

void breathing_cube(vector<Spring> &spring, float T) {
    for (int i = 0; i < NUM_OF_SPRINGS; i++) {
        spring[i].l0 = spring[i].a + spring[i].b * sin(OMEGA * T + spring[i].c);
    }
}

void breathing_cube_gpu(Spring* springs, float T, Gpu_info* info) {
    info->err = clEnqueueWriteBuffer(info->commands, 
                                     info->input_bc, 
                                     CL_TRUE, 
                                     0, 
                                     sizeof(Spring) * NUM_OF_SPRINGS, 
                                     springs, 
                                     0, 
                                     NULL, 
                                     NULL
                                    );

    //cout << "err after write buffer: " << info->err << endl;


    info->err = 0;
    info->err  = clSetKernelArg(info->kernel_bc, 
                                0, 
                                sizeof(cl_mem), 
                                &(info->input_bc)
                               );
    //cout << "err after set args: " << info->err << endl;
    if (info->err != 0) exit(1);

    float w = OMEGA; 
    info->err |= clSetKernelArg(info->kernel_bc, 
                                1, 
                                sizeof(cl_float), 
                                &w
                               );
    info->err |= clSetKernelArg(info->kernel_bc, 
                                2, 
                                sizeof(cl_float), 
                                &T
                               );
    int num_springs = NUM_OF_SPRINGS;
    info->err |= clSetKernelArg(info->kernel_bc, 
                                3, 
                                sizeof(cl_int), 
                                &num_springs
                               );
    //cout << "err after set args: " << info->err << endl;
    if (info->err != 0) exit(1);

    

    info->err = clGetKernelWorkGroupInfo(info->kernel_bc, 
                                         info->device_id, 
                                         CL_KERNEL_WORK_GROUP_SIZE, 
                                         sizeof(info->local_bc), 
                                         &(info->local_bc), 
                                         NULL
                                        );
    if (info->err != 0) {
        cout << "kernel workgroup info err\n";
        exit(1);
    }

    //cout << "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: " << CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS << endl;
    //cout << "CL_KERNEL_WORK_GROUP_SIZE: " << CL_KERNEL_WORK_GROUP_SIZE << endl;
    //cout << "sizeof(Spring) :" << sizeof(Spring) <<endl;
    //cout << "local size: " << info->local_bc << endl;
    //cout << "num springs (global): " << NUM_OF_SPRINGS << endl;

    info->global_bc = NUM_OF_SPRINGS;
    info->err = clEnqueueNDRangeKernel(info->commands, 
                                       info->kernel_bc, 
                                       1, 
                                       NULL, 
                                       &(info->global_bc), 
                                       NULL, 
                                       0, 
                                       NULL, 
                                       NULL
                                      );
    if (info->err != 0) {
        cout << "enqueue kernel err: " << info->err << endl;
        exit(1);
    }


    clFinish(info->commands);
    info->err = clEnqueueReadBuffer(info->commands, 
                                    info->input_bc, 
                                    CL_TRUE, 
                                    0, 
                                    sizeof(float) * NUM_OF_SPRINGS, 
                                    springs, 
                                    0, 
                                    NULL, 
                                    NULL 
                                   );
    if (info->err != 0) {
        cout << "read buffer error\n";
        exit(1);
    }

}

/* 
 * returns by value a vector of size 3 containing the center of mass of the 
 * masses passed to this function
 */
vector<double> calculate_center_of_mass(vector<Mass> &masses) {
    double total_mass = 0.0;
    double pos_x = 0.0;
    double pos_y = 0.0;
    double pos_z = 0.0;

    for (int i=0; i<masses.size(); i++) {
        total_mass += masses[i].m;
        pos_x += masses[i].m * masses[i].p[0];
        pos_y += masses[i].m * masses[i].p[1];
        pos_z += masses[i].m * masses[i].p[2];
    }

    pos_x /= total_mass;
    pos_y /= total_mass;
    pos_z /= total_mass;

    vector<double> com;
    com.emplace_back(pos_x);
    com.emplace_back(pos_y);
    com.emplace_back(pos_z);

    return com;
}

