#include <iostream>
#include <fstream>
#include <thread>

#include "simulate.hpp"
#include "geneticalgorithm.hpp"

using namespace std;

void simulation_loop(Cube &individual, bool opengl) {
    // initialize files
    ofstream energy_file;
    energy_file.open(ENERGY_TXT);
    ofstream opengl_file;
    opengl_file.open(OPENGL_TXT);

    // declare variables
    double T = 0.0;
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

    // write to file for opengl
    if (opengl) {
        write_to_opengl_file(mass, opengl_file);
    }

    int num_iterations = NUM_OF_ITERATIONS;
    if (opengl) num_iterations *= 3;

    // simulation loop
    thread threads[NUM_OF_MASSES];
    for (int iteration = 0; iteration < num_iterations; iteration++) {
        // initialize force vector
        vector<vector<double>> force(NUM_OF_MASSES, vector<double>(DIMENSIONS));

        // breathing cube
        breathing_cube(spring, T);
        // calculate force on each spring
        calculate_force(mass, spring, force);  // Potential threading race condition in this function
//        add_external_force(mass, spring, force);

        // these two functions can be easily threaded:
        for (int i=0; i<NUM_OF_MASSES; i++) {
            // add friction and normal forces:
            //add_ground_force(mass, spring, force);
            threads[i] = thread(add_ground_force, 
                                ref(mass[i]), ref(spring), ref(force[i]));
        }
        
        for (int i=0; i<NUM_OF_MASSES; i++) {
            // wait for add_ground_force[i] to finish:
            threads[i].join();
        }
        for (int i=0; i<NUM_OF_MASSES; i++) {

            // update position of cube:
            //update_position(mass, spring, force);
            threads[i] = thread(update_position, 
                                ref(mass[i]), ref(spring), ref(force[i]));
        }
        // wait for all threads:
        
        for (int i=0; i<NUM_OF_MASSES; i++) {
            threads[i].join();
        }

        // calculate energy
        //kinetic_energy.emplace_back(calculate_kinetic_energy(mass, spring));
        //potential_energy.emplace_back(calculate_potential_energy(mass, spring));

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
    double dist_travelled_z = ending_com[2] - starting_com[2];

    // assign fitness equal to distance travelled in the positive x direction
    // TODO: maybe change later to be a function of dist_travelled as well?
    //individual.fitness = dist_travelled - dist_travelled_z;
    individual.fitness = dist_travelled_x; // - abs(dist_travelled_z);


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
    }
    opengl_file.close();
    energy_file.close();
}

double dist(vector<double> a, vector<double> b) {
    return sqrt(pow(b[0]-a[0],2) + pow(b[1]-a[1],2) + pow(b[2]-a[2],2));
}

void calculate_force(vector<Mass> &mass, 
                     vector<Spring> &spring, 
                     vector<vector<double>> &force) {

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
        // TODO: (warning) this would be a race condition for threading
        // may introduce invalid values
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

void add_external_force(vector<Mass> &mass, 
                        vector<Spring> &spring, 
                        vector<vector<double>> &force) {
    // this is where i would calculate the external forces ... if there were any
}

void add_ground_force(Mass &mass, 
                      vector<Spring> &spring, 
                      vector<double> &force) {

    // un-cast constness:
    /*
    Mass &mass = const_cast<Mass &>(_mass);
    vector<Spring> &spring = const_cast<vector<Spring> &>(_spring);
    vector<double> &force = const_cast<vector<double> &>(_force);
    */

    //for (int i = 0; i < NUM_OF_MASSES; i++) {
        // if mass is "under" ground
        if (mass.p[2] <= 0) {

            // calculate force due to friction:
            /**/

            double force_horizontal = sqrt(pow(force[0], 2) + pow(force[1], 2));
            // get angle of horizontal force:
            double cos_theta = force[0]/force_horizontal;  // F_x/F
            double sin_theta = force[1]/force_horizontal;  // F_y/F

            if (force[2] < 0) {
                if (force_horizontal < (-1) * force[2] * U_S) {
                    force[0] = 0;
                    force[1] = 0;
                    //cout << "STATIC\n";
                }
                else {
                    // Add force of kineteic friction in the opposite direction
                    // of the force. We know that force[i][2] is always negative
                    // here, so we add it to force_horizontal (because we want
                    // to actually subtract it):
                    double force_kinetic_friction = U_K * force[2];
                    //cout << "Fk = " << force_kinetic_friction << endl;
                    //cout << "Fh = " << force_horizontal << endl;
                    force_horizontal += force_kinetic_friction;
                    //cout << "Fh' = " << force_horizontal << endl;

                    // get back the x and y components of the horizontal force:
                    //cout << "Before: Fx = " << force[i][0] << ", Fy = " << force[i][1] << endl;
                    force[0] = force_horizontal*cos_theta;
                    force[1] = force_horizontal*sin_theta;
                    //cout << "After: Fx = " << force[i][0] << ", Fy = " << force[i][1] << endl;

                    
                    //force[i][0] += force[i][2] * U_K * (force[i][0] / force_horizontal);
                    //force[i][1] += force[i][2] * U_K * (force[i][1] / force_horizontal);
                    
                    
                }
            }/**/  // END FRICTION

            // Apply restorative force
            force[2] += K_GROUND * abs(mass.p[2]);
        }
    //}
}

void update_position(Mass &mass, 
                     vector<Spring> &spring, 
                     vector<double> &force) {

    // un-cast constness:
    /*
    Mass &mass = const_cast<Mass &>(_mass);
    vector<Spring> &spring = const_cast<vector<Spring> &>(_spring);
    vector<double> &force = const_cast<vector<double> &>(_force);
    */

    //for (int i = 0; i < NUM_OF_MASSES; i++) {
        // acceleration, velocity, position calculation
        for (int j = 0; j < DIMENSIONS; j++) {
            mass.a[j] = force[j] / mass.m;
            mass.v[j] += mass.a[j] * DT;
            mass.v[j] = mass.v[j] * V_DAMP_CONST; // velocity dampening
            mass.p[j] += mass.v[j] * DT;
        }
    //}
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

void breathing_cube(vector<Spring> &spring, double T) {
    for (int i = 0; i < NUM_OF_SPRINGS; i++) {
        spring[i].l0 = spring[i].a + spring[i].b * sin(OMEGA * T + spring[i].c);
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

