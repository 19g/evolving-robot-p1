#include <iostream>
#include <fstream>
#include <ctime>
#include "main.hpp"

using namespace std;

int main() {

    // initialize files
    ofstream energy_file;
    energy_file.open(ENERGY_TXT);
    ofstream opengl_file;
    opengl_file.open(OPENGL_TXT);

    // declare variables
    double T = 0.0;
    vector<Mass> mass;
    vector<Spring> spring;
    vector<double> kinetic_energy;
    vector<double> potential_energy;

    // create initial cube
    initialize_cube(mass, spring);
    cout << mass[0].p[2] << "\n";


//    for (int i = 0; i < NUM_OF_MASSES; i++) {
//        cout << mass[i].p[0] << "\t";
//        cout << mass[i].p[1] << "\t";
//        cout << mass[i].p[2] << "\t\n";
//    }

    // begin timer
    clock_t begin = clock();
    // simulation loop
    for (int iteration = 0; iteration < NUM_OF_ITERATIONS; iteration++) {
        // initialize force vector
        vector<vector<double>> force(NUM_OF_MASSES, vector<double>(DIMENSIONS));

        // breathing cube
//        breathing_cube(spring, T);
        // calculate force on each spring
        calculate_force(mass, spring, force);
//        add_external_force(mass, spring, force);
        add_ground_force(mass, spring, force);

        // update position of cube
        update_position(mass, spring, force);
        // calculate energy
        kinetic_energy.emplace_back(calculate_kinetic_energy(mass, spring));
        potential_energy.emplace_back(calculate_potential_energy(mass, spring));

        cout << "T: " << T << "\n";
//        for (int i = 0; i < NUM_OF_MASSES; i++) {
//            print_mass(mass, 0);
//        }
        cout << "e: " << kinetic_energy[iteration] + potential_energy[iteration] << "\n\n";

        // write to file for opengl
        write_to_opengl_file(mass, opengl_file);

        // update time
        T += DT;
    }

    // end timer
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    double iters_per_sec = NUM_OF_ITERATIONS / elapsed_secs;
    cout << "iter/sec: " << iters_per_sec << "\n";


    // write energy to file
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

    opengl_file.close();
    energy_file.close();
}

void initialize_cube(vector<Mass> &mass, vector<Spring> &spring) {
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

    spring.emplace_back(Spring {K_SPRING, dist(mass[0].p, mass[1].p), 0, 1});
    spring.emplace_back(Spring {K_SPRING, dist(mass[0].p, mass[2].p), 0, 2});
    spring.emplace_back(Spring {K_SPRING, dist(mass[0].p, mass[3].p), 0, 3});
    spring.emplace_back(Spring {K_SPRING, dist(mass[0].p, mass[4].p), 0, 4});
    spring.emplace_back(Spring {K_SPRING, dist(mass[0].p, mass[5].p), 0, 5});
    spring.emplace_back(Spring {K_SPRING, dist(mass[0].p, mass[6].p), 0, 6});
    spring.emplace_back(Spring {K_SPRING, dist(mass[0].p, mass[7].p), 0, 7});
    spring.emplace_back(Spring {K_SPRING, dist(mass[1].p, mass[2].p), 1, 2});
    spring.emplace_back(Spring {K_SPRING, dist(mass[1].p, mass[3].p), 1, 3});
    spring.emplace_back(Spring {K_SPRING, dist(mass[1].p, mass[4].p), 1, 4});
    spring.emplace_back(Spring {K_SPRING, dist(mass[1].p, mass[5].p), 1, 5});
    spring.emplace_back(Spring {K_SPRING, dist(mass[1].p, mass[6].p), 1, 6});
    spring.emplace_back(Spring {K_SPRING, dist(mass[1].p, mass[7].p), 1, 7});
    spring.emplace_back(Spring {K_SPRING, dist(mass[2].p, mass[3].p), 2, 3});
    spring.emplace_back(Spring {K_SPRING, dist(mass[2].p, mass[4].p), 2, 4});
    spring.emplace_back(Spring {K_SPRING, dist(mass[2].p, mass[5].p), 2, 5});
    spring.emplace_back(Spring {K_SPRING, dist(mass[2].p, mass[6].p), 2, 6});
    spring.emplace_back(Spring {K_SPRING, dist(mass[2].p, mass[7].p), 2, 7});
    spring.emplace_back(Spring {K_SPRING, dist(mass[3].p, mass[4].p), 3, 4});
    spring.emplace_back(Spring {K_SPRING, dist(mass[3].p, mass[5].p), 3, 5});
    spring.emplace_back(Spring {K_SPRING, dist(mass[3].p, mass[6].p), 3, 6});
    spring.emplace_back(Spring {K_SPRING, dist(mass[3].p, mass[7].p), 3, 7});
    spring.emplace_back(Spring {K_SPRING, dist(mass[4].p, mass[5].p), 4, 5});
    spring.emplace_back(Spring {K_SPRING, dist(mass[4].p, mass[6].p), 4, 6});
    spring.emplace_back(Spring {K_SPRING, dist(mass[4].p, mass[7].p), 4, 7});
    spring.emplace_back(Spring {K_SPRING, dist(mass[5].p, mass[6].p), 5, 6});
    spring.emplace_back(Spring {K_SPRING, dist(mass[5].p, mass[7].p), 5, 7});
    spring.emplace_back(Spring {K_SPRING, dist(mass[6].p, mass[7].p), 6, 7});

}

double dist(vector<double> a, vector<double> b) {
    return sqrt(pow(b[0]-a[0],2) + pow(b[1]-a[1],2) + pow(b[2]-a[2],2));
}

void calculate_force(vector<Mass> &mass, vector<Spring> &spring, vector<vector<double>> &force) {
    // force due to spring
    for (int i = 0; i < NUM_OF_SPRINGS; i++) {
        // calculate force vector for spring
        double length = dist(mass[spring[i].m1].p, mass[spring[i].m2].p);

        double forceNormalized = K_SPRING * (length - spring[i].l0);
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
        // if "under" ground, then apply restorative force
        if (mass[i].p[2] < 0) {
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
        spring_potential_energy += 0.5 * K_SPRING * pow(spring[i].l0 - dist(mass[spring[i].m2].p, mass[spring[i].m1].p), 2);
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

void print_mass (vector<Mass> &mass, int i) {
    cout << "p: " << mass[i].p[0] << " " << mass[i].p[1] << " " << mass[i].p[2] << "\n";
    cout << "v: " << mass[i].v[0] << " " << mass[i].v[1] << " " << mass[i].v[2] << "\n";
    cout << "a: " << mass[i].a[0] << " " << mass[i].a[1] << " " << mass[i].a[2] << "\n";
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
        spring[i].l0 = 0.9 * L0_SIDE * cos(5 * T) + L0_SIDE;
    }
}
