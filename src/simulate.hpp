#ifndef EVOLVING_ROBOT_SIMULATE_HPP
#define EVOLVING_ROBOT_SIMULATE_HPP

using namespace std;
#include <vector>
#include <cmath>

// cube variables
#define WEIGHT_PER_MASS 0.1
#define L0_SIDE 0.1
//#define K_SPRING 1000.0
//#define INITIAL_HEIGHT 2.0
#define NUM_OF_MASSES 8
#define NUM_OF_SPRINGS 28
#define DIMENSIONS 3

// world variables
#define G 9.807
#define K_GROUND 100000.0
#define DT 0.0001
#define V_DAMP_CONST 0.999999 //0.999999
#define NUM_OF_ITERATIONS 20000

// opengl
#define NUM_OPEN_GL_CUBE_PTS 36

// file names
#define ENERGY_TXT "energy.txt"
#define OPENGL_TXT "cubes.txt"

struct Mass {
    double m; // mass in kg
    vector<double> p; // position vector in m
    vector<double> v; // velocity vector in m/s
    vector<double> a; // acceleration vector in m/s^2
};

struct Spring {
    double k; // spring constant in N/m
    double l0; // original rest length in m
    int m1; // index of first mass object in spring
    int m2; // index of second mass object in spring
    double a; // a + b*sin(wt+c) + d*sin(2wt+e)
    double b; // a + b*sin(wt+c) + d*sin(2wt+e)
    double c; // a + b*sin(wt+c) + d*sin(2wt+e)
    double d; // a + b*sin(wt+c) + d*sin(2wt+e)
    double e; // a + b*sin(wt+c) + d*sin(2wt+e)
};

// perform the physics simulation
double simulation_loop(vector<Mass> &, vector<Spring> &);
// initialize cube with masses and springs
void initialize_cube(vector<Mass> &, vector<Spring> &);
// calculate distance between two points
double dist(vector<double>, vector<double>);
// calculate force from springs, gravity
void calculate_force(vector<Mass> &, vector<Spring> &, vector<vector<double>> &);
// add external forces
void add_external_force(vector<Mass> &, vector<Spring> &, vector<vector<double>> &);
// add force due to ground
void add_ground_force(vector<Mass> &, vector<Spring> &, vector<vector<double>> &);
// update position of masses due to force
void update_position(vector<Mass> &, vector<Spring> &, vector<vector<double>> &);
// calculate total energy of cube
double calculate_potential_energy(vector<Mass> &, vector<Spring> &);
// calculate kinetic energy of cube
double calculate_kinetic_energy(vector<Mass> &, vector<Spring> &);
// print mass object
void print_mass(Mass &);
// write to opengl file
void write_to_opengl_file(vector<Mass> &, ofstream &);
// breathing cube function
void breathing_cube(vector<Spring> &, double);




#endif //EVOLVING_ROBOT_SIMULATE_HPP
