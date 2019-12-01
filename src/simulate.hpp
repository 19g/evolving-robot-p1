#ifndef EVOLVING_ROBOT_SIMULATE_HPP
#define EVOLVING_ROBOT_SIMULATE_HPP

#include <vector>
#include <cmath>

// cube variables
#define WEIGHT_PER_MASS 0.1
#define L0_SIDE 0.1
//#define K_SPRING 10000.0
//#define INITIAL_HEIGHT 2.0
#define NUM_OF_MASSES 8
#define NUM_OF_SPRINGS 28
#define DIMENSIONS 3

// world variables
#define G 9.807
#define OMEGA (4*PI)
#define K_GROUND 100000.0
#define DT 0.0001
#define V_DAMP_CONST 0.999 //0.999999
#define NUM_OF_ITERATIONS 15000
// 15000 = 3 cycles (for omega = 4pi and DT = 0.0001)
#define U_S 1.0
#define U_K 0.8
#define Z_PENALTY 5

// opengl
#define NUM_OPEN_GL_CUBE_PTS 36

// file names
#define ENERGY_TXT "energy.txt"

struct Mass {
    double m; // mass in kg
    std::vector<double> p; // position std::vector in m
    std::vector<double> v; // velocity std::vector in m/s
    std::vector<double> a; // acceleration std::vector in m/s^2
};

struct Spring {
    double k; // spring constant in N/m
    double l0; // original rest length in m
    int m1; // index of first mass object in spring
    int m2; // index of second mass object in spring
    double a; // a + b*sin(wt+c) + d*sin(2wt+e)
    double b; // a + b*sin(wt+c) + d*sin(2wt+e)
    double c; // a + b*sin(wt+c) + d*sin(2wt+e)
};

struct Cube {
    std::vector<Mass> mass;
    std::vector<Spring> spring;
    double fitness;
    double distance;
};

// perform the physics simulation
void simulation_loop(Cube &, int, bool);
// initialize cube with masses and springs
void initialize_cube(std::vector<Mass> &, std::vector<Spring> &);
// calculate distance between two 3D points
double dist(std::vector<double>, std::vector<double>);
// calculate force from springs, gravity
void calculate_force(std::vector<Mass> &, std::vector<Spring> &, std::vector<std::vector<double>> &);
// add external forces
void add_external_force(std::vector<Mass> &, std::vector<Spring> &, std::vector<std::vector<double>> &);
// add force due to ground
void add_ground_force(std::vector<Mass> &, std::vector<Spring> &, std::vector<std::vector<double>> &);
// update position of masses due to force
void update_position(std::vector<Mass> &, std::vector<Spring> &, std::vector<std::vector<double>> &);
// calculate total energy of cube
double calculate_potential_energy(std::vector<Mass> &, std::vector<Spring> &);
// calculate kinetic energy of cube
double calculate_kinetic_energy(std::vector<Mass> &, std::vector<Spring> &);
// print mass object
void print_mass(Mass &);
// write to opengl file
void write_to_opengl_file(std::vector<Mass> &, std::ofstream &);
// breathing cube function
void breathing_cube(std::vector<Spring> &, double);

// calculate center of mass of a std::vector of masses:
std::vector<double> calculate_center_of_mass(std::vector<Mass> &);


#endif //EVOLVING_ROBOT_SIMULATE_HPP

