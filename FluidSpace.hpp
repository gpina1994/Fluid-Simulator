#ifndef FLUIDSPACE_HPP
#define FLUIDSPACE_HPP
#include <vector>
#include <list>
#include <iostream>
#include <SFML/System.hpp>

/*
sf::Clock collideAvgTime;
double collideTotalTime = 0;
int collideCallCount = 0;

sf::Clock streamAvgTime;
double streamTotalTime = 0;
int streamCallCount = 0;
*/

struct Velocity {
    double x, y;
public:
    Velocity(): x(0.0), y(0.0) {};
};

class Lattices;

class FluidSpace {
    /**
      this class contains everything related to the physics of the
      simulation and nothing related to the external library used
      to create an interface with the user.
    **/
    std::list<Lattices> latticesList; //holds the 2 grids
    int width;
    int height;
    std::vector< double > densityVect;
    std::vector< Velocity > velocityVect; // may be better to just make a velocity struct
    void streamE();
    void streamN();
    void streamW();
    void streamS();
    void streamNE();
    void streamNW();
    void streamSE();
    void streamSW();
    void stream0();
    bool isValidCoords(int x, int y);
    int xyToIndex(int x, int y);
public:
    FluidSpace();
    FluidSpace(int w, int h);
    void stream();
    void collide();
    int getWidth();
    int getHeight();
    std::vector<double>& getDensity();
    std::vector< Velocity > getVelocity();
    void setDensity(int x, int y, double val);
    void incrDensity(int x, int y);
    void incrDensityPx(int x, int y);
    void applyWind(int x, int y, double directionRad);
    void applyWindPx(int x, int y, double directionRad);
    void update(); // tbi: delta time argument
    void printLattice(int x, int y);
    void printLattice(int index);
    void printTotalMass();
};

class Lattices {
    std::vector<double> eqNumD; // equilibrium number densities
    std::vector<double> v; // velocity vector
    double latticeSpeed;
public:
    int width;
    int height;
    std::vector<double>
        nE, nN, nW, nS, nNE, nNW, nSW, nSE, n0;
    std::vector<double>
        eqE, eqN, eqW, eqS, eqNE, eqNW, eqSW, eqSE, eq0;
    std::vector<double> densities;
    Lattices() {};
    Lattices(int w, int h);
    void resetBorder();
    void resetToZeroes();
    double getDensity(int x, int y);
    double getDensity(int index);
    //void updateDensities();
    std::vector<double>& updateDensities();
    void updateEqDensities();
    std::vector<double>& getVelocity(int x, int y);
    std::vector<double>& getVelocity(int index);
    void print(int x, int y);
    void print(int index);
    void printEq(int index);
};

#endif // FLUIDSPACE_HPP
