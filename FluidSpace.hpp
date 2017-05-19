#ifndef FLUIDSPACE_HPP
#define FLUIDSPACE_HPP
#include <vector>
#include <list>
#include <iostream>
#include <SFML/System.hpp>

/*
sf::Clock collideAvgTime;
float collideTotalTime = 0;
int collideCallCount = 0;

sf::Clock streamAvgTime;
float streamTotalTime = 0;
int streamCallCount = 0;
*/

struct Velocity {
    float x, y;
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
    std::vector< float > densityVect;
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
public:
    FluidSpace();
    FluidSpace(int w, int h);
    void stream();
    void collide();
    int getWidth();
    int getHeight();
    std::vector<float>& getDensity();
    std::vector< Velocity > getVelocity();
    void setDensity(int x, int y, float val);
    void incrDensity(int x, int y);
    void incrDensityPx(int x, int y);
    void update(); // tbi: delta time argument
    void printLattice(int x, int y);
};

class Lattices {
    std::vector<float> eqNumD; // equilibrium number densities
    std::vector<float> v; // velocity vector
public:
    int width;
    int height;
    std::vector<float>
        nE, nN, nW, nS, nNE, nNW, nSW, nSE, n0;
    std::vector<float>
        eqE, eqN, eqW, eqS, eqNE, eqNW, eqSW, eqSE, eq0;
    std::vector<float> densities;
    Lattices() {};
    Lattices(int w, int h);
    void resetBorder();
    void resetToZeroes();
    float getDensity(int x, int y);
    float getDensity(int index);
    //void updateDensities();
    std::vector<float>& updateDensities();
    void updateEqDensities();
    std::vector<float>& getVelocity(int x, int y);
    std::vector<float>& getVelocity(int index);
    std::vector<float> getEquilibNumDensities(int index);
    void print(int x, int y);
};

#endif // FLUIDSPACE_HPP
