#include "FluidSpace.hpp"
#include <stdlib.h> // rand, srand
#include <time.h> // seeding rand
#include <math.h> // isnan()

FluidSpace::FluidSpace()
{
    // nothing to do here
}

FluidSpace::FluidSpace(int w, int h)
{
    width = w;
    height = h;
    int vectSize = (w+2)*(h+2);
    densityVect.resize(vectSize);
    velocityVect.resize(vectSize);

    // initialize lattices
    Lattices lats(w, h);
    //latticesVect.resize(2, lats);
    latticesList.push_front(lats);
    latticesList.push_front(lats);
    std::cout << "nE lattice @ (100,100): " << latticesList.front().nE[(100+1)*(width+2)+(100+1)] << std::endl;
    // dummy functionality: initialize density array to 1's.
    for (int i=0; i<width+2; ++i) {
        for (int j=0; j<height+2; ++j) {
            densityVect[i*(height+2) + j] = 1.0;
        }
    }
}

std::vector<float>& FluidSpace::getDensity()
{
    return latticesList.front().updateDensities();
    //return latticesList.front().densities;
}

int FluidSpace::getWidth()
{
    return width;
}

int FluidSpace::getHeight()
{
    return height;
}

void FluidSpace::incrDensity(int x, int y)
{
    for (int i=-3; i<4; ++i) {
        for (int j=-3; j<4; ++j) {
            incrDensityPx(x+i,y+j);
        }
    }
}

void FluidSpace::incrDensityPx(int x, int y)
{
    /// check if position is valid (should just be a function)
    /// also, this class isn't the appropriate place to check that
    if (x < 0 || width < x || y < 0 || height < y) {
        return;
    }
    int index = (y+1)*(width+2)+(x+1);
    /// get original density
    float orig = latticesList.front().getDensity(index);
    //float incrAmt = (10.0 - orig)*0.60;
    float incrRatio = (orig+60)/orig; /// increment density by 60
    latticesList.front().n0[index] *= incrRatio;
    latticesList.front().nE[index] *= incrRatio;
    latticesList.front().nN[index] *= incrRatio;
    latticesList.front().nW[index] *= incrRatio;
    latticesList.front().nS[index] *= incrRatio;
    latticesList.front().nNE[index] *= incrRatio;
    latticesList.front().nNW[index] *= incrRatio;
    latticesList.front().nSW[index] *= incrRatio;
    latticesList.front().nSE[index] *= incrRatio;
}

void FluidSpace::collide()
{
    /// called through update()
    std::vector<float> equilibDensities;
    // update equilibrium densities
    latticesList.front().updateEqDensities();
    int ind = width+1;
    for (int y=0; y<height; ++y) {
        ind += 2;
        for (int x=0; x<width; ++x) {
            ++ind;
            //equilibDensities = latticesList.front().getEquilibNumDensities(ind);
            latticesList.front().n0[ind] = latticesList.front().n0[ind] + 1.8*(latticesList.front().eq0[ind] - latticesList.front().n0[ind]);
            latticesList.front().nE[ind] = latticesList.front().nE[ind] + 1.8*(latticesList.front().eqE[ind] - latticesList.front().nE[ind]);
            latticesList.front().nN[ind] = latticesList.front().nN[ind] + 1.8*(latticesList.front().eqN[ind] - latticesList.front().nN[ind]);
            latticesList.front().nW[ind] = latticesList.front().nW[ind] + 1.8*(latticesList.front().eqW[ind] - latticesList.front().nW[ind]);
            latticesList.front().nS[ind] = latticesList.front().nS[ind] + 1.8*(latticesList.front().eqS[ind] - latticesList.front().nS[ind]);
            latticesList.front().nNE[ind] = latticesList.front().nNE[ind] + 1.8*(latticesList.front().eqNE[ind] - latticesList.front().nNE[ind]);
            latticesList.front().nNW[ind] = latticesList.front().nNW[ind] + 1.8*(latticesList.front().eqNW[ind] - latticesList.front().nNW[ind]);
            latticesList.front().nSW[ind] = latticesList.front().nSW[ind] + 1.8*(latticesList.front().eqSW[ind] - latticesList.front().nSW[ind]);
            latticesList.front().nSE[ind] = latticesList.front().nSE[ind] + 1.8*(latticesList.front().eqSE[ind] - latticesList.front().nSE[ind]);
        }
    }
}

void FluidSpace::stream()
{ /// should be called through update()
    latticesList.back().resetToZeroes();
    streamE();
    streamN();
    streamW();
    streamS();
    streamNE();
    streamNW();
    streamSW();
    streamSE();
    stream0();
    //latticesList.back() = latticesList.front();
    latticesList.back().resetBorder();
    latticesList.reverse();
    return;
}

void FluidSpace::streamE()
{
    int offset = 1;
    for (int i=1; i<height+1; ++i) {
        for (int j=0; j<width+1; ++j) {
            latticesList.back().nE[j+i*(width+2) + offset] = latticesList.front().nE[j+i*(width+2)];
        }
    }
    return;
}

void FluidSpace::streamN()
{
    int offset = 0-width-2;
    for (int i=1; i<height+2; ++i) {
        for (int j=1; j<width+1; ++j) {
            latticesList.back().nN[j+i*(width+2) + offset] = latticesList.front().nN[j+i*(width+2)];
        }
    }
    return;
}

void FluidSpace::streamW()
{
    int offset = -1;
    for (int i=1; i<height+1; ++i) {
        for (int j=1; j<width+2; ++j) {
            latticesList.back().nW[j+i*(width+2) + offset] = latticesList.front().nW[j+i*(width+2)];
        }
    }
    return;
}

void FluidSpace::streamS()
{
    int offset = width+2;
    for (int i=0; i<height+1; ++i) {
        for (int j=1; j<width+1; ++j) {
            latticesList.back().nS[j+i*(width+2) + offset] = latticesList.front().nS[j+i*(width+2)];
        }
    }
    return;
}

void FluidSpace::streamNE()
{
    int offset = -width-1;
    for (int i=1; i<height+2; ++i) {
        for (int j=0; j<width+1; ++j) {
            latticesList.back().nNE[j+i*(width+2) + offset] = latticesList.front().nNE[j+i*(width+2)];
        }
    }
    return;
}

void FluidSpace::streamNW()
{
    int offset = -width-3;
    for (int i=1; i<height+2; ++i) {
        for (int j=1; j<width+2; ++j) {
            latticesList.back().nNW[j+i*(width+2) + offset] = latticesList.front().nNW[j+i*(width+2)];
        }
    }
    return;
}

void FluidSpace::streamSW()
{
    int offset = width+1;
    for (int i=0; i<height+1; ++i) {
        for (int j=1; j<width+2; ++j) {
            latticesList.back().nSW[j+i*(width+2) + offset] = latticesList.front().nSW[j+i*(width+2)];
        }
    }
    return;
}

void FluidSpace::streamSE()
{
    int offset = width+3;
    for (int i=0; i<height+1; ++i) {
        for (int j=0; j<width+1; ++j) {
            latticesList.back().nSE[j+i*(width+2) + offset] = latticesList.front().nSE[j+i*(width+2)];
        }
    }
    return;
}

void FluidSpace::stream0()
{
    for (int i=1; i<height+1; ++i) {
        for (int j=1; j<width+1; ++j) {
            latticesList.back().n0[j+i*(width+2)] = latticesList.front().n0[j+i*(width+2)];
        }
    }
}

void FluidSpace::update()
{
    /// step 1: stream
    stream();

    ///step 2: collide
    collide();

    return;
}

void FluidSpace::printLattice(int x, int y)
{
    latticesList.front().print(x,y);
}

Lattices::Lattices(int w, int h)
{
    width = w;
    height = h;
    int vectSize = (w+2)*(h+2); // Grid is buffered on each edge by 1 lattice
    n0.resize(vectSize, (float)4/9);
    nE.resize(vectSize, (float)1/9);
    nN.resize(vectSize, (float)1/9);
    nW.resize(vectSize, (float)1/9);
    nS.resize(vectSize, (float)1/9);
    nNE.resize(vectSize, (float)1/36);
    nNW.resize(vectSize, (float)1/36);
    nSW.resize(vectSize, (float)1/36);
    nSE.resize(vectSize, (float)1/36);

    eq0.resize(vectSize);
    eqE.resize(vectSize);
    eqN.resize(vectSize);
    eqW.resize(vectSize);
    eqS.resize(vectSize);
    eqNE.resize(vectSize);
    eqNW.resize(vectSize);
    eqSW.resize(vectSize);
    eqSE.resize(vectSize);

    eqNumD.resize(9);
    v.resize(2);
    densities.resize(w*h);
}

void Lattices::resetBorder()
{
    //set top edge
    for (int i=0; i<width+2; ++i) {
        n0[i] = (float)4/9;
        nE[i] = (float)1/9;
        nN[i] = (float)1/9;
        nW[i] = (float)1/9;
        nS[i] = (float)1/9;
        nNE[i] = (float)1/36;
        nNW[i] = (float)1/36;
        nSW[i] = (float)1/36;
        nSE[i] = (float)1/36;
    }
    //set left edge
    for (int i=0; i<(width+2)*(height+2); i += width+2) {
        n0[i] = (float)4/9;
        nE[i] = (float)1/9;
        nN[i] = (float)1/9;
        nW[i] = (float)1/9;
        nS[i] = (float)1/9;
        nNE[i] = (float)1/36;
        nNW[i] = (float)1/36;
        nSW[i] = (float)1/36;
        nSE[i] = (float)1/36;
    }
    //set right edge
    for (int i=width+1; i<(width+2)*(height+2); i += width+2) {
        n0[i] = (float)4/9;
        nE[i] = (float)1/9;
        nN[i] = (float)1/9;
        nW[i] = (float)1/9;
        nS[i] = (float)1/9;
        nNE[i] = (float)1/36;
        nNW[i] = (float)1/36;
        nSW[i] = (float)1/36;
        nSE[i] = (float)1/36;
    }
    //set bottom edge
    for (int i=(width+2)*(height+1); i<(width+2)*(height+2); ++i) {
        n0[i] = (float)4/9;
        nE[i] = (float)1/9;
        nN[i] = (float)1/9;
        nW[i] = (float)1/9;
        nS[i] = (float)1/9;
        nNE[i] = (float)1/36;
        nNW[i] = (float)1/36;
        nSW[i] = (float)1/36;
        nSE[i] = (float)1/36;
    }
}

void Lattices::resetToZeroes()
{
    for (int i=0; i<(width+2)*(height+2); ++i) {
        n0[i] = 0.0;
        nE[i] = 0.0;
        nN[i] = 0.0;
        nW[i] = 0.0;
        nS[i] = 0.0;
        nNE[i] = 0.0;
        nNW[i] = 0.0;
        nSW[i] = 0.0;
        nSE[i] = 0.0;
    }
}

inline float Lattices::getDensity(int x, int y)
{
    int index = (y+1)*(width+2)+(x+1);
    float density = 0.0;
    density += n0[index];
    density += nE[index];
    density += nN[index];
    density += nW[index];
    density += nS[index];
    density += nNE[index];
    density += nNW[index];
    density += nSW[index];
    density += nSE[index];
    return density;
}

inline float Lattices::getDensity(int index)
{
    return n0[index] + nE[index] + nN[index] + nW[index] + nS[index] +\
           nNE[index] + nNW[index] + nSW[index] + nSE[index];
}

/*
void Lattices::updateDensities()
{
    int ind;
    for (int i=0; i<height; ++i) {
        for (int j=0; j<width; ++j) {
            ind = (i+1)*(width+2)+(j+1);
            densities[i*width+j] = n0[ind] + nE[ind] + nN[ind] + nW[ind] + nS[ind] \
                                    + nNE[ind] + nNW[ind] + nSW[ind] + nSE[ind];
        }
    }
}
*/

std::vector<float>& Lattices::updateDensities()
{
    int ind;
    for (int i=0; i<height; ++i) {
        for (int j=0; j<width; ++j) {
            ind = (i+1)*(width+2)+(j+1);
            densities[i*width+j] = n0[ind] + nE[ind] + nN[ind] + nW[ind] + nS[ind] \
                                    + nNE[ind] + nNW[ind] + nSW[ind] + nSE[ind];
        }
    }
    return densities;
}

std::vector<float>& Lattices::getVelocity(int x, int y)
{
    float vx;
    float vy;
    int index = (y+1)*(width+2)+(x+1);
    return getVelocity(index);
}

inline std::vector<float>& Lattices::getVelocity(int index)
{
    v[0] = (nE[index] - nW[index] + 0.7071*(nNE[index] + nSE[index]) - 0.7071*(nNW[index] + nSW[index]));
    v[1] = (nS[index] - nN[index] + 0.7071*(nSW[index] + nSE[index]) - 0.7071*(nNE[index] + nNW[index]));
    // cap the velocity at 1/12
    v[0] = v[0]/(6*fabs(v[0])+2);
    v[1] = v[1]/(6*fabs(v[1])+2);
    return v;
}

std::vector<float> Lattices::getEquilibNumDensities(int index)
{
    //int index = (y+1)*(width+2)+(x+1); // unused
    float density = getDensity(index);
    v = getVelocity(index);
    float vDotv = v[0]*v[0]+v[1]*v[1]; // v dotted with itself

    eqNumD[0] = density * 4.0/9.0 * (1 - (3.0/2.0)*vDotv);
    eqNumD[1] = density * 1.0/9.0 * (1 + 3*v[0] + (9.0/2.0)*v[0]*v[0] - (3.0/2.0)*vDotv);
    eqNumD[2] = density * 1.0/9.0 * (1 + 3*(0-v[1]) + (9.0/2.0)*v[1]*v[1] - (3.0/2.0)*vDotv);
    eqNumD[3] = density * 1.0/9.0 * (1 + 3*(0-v[0]) + (9.0/2.0)*v[0]*v[0] - (3.0/2.0)*vDotv);
    eqNumD[4] = density * 1.0/9.0 * (1 + 3*v[1] + (9.0/2.0)*v[1]*v[1] - (3.0/2.0)*vDotv);
    eqNumD[5] = density * 1.0/36.0 * (1 + 3*(v[0]-v[1]) + (9.0/2.0)*(v[0]*v[0]-v[1]*v[1]) - (3.0/2.0)*vDotv);
    eqNumD[6] = density * 1.0/36.0 * (1 + 3*(0-v[0]-v[1]) + (9.0/2.0)*(0-v[0]-v[1])*(0-v[0]-v[1]) - (3.0/2.0)*vDotv);
    eqNumD[7] = density * 1.0/36.0 * (1 + 3*(v[1]-v[0]) + (9.0/2.0)*(v[1]*v[1]-v[0]*v[0]) - (3.0/2.0)*vDotv);
    eqNumD[8] = density * 1.0/36.0 * (1 + 3*(v[0]+v[1]) + (9.0/2.0)*(v[0]+v[1])*(v[0]-v[1]) - (3.0/2.0)*vDotv);

    return eqNumD;
}

void Lattices::updateEqDensities()
{
    float nineHalfs = 9.0/2.0;
    float oneNinth = 1.0/9.0;
    float oneThirtySixth = 1.0/36.0;
    float density, vDotv, vDotvTimes3over2, vx2, vy2;
    int index = width+1;
    for (int i=0; i<height; ++i) {
        index += 2;
        for (int j=0; j<width; ++j) {
            ++index;
            density = getDensity(index);
            v = getVelocity(index);
            vDotv = v[0]*v[0]+v[1]*v[1]; // v dotted with itself
            vDotvTimes3over2 = vDotv * 3.0/2.0;
            vx2 = v[0]*v[0];
            vy2 = v[1]*v[1];

            eq0[index] = density * 4.0/9.0 * (1 - vDotvTimes3over2);
            eqE[index] = density * oneNinth * (1 + 3*v[0] + (nineHalfs)*vx2 - vDotvTimes3over2);
            eqN[index] = density * oneNinth * (1 + 3*(0-v[1]) + (nineHalfs)*vy2 - vDotvTimes3over2);
            eqW[index] = density * oneNinth * (1 + 3*(0-v[0]) + (nineHalfs)*vx2 - vDotvTimes3over2);
            eqS[index] = density * oneNinth * (1 + 3*v[1] + (nineHalfs)*vy2 - vDotvTimes3over2);
            eqNE[index] = density * oneThirtySixth * (1 + 3*(v[0]-v[1]) + (nineHalfs)*(v[0]-v[1])*(v[0]-v[1]) - vDotvTimes3over2);
            eqNW[index] = density * oneThirtySixth * (1 + 3*(0-v[0]-v[1]) + (nineHalfs)*(0-v[0]-v[1])*(0-v[0]-v[1]) - vDotvTimes3over2);
            eqSW[index] = density * oneThirtySixth * (1 + 3*(v[1]-v[0]) + (nineHalfs)*(v[1]-v[0])*(v[1]-v[0]) - vDotvTimes3over2);
            eqSE[index] = density * oneThirtySixth * (1 + 3*(v[0]+v[1]) + (nineHalfs)*(v[0]+v[1])*(v[0]+v[1]) - vDotvTimes3over2);
        }
    }
}

void Lattices::print(int x, int y)
{
    int index = (y+1)*(width+2)+(x+1);
    std::vector<float> v = getVelocity(x,y);
    std::cout << "index:" << index << std::endl;
    std::cout << "Lattice:" << std::endl
              << nNW[index] << "," << nN[index] << "," << nNE[index] << std::endl
              << nW[index]  << "," << n0[index] << "," << nE[index]  << std::endl
              << nSW[index] << "," << nS[index] << "," << nSE[index] << std::endl
              << "Velocity: [" << v[0] << ", " << v[1] << "]" << std::endl
              << "Density: " << getDensity(x,y) << std::endl;
}
