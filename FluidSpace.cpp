#include "FluidSpace.hpp"
#include <stdlib.h> // rand, srand
#include <time.h> // seeding rand
#include <math.h> // isnan()

/**

Configuration details
    the fluid kinematic viscosity is proportional to the lattice speed (Lattices::latticeSpeed)
                                              and to 2*relaxation - 1
**/

namespace worldParams
{
    const double latticeSpeed = 1.0/sqrt(3); /// apparently the lattice speed *has* to be 1/sqrt(3) for D2Q9 lattices
                                            /// apparently this is also sqrt(3 * gas constant * temperature)
                                            /// but it's 1/sqrt(3) because our units are special ???
    const double relaxation = 0.6; /// relaxation of 1.0 gives a viscosity of 1/6; relaxation should be between 0.5 and inf
}

FluidSpace::FluidSpace()
{
    // nothing to do here
}

FluidSpace::FluidSpace(int w, int h)
{
    width = w;
    height = h;
    int vectSize = (w+2)*(h+2); // +2's account for the one cell thick border around the grid
    densityVect.resize(vectSize);
    velocityVect.resize(vectSize);

    // initialize lattices
    Lattices lats(w, h);
    latticesList.push_front(lats);
    latticesList.push_front(lats);

    // test code:
    std::cout << "nE lattice @ (100,100): " << latticesList.front().nE[(100+1)*(width+2)+(100+1)] << std::endl;

    for (int i=0; i<width+2; ++i) {
        for (int j=0; j<height+2; ++j) {
            densityVect[i*(height+2) + j] = 1.0;
        }
    }
}

std::vector<double>& FluidSpace::getDensity()
{
    return latticesList.front().updateDensities();
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
    /**
      Increments the density at a square of points centered at (x,y) on
      the cell grid
    **/
    //incrDensityPx(x, y);

    double radius = 4;
    for (int i=-4; i<5; ++i) {
        for (int j=-4; j<5; ++j) {
            if (i*i + j*j > radius*radius) continue;
            incrDensityPx(x+i,y+j);
        }
    }
}

void FluidSpace::incrDensityPx(int x, int y)
{
    /// check if position is valid (should just be a function)
    /// also, this class isn't the appropriate place to check that
    /// (since the screen may be resizable in future versions)
    if (!isValidCoords(x, y)) {
        return;
    }
    int index = xyToIndex(x, y);
    /// get original density
    //double incrAmt = (10.0 - orig)*0.60;

    double incrAmt = .05;
    latticesList.front().n0[index] += (4.0/9.0)*incrAmt;
    latticesList.front().nE[index] += (1.0/9.0)*incrAmt;
    latticesList.front().nN[index] += (1.0/9.0)*incrAmt;
    latticesList.front().nW[index] += (1.0/9.0)*incrAmt;
    latticesList.front().nS[index] += (1.0/9.0)*incrAmt;
    latticesList.front().nNE[index] += (1.0/36.0)*incrAmt;
    latticesList.front().nNW[index] += (1.0/36.0)*incrAmt;
    latticesList.front().nSW[index] += (1.0/36.0)*incrAmt;
    latticesList.front().nSE[index] += (1.0/36.0)*incrAmt;
}

void FluidSpace::applyWind(int x, int y, double directionRad)
{
    //applyWindPx(x, y, directionRad);
    double radius = 6.0;
    for (int i=-6; i<7; ++i) {
        for (int j=-6; j<7; ++j) {
            if (i*i + j*j > radius*radius) continue;
            applyWindPx(x+i, y+j, directionRad);
        }
    }
}

void FluidSpace::applyWindPx(int x, int y, double directionRad)
{
    if (!isValidCoords(x, y)) {
        return;
    }
    const int index = xyToIndex(x, y);
    Lattices& currentLattice = latticesList.front();
    const double density = currentLattice.getDensity(index);

    /// form vector in direction 'directionRad' of magnitude 'latticeSpeed' (speed of sound)
    double windX = cos(directionRad)*worldParams::latticeSpeed;
    double windY = sin(directionRad)*worldParams::latticeSpeed;

//    std::cout << "windX: " << windX << std::endl
//        << "windY: " << windY << std::endl;

    const double vDotv = windX * windX + windY * windY;


//    std::cout << "vDotv: " << vDotv << std::endl;

    const double latticeSpeed = worldParams::latticeSpeed;
    const double latticeSpeed2 = latticeSpeed * latticeSpeed;
    const double latticeSpeed4 = latticeSpeed2 * latticeSpeed2;

//    std::cout << "initial density: " << density << std::endl;
//    std::cout << "initial lattice:" << std::endl;
//    printLattice(index);

    currentLattice.n0[index] = density * 4.0/9.0 * (1 - (1/2.0)*vDotv / latticeSpeed2);
    currentLattice.nE[index] = density * (1.0/9.0) * (1 + windX/latticeSpeed2 + (1/2.0)*(windX*windX)/latticeSpeed4 - (1/2.0)*vDotv / latticeSpeed2);
    currentLattice.nN[index] = density * (1.0/9.0) * (1 + (0-windY)/latticeSpeed2 + (1/2.0)*(windY*windY)/latticeSpeed4 - (1/2.0)*vDotv / latticeSpeed2);
    currentLattice.nW[index] = density * (1.0/9.0) * (1 + (0-windX)/latticeSpeed2 + (1/2.0)*(windX*windX)/latticeSpeed4 - (1/2.0)*vDotv / latticeSpeed2);
    currentLattice.nS[index] = density * (1.0/9.0) * (1 + windY/latticeSpeed2 + (1/2.0)*(windY*windY)/latticeSpeed4 - (1/2.0)*vDotv / latticeSpeed2);
    currentLattice.nNE[index] = density * (1.0/36.0) * (1 + (windX-windY)/latticeSpeed2 + (1/2.0)*(windX-windY)*(windX-windY)/latticeSpeed4 - (1/2.0)*vDotv / latticeSpeed2);
    currentLattice.nNW[index] = density * (1.0/36.0) * (1 + (0-windX-windY)/latticeSpeed2 + (1/2.0)*(0-windX-windY)*(0-windX-windY)/latticeSpeed4 - (1/2.0)*vDotv / latticeSpeed2);
    currentLattice.nSW[index] = density * (1.0/36.0) * (1 + (windY-windX)/latticeSpeed2 + (1/2.0)*(windY-windX)*(windY-windX)/latticeSpeed4 - (1/2.0)*vDotv / latticeSpeed2);
    currentLattice.nSE[index] = density * (1.0/36.0) * (1 + (windX+windY)/latticeSpeed2 + (1/2.0)*(windX+windY)*(windX+windY)/latticeSpeed4 - (1/2.0)*vDotv / latticeSpeed2);

//    std::cout << "final density: " << currentLattice.getDensity(index) << std::endl;
//    std::cout << "final lattice:" << std::endl;
//    printLattice(index);
}

bool FluidSpace::isValidCoords(int x, int y)
{
    return !(x < 0 || width < x || y < 0 || height < y);
}

void FluidSpace::update()
{
    /// step 1: stream
    stream();

    //std::cout << "after streaming: ";
    //printTotalMass();

    ///step 2: collide
    collide();

    //std::cout << "after collision: ";
    //printTotalMass();

    return;
}

void FluidSpace::collide()
{
    /// called through update()
    /**
      This method simulates interactions between particles within a cell: every
      number stored in each lattice is brought closer to the computed equilibrium
      density for each space in the lattice.
    **/

    // update equilibrium densities
    latticesList.front().updateEqDensities();

    Lattices& currentLatticeGrid = latticesList.front();
    int ind = width+1;
    const double w = 1.0 / worldParams::relaxation;
    for (int y=0; y<height; ++y) {
        ind += 2;
        for (int x=0; x<width; ++x) {
            ++ind;

            double density = currentLatticeGrid.getDensity(ind);

            if (density > 2)
            {
                int wow = 0xfaded;
            }

            if (ind == -1)
            {
                std::cout << "lattice before collisions:" << std::endl;
                currentLatticeGrid.print(ind);
            }

            //currentLatticeGrid.n0[ind] = (1.0 - w)*currentLatticeGrid.n0[ind] + w*currentLatticeGrid.eq0[ind];
            currentLatticeGrid.n0[ind] = currentLatticeGrid.n0[ind] + (1/worldParams::relaxation)*(currentLatticeGrid.eq0[ind] - currentLatticeGrid.n0[ind]);
            currentLatticeGrid.nE[ind] = (1.0 - w)*currentLatticeGrid.nE[ind] + w*currentLatticeGrid.eqE[ind];
            currentLatticeGrid.nN[ind] = (1.0 - w)*currentLatticeGrid.nN[ind] + w*currentLatticeGrid.eqN[ind];
            currentLatticeGrid.nW[ind] = (1.0 - w)*currentLatticeGrid.nW[ind] + w*currentLatticeGrid.eqW[ind];
            currentLatticeGrid.nS[ind] = (1.0 - w)*currentLatticeGrid.nS[ind] + w*currentLatticeGrid.eqS[ind];
            currentLatticeGrid.nNE[ind] = (1.0 - w)*currentLatticeGrid.nNE[ind] + w*currentLatticeGrid.eqNE[ind];
            currentLatticeGrid.nNW[ind] = (1.0 - w)*currentLatticeGrid.nNW[ind] + w*currentLatticeGrid.eqNW[ind];
            currentLatticeGrid.nSW[ind] = (1.0 - w)*currentLatticeGrid.nSW[ind] + w*currentLatticeGrid.eqSW[ind];
            currentLatticeGrid.nSE[ind] = (1.0 - w)*currentLatticeGrid.nSE[ind] + w*currentLatticeGrid.eqSE[ind];

            if (ind == -1)
            {
                std::cout << "equilibrium lattice:" << std::endl;
                currentLatticeGrid.printEq(ind);
                std::cout << "lattice after collisions:" << std::endl;
                currentLatticeGrid.print(ind);
            }
        }
    }
}

void FluidSpace::printTotalMass()
{
    double totalMass = 0.0;
    Lattices& currentLatticeGrid = latticesList.front();

    int ind = width+1;
    for (int y=0; y<height; ++y) {
        ind += 2;
        for (int x=0; x<width; ++x) {
            ++ind;
            totalMass += currentLatticeGrid.getDensity(ind);
        }
    }

    std::cout << "total mass: " << totalMass << std::endl;
}

void FluidSpace::stream()
{ /// should be called through update()
    /**
      This method moves the densities stored in each lattice towards the direction
      they are moving (determined by what space they occupy in the lattice)
    **/

    streamE();
    streamN();
    streamW();
    streamS();
    streamNE();
    streamNW();
    streamSW();
    streamSE();
    stream0();

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

void FluidSpace::printLattice(int x, int y)
{
    latticesList.front().print(x,y);
}

void FluidSpace::printLattice(int index)
{
    latticesList.front().print(index);
}

inline int FluidSpace::xyToIndex(int x, int y)
{
    return (y+1)*(width+2)+(x+1);
}

Lattices::Lattices(int w, int h)
{
    latticeSpeed = worldParams::latticeSpeed;
    width = w;
    height = h;
    int vectSize = (w+2)*(h+2); // Grid is buffered on each edge by 1 lattice
    n0.resize(vectSize, (double)4/9);
    nE.resize(vectSize, (double)1/9);
    nN.resize(vectSize, (double)1/9);
    nW.resize(vectSize, (double)1/9);
    nS.resize(vectSize, (double)1/9);
    nNE.resize(vectSize, (double)1/36);
    nNW.resize(vectSize, (double)1/36);
    nSW.resize(vectSize, (double)1/36);
    nSE.resize(vectSize, (double)1/36);

    eq0.resize(vectSize);
    eqE.resize(vectSize);
    eqN.resize(vectSize);
    eqW.resize(vectSize);
    eqS.resize(vectSize);
    eqNE.resize(vectSize);
    eqNW.resize(vectSize);
    eqSW.resize(vectSize);
    eqSE.resize(vectSize);


    int index = (150+1)*(width+2)+(150+1);
    n0[index] = 400/9.;
    nE[index] = 100/9.;
    nN[index] = 100/9.;
    nW[index] = 100/9.;
    nS[index] = 100/9.;
    nNE[index] = 100/36.;
    nNW[index] = 100/36.;
    nSW[index] = 100/36.;
    nSE[index] = 100/36.;

    eqNumD.resize(9);
    v.resize(2);
    densities.resize(w*h);
}

void Lattices::resetBorder()
{
    /**
      Called every iteration to keep the outside border constant. The values are
      equilibrium values for 0 velocity cells of density 1.
    **/
    //set top edge
    for (int i=0; i<width+2; ++i) {
        n0[i] = (double)4/9;
        nE[i] = (double)1/9;
        nN[i] = (double)1/9;
        nW[i] = (double)1/9;
        nS[i] = (double)1/9;
        nNE[i] = (double)1/36;
        nNW[i] = (double)1/36;
        nSW[i] = (double)1/36;
        nSE[i] = (double)1/36;
    }
    //set left edge
    for (int i=0; i<(width+2)*(height+2); i += width+2) {
        n0[i] = (double)4/9;
        nE[i] = (double)1/9;
        nN[i] = (double)1/9;
        nW[i] = (double)1/9;
        nS[i] = (double)1/9;
        nNE[i] = (double)1/36;
        nNW[i] = (double)1/36;
        nSW[i] = (double)1/36;
        nSE[i] = (double)1/36;
    }
    //set right edge
    for (int i=width+1; i<(width+2)*(height+2); i += width+2) {
        n0[i] = (double)4/9;
        nE[i] = (double)1/9;
        nN[i] = (double)1/9;
        nW[i] = (double)1/9;
        nS[i] = (double)1/9;
        nNE[i] = (double)1/36;
        nNW[i] = (double)1/36;
        nSW[i] = (double)1/36;
        nSE[i] = (double)1/36;
    }
    //set bottom edge
    for (int i=(width+2)*(height+1); i<(width+2)*(height+2); ++i) {
        n0[i] = (double)4/9;
        nE[i] = (double)1/9;
        nN[i] = (double)1/9;
        nW[i] = (double)1/9;
        nS[i] = (double)1/9;
        nNE[i] = (double)1/36;
        nNW[i] = (double)1/36;
        nSW[i] = (double)1/36;
        nSE[i] = (double)1/36;
    }
}

void Lattices::resetToZeroes()
{
    /**
      Resets each lattice to 0's for immediate overwriting by the stream methods.
    **/
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

inline double Lattices::getDensity(int x, int y)
{
    int index = (y+1)*(width+2)+(x+1);
    return getDensity(index);
}

inline double Lattices::getDensity(int index)
{
    return n0[index] + nE[index] + nN[index] + nW[index] + nS[index] +\
           nNE[index] + nNW[index] + nSW[index] + nSE[index];
}

std::vector<double>& Lattices::updateDensities()
{
    int ind;
    for (int i=0; i<height; ++i) {
        for (int j=0; j<width; ++j) {
            ind = (i+1)*(width+2)+(j+1);
            densities[i*width+j] = getDensity(ind);
        }
    }
    return densities;
}

std::vector<double>& Lattices::getVelocity(int x, int y)
{
    int index = (y+1)*(width+2)+(x+1);
    return getVelocity(index);
}

inline std::vector<double>& Lattices::getVelocity(int index)
{
    double density = getDensity(index);

    if (density != 0.0)
    {
//        v[0] = (1/9.)*(nE[index] - nW[index]) + (1/36.)*(nNE[index] + nSE[index] - nNW[index] - nSW[index]);
//        v[1] = (1/9.)*(nS[index] - nN[index]) + (1/36.)*(nSW[index] + nSE[index] - nNE[index] - nNW[index]);
        v[0] = (nE[index] - nW[index]) + (nNE[index] + nSE[index] - nNW[index] - nSW[index]);
        v[1] = (nS[index] - nN[index]) + (nSW[index] + nSE[index] - nNE[index] - nNW[index]);
        v[0] /= density;
        v[1] /= density;

        /// capping to the lattice speed here is contrary to the algorithm (and therefore probably wrong from the standpoint of
        /// simulation accuracy), but it is being done for stability and also makes *some* sense because the velocity of the
        /// fluid *cant* exceed the speed of sound (1/sqrt(3)); without this, the max calculable velocity is 1.
        const double speed2 = v[0]*v[0] + v[1]*v[1];
        const double latticeSpeed2 = worldParams::latticeSpeed*worldParams::latticeSpeed;
        if (speed2 > latticeSpeed2)
        {
            v[0] *= worldParams::latticeSpeed / sqrt(speed2);
            v[1] *= worldParams::latticeSpeed / sqrt(speed2);
        }
    }
    else
    {
        v[0] = 0.0;
        v[1] = 0.0;
    }
    return v;
}

void Lattices::updateEqDensities()
{
    /**
      Updates the equilibrium densities for every lattice based on the
      computed velocity of each lattice
    **/
    const double latticeSpeed2 = latticeSpeed * latticeSpeed;
    const double latticeSpeed4 = latticeSpeed2 * latticeSpeed2;
    int index = width+1;
    for (int i=0; i<height; ++i) {
        index += 2;
        for (int j=0; j<width; ++j) {
            ++index;

            const double density = getDensity(index);
            v = getVelocity(index);
            const double vx = v[0];
            const double vy = v[1];
            const double vDotv = vx*vx + vy*vy; // v dotted with itself
            const double vx2 = vx*vx;
            const double vy2 = vy*vy;

            eq0[index]  = density*4/9. *(1                                                                        -  (1/2.)*vDotv/latticeSpeed2);
            eqE[index]  = density*1/9. *(1  +  vx      /latticeSpeed2  +  (1/2.)*vx2              /latticeSpeed4  -  (1/2.)*vDotv/latticeSpeed2);
            eqN[index]  = density*1/9. *(1  +  -vy     /latticeSpeed2  +  (1/2.)*vy2              /latticeSpeed4  -  (1/2.)*vDotv/latticeSpeed2);
            eqW[index]  = density*1/9. *(1  +  -vx     /latticeSpeed2  +  (1/2.)*vx2              /latticeSpeed4  -  (1/2.)*vDotv/latticeSpeed2);
            eqS[index]  = density*1/9. *(1  +  vy      /latticeSpeed2  +  (1/2.)*vy2              /latticeSpeed4  -  (1/2.)*vDotv/latticeSpeed2);
            eqNE[index] = density*1/36.*(1  +  (vx-vy) /latticeSpeed2  +  (1/2.)*(vx-vy)*(vx-vy)  /latticeSpeed4  -  (1/2.)*vDotv/latticeSpeed2);
            eqNW[index] = density*1/36.*(1  +  (-vx-vy)/latticeSpeed2  +  (1/2.)*(-vx-vy)*(-vx-vy)/latticeSpeed4  -  (1/2.)*vDotv/latticeSpeed2);
            eqSW[index] = density*1/36.*(1  +  (vy-vx) /latticeSpeed2  +  (1/2.)*(vy-vx)*(vy-vx)  /latticeSpeed4  -  (1/2.)*vDotv/latticeSpeed2);
            eqSE[index] = density*1/36.*(1  +  (vx+vy) /latticeSpeed2  +  (1/2.)*(vx+vy)*(vx+vy)  /latticeSpeed4  -  (1/2.)*vDotv/latticeSpeed2);
        }
    }
}

void Lattices::print(int x, int y)
{
    int index = (y+1)*(width+2)+(x+1);
    print(index);
}

void Lattices::print(int index)
{
    std::vector<double> v = getVelocity(index);
    std::cout << "index:" << index << std::endl;
    std::cout << "Lattice:" << std::endl
              << nNW[index] << "," << nN[index] << "," << nNE[index] << std::endl
              << nW[index]  << "," << n0[index] << "," << nE[index]  << std::endl
              << nSW[index] << "," << nS[index] << "," << nSE[index] << std::endl
              << "Velocity: [" << v[0] << ", " << v[1] << "]" << std::endl
              << "Density: " << getDensity(index) << std::endl;
}

void Lattices::printEq(int index)
{
    std::vector<double> v = getVelocity(index);
    std::cout << "index:" << index << std::endl;
    std::cout << "Lattice:" << std::endl
              << eqNW[index] << "," << eqN[index] << "," << eqNE[index] << std::endl
              << eqW[index]  << "," << eq0[index] << "," << eqE[index]  << std::endl
              << eqSW[index] << "," << eqS[index] << "," << eqSE[index] << std::endl
              << "Velocity: [" << v[0] << ", " << v[1] << "]" << std::endl
              << "Density: " << getDensity(index) << std::endl;
}
