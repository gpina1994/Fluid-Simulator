#ifndef FLUIDSPACEIMAGE_HPP
#define FLUIDSPACEIMAGE_HPP

#include <SFML/Graphics.hpp>
#include "FluidSpace.hpp"

//class FluidSpace;

class FluidSpaceImage {
/* converts data from FluidSpace class in the form of 2d vectors
   to image data in the form of sf::Image's
*/
    //FluidSpace& fs;
    sf::Color densityToColor(double density);
    std::vector<sf::Uint8> img;
    std::vector<double> densityVect;
public:
    FluidSpaceImage() {};
    FluidSpaceImage(int w, int h);
    //void setFluidSpace(const FluidSpace& fs);
    sf::Image getDensityImage(void); // dummy function
    sf::Image getDensityImage(FluidSpace& fs);
};

#endif // FLUIDSPACEIMAGE_HPP
