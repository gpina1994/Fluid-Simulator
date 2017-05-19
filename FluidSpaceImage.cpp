#include "FluidSpaceImage.hpp"
#include <vector>
#include <math.h>

/*
void FluidSpaceImage::setFluidSpace(const FluidSpace& fluSpace)
{
    // to be implemented
}
*/

FluidSpaceImage::FluidSpaceImage(int w, int h)
{
    img.resize(w*h*4);
    densityVect.resize((w+2)*(h+2));
}

sf::Image FluidSpaceImage::getDensityImage(FluidSpace& fs)
{
    /**
    This function converts data from the FluidSpace object to an image
    usable by SFML to update the screen.
    **/
    int width = fs.getWidth();
    int height = fs.getHeight();
    densityVect = fs.getDensity();

    sf::Color color;
    int imgIndex;
    for (int i=0; i<height; ++i) {
        for (int j=0; j<width; ++j) {
            color = densityToColor(densityVect[i*width+j]);
            imgIndex = i*width*4 + j*4;
            img[imgIndex] = color.r;
            img[imgIndex + 1] = color.g;
            img[imgIndex + 2] = color.b;
            img[imgIndex + 3] = color.a;
        }
    }

    sf::Image outputImg;
    outputImg.create(width, height, &img[0]);
    return outputImg;
}

inline sf::Color FluidSpaceImage::densityToColor(float density)
{
    /// takes densities b/w 0 and infinity and returns a color
    // float g = 16*sqrt(density*25.5);
    //float g = (255.0*density)/(density+5.0);
    return sf::Color(0,(int)(255.0*density)/(density+40.0),0,255);
}





