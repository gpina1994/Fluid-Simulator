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

inline sf::Color FluidSpaceImage::densityToColor(double density)
{
    /// takes densities b/w 0 and infinity and returns a color
    // double g = 16*sqrt(density*25.5);
    //double g = (255.0*density)/(density+5.0);
    const double maxDisplayedDensity = 1.0005;
    const double minDisplayedDensity = 0.999;
    const double neutralDensity = 1.0;

    if (density > maxDisplayedDensity) density = maxDisplayedDensity;
    if (density < minDisplayedDensity) density = minDisplayedDensity;

    //double colorIndex = density/maxDisplayedDensity;

    double green = (density - neutralDensity)/(maxDisplayedDensity - neutralDensity);
    double blue = 1.0 - (density - minDisplayedDensity)/(neutralDensity - minDisplayedDensity);

    if (blue < 0.0) blue = 0.0;
    if (green < 0.0) green = 0.0;
    if (blue > 1.0) blue = 1.0;
    if (green > 1.0) green = 1.0;

    green *= 255.0;
    blue *= 255.0;

    return sf::Color(0,green,blue,255);
}





