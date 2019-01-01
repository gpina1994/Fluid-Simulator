#include <SFML/Graphics.hpp>
#include "FluidSpaceImage.hpp"
//#include "FluidSpace.hpp"
#include <iostream>
#include <time.h>

/*
    0.0.1: This version generates a window and tracks (prints) where the mouse is clicked and unclicked
    0.0.2: In addition, there is now an output image coming from a class that will interface with the FluidSpace class
    0.0.3: The FluidSpace class is now generating its own dummy data that is being turned into an image by the
           FluidSpaceImage class.
           There is also now a display of the average iterations/second of the main loop after the program window is
           closed.
    0.0.4: The FluidSpace class is generating data from user input now. Streaming has been implemented, but not
           collisions. The enter key advances the simulation by one step.
    0.0.5: The Lattices class now computes velocity. The collide function is implemented. There was a bug with NaN
           values showing up in the lattices, which was resolved by dividing the output of getVelocity(). The density
           coloring function has been modified to take an input domain containing [0,inf). Basically the simulation
           works but needs optimization and tuning.
    0.1.0: Optimizing. Reduced size to 400x400. Speed of 0.03 s/update is typical (from 0.7 s/update). The cause of the
           NaN bug is still unknown, but fixing it is a goal for the next version. Dividing the velocity vector
           alleviates it, but does not remove it. Some tuning modifications. Also, density used to not propagate evenly
           in every direction--this is now fixed.
    0.2.0: Trying to figure out cause of NaN bug. Also, air propagates in an expanding square shape--this needs to be
           fixed.
    0.3.0: NaN bug fixed by dividing by density + adding a fudge factor. The shape of air propagation is not a
           problem--it's actually circular. Window can now be resized.
*/

/**
To do:
    Add in solid objects (draw mode)
    Add in pausing/unpausing
    Possibly change the updateEqDensities() method
    Add in some sort of "streamer"
    Fix bug on the west edge
**/

void handleMouse(const sf::RenderWindow&, FluidSpace& fs);
void displayFs(FluidSpace& fs, FluidSpaceImage& fsi, sf::RenderWindow& w);

int main()
{
    /// Initialize objects
    const int WINDOW_WIDTH = 900;
    const int WINDOW_HEIGHT = 900;
    const int FLUIDSPACE_WIDTH = 300;
    const int FLUIDSPACE_HEIGHT = 300;
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "SFML works!", sf::Style::Default);
    FluidSpaceImage fsi(FLUIDSPACE_WIDTH, FLUIDSPACE_HEIGHT);
    FluidSpace fs(FLUIDSPACE_WIDTH, FLUIDSPACE_HEIGHT);
    sf::Clock updateClk;
    bool drawMode = false; // drawmode: draw solid objects to interact with the fluid environment (UNIMPLEMENTED)
    bool paused = true;
    bool stepNextFrame = false;

    /// "Game loop"
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            sf::Vector2i pos; // mouse position
            switch (event.type)
            {
            case sf::Event::Closed:
            {
                window.close();
                break;
            }
            case sf::Event::MouseButtonPressed:
            {
                pos = sf::Mouse::getPosition(window);
                std::cout << "Mouse Button Pressed at " << pos.x << ", " << pos.y << std::endl;
                break;
            }
            case sf::Event::MouseButtonReleased:
            {
                pos = sf::Mouse::getPosition(window);
                std::cout << "Mouse Button Released at " << pos.x << "," << pos.y << std::endl;
                break;
            }
            case sf::Event::KeyPressed:
            {
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Return)) {
                    std::cout << "return key pressed" << std::endl;
                    paused = !paused;
//                    std::cout<<"updating"<<std::endl;
//                    updateClk.restart();
//                    fs.update();
//                    std::cout << "updated in " << updateClk.getElapsedTime().asSeconds() << "seconds." << std::endl;
//                    break;
                }
                else if(sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) {
                    stepNextFrame = true;
                }
                break;
            }
            default:
                break;
            }
        }

        if (!paused)
        {
            fs.update();
        }
        else if (stepNextFrame)
        {
            fs.update();
            stepNextFrame = false;
        }
        handleMouse(window, fs);
        displayFs(fs, fsi, window);
    }

    return 0;
}

void displayFs(FluidSpace& fs, FluidSpaceImage& fsi, sf::RenderWindow& w)
{
    /**
    Draws FluidSpace data onto the window
    **/

    // create sf::Sprite from FluidSpace
    sf::Sprite outSprite;
    sf::Texture outTexture;
    outTexture.loadFromImage(fsi.getDensityImage(fs));
    outSprite.setTexture(outTexture);

    //scale sprite to screen size
    sf::Vector2u wSize = w.getSize();
    double scaleFactorX = wSize.x/(double)fs.getWidth();
    double scaleFactorY = wSize.y/(double)fs.getHeight();
    outSprite.setScale(scaleFactorX, scaleFactorY);

    // draw sprite
    w.clear();
    w.draw(outSprite);
    w.display();
}

void handleMouse(const sf::RenderWindow& w, FluidSpace& fs)
{
    /**
    This function updates data on the density field through the FluidSpace object
    according to where the user has clicked.
    **/
    sf::Vector2i pos = sf::Mouse::getPosition(w);
    sf::Vector2u wSize = w.getSize();
    double scaleFactorX = (double)fs.getWidth()/wSize.x;
    double scaleFactorY = (double)fs.getHeight()/wSize.y;
    int scaledPosX, scaledPosY;
    scaledPosX = pos.x*scaleFactorX;
    scaledPosY = pos.y*scaleFactorY;
    // first check if mouse is even in the window
    if (0 < pos.x && pos.x < wSize.x &&
        0 < pos.y && pos.y < wSize.y) {
        // left mouse button -> increment density on clicked area
        if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
            fs.incrDensity(scaledPosX, scaledPosY);
        // right mouse button -> show lattice info
        } else if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
            fs.applyWind(scaledPosX, scaledPosY, 0.0);
            //fs.printLattice(scaledPosX, scaledPosY);
        }
    }
}
