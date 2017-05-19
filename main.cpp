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
*/

void handleMouse(const sf::RenderWindow&, FluidSpace& fs);

int main()
{
    /// Initialize objects
    sf::RenderWindow window(sf::VideoMode(400,400), "SFML works!", sf::Style::Titlebar|sf::Style::Close);
    FluidSpaceImage fsi(400,400);
    sf::Image output;
    sf::Texture outTexture;
    sf::Sprite outSprite;
    srand(time(NULL));
    FluidSpace fs(400,400);
    sf::Clock clk;
    sf::Clock updateClk;
    int iterations = 0;

    /// "Game loop"
    while (window.isOpen())
    {
        ++iterations;
        sf::Event event;
        bool mouseIsPressed = false; // unused
        bool enterIsPressed = false; // unused
        bool updateScreen = false;
        while (window.pollEvent(event))
        {
            sf::Vector2i pos; // mouse position
            switch (event.type)
            {
            case sf::Event::Closed:
                window.close();
                break;
            case sf::Event::MouseButtonPressed:
                mouseIsPressed = true;
                pos = sf::Mouse::getPosition(window);
                std::cout << "Mouse Button Pressed at " << pos.x << ", " << pos.y << std::endl;
                break;
            case sf::Event::MouseButtonReleased:
                mouseIsPressed = false;
                pos = sf::Mouse::getPosition(window);
                std::cout << "Mouse Button Released at " << pos.x << "," << pos.y << std::endl;
                break;
            case sf::Event::KeyPressed:
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Return)) {
                    std::cout<<"updating"<<std::endl;
                    updateClk.restart();
                    fs.update();
                    std::cout << "updated in " << updateClk.getElapsedTime().asSeconds() << "seconds." << std::endl;
                    updateScreen = true;
                    break;
                }
                break;
            default:
                break;
            }
            if (updateScreen) {
                updateScreen = false;
                break;
            }
        }

        handleMouse(window, fs);

        //output = fsi.getDensityImage();
        output = fsi.getDensityImage(fs);
        outTexture.loadFromImage(output);
        outSprite.setTexture(outTexture);

        window.clear();
        window.draw(outSprite);
        window.display();
    }

    std::cout << "Average iterations per second: " <<  (float) iterations/clk.getElapsedTime().asSeconds()
              << std::endl;
/*
    std::cout << "=======Average times=======" << std::endl
              << "stream(): " << streamTotalTime/streamCallCount << std::endl
              << "collide(): " << collideTotalTime/collideCallCount << std::endl;
*/
    return 0;
}

void handleMouse(const sf::RenderWindow& w, FluidSpace& fs)
{
    /**
    This function updates data on the density field through the FluidSpace object
    according to where the user has clicked.
    **/
    sf::Vector2i pos = sf::Mouse::getPosition(w);
    sf::Vector2u wSize = w.getSize();
    //if the mouse button is pressed in the window, draw where it is pressed
    if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
        if (0 < pos.x && pos.x < wSize.x &&
            0 < pos.y && pos.y < wSize.y) {
            // draw
            fs.incrDensity(pos.x, pos.y);
        }
    } else if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
        if (0 < pos.x && pos.x < wSize.x &&
            0 < pos.y && pos.y < wSize.y) {
            fs.printLattice(pos.x,pos.y);

        }
    }
}
