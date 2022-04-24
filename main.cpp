#include <SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>
#include <sstream>
#include "src/Graphics.h"
#include "src/Grid.h"
#include "src/Simulation.h"
#include "src/Constants.h"
// dU/dt = a*dU/dx


int main()
{
    double time=0;
    bool paused=false;
    bool shift =false;
    double upperbound[9] ={0.02,1,1,1,1,1,1,1e11,100000};
    int currentmode = 7;
    sf::RenderWindow window(sf::VideoMode(700, 700), "wave", sf::Style::Default, sf::ContextSettings(32));
    window.setActive(true);
    window.setFramerateLimit(60);
    //Grid grid(70,70);
    SphericalGrid grid(30,30,1,1.497131e10,2.28e11,M_PI_2);
    InitialConditions(grid);
//grid.Fill(10);
    sf::View w;
    w = window.getDefaultView();
    w.setCenter(0, 0);
    window.setView(w);
    sf::Font f;
    f.loadFromFile("arial.ttf");

    sf::Text t;
    t.setFont(f);
    t.setCharacterSize(12);
    t.setFillColor(sf::Color::White);
    t.setOutlineColor(sf::Color::Black);
    t.setOutlineThickness(1);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            else if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Space)
                paused = !paused;
            else if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::LShift)
                shift = true;
            else if (event.type == sf::Event::KeyReleased && event.key.code == sf::Keyboard::LShift)
                shift = false;
            else if (event.type == sf::Event::MouseWheelMoved && !shift) {
                int delta = event.mouseWheel.delta;
                while(delta>0) { upperbound[currentmode] *= 1.1; delta--;}
                while(delta<0)
                {upperbound[currentmode] /= 1.1; delta++;}
            }
            else if (event.type == sf::Event::MouseWheelMoved && shift) {
                currentmode += event.mouseWheel.delta;
                currentmode = clamp(currentmode,0, 8);
            }
        }


        window.clear(sf::Color::Black);
        if (!paused) {
            for (int i = 0; i < 1; i++)
            {
                RKIntegrator(grid, DT,time);
            }

        }
        double sum = 0;
        double vel = 0;

        show(grid, window,t,upperbound[currentmode],currentmode);

        window.display();
    }

    return 0;
}