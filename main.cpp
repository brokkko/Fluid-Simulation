#include <SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>
#include <sstream>
#include "src/Graphics.h"
#include "src/Grid.h"
#include "src/Simulation.h"
// dU/dt = a*dU/dx


int main()
{
    sf::RenderWindow window(sf::VideoMode(700, 700), "wave", sf::Style::Default, sf::ContextSettings(32));
    window.setActive(true);
    window.setFramerateLimit(60);
    //InitGL();
    Grid grid(70,70);
    InitialConditions(grid);
    //grid.Fill(0.01);
    /*for (int y = 45; y < 55; y++){
        for (int x = 20; x < 40; x++){
            grid.mesh[x][y].rho=5;
            grid.mesh[x][y].vx=10;
            grid.mesh[x][y].vy=0;
            grid.mesh[x][y].Bz=0;
            grid.mesh[x][y].E = 3000;
        }
    }*/

    double h = 0.001;
    sf::View w;
    w = window.getDefaultView();
    w.setCenter(0, 0);
    window.setView(w);
    sf::Font f;
    f.loadFromFile("arial.ttf");

    sf::Text t;
    t.setFont(f);
    t.setCharacterSize(12);
    t.setFillColor(sf::Color::Magenta);

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }


        window.clear(sf::Color::Black);
        for(int i=0;i<1;i++)
            RKIntegrator(grid, h);
        double sum = 0;
        double vel = 0;
        for (int y = 0; y < grid.sizeY; y++)
        {
            for (int x = 0; x < grid.sizeX; x++){
                sum += grid.mesh[x][y].E;
               // vel+=grid.mesh[x][y].E);

            }
        }
        printf("E = %0.10e\n",sum);
        // printf("vel = %0.10e\n",vel);

        show(grid, window,t);

        window.display();
    }

    return 0;
}