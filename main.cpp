
#include <SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>
#include <sstream>
#include "Graphics.h"
#include "Simulation.h"
// dU/dt = a*dU/dx


int main()
{
    sf::RenderWindow window(sf::VideoMode(700, 700), "wave", sf::Style::Default, sf::ContextSettings(32));
    window.setActive(true);
    window.setFramerateLimit(60);
    //InitGL();
    Grid grid(70,70);
    grid.Fill(0.01);
    for (int y = 45; y < 55; y++){
        for (int x = 20; x < 40; x++){
          // grid.mesh[x][y].rho = 5;
        }
    }
    //grid.mesh[5][5].val = 10e-300;
    //grid.mesh[10].val = 50;
    //grid.mesh[11].val = 50;
    // dataArray[20] = 50;
    //dataArray[80] = 50;
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
        //solve(dataArray, arraySize, h);
        for(int i=0;i<5;i++)
            RKIntegrator(grid, h);
        double sum = 0;
        double vel = 0;
        for (int y = 0; y < grid.sizeY; y++)
        {
            for (int x = 0; x < grid.sizeX; x++){
                sum += grid.mesh[x][y].rho;
                vel+=std::pow( grid.mesh[x][y].vx,2)+std::pow( grid.mesh[x][y].vy,2);
            }
        }
       // printf("p = %0.10e\n",sum);
       // printf("vel = %0.10e\n",vel);
        // std::cout << dataArray[5] << "\n";
        //grid.Print();

        show(grid, window,t);

        window.display();
    }

    return 0;
}

