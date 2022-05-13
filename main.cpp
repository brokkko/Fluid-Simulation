#include <SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>
#include <sstream>
#include "src/Graphics.h"
#include "src/Grid.h"
#include "src/Simulation.h"
#include "src/Constants.h"
#include "src/ReadData.h"
// dU/dt = a*dU/dx

double* getDensity(){
    double *data = nullptr; int dataSize = 0;
    ReadData reader("../data/bnd.nc");
    reader.readData("D", &data, &dataSize);
    auto *density = new double [180];
    int index = 0;
    for(int i=0; i<180; i++){
        density[index++] = data[30+60*i];
    }
    free(data);
    return density;
}

double* getVelocity(){
    double *data = nullptr; int dataSize = 0;
    ReadData reader("../data/bnd.nc");
    reader.readData("V1", &data, &dataSize);
    auto *velocity = new double [180];
    int index = 0;
    for(int i=0; i<180; i++){
        velocity[index++] = data[30+60*i];
    }
    free(data);
    return velocity;
}


int main()
{
    double *densities= getDensity();
    double *vels= getVelocity();


    double time=0;
    bool paused=false;
    bool shift =false;
    double upperbound[9] ={small_rho*20,1000000,1,1,1,1,1,1,1};
    int currentmode = 0;
    sf::RenderWindow window(sf::VideoMode(700, 700), "wave", sf::Style::Default, sf::ContextSettings(32));
    window.setActive(true);
    window.setFramerateLimit(60);
    //Grid grid(70,70);
    SphericalGrid grid(SIZE_R,SIZE_PH,SIZE_TH,MIN_RADIUS,MAX_RADIUS,M_PI/3);
    Simulation simulation(grid);
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

        double dt = DT;
        window.clear(sf::Color::Black);
        if (!paused) {
            for (int i = 0; i < 1; i++)
            {
                simulation.RKIntegrator(dt,time);
                ApplyBoundaryConditions(grid,time,densities,vels);

            }

        }
        double sum = 0;
        double vel = 0;

        show(grid, window,t,upperbound[currentmode],currentmode, dt);

        window.display();
        //paused=true;
    }

    return 0;
}