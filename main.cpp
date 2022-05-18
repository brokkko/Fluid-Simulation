#include <SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>
#include "src/Graphics.h"
#include "src/Simulation.h"
#include "src/Constants.h"
#include "src/ReadData.h"

double* getDensity(){
    double *data = nullptr; int dataSize = 0;
    ReadData reader("../data/bnd.nc");
    reader.readData("D", &data, &dataSize);
    return data;
}

double* getVelocity(){
    double *data = nullptr; int dataSize = 0;
    ReadData reader("../data/bnd.nc");
    reader.readData("V1", &data, &dataSize);
    return data;
}

double* getTemperature(){
    double *data = nullptr; int dataSize = 0;
    ReadData reader("../data/bnd.nc");
    reader.readData("T", &data, &dataSize);
    return data;
}

double* getMagneticField(){
    double *data = nullptr; int dataSize = 0;
    ReadData reader("../data/bnd.nc");
    reader.readData("B1", &data, &dataSize);
    return data;
}


int main(){

    double *densities = getDensity();
    double *vels = getVelocity();
    double *temperature = getTemperature();
    double *magneticField = getMagneticField();

    double time=0;
    bool paused=false;
    bool shift =false;
    double upperbound[9] ={small_rho*20,1000000,1,1,1,1,1,1,1};
    int currentmode = 0;
    sf::RenderWindow window(sf::VideoMode(1300, 700), "wave", sf::Style::Default, sf::ContextSettings(32));
    window.setActive(true);
    window.setFramerateLimit(60);

    SphericalGrid grid(SIZE_R,SIZE_PH,SIZE_TH,MIN_RADIUS,MAX_RADIUS,M_PI/3);
    Simulation simulation(grid);
    InitialConditions(grid);

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
                ApplyBoundaryConditions(grid,time,densities, vels, temperature, magneticField);

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