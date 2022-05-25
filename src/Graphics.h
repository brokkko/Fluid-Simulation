//
// Created by alex on 18.03.2022.
//

#ifndef FLUID_SIMULATION_GRAPHICS_H
#define FLUID_SIMULATION_GRAPHICS_H
#include <SFML/Graphics.hpp>
#include "Grid.h"
#include "SphericalGrid.h"
#include <sstream>
#include <iostream>
#include <cmath>



sf::Color toColor(double val,double min,double max);
void show(SphericalGrid& grid, sf::RenderWindow& window,sf::Text& t,double upperbound,int mode, double dt);
//void show(Grid& grid, sf::RenderWindow& window,sf::Text& t, double upperbound,int mode);
#endif //FLUID_SIMULATION_GRAPHICS_H