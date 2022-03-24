//#include <SFML/Window.hpp>
//#include <GL/gl.h>
//#include <GL/glu.h>
//#include "src/ReadData.h"
//#include <iostream>
//#include <SFML/Graphics.hpp>
//#include <algorithm>
//#include <cmath>
//
//#include "src/Grid.h"
//#include "src/Simulation.h"
//
//struct Vector3D{
//    double x, y, z;
//};
//
//Vector3D Spherical(double r, double theta, double phi)
//{
//    Vector3D pt = {0, 0, 0};
//    double snt = (float)std::sin(theta);
//    double cnt = (float)std::cos(theta);
//    double snp = (float)std::sin(phi);
//    double cnp = (float)std::cos(phi);
//    pt.x = r * snt * cnp;
//    pt.y = r * cnt;
//    pt.z = -r * snt * snp;
//    return pt;
//}
//
//sf::Color toColor(double val,double min,double max)
//{
//
//    double mid = (min+max)/2;
//    if (std::isnan(val)){
//        return sf::Color::Magenta;
//    }
//    double r=std::min(1.0,1-(mid-val)/(max-mid));
//    double g=1-std::abs(val-mid)/(max-mid);
//    double b=std::min(1.0,1-(val-mid)/(max-mid));
//    double n=(val-min)/(max-min);
////    double r=-5*n*n+7.5*n+1-5.0/1.777777777;
////    double g=-20*n*n+20*n-4;
////    double b=-5*n*n+2.5*n+0.6875;
//    sf::Color c((int)std::max(0.0,r*255),(int)std::max(0.0,g*255),(int)std::max(0.0,b*255));
//    return c;
//}
//
//
//void cameraControl(sf::Window *window, bool *firstTouch, double *ax, double *ay, sf::Vector2i *lastPosition, double *dist){
//    if (sf::Mouse::isButtonPressed(sf::Mouse::Left)){
//        window->setMouseCursorVisible(false);
//        sf::Vector3<int> pos;
//        pos.x = sf::Mouse::getPosition().x - lastPosition->x;
//        pos.y = sf::Mouse::getPosition().y - lastPosition->y;
//
//        lastPosition->x = sf::Mouse::getPosition().x;
//        lastPosition->y = sf::Mouse::getPosition().y;
//
//        if(*firstTouch){
//            pos.x = 0; pos.y = 0;
//            *firstTouch = false;
//        }
//
//        *ax += -(double) pos.x / 300;
//        *ay += -(double) pos.y / 300;
//
//
//    } else{
//        *firstTouch = true;
//        window->setMouseCursorVisible(true);
//    }
//
//    gluLookAt(0 + *dist * sin(*ax)* cos(*ay), 0 + *dist * sin(*ax)*sin(*ay), 0 + *dist * cos(*ax),
//              0, 0, 0, 0, 1, 0);
//}
//
//
//void drawSphere(Grid &grid){
//    glPushMatrix();
//    glTranslatef(0, 0, 0);
//    GLUquadric* gluq = gluNewQuadric();
//    glColor3f(255.0, 242.0, 0.0);
//    gluSphere(gluq, 4, 10, 10);
//    glPopMatrix();
//
//    glPointSize(2);
//    glBegin(GL_POINTS);
//
//    int maxLongitude = 180/3; // (1) количество секторов на круге
//    int maxLatitude = 60/3; // (2)
//    int numOfSectors = 60; // (3)
//    Vector3D dataHelper [numOfSectors+1][maxLatitude+1][maxLongitude+1];
//    double radius = 6;
//    int tmp=0;
//    while(tmp <= numOfSectors) {
//        glColor3f(255.0, 255.0, 255.0);
//        for (int i = 0; i < maxLatitude; i++) {
//            for (int j = 0; j < maxLongitude; j++) {
//                double phi = ((j)*2*M_PI) / maxLongitude;
//                double theta = ((i)*M_PI) / maxLatitude;
//                Vector3D coord = Spherical(radius, theta, phi);
//                //glVertex3d(coord.x, coord.y, coord.z);
//                dataHelper[tmp][i][j] = {coord.x, coord.y, coord.z};
//            }
//        }
//        radius+=4;
//        tmp++;
//    }
//
//    glEnd();
//
//    glBegin(GL_QUADS);
//    glLineWidth(2);
//    glColor3f(255.0, 255.0, 255.0);
//    for(int i=0; i<numOfSectors; i++){
//        //for(int y=0; y<maxLatitude; y++) {
//            int y = maxLatitude/2;
//            for (int x = 0; x < maxLongitude-1; x++) {
//                sf::Color col = toColor(grid.mesh[i][x].rho, 0, 5);
//                glColor3f(col.r/255.0, col.g/255.0, col.b/255.0);
//
//                glVertex3d(dataHelper[i][y][x].x, dataHelper[i][y][x].y, dataHelper[i][y][x].z);
//                glVertex3d(dataHelper[i+1][y][x].x, dataHelper[i+1][y][x].y, dataHelper[i+1][y][x].z);
//                glVertex3d(dataHelper[i+1][y][x+1].x, dataHelper[i+1][y][x+1].y, dataHelper[i+1][y][x+1].z);
//                glVertex3d(dataHelper[i][y][x + 1].x, dataHelper[i][y][x + 1].y, dataHelper[i][y][x + 1].z);
//
//            }
//            glVertex3d(dataHelper[i][y][maxLongitude-1].x, dataHelper[i][y][maxLongitude-1].y, dataHelper[i][y][maxLongitude-1].z);
//            glVertex3d(dataHelper[i+1][y][maxLongitude-1].x, dataHelper[i+1][y][maxLongitude-1].y, dataHelper[i+1][y][maxLongitude-1].z);
//            glVertex3d(dataHelper[i+1][y][0].x, dataHelper[i+1][y][0].y, dataHelper[i+1][y][0].z);
//            glVertex3d(dataHelper[i][y][0].x, dataHelper[i][y][0].y, dataHelper[i][y][0].z);
//
//        //}
//
//    }
//    glEnd();
//
//}
//
//
//void InitGL()
//{
//    sf::Window window(sf::VideoMode(800, 600), "OpenGL", sf::Style::Default, sf::ContextSettings(32));
//
//    window.setVerticalSyncEnabled(true);
//
//    // activate the window
//    window.setActive(true);
//
//    // load resources, initialize the OpenGL states, ...
//    glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // устанавливаем фоновый цвет на черный
//
//    glClearDepth(1.0);
//    glDepthFunc(GL_LESS);
//    glEnable(GL_DEPTH_TEST); // включаем тест глубины
//    glShadeModel(GL_SMOOTH);
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    gluPerspective(45.0f, (float)800 / (float)600, 0.1f, 1000.0f); // настраиваем трехмерную перспективу
//    glMatrixMode(GL_MODELVIEW);
//
//    // draw...
//    Grid grid(60,60); //TODO: grid size -> lat and lon
//    grid.Fill(0.01);
//    //grid.mesh[5][5].val = 10e-300;
//    //grid.mesh[10].val = 50;
//    //grid.mesh[11].val = 50;
//    // dataArray[20] = 50;
//    //dataArray[80] = 50;
//    double h = 0.0005;
//
//    double ax = 0, ay = 0;
//    bool firstTouch;
//    sf::Vector2i lastPosition = {0, 0};
//    double dist = 500;
//
//    while (window.isOpen())
//    {
//        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//        sf::Event event;
//        while (window.pollEvent(event))
//        {
//            if (event.type == sf::Event::Closed)
//                window.close();
//            if(event.type == sf::Event::MouseWheelMoved)
//            {
//                // display number of ticks mouse wheel has moved
//                dist -= event.mouseWheel.delta*2;
//            }
//        }
//
//        glPushMatrix();
//
//        cameraControl(&window, &firstTouch, &ax, &ay, &lastPosition, &dist);
//        drawSphere(grid);
//
//        for(int i=0;i<5;i++)
//            RKIntegrator(grid, h);
//        double sum = 0;
//        double vel = 0;
//        for (int y = 0; y < grid.sizeY; y++)
//        {
//            for (int x = 0; x < grid.sizeX; x++){
//                sum += grid.mesh[x][y].rho;
//                vel+=std::pow( grid.mesh[x][y].vx,2)+std::pow( grid.mesh[x][y].vy,2);
//            }
//        }
//
//        glPopMatrix();
//        // end the current frame (internally swaps the front and back buffers)
//
//
//        window.display();
//    }
//
//}
//
//
//
//int main()
//{
//    InitGL();
//
//    return 0;
//}
//
//
////int main(){
////    InitGL();
//////    ReadData reader("../data/bnd.nc");
//////    double *data = nullptr; int dataSize = 0;
//////    reader.readData("D", &data, &dataSize);
//////    for(int i=0; i<dataSize; i++){
//////        std::cout << data[i]*10000 << " ";
//////    }
////
////    return 0;
////}




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
    window.setFramerateLimit(20);
    //InitGL();
    Grid grid(70,70);
    grid.Fill(0.01);
//    for (int y = 45; y < 55; y++){
//        for (int x = 20; x < 40; x++){
//            grid.mesh[x][y].rho=5;
//            grid.mesh[x][y].vx=10;
//            grid.mesh[x][y].vy=0;
//            grid.mesh[x][y].Bz=0;
//        }
//    }
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
        for(int i=0;i<1;i++)
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