#include <SFML/Window.hpp>
#include <GL/gl.h>
#include <GL/glu.h>
#include "src/ReadData.h"
#include <iostream>
#include <SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>

struct Vector3D{
    double x, y, z;
};

Vector3D Spherical(double r, double theta, double phi)
{
    Vector3D pt = {0, 0, 0};
    double snt = (float)std::sin(theta);// * M_PI / 180);
    double cnt = (float)std::cos(theta);// * M_PI / 180);
    double snp = (float)std::sin(phi);// * M_PI / 180);
    double cnp = (float)std::cos(phi);// * M_PI / 180);
    pt.x = r * snt * cnp;
    pt.y = r * cnt;
    pt.z = -r * snt * snp;
    return pt;
}

double mapFunction(int n, double start1, double stop1, double start2, double stop2) {
    return ((n-start1)/(stop1-start1))*(stop2-start2)+start2;
};

sf::Color toColor(double val,double min,double max)
{

    double mid = (min+max)/2;
    if (std::isnan(val)){
        return sf::Color::Magenta;
    }
    // double r=std::min(1.0,1-(mid-val)/(max-mid));
    //double g=1-std::abs(val-mid)/(max-mid);
    // double b=std::min(1.0,1-(val-mid)/(max-mid));
    double n=(val-min)/(max-min);
    double r=-5*n*n+7.5*n+1-5.0/1.777777777;
    double g=-20*n*n+20*n-4;
    double b=-5*n*n+2.5*n+0.6875;
    sf::Color c((int)std::max(0.0,r*255),(int)std::max(0.0,g*255),(int)std::max(0.0,b*255));
    return c;
}


void cameraControl(sf::Window *window, bool *firstTouch, double *ax, double *ay, sf::Vector2i *lastPosition, double *dist){
    if (sf::Mouse::isButtonPressed(sf::Mouse::Left)){
        window->setMouseCursorVisible(false);
        sf::Vector3<int> pos;
        pos.x = sf::Mouse::getPosition().x - lastPosition->x;
        pos.y = sf::Mouse::getPosition().y - lastPosition->y;

        lastPosition->x = sf::Mouse::getPosition().x;
        lastPosition->y = sf::Mouse::getPosition().y;

        if(*firstTouch){
            pos.x = 0; pos.y = 0;
            *firstTouch = false;
        }

        *ax += -(double) pos.x / 300;
        *ay += -(double) pos.y / 300;


    } else{
        *firstTouch = true;
        window->setMouseCursorVisible(true);
    }

    gluLookAt(0 + *dist * sin(*ax)* cos(*ay), 0 + *dist * sin(*ax)*sin(*ay), 0 + *dist * cos(*ax),
              0, 0, 0, 0, 1, 0);
}

// dU/dt = a*dU/dx


struct Cell
{
    double val;
    double vx;
    double vy;
    double B;

    Cell operator+(const Cell r) const
    {
        return Cell{ val + r.val ,vx + r.vx,vy + r.vy,B+r.B};
    }
    Cell operator-(const Cell r) const
    {
        return Cell{ val - r.val ,vx - r.vx,vy - r.vy,B-r.B};
    }
    Cell operator*(const Cell r) const
    {
        return Cell{ val * r.val ,vx * r.vx,vy * r.vy,B*r.B};
    }
    Cell operator/(const Cell r) const
    {
        return Cell{ val / r.val ,vx / r.vx,vy / r.vy,B/r.B};
    }
    Cell operator*(const double r) const
    {
        return Cell{ val * r ,vx*r,vy*r,B*r};
    }
    Cell operator/(const double r) const
    {
        return Cell{ val / r ,vx/r,vy/r,B/r};
    }
    friend Cell operator*(const double l, const Cell r)
    {
        return { l * r.val ,l*r.vx,l*r.vy,l*r.B};
    }
};


struct Grid
{
    Cell** mesh;
    Cell** fluxMesh;
    unsigned int sizeX;
    unsigned int sizeY;
    Grid(unsigned int sizeX,unsigned int sizeY)
    {
        this->sizeX = sizeX;
        this->sizeY = sizeY;
        mesh = new Cell*[sizeX+1];
        for (int i=0; i<sizeX+1; i++){
            mesh[i] = new Cell[sizeY+1];
        }
        for (int i=0; i<sizeX+1; i++){
            mesh[i][0] = {0,0,0,0};
        }
        for (int i=0; i<sizeY+1; i++){
            mesh[0][i] = {0,0,0,0};
        }
        fluxMesh= nullptr;
    }
    void fluxMeshInit()
    {
        fluxMesh = new Cell*[sizeX+1];
        for (int i=0; i<sizeX+1; i++){
            fluxMesh[i] = new Cell[sizeY+1];
        }
        for (int i=0; i<sizeX+1; i++){
            fluxMesh[i][0] = {0,0,0,0};
        }
        for (int i=0; i<sizeY+1; i++){
            fluxMesh[0][i] = {0,0,0,0};
        }
    }
    void Fill(double v)
    {
        for (int i = 0; i < sizeX + 1; i++){
            for(int j=0; j < sizeY + 1; j++){
                mesh[i][j] = {v,0,0,-20};
                //if (j<50 && i < 50) {mesh[i][j].vx=0;mesh[i][j].vy=30;}
                //if (j>50 && i < 50) {mesh[i][j].vx=30;mesh[i][j].vy=0;}
                // if (j>50 && i > 50) {mesh[i][j].vx=0;mesh[i][j].vy=-30;}
                // if (j<50 && i > 50) {mesh[i][j].vx=-30;mesh[i][j].vy=0;}
            }
        }
    }
    void Print()
    {
        for (int i = 0; i < sizeY + 1; i++){
            for(int j=0; j < sizeX + 1; j++){
                std::cout << mesh[j][i].val<< " ";
            }
            std::cout <<"\n";
        }
    }
    ~Grid() {
        for (int i = 0; i < sizeX + 1; i++) {
            delete[] mesh[i];
        }
        delete[] mesh;
        if (fluxMesh) {
            for (int i = 0; i < sizeX + 1; i++) {
                delete[] fluxMesh[i];
            }
            delete[] fluxMesh;
        }
    }
};

// one variable!!
// ospre flux limiter
Cell SlopeLim(Cell r)
{
    /* return {std::max(0.0,std::max(std::min(2*r.val,1.0),std::min(r.val,2.0))),
             std::max(0.0,std::max(std::min(2*r.vx,1.0),std::min(r.vx,2.0))),
             std::max(0.0,std::max(std::min(2*r.vy,1.0),std::min(r.vy,2.0)))};*/
    return {std::max(0.0, std::min(1.0, r.val)),
            std::max(0.0, std::min(1.0, r.vx)),
            std::max(0.0, std::min(1.0, r.vy)),
            std::max(0.0, std::min(1.0, r.B))};
    /*if (r.val > 0)
        return Cell{ 1.5 * (r.val * r.val + r.val) / (r.val * r.val + r.val + 1) };
    return Cell{ 0 };*/


}

double nonZeroDouble(double val)
{
    if (std::abs(val) < __DBL_EPSILON__)
    {
        if (val < 0)
            return -__DBL_EPSILON__;
        else
            return __DBL_EPSILON__;
    }
    return val;
}

//one vastd::absble !!
Cell nonZeroDenom(Cell denom)
{
    return {nonZeroDouble(denom.val),nonZeroDouble(denom.vx),nonZeroDouble(denom.vy),nonZeroDouble(denom.B)};
}

Cell F(Cell dudx,Cell dudy,Cell U)
{
    // dU/dt = U(dv/dx+dv/dy)
    //dp/dt= -(p * vx/dx + p/dx * vx + p*vy/dy + p/dy * vy)
    //dp/dt= -(p/dx * v + p * v/dx + p/dy * v + p * v/dy)
    //dvx/dt = -(vx*dvx/dx + vy*dvy/dy) - dp/dx/p
    //dvy/dt = -(vx*dvx/dx + vy*dvy/dy) - dp/dy/p

    double dvaldt=   dudx.val * U.vx + U.val * dudx.vx +
                     dudy.val * U.vy + U.val * dudy.vy;
    //vx,vy = pvx,pvy
    //double dvaldt = dudx.vx + dudy.vy;
    //p/dx*vx*vx + 2*p*vx* vx/dx + p/dx + p/dy*vy*vx + p*vy/dy*vx +p*vy/*vx/dy

    // ∂u/∂t = u ∂u/∂x + v ∂u/∂y
    double dpvxdt = (U.vx*dudx.vx + U.vy * dudy.vx) + dudx.val*10 /*+ U.B*(dudx.B + dudy.B)*/;
    double dpvydt = (U.vx*dudx.vy + U.vy * dudy.vy) + dudy.val*10/* + U.B*(dudx.B + dudy.B)*/;
    //dpvxdt=0;
    //dpvydt=0;
    //double dpvxdt = dudx.val * U.vx * U.vx + 2 * U.val * U.vx * dudx.vx + dudx.val*20000 + dudy.val * U.vx*U.vy + U.val *U.vy * dudy.vx + U.val * U.vx * dudy.vy ;
    //double dpvydt = dudy.val * U.vy * U.vy + 2 * U.val * U.vy * dudy.vy + dudy.val*20000 + dudx.val * U.vx*U.vy + U.val *U.vy * dudx.vx + U.val * U.vx * dudx.vy ;

    //double dpvxdt = dudx.val*10;
    //double dpvydt = dudy.val*10;
    // double dpvxdt=dudx.val;


    return Cell{/*(dudx.val*dudx.vx+dudy.val*dudy.vy)*/dvaldt, dpvxdt,dpvydt,0};
}


void CalculateFlux(Grid& out, Grid& in)
{
    //|  E  |  P  |  W  |
    //      e     p
    //in.Print();
    for(int y = 0; y < in.sizeY; y++){
        unsigned int y_1= (y-1+in.sizeY)%(int)in.sizeY;
        unsigned int y1= (y+1+in.sizeY)%(int)in.sizeY;
        unsigned int y_2= (y-2+in.sizeY)%(int)in.sizeY;
        for (int x = 0; x < in.sizeX; x++) // 2->n ??
        {
            //r_i= (u{i} - u{i-1}) / (u{i+1}-u{i})   r_P = (P - E) / (W - P)
            unsigned int x_1= (x-1+in.sizeX)%(int)in.sizeX;
            unsigned int x1= (x+1+in.sizeX)%(int)in.sizeX;
            unsigned int x_2= (x-2+in.sizeX)%(int)in.sizeX;



            auto rx_current = (in.mesh[x][y]     - in.mesh[x_1][y]) / nonZeroDenom( in.mesh[x1][y]- in.mesh[x][y]);
            auto rx_prev =    (in.mesh[x_1][y] - in.mesh[x_2][y]) / nonZeroDenom(in.mesh[x][y] - in.mesh[x_1][y]);




            auto ry_current = (in.mesh[x][y]   - in.mesh[x][y_1]) / nonZeroDenom(in.mesh[x][y1] - in.mesh[x][y]);
            auto ry_prev =    (in.mesh[x][y_1] - in.mesh[x][y_2]) / nonZeroDenom(in.mesh[x][y] - in.mesh[x][y_1]);
            //uL{i-0.5} = u{i-1} + 0.5 * phi(r{i-1}) * (u{i} - u{i - 1})
            //uL{P} = E + 0.5 * phi * (P - E)
            //uR{i-0.5} = u{i} - 0.5 * phi(r{i}) * (u{i+1} - u{i})
            //uR{P} = P - 0.5 * phi * (W - P)
            //phi - slope limiter
            auto uL_x = in.mesh[x_1][y] + 0.5 * SlopeLim(rx_prev) * (in.mesh[x][y] - in.mesh[x_1][y]);
            auto uR_x = in.mesh[x][y]   - 0.5 * SlopeLim(rx_current)   * (in.mesh[x1][y] - in.mesh[x][y]);

            auto uL_y = in.mesh[x][y_1] + 0.5 * SlopeLim(ry_prev) * (in.mesh[x][y] - in.mesh[x][y_1]);
            auto uR_y = in.mesh[x][y]   - 0.5 * SlopeLim(ry_current)   * (in.mesh[x][y1] - in.mesh[x][y]);

            //eigenvalues? -> max speed a
            double a = 50;
            // F{i-0.5} = 0.5 * (F(uR{i-0.5}) + F(uL{i-0.5}) - a * (uR{i-0.5} - uL{i-0.5}))
            out.mesh[x][y] =  (  0.5 * (F(uR_x, {0,0,0,0},in.mesh[x][y]) + F(uL_x, {0,0,0,0},in.mesh[x][y]) - a * (uR_x - uL_x)));
            out.fluxMesh[x][y] = (0.5 * (F({0,0,0,0}, uR_y,in.mesh[x][y]) + F({0,0,0,0}, uL_y,in.mesh[x][y]) - a * (uR_y - uL_y)));

        }

    }

}

Cell S(int x,int y,Cell val)
{
    if (x>3 && x<7 && y>0 && y< 60)
        return {500,300,0,0};
    else return {0,0,0,0};
}

void RKIntegrator(Grid& grid, double dt)
{
    double dx = 1;
    //f=f(grid)
    //k1 = -dt/dx * (f(u){i+0.5} - f(u){i-0.5}
    //k2=f(grid+0*5*k1*dt)
    //grid=grid+k2*dt
    Grid flux(grid.sizeX,grid.sizeY);
    flux.fluxMeshInit();
    CalculateFlux(flux,grid);
    //flux.Print();
    Grid k(grid.sizeX,grid.sizeY);
    for (int y = 0; y < grid.sizeY; y++)
    {
        unsigned int y1= (y+1+grid.sizeY)%(int)grid.sizeY;
        for (int x = 1; x < grid.sizeX; x++)
        {
            unsigned int x1= (x+1+grid.sizeX)%(int)grid.sizeX;
            k.mesh[x][y] = grid.mesh[x][y] - 0.5 * dt / dx * (flux.mesh[x1][y]+ flux.fluxMesh[x][y1] - flux.mesh[x][y] - flux.fluxMesh[x][y]) +0.5*dt * S(x,y,grid.mesh[x][y]);
        }
        //grid.mesh[grid.sizeX-1][y]={0};
    }
    CalculateFlux(flux, k);
    for (int y = 0; y < grid.sizeY; y++)
    {
        unsigned int y1= (y+1+grid.sizeY)%(int)grid.sizeY;
        for (int x = 1; x < grid.sizeX; x++)
        {
            unsigned int x1= (x+1+grid.sizeX)%(int)grid.sizeX;
            grid.mesh[x][y] = grid.mesh[x][y] - dt / dx * (flux.mesh[x1][y]+ flux.fluxMesh[x][y1] - flux.mesh[x][y] - flux.fluxMesh[x][y]) + dt * S(x,y,grid.mesh[x][y]);
        }
        //grid.mesh[grid.sizeX][y]=grid.mesh[grid.sizeX-1][y];
        //grid.mesh[grid.sizeX-1][y]={0};
    }
}

void drawSphere(Grid &grid){
    glPushMatrix();
    glTranslatef(0, 0, 0);
    GLUquadric* gluq = gluNewQuadric();
    glColor3f(255.0, 242.0, 0.0);
    gluSphere(gluq, 4, 10, 10);
    glPopMatrix();

    glPointSize(2);
    glBegin(GL_POINTS);

    int maxLongitude = 180/3; // (1) количество секторов на круге
    int maxLatitude = 60/3; // (2)
    int numOfSectors = 60; // (3)
    Vector3D dataHelper [numOfSectors+1][maxLatitude+1][maxLongitude+1];
    double radius = 6;
    int tmp=0;
    while(tmp <= numOfSectors) {
        glColor3f(255.0, 255.0, 255.0);
        for (int i = 0; i < maxLatitude; i++) {
            for (int j = 0; j < maxLongitude; j++) {
                double phi = ((j)*2*M_PI) / maxLongitude;
                double theta = ((i)*M_PI) / maxLatitude;
                Vector3D coord = Spherical(radius, theta, phi);
                //glVertex3d(coord.x, coord.y, coord.z);
                dataHelper[tmp][i][j] = {coord.x, coord.y, coord.z};
            }
        }
        radius+=4;
        tmp++;
    }

    glEnd();

    glBegin(GL_QUADS);
    glLineWidth(2);
    glColor3f(255.0, 255.0, 255.0);
    for(int i=0; i<numOfSectors; i++){
        //for(int y=0; y<maxLatitude; y++) {
            int y = maxLatitude/2;
            for (int x = 0; x < maxLongitude-1; x++) {
                sf::Color col = toColor(grid.mesh[i][x].val, 0, 50);
                glColor3f(col.r/255.0, col.g/255.0, col.b/255.0);

                glVertex3d(dataHelper[i][y][x].x, dataHelper[i][y][x].y, dataHelper[i][y][x].z);
                glVertex3d(dataHelper[i+1][y][x].x, dataHelper[i+1][y][x].y, dataHelper[i+1][y][x].z);
                glVertex3d(dataHelper[i+1][y][x+1].x, dataHelper[i+1][y][x+1].y, dataHelper[i+1][y][x+1].z);
                glVertex3d(dataHelper[i][y][x + 1].x, dataHelper[i][y][x + 1].y, dataHelper[i][y][x + 1].z);

            }
            glVertex3d(dataHelper[i][y][maxLongitude-1].x, dataHelper[i][y][maxLongitude-1].y, dataHelper[i][y][maxLongitude-1].z);
            glVertex3d(dataHelper[i+1][y][maxLongitude-1].x, dataHelper[i+1][y][maxLongitude-1].y, dataHelper[i+1][y][maxLongitude-1].z);
            glVertex3d(dataHelper[i+1][y][0].x, dataHelper[i+1][y][0].y, dataHelper[i+1][y][0].z);
            glVertex3d(dataHelper[i][y][0].x, dataHelper[i][y][0].y, dataHelper[i][y][0].z);

        //}

    }
    glEnd();

}


void InitGL()
{
    sf::Window window(sf::VideoMode(800, 600), "OpenGL", sf::Style::Default, sf::ContextSettings(32));

    window.setVerticalSyncEnabled(true);

    // activate the window
    window.setActive(true);

    // load resources, initialize the OpenGL states, ...
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // устанавливаем фоновый цвет на черный

    glClearDepth(1.0);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST); // включаем тест глубины
    glShadeModel(GL_SMOOTH);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, (float)800 / (float)600, 0.1f, 1000.0f); // настраиваем трехмерную перспективу
    glMatrixMode(GL_MODELVIEW);

    // draw...
    Grid grid(60,60); //TODO: grid size -> lat and lon
    grid.Fill(1);
    //grid.mesh[5][5].val = 10e-300;
    //grid.mesh[10].val = 50;
    //grid.mesh[11].val = 50;
    // dataArray[20] = 50;
    //dataArray[80] = 50;
    double h = 0.005;

    double ax = 0, ay = 0;
    bool firstTouch;
    sf::Vector2i lastPosition = {0, 0};
    double dist = 500;

    while (window.isOpen())
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
            if(event.type == sf::Event::MouseWheelMoved)
            {
                // display number of ticks mouse wheel has moved
                dist -= event.mouseWheel.delta;
            }
        }

        glPushMatrix();

        cameraControl(&window, &firstTouch, &ax, &ay, &lastPosition, &dist);
        drawSphere(grid);

        //solve(dataArray, arraySize, h);
        for(int i=0;i<1;i++)
            RKIntegrator(grid, h);
        double sum = 0;
        double vel = 0;
        for (int y = 0; y < grid.sizeY; y++)
        {
            for (int x = 0; x < grid.sizeX; x++){
                sum += grid.mesh[x][y].val;
                vel+=std::pow( grid.mesh[x][y].vx,2)+std::pow( grid.mesh[x][y].vy,2);
            }
        }
        // printf("p = %0.10e\n",sum);
        // printf("vel = %0.10e\n",vel);
        // std::cout << dataArray[5] << "\n";
        //grid.Print();

        //show(grid, window,t);

        glPopMatrix();
        // end the current frame (internally swaps the front and back buffers)


        window.display();
    }

}



int main()
{
    InitGL();

    return 0;
}


//int main(){
//    InitGL();
////    ReadData reader("../data/bnd.nc");
////    double *data = nullptr; int dataSize = 0;
////    reader.readData("D", &data, &dataSize);
////    for(int i=0; i<dataSize; i++){
////        std::cout << data[i]*10000 << " ";
////    }
//
//    return 0;
//}
