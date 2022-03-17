#include <iostream>
#include <SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>
#include <sstream>

// dU/dt = a*dU/dx


struct Cell
{
    double rho;
    double vx;
    double vy;
    double B;

    Cell operator+(const Cell r) const
    {
        return Cell{ rho + r.rho ,vx + r.vx,vy + r.vy,B+r.B};
    }
    Cell operator-(const Cell r) const
    {
        return Cell{ rho - r.rho ,vx - r.vx,vy - r.vy,B-r.B};
    }
    Cell operator*(const Cell r) const
    {
        return Cell{ rho * r.rho ,vx * r.vx,vy * r.vy,B*r.B};
    }
    Cell operator/(const Cell r) const
    {
        return Cell{ rho / r.rho ,vx / r.vx,vy / r.vy,B/r.B};
    }
    Cell operator*(const double r) const
    {
        return Cell{ rho * r ,vx*r,vy*r,B*r};
    }
    Cell operator/(const double r) const
    {
        return Cell{ rho / r ,vx/r,vy/r,B/r};
    }
    friend Cell operator*(const double l, const Cell r)
    {
        return { l * r.rho ,l*r.vx,l*r.vy,l*r.B};
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
        for (int i = 0; i < sizeX ; i++){
            for(int j=0; j < sizeY ; j++){
                mesh[i][j] = {v,0,0,50};
            }
        }
    }
    void Print()
    {
        for (int i = 0; i < sizeY + 1; i++){
            for(int j=0; j < sizeX + 1; j++){
                std::cout << mesh[j][i].rho<< " ";
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
   /* return {std::max(0.0,std::max(std::min(2*r.rho,1.0),std::min(r.rho,2.0))),
            std::max(0.0,std::max(std::min(2*r.vx,1.0),std::min(r.vx,2.0))),
            std::max(0.0,std::max(std::min(2*r.vy,1.0),std::min(r.vy,2.0))),
                std::max(0.0,std::max(std::min(2*r.B,1.0),std::min(r.B,2.0)))};*/
    return {std::max(0.0, std::min(1.0, r.rho)),
            std::max(0.0, std::min(1.0, r.vx)),
            std::max(0.0, std::min(1.0, r.vy)),
             std::max(0.0, std::min(1.0, r.B))};

       /* return Cell{std::max(0.0, 1.5 * (r.rho * r.rho + r.rho) / (r.rho * r.rho + r.rho + 1)),
                    std::max(0.0, 1.5 * (r.vx * r.vx + r.vx) / (r.vx * r.vx + r.vx + 1)),
                             std::max(0.0, 1.5 * (r.vy * r.vy + r.vy) / (r.vy * r.vy + r.vy + 1)),
                             std::max(0.0, 1.5 * (r.B * r.B + r.B) / (r.B * r.B + r.B + 1))};*/



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
    return {nonZeroDouble(denom.rho),nonZeroDouble(denom.vx),nonZeroDouble(denom.vy),nonZeroDouble(denom.B)};
}

Cell F(Cell DX,Cell DY,Cell U)
{
  // double dvaldt=   dudx.rho * U.vx + U.rho * dudx.vx +
                  //  dudy.rho * U.vy + U.rho * dudy.vy;
//dp/dt = p*Vx/dx + p/dx * Vx + p*Vy/dy + p/dy*Vy
    double drhodt = U.rho * DX.vx + DX.rho * U.vx + U.rho * DY.vy + DY.rho * U.vy;
    //p * Vx/dt = p/dx*Vx*Vx + 2*p*Vx/dx*Vx +
    //			p/dy*Vx*Vy + p*Vx/dy*Vy   + p*Vx*Vy/dy
    //			- P/dx
    //			-2*Bx*Bx/dx
    //			-Bx/dy*By - Bx*By/dy
    //			-p/dt * Vx
    double Pdx = 10 * DX.rho/U.rho;
    double dvxdt = DX.rho*U.vx*U.vx + 2*U.rho*DX.vx*U.vx +
                   DY.rho*U.vx*U.vy + U.rho * DY.vx * U.vy + U.rho * U.vx * DY.vy
                    +Pdx
                    -drhodt * U.vx;
                     - (U.vx*U.B)+ (U.vy*U.B);
    //double pdvxdt = U.vx * DX.vx + U.vy * DY.vx + Pdx;
    //p * Vy/dt = p/dx*Vx*Vy + p*Vx/dx*Vy   + p*Vx*Vy/dx
    //			p/dy*Vy*Vy + 2*p*Vy/dy*Vy +
    //			-P/dy
    //			-Bx/dx*By - Bx*By/dx
    //			-2*By*By/dy
    //			-p/dt * Vy
    double Pdy = 10 * DY.rho/U.rho;
    double dvydt = DX.rho*U.vx*U.vy + U.rho * DX.vx * U.vy + U.rho * U.vx * DX.vy +
                    DY.rho*U.vy*U.vy + 2*U.rho*DY.vy*U.vy +
                    +Pdy
                    -drhodt * U.vy;
                     + (U.vy*U.B)-(U.vx*U.B);
    //double pdvydt =U.rho *  (U.vx * DX.vy + U.vy * DY.vy) + Pdy;
    //pdvxdt = (U.vx*DX.vx + U.vy * DY.vx) + DX.rho*10 /*+ U.B*(dudx.B + dudy.B)*/;
    //pdvydt = (U.vx*DX.vy + U.vy * DY.vy) + DY.rho*10/* + U.B*(dudx.B + dudy.B)*/;
   return Cell{drhodt, dvxdt ,dvydt ,0};
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
            double a = 300;
            // F{i-0.5} = 0.5 * (F(uR{i-0.5}) + F(uL{i-0.5}) - a * (uR{i-0.5} - uL{i-0.5}))
            out.mesh[x][y] =     (0.5 * (F(uR_x, {0,0,0,0},in.mesh[x][y]) + F(uL_x, {0,0,0,0},in.mesh[x][y]) - a * (uR_x - uL_x)));
            out.fluxMesh[x][y] = (0.5 * (F({0,0,0,0}, uR_y,in.mesh[x][y]) + F({0,0,0,0}, uL_y,in.mesh[x][y]) - a * (uR_y - uL_y)));

        }

    }

}

Cell S(int x,int y,Cell val)
{
    if (x>3 && x<7 && y>32 && y< 58)
            return {0,0,0,0};
    else return {0,0,0,0};
}

void ApplyBoundaryConditions(Grid& grid)
{
    for (int x=3;x<7;x++)
    {
        for (int y=32; y< 58;y++)
        {
            grid.mesh[x][y].rho=10;
            grid.mesh[x][y].vx=10;
            grid.mesh[x][y].vy=0;
        }
    }
    for (int y=0;y<grid.sizeY;y++) {
        grid.mesh[0][y].rho=0.01;
        grid.mesh[0][y].vx=0;
        grid.mesh[0][y].vy=0;
        grid.mesh[grid.sizeX-1][y].rho=0.01;
        grid.mesh[grid.sizeX-1][y].vx=0;
        grid.mesh[grid.sizeX-1][y].vy=0;
    }
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
        for (int x = 0; x < grid.sizeX; x++)
        {
            unsigned int x1= (x+1+grid.sizeX)%(int)grid.sizeX;
            k.mesh[x][y] = grid.mesh[x][y] - 0.5 * dt / dx * (flux.mesh[x1][y]+ flux.fluxMesh[x][y1] - flux.mesh[x][y] - flux.fluxMesh[x][y]) +0.5*dt * S(x,y,grid.mesh[x][y]);
        }
    }
    ApplyBoundaryConditions(k);
    CalculateFlux(flux, k);
    for (int y = 0; y < grid.sizeY; y++)
    {
        unsigned int y1= (y+1+grid.sizeY)%(int)grid.sizeY;
        for (int x = 0; x < grid.sizeX; x++)
        {
            unsigned int x1= (x+1+grid.sizeX)%(int)grid.sizeX;
            grid.mesh[x][y] = grid.mesh[x][y] - dt / dx * (flux.mesh[x1][y]+ flux.fluxMesh[x][y1] - flux.mesh[x][y] - flux.fluxMesh[x][y]) + dt * S(x,y,grid.mesh[x][y]);
        }
    }
    ApplyBoundaryConditions(grid);

}

sf::Color toColor(double val,double min,double max)
{

    double mid = (min+max)/2;
    if (std::isnan(val)){
        return sf::Color::Magenta;
    }
    double r=std::min(1.0,1-(mid-val)/(max-mid));
   double g=1-std::abs(val-mid)/(max-mid);
    double b=std::min(1.0,1-(val-mid)/(max-mid));
    double n=(val-min)/(max-min);
   // double r=-5*n*n+7.5*n+1-5.0/1.777777777;
   //double g=-20*n*n+20*n-4;
  //double b=-5*n*n+2.5*n+0.6875;
    sf::Color c((int)std::max(0.0,r*255),(int)std::max(0.0,g*255),(int)std::max(0.0,b*255));
    return c;
}

void show(Grid& grid, sf::RenderWindow& window,sf::Text& t) {
    //sf::Vertex* points=new sf::Vertex[grid.sizeX*grid.sizeY];
    sf::Vertex varr[grid.sizeX];
    sf::Vertex varry[grid.sizeY];
    sf::Vertex l[2*((grid.sizeX+1)/5)*((grid.sizeY+1)/5+2)];
    double graphm=1.0/3;
    int graph_h=100;

    auto mPos=sf::Mouse::getPosition()- window.getPosition();

    unsigned int windowsizeX = window.getSize().x;
    unsigned int windowsizeY = window.getSize().y;
    int segmentX = windowsizeX / grid.sizeX;
    int segmentY = (windowsizeY-graph_h) / grid.sizeY;
    int mposx=std::clamp(mPos.x/segmentX,0,(int)grid.sizeX);
    int mposy=std::clamp(mPos.y/segmentY,0,(int)grid.sizeY);

    //std::cout<<grid.mesh[mposx][mposy].val<<"\n";
    sf::RectangleShape r;
    r.setSize(sf::Vector2f{float(segmentX),float(segmentY)});
    for (int y = 0; y < grid.sizeY+1; y++){
        for (int x = 0; x < grid.sizeX+1; x++){
          /*  points[x+y*grid.sizeX] = sf::Vertex(sf::Vector2f(  float(x* segmentX)- (float)windowsizeX ,
                                                    float(y* segmentY)- (float)windowsizeY),
                                     toColor(grid.mesh[x][y].val,-300,300));*/
          double radius=5;
         // if (std::sqrt(std::pow((x-mposx),2)+std::pow((y-mposy),2))<30) radius=10e-324;
         auto displayvar=grid.mesh[x][y].rho;
          r.setFillColor(toColor(displayvar,-radius*0,radius));
          r.setPosition(sf::Vector2f(  float(x* segmentX)- (float)windowsizeX/2, float(y* segmentY)- (float)windowsizeY/2));
          window.draw(r);
          if (x%5 ==0 && y % 5 ==0)
          {
              auto vec=sf::Vector2f(  float(x* segmentX)- (float)windowsizeX/2, float(y* segmentY)- (float)windowsizeY/2);
              sf::Vector2f dir={(float)grid.mesh[x][y].vx*0.9f,(float)grid.mesh[x][y].vy*0.9f};
              l[2*(x/5+y/5*grid.sizeX/5)] = sf::Vertex(vec,sf::Color::Black);
              l[2*(x/5+y/5*grid.sizeX/5)+1]=  sf::Vertex(vec+dir,sf::Color::Black);

          }
          if (y==mposy)
             varr[x] = sf::Vertex(sf::Vector2f(float(x* segmentX)/2- (float)windowsizeX/2,(float)windowsizeY/2-graph_h/2-grid.mesh[x][y].rho*graphm ),sf::Color::Green);
        }

        if (y==mposy)
        {
            window.draw(varr, grid.sizeX, sf::LinesStrip);
        }
        varry[y] = sf::Vertex(sf::Vector2f(float(y* segmentX)/2,(float)windowsizeY/2-graph_h/2-grid.mesh[mposx][y].rho*graphm ),sf::Color::Green);
    }
    window.draw(l, 2*((grid.sizeX)/5)*((grid.sizeY)/5), sf::Lines);
    window.draw(varry, grid.sizeY, sf::LinesStrip);

    //auto v=grid.mesh[mposx][mposy].val;
    //std::string varAsString = std::to_string(v);
    //std::cout << varAsString << std::endl;

    //sf::Text* t = new sf::Text();
    //t->setFont(f);
    //t->setCharacterSize(12);
    //t->setFillColor(sf::Color::Black);
    std::stringstream ss;
    ss<<grid.mesh[mposx][mposy].rho;
    //ss<<varAsString;
    t.setString(ss.str());

    //t->setString("a34434434434344334\0\0\0\0\0\0\n");

    t.setPosition(window.mapPixelToCoords( {mPos.x,mPos.y}));
    window.draw(t);


   /* for (int i = 0; i < grid.sizeX; i++)
    {

    }*/

   // window.draw(points, grid.sizeX*grid.sizeY, sf::Points);
    //delete[] points;
    //delete t;
}


int main()
{
    sf::RenderWindow window(sf::VideoMode(700, 700), "wave", sf::Style::Default, sf::ContextSettings(32));
    window.setActive(true);
    window.setFramerateLimit(60);
    //InitGL();
    Grid grid(70,70);
    grid.Fill(0.01);
    for (int y = 90; y < 99; y++){
        for (int x = 30; x < 60; x++){
           //grid.mesh[x][y].rho = 50;
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

