#include <iostream>
#include <SFML/Graphics.hpp>
#include <cmath>

// dU/dt = a*dU/dx

struct Cell
{
    double val;

    Cell operator+(const Cell r) const
    {
        return Cell{ val + r.val };
    }
    Cell operator-(const Cell r) const
    {
        return Cell{ val - r.val };
    }
    Cell operator*(const Cell r) const
    {
        return Cell{ val * r.val };
    }
    Cell operator/(const Cell r) const
    {
        return Cell{ val / r.val };
    }
    Cell operator*(const double r) const
    {
        return Cell{ val * r };
    }
    Cell operator/(const double r) const
    {
        return Cell{ val / r };
    }
    friend Cell operator*(const double l, const Cell r)
    {
        return { l * r.val };
    }
};


struct Grid
{
    Cell** mesh;
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
            mesh[i][0] = {0};
        }
        for (int i=0; i<sizeY+1; i++){
            mesh[0][i] = {0};
        }
    }
    void Fill(double v)
    {
        for (int i = 0; i < sizeX + 1; i++){
            for(int j=0; j < sizeY + 1; j++){
                mesh[i][j] = {v};
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
    ~Grid()
    {
        for (int i=0; i<sizeX+1; i++){
            delete[] mesh[i];
        }
        delete[] mesh;
    }
};

// one variable!!
// ospre flux limiter
Cell SlopeLim(Cell r)
{
    if (r.val > 0)
        return Cell{ 1.5 * (r.val * r.val + r.val) / (r.val * r.val + r.val + 1) };
    return Cell{ 0 };
}

Cell F(Cell dudx,Cell dudy)
{
    // dU/dt = a*(dU/dx+dU/dy)
    double a = -10;
    return (dudy) * a;
}

//one vastd::absble !!
Cell nonZeroDenom(Cell denom)
{
    if (std::abs(denom.val) < __DBL_EPSILON__)
    {
        if (denom.val < 0)
            return {-__DBL_EPSILON__};
        else
            return {__DBL_EPSILON__};
    }
    return denom;
}

void CalculateFlux(Grid& out, Grid& in)
{
    //|  E  |  P  |  W  |
    //      e     p
    //in.Print();
    for(int y = 2; y < in.sizeY; y++){
        for (int x = 2; x < in.sizeX; x++) // 2->n ??
        {
            //r_i= (u{i} - u{i-1}) / (u{i+1}-u{i})   r_P = (P - E) / (W - P)
            auto rx_current = (in.mesh[x][y]     - in.mesh[x - 1][y]) / nonZeroDenom(in.mesh[x + 1][y] - in.mesh[x][y]);
            auto rx_prev =    (in.mesh[x - 1][y] - in.mesh[x - 2][y]) / nonZeroDenom(in.mesh[x][y] - in.mesh[x - 1][y]);

            auto ry_current = (in.mesh[x][y]     - in.mesh[x][y - 1]) / nonZeroDenom(in.mesh[x][y + 1] - in.mesh[x][y]);
            auto ry_prev =    (in.mesh[x][y - 1] - in.mesh[x][y - 2]) / nonZeroDenom(in.mesh[x][y] - in.mesh[x][y - 1]);
            //uL{i-0.5} = u{i-1} + 0.5 * phi(r{i-1}) * (u{i} - u{i - 1})
            //uL{P} = E + 0.5 * phi * (P - E)
            //uR{i-0.5} = u{i} - 0.5 * phi(r{i}) * (u{i+1} - u{i})
            //uR{P} = P - 0.5 * phi * (W - P)
            //phi - slope limiter
            auto uL_x = in.mesh[x - 1][y] + 0.5 * SlopeLim(rx_prev) * (in.mesh[x][y] - in.mesh[x - 1][y]);
            auto uR_x = in.mesh[x][y]     - 0.5 * SlopeLim(rx_current)   * (in.mesh[x + 1][y] - in.mesh[x][y]);

            auto uL_y = in.mesh[x][y - 1] + 0.5 * SlopeLim(ry_prev) * (in.mesh[x][y] - in.mesh[x][y - 1]);
            auto uR_y = in.mesh[x][y]     - 0.5 * SlopeLim(ry_current)   * (in.mesh[x][y + 1] - in.mesh[x][y]);

            //eigenvalues? -> max speed a
            double a = 10;
            // F{i-0.5} = 0.5 * (F(uR{i-0.5}) + F(uL{i-0.5}) - a * (uR{i-0.5} - uL{i-0.5}))
            out.mesh[x][y] = 0.5 * (F(uR_x, uR_y) + F(uL_x, uL_y) - a * ((uR_x + uR_y) - (uL_x + uL_y)));
        }
        out.mesh[0][y] = { 0 };
        out.mesh[1][y] = { 0 };
        for (int x = 0; x < 2; x++) {
            out.mesh[x][0] = { 0 };
            out.mesh[x][1] = { 0 };
        }
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
    CalculateFlux(flux,grid);
    //flux.Print();
    Grid k(grid.sizeX,grid.sizeY);
    for (int y = 1; y < grid.sizeY; y++)
    {
        for (int x = 1; x < grid.sizeX; x++)
        {
            k.mesh[x][y] = grid.mesh[x][y] - 0.5 * dt / dx * (flux.mesh[x + 1][y] + flux.mesh[x][y + 1] - flux.mesh[x][y] * 2);
        }
    }
    CalculateFlux(flux, k);
    for (int y = 1; y < grid.sizeY; y++)
    {
        for (int x = 1; x < grid.sizeX; x++)
        {
            grid.mesh[x][y] = grid.mesh[x][y] - dt / dx * (flux.mesh[x + 1][y] + flux.mesh[x][y + 1] - flux.mesh[x][y] * 2);
        }
    }
}

sf::Color toColor(double val,double min,double max)
{

    double mid = (min+max)/2;
    if (std::isnan(val)){
        return sf::Color::Black;
    }
    double r=std::min(1.0,1-(mid-val)/(max-mid));
    double g=1-std::abs(val-mid)/(max-mid);
    double b=std::min(1.0,1-(val-mid)/(max-mid));
    double n=(val-min)/(max-min);
    //double r=-5*n*n+7.5*n+1-5.0/1.777777777;
    //double g=-20*n*n+20*n-4;
    //double b=-5*n*n+2.5*n+0.6875;
    sf::Color c((int)std::max(0.0,r*255),(int)std::max(0.0,g*255),(int)std::max(0.0,b*255));
    return c;
}

void show(Grid& grid, sf::RenderWindow& window) {
    //sf::Vertex* points=new sf::Vertex[grid.sizeX*grid.sizeY];

    unsigned int windowsizeX = window.getSize().x/2;
    unsigned int windowsizeY = window.getSize().y/2;
    unsigned int segmentX = window.getSize().x / grid.sizeX;
    unsigned int segmentY = window.getSize().y / grid.sizeY;
    sf::RectangleShape r;
    r.setSize(sf::Vector2f{float(segmentX),float(segmentY)});
    for (int y = 0; y < grid.sizeY; y++){
        for (int x = 0; x < grid.sizeX; x++){
          /*  points[x+y*grid.sizeX] = sf::Vertex(sf::Vector2f(  float(x* segmentX)- (float)windowsizeX ,
                                                    float(y* segmentY)- (float)windowsizeY),
                                     toColor(grid.mesh[x][y].val,-300,300));*/
          r.setFillColor(toColor(grid.mesh[x][y].val,-300,300));
          r.setPosition(sf::Vector2f(  float(x* segmentX)- (float)windowsizeX, float(y* segmentY)- (float)windowsizeY));
          window.draw(r);
        }
    }

   // window.draw(points, grid.sizeX*grid.sizeY, sf::Points);
    //delete[] points;
}


int main()
{
    sf::RenderWindow window(sf::VideoMode(1000, 600), "wave");

    Grid grid(100,100);
    grid.Fill(0);
    for (int y = 10; y < 30; y++){
        for (int x = 10; x < 30; x++){
            grid.mesh[x][y].val = 100;
        }
    }
    //grid.mesh[10].val = 50;
    //grid.mesh[11].val = 50;
    // dataArray[20] = 50;
    //dataArray[80] = 50;
    double h = 0.01;
    sf::View w;
    w = window.getDefaultView();
    w.setCenter(0, 0);
    window.setView(w);
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
        RKIntegrator(grid, h);
        double sum = 0;
        for (int y = 0; y < grid.sizeY; y++)
        {
            for (int x = 0; x < grid.sizeX; x++){
                //sum += abs(grid.mesh[x+1].val - grid.mesh[i].val);
            }
        }
        // printf("Total Variation = %0.10e\n",sum);
        // std::cout << dataArray[5] << "\n";
        //grid.Print();
        show(grid, window);

        window.display();
    }

    return 0;
}

