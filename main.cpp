#include <iostream>
#include <SFML/Graphics.hpp>

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
    Cell* mesh;
    unsigned int size;
    Grid(unsigned int size)
    {
        this->size = size;
        mesh = new Cell[size+1];  
        mesh[0] = { 0 };
        mesh[size] = { 0 };
    }
    void Fill(double v)
    {
        for (int i = 0; i < size + 1; i++)
            mesh[i].val = v;
    }
    void Print()
    {
        for (int i = 0; i < size + 1; i++)
            std::cout << mesh[i].val<< " ";
        std::cout <<"\n";
    }
    ~Grid()
    {
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

Cell F(Cell u)
{
    double a = 50;
    return u * a;
}

//one variable !!
Cell nonZeroDenom(Cell denom)
{
    if (abs(denom.val) < __DBL_EPSILON__)
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
    for (int i = 2; i < in.size; i++) // 2->n ??
    {   
        //r_i= (u{i} - u{i-1}) / (u{i+1}-u{i})   r_P = (P - E) / (W - P)
        auto ri =   (in.mesh[i]     - in.mesh[i - 1]) / nonZeroDenom(in.mesh[i + 1] - in.mesh[i]);
        auto ri_1 = (in.mesh[i - 1] - in.mesh[i - 2]) / nonZeroDenom(in.mesh[i] - in.mesh[i - 1]);

        //uL{i-0.5} = u{i-1} + 0.5 * phi(r{i-1}) * (u{i} - u{i - 1})
        //uL{P} = E + 0.5 * phi * (P - E)
        //uR{i-0.5} = u{i} - 0.5 * phi(r{i}) * (u{i+1} - u{i})
        //uR{P} = P - 0.5 * phi * (W - P)
        //phi - slope limiter
        auto uL = in.mesh[i - 1] + 0.5 * SlopeLim(ri_1) * (in.mesh[i] - in.mesh[i - 1]);
        auto uR = in.mesh[i]     - 0.5 * SlopeLim(ri)   * (in.mesh[i + 1] - in.mesh[i]);

        //eigenvalues? -> max speed a
        double a = 50;
        // F{i-0.5} = 0.5 * (F(uR{i-0.5}) + F(uL{i-0.5}) - a * (uR{i-0.5} - uL{i-0.5}))
        out.mesh[i] = 0.5 * (F(uR) + F(uL) - a * (uR - uL));
    }
    out.mesh[0] = { 0 };
    out.mesh[1] = { 0 };
}


void RKIntegrator(Grid& grid, double dt)
{
    double dx = 1;
    //f=f(grid)
    //k1 = -dt/dx * (f(u){i+0.5} - f(u){i-0.5}
    //k2=f(grid+0*5*k1*dt)
    //grid=grid+k2*dt
    Grid flux(grid.size);
    CalculateFlux(flux,grid);
    Grid k(grid.size);
    for (int i = 1; i < grid.size; i++)
    {
        k.mesh[i] = grid.mesh[i] - 0.5 * dt / dx * (flux.mesh[i + 1] - flux.mesh[i]);
    }
    CalculateFlux(flux, k);
    for (int i = 1; i < grid.size; i++)
    {
        grid.mesh[i] = grid.mesh[i] - dt / dx * (flux.mesh[i + 1] - flux.mesh[i]);
    }
}


void show(Grid& grid, sf::RenderWindow& window) {
    sf::Vertex* line=new sf::Vertex[grid.size];

    int windowsize = window.getSize().x/2;
    int segment = window.getSize().x / grid.size;
    for (int i = 0; i < grid.size; i++)
    {
        line[i] = sf::Vertex(sf::Vector2f(double(i* segment)- windowsize, -grid.mesh[i].val));
       
    }
    window.draw(line, grid.size, sf::LineStrip);
    delete[] line;
}


int main()
{
    sf::RenderWindow window(sf::VideoMode(1000, 600), "wave");
    
    Grid grid(1000);
    grid.Fill(0);
    for (int i = 10; i < 300; i++) grid.mesh[i].val = i % 50 > 25 ? 200 : -200;
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
        for (int i = 0; i < grid.size; i++)
        {
            sum += abs(grid.mesh[i+1].val - grid.mesh[i].val);
        }
        printf("Total Variation = %0.10e\n",sum);
       // std::cout << dataArray[5] << "\n";
        show(grid, window);

        window.display();
    }

    return 0;
}

