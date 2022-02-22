#include "main.h"

// dU/dt = a*dU/dx

double f(double i)
{
    double a = 1;
    return a * i;

}

double h(double U1, double U2) {
    return 0.5 * (f(U1) + f(U2));
}

double slopeU(double* arr, int i)
{
    double j_12 = arr[i] - arr[i - 1];
    double j12 = arr[i + 1] - arr[i];
    double j_12sign = signbit(j_12) ? -1 : 1;
    return j_12sign * std::max(0.0, std::min(abs(j_12), j_12sign * j12));

}


void solve(double* arr, int size, double dt)
{
    double dx = 1;
    double* temp = new double[size];
    double* res = new double[size];
    //double tt[100];
    temp[size - 1] = 0;
    res[size-1] = 0;
    res[0] = 0;
    for (int i = 0; i < size - 1; i++)
    {
        double curr = arr[i];



        double v_12 = curr - 0.5 * dt / dx * (
            f(curr + 0.5 * (curr + 0.5 * arr[i + 1])) -
            f(curr + 0.5 * (curr - 0.5 * arr[i + 1]))
            );
        temp[i] = v_12;

    }
   // memcpy(tt, temp, size * sizeof(double));
    for (int i = 1; i < size - 1; i++)
    {
        double h1 = h(temp[i + 1] - 0.5 * slopeU(arr, i + 1), temp[i] + 0.5 * slopeU(arr, i));
        double h2 = h(temp[i] - 0.5 * slopeU(arr, i), temp[i - 1] + 0.5 * slopeU(arr, i - 1));
        double v = arr[i] - dt / dx * (h1 - h2);
        res[i] = v;
    }
    memcpy(arr, res, size * sizeof(double));
    delete temp;
    delete res;
}


void show(double* arr, int size, sf::RenderWindow& window) {
    sf::Vertex line[2];

    for (int i = 0; i < size-1; i++)
    {
        line[0] = sf::Vertex(sf::Vector2f(double(i*1)-500, arr[i]));    
        line[1] = sf::Vertex(sf::Vector2f(double((i+1) * 1-500), arr[i+1]));
        window.draw(line, 2, sf::Lines);
    }
    
}


int main()
{
    sf::RenderWindow window(sf::VideoMode(1000, 600), "bike");
    
    int arraySize = 1000;
    double* dataArray = new double [arraySize];
    
    for(int i = 0; i<arraySize; i++){
        dataArray[i] = 0.0;
    }
    dataArray[20] = 50;
    dataArray[80] = 50;
    double h = 0.1;
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
        solve(dataArray, arraySize, h);
        std::cout << dataArray[5] << "\n";
        show(dataArray, arraySize, window);

        window.display();
    }

    return 0;
}

