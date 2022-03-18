#include "Graphics.h"
#define COLOR_SCHEME 0
sf::Color toColor(double val,double min,double max)
{

    double mid = (min+max)/2;
    if (std::isnan(val)){
        return sf::Color::Magenta;
    }
#if COLOR_SCHEME == 0
    double r = std::min(1.0, 1 - (mid - val) / (max - mid));
    double g = 1 - std::abs(val - mid) / (max - mid);
    double b = std::min(1.0, 1 - (val - mid) / (max - mid));
#else
    double n = (val - min) / (max - min);

        double r=-5*n*n+7.5*n+1-5.0/1.777777777;
        double g=-20*n*n+20*n-4;
        double b=-5*n*n+2.5*n+0.6875;
#endif
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

    sf::RectangleShape r;
    r.setSize(sf::Vector2f{float(segmentX),float(segmentY)});
    for (int y = 0; y < grid.sizeY+1; y++){
        for (int x = 0; x < grid.sizeX+1; x++){
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

    std::stringstream ss;
    ss<<grid.mesh[mposx][mposy].rho;
    t.setString(ss.str());
    t.setPosition(window.mapPixelToCoords( {mPos.x,mPos.y}));
    window.draw(t);
}