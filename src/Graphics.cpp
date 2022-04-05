#include "Graphics.h"
#include "Constants.h"

sf::Color toColor(double val,double min,double max)
{

    double mid = (min+max)/2;
    if (std::isnan(val)){
        return sf::Color::Magenta;
    }
    if (val<min){
        return sf::Color::Green;
    }
    if (val>max){
        return sf::Color::Yellow;
    }
#if COLOR_SCHEME == 0
    double r = std::min(1.0, 1 - (mid - val) / (max - mid));
    double g = 1 - std::abs(val - mid) / (max - mid);
    double b = std::min(1.0, 1 - (val - mid) / (max - mid));
#elif COLOR_SCHEME == 1
    double n = (val - min) / (max - min);

    double r=-5*n*n+7.5*n+1-5.0/1.777777777;
    double g=-20*n*n+20*n-4;
    double b=-5*n*n+2.5*n+0.6875;
#else
    double seglen = (max - min)/7;
    int segnum = (int)((val-min) / seglen);
    int segnum1 = std::min(segnum+1,6);
    unsigned char rcomp[7] = {252,254,252,25,22,5,179};
    unsigned char gcomp[7] = {0,173,245,225,166,12,15};
    unsigned char bcomp[7] = {26,41,52,39,252,193,248};
    double r = std::lerp(rcomp[6-segnum],rcomp[6-segnum1],((val-min)-segnum*seglen)/seglen)/255;
    double g = std::lerp(gcomp[6-segnum],gcomp[6-segnum1],((val-min)-segnum*seglen)/seglen)/255;
    double b = std::lerp(bcomp[6-segnum],bcomp[6-segnum1],((val-min)-segnum*seglen)/seglen)/255;
#endif

    sf::Color c((int)std::max(0.0,r*255),(int)std::max(0.0,g*255),(int)std::max(0.0,b*255));
    return c;
}

void show(Grid& grid, sf::RenderWindow& window,sf::Text& t,double upperbound,int mode) {
    //sf::Vertex* points=new sf::Vertex[grid.sizeX*grid.sizeY];
    sf::Vertex varr[grid.sizeX];
    sf::Vertex varry[grid.sizeY];
    sf::Vertex l[2*((grid.sizeX+1)/5)*((grid.sizeY+1)/5+2)];
    double graphm=1.0/2000;
    int graph_h=100;
    int graph_offset=50;

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
            double radius=upperbound;

            Cell U =grid.mesh[x][y];
            double p = gamma * (U.E - 1.0/2*U.rho* (U.vx*U.vx + U.vy*U.vy + U.vz*U.vz)
                                - 1.0/(2*mu)*(U.Bx*U.Bx + U.By*U.By + U.Bz*U.Bz));
            double P = p + (U.Bx * U.Bx + U.By * U.By + U.Bz * U.Bz) / (2 * mu);
            double displayvar=0;
            if (mode<8)
                displayvar=reinterpret_cast<double*>(&U)[mode];
            else displayvar=P;
            int lowerbound=0;
            if (mode >0 && mode < 7) lowerbound = 1;
            r.setFillColor(toColor(displayvar,-radius*lowerbound,radius));
            r.setPosition(sf::Vector2f(  float(x* segmentX)- (float)windowsizeX/2, float(y* segmentY)- (float)windowsizeY/2));
            window.draw(r);
            if (x%5 ==0 && y % 5 ==0)
            {
                auto vec=sf::Vector2f(  float(x* segmentX)- (float)windowsizeX/2, float(y* segmentY)- (float)windowsizeY/2);
                sf::Vector2f dir={(float)grid.mesh[x][y].vx*ARROW_LEN_MULT,(float)grid.mesh[x][y].vy*ARROW_LEN_MULT};
                l[2*(x/5+y/5*grid.sizeX/5)] = sf::Vertex(vec,sf::Color::Black);
                l[2*(x/5+y/5*grid.sizeX/5)+1]=  sf::Vertex(vec+dir,sf::Color::Black);

            }
            if (y==mposy)
                varr[x] = sf::Vertex(sf::Vector2f(float(x* segmentX)/2- (float)windowsizeX/2,(float)windowsizeY/2+graph_offset-graph_h/2-grid.mesh[x][y].E*graphm ),sf::Color::Magenta);
        }

        if (y==mposy)
        {

        }
        varry[y] = sf::Vertex(sf::Vector2f(float(y* segmentX)/2,(float)windowsizeY/2+graph_offset-graph_h/2-grid.mesh[mposx][y].E*graphm ),sf::Color::Magenta);
    }
    window.draw(varr, grid.sizeX, sf::LinesStrip);
    window.draw(l, 2*((grid.sizeX)/5)*((grid.sizeY)/5), sf::Lines);
    window.draw(varry, grid.sizeY, sf::LinesStrip);

    Cell U =grid.mesh[mposx][mposy];
    double p = gamma * (U.E - 1.0/2*U.rho* (U.vx*U.vx + U.vy*U.vy + U.vz*U.vz)
                        - 1.0/(2*mu)*(U.Bx*U.Bx + U.By*U.By + U.Bz*U.Bz));
    double P = p + (U.Bx * U.Bx + U.By * U.By + U.Bz * U.Bz) / (2 * mu);


    double sum=0;
    for (int y = 0; y < grid.sizeY; y++)
    {
        for (int x = 0; x < grid.sizeX; x++){
            if (mode<8)
                sum +=reinterpret_cast<double*>(&grid.mesh[x][y])[mode];
            else sum +=P;
        }
    }


    std::stringstream ss;
    std::stringstream ss2;
    if (mode<8)
        ss<<reinterpret_cast<double*>(&U)[mode];
    else ss<<P;

    t.setString(ss.str());
    t.setPosition(window.mapPixelToCoords( {mPos.x,mPos.y}));
    window.draw(t);

    std::string names[] ={"RHO","Vx","Vy","Vz","Bx","By","Bz","E","P"};
    ss2<<"upperlimit: "<<upperbound << " mode: "<<names[mode]<<" sum: "<<sum;
    t.setString(ss2.str());
    t.setPosition(window.mapPixelToCoords( {0,(int)windowsizeY-20}));
    window.draw(t);
}