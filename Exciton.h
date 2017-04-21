#ifndef EXCITON_H
#define EXCITON_H

#include "Utils.h"
#include "Object.h"
#include "Event.h"
#include <string>

using namespace std;

class Exciton : public Object{
    public:
        static const string name;
        static double lifetime;
        static double R_hop;
        static int FRET_cutoff;
        Exciton(const double time,const int tag_num,const Coords& start_coords) : Object(time,tag_num,start_coords){}
        string getName(){return name;}
        double getLifetime(){return lifetime;}
    private:

};

class Exciton_Creation : public Event{
    public:
        void calculateEvent(const Coords& dest_coords,const double distance,const double E_delta,const int temperature, const double prefactor){
            setDestCoords(dest_coords);
            // No target object
            setWaitTime(-1*prefactor*log(rand01()));
        }
        string getName(){return name;}
    private:
        static const string name;

};

class Exciton_Hop : public Event{
    public:
        void calculateEvent(const Coords& dest_coords,const double distance,const double E_delta,const int temperature, const double prefactor){
            setDestCoords(dest_coords);
            // No target object
            setWaitTime((-1/(Exciton::R_hop*intpow(1/distance,6)*exp(-E_delta/(K_b*temperature))))*log(rand01()));
        }
        string getName(){return name;}
    private:
        static const string name;
};

class Exciton_Recombine : public Event{
    public:
        void calculateEvent(const Coords& dest_coords,const double distance,const double E_delta,const int temperature, const double prefactor){
            setDestCoords(dest_coords);
            // No target object
            setWaitTime(-1*Exciton::lifetime*log(rand01()));
        }
        string getName(){return name;}
    private:
        static const string name;
};

#endif // EXCITON_H
