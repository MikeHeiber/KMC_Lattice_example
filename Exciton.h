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
        static const double lifetime;
        static const double R_hop;
        static const int FRET_cutoff;
        Exciton();
        string getName(){return name;}
        double getLifetime(){return lifetime;}
    private:

};

class Exciton_Creation : public Event{
    public:
        void calculateEvent(const double generation_rate){

        }
        bool executeEvent(){

        }
        string getName(){return name;}
    private:
        static const string name;

};

class Exciton_Hop : public Event{
    public:
        void calculateEvent(const Coords& dest_coords,const double distance,const double E_delta,const int temperature){
            list<unique_ptr<Object>>::iterator it = getObjectIt();
            setDestCoords(dest_coords);
            setWaitTime((-1/(Exciton::R_hop*intpow(1/distance,6)*exp(-E_delta/(K_b*temperature))))*log(rand01()));
        }
        bool executeEvent(){

        }
        string getName(){return name;}
    private:
        static const string name;
};

class Exciton_Recombine : public Event{
    public:
        void calculateEvent(){
            list<unique_ptr<Object>>::iterator it = getObjectIt();
            setDestCoords((*it)->getCoords());
            // No target object
            setWaitTime(-1*Exciton::lifetime*log(rand01()));
        }
        bool executeEvent(){

        }
        string getName(){return name;}
    private:
        static const string name;

};

// Initialize names
const string Exciton::name = "Exciton";
const string Exciton_Creation::name = "Exciton Creation";
const string Exciton_Hop::name = "Exciton Hop";
const string Exciton_Recombine::name = "Exciton Recombine";

#endif // EXCITON_H
