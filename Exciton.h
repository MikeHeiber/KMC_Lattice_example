// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

#ifndef EXCITON_H
#define EXCITON_H

#include "KMC_Lattice/Utils.h"
#include "KMC_Lattice/Object.h"
#include "KMC_Lattice/Event.h"
#include <string>

using namespace std;

class Exciton : public Object{
    public:
        static const string name;
        Exciton(const double time,const int tag_num,const Coords& start_coords) : Object(time,tag_num,start_coords){}
        string getName() const{return name;}
    private:
};

class Exciton_Creation : public Event{
    public:
        static const string name;
        string getName() const{return name;}
    private:


};

class Exciton_Hop : public Event{
    using Event::calculateExecutionTime;
    public:
        static const string name;
        void calculateExecutionTime(const double prefactor,const double distance,const double E_delta,const double temp,const double current_time){
            double rate = prefactor*intpow(1/distance,6);
            if(E_delta>0){
                rate *= exp(-E_delta/(K_b*temp));
            }
            calculateExecutionTime(rate,current_time);
        }
        string getName() const{return name;}
    private:

};

class Exciton_Recombination : public Event{
    public:
        static const string name;
        string getName() const{return name;}
    private:
};

#endif // EXCITON_H
