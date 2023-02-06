#pragma once

#include "ATK.hpp"

namespace AudioTK
{
    struct AttackReleaseHysteresis : public ATKFilter
    {
        ATK::AttackReleaseHysteresisFilter<DspFloatType> * filter;

        AttackReleaseHysteresis()
        {
            filter = new ATK::AttackReleaseHysteresisFilter<DspFloatType>(2);
            this->setFilter(filter);
        }
        ~AttackReleaseHysteresis() {
            if(filter) delete filter;
        }
        enum {
            PORT_ATTACK,
            PORT_ATTACKDB,
            PORT_RELEASE,
            PORT_RELEASEDB,           
        };
        void printPorts() {
            printf("PORTS\n");
            printf("PORT_ATTACK=0\n");
            printf("PORT_ATTACKDB=1\n");
            printf("PORT_RELEASE=2\n");
            printf("PORT_RELEASEDB=3\n");
        }
        void setPort(int port, DspFloatType value) {
            switch(port) {
                case PORT_ATTACK: filter->set_attack_hysteresis(value); break;
                case PORT_RELEASE: filter->set_release_hysteresis(value); break;
                case PORT_ATTACKDB: filter->set_attack_hysteresis_db(value); break;
                case PORT_RELEASEDB: filter->set_release_hysteresis_db(value); break;            
            }
        }
    };
}