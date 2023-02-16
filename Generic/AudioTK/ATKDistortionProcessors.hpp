#pragma once

#include "ATK.hpp"

namespace AudioTK
{

    struct DiodeClipper : public ATKMonoFilter
    {
        ATK::DiodeClipperFilter<DspFloatType> * filter;

        DiodeClipper() : ATKMonoFilter()
        {
            filter = new ATK::DiodeClipperFilter<DspFloatType>;
            this->setFilter(filter);
        }
        ~DiodeClipper() {
            if(filter) delete filter;
        }
    };

    struct MonoHalfTanhShaper : public ATKMonoFilter
    {
        ATK::HalfTanhShaperFilter<DspFloatType> * filter;

        MonoHalfTanhShaper() : ATKMonoFilter()
        {
            filter = new ATK::HalfTanhShaperFilter<DspFloatType>(1);
            this->setFilter(filter);
        }
        ~MonoHalfTanhShaper() {
            if(filter) delete filter;
        }
        enum {
            PORT_COEFF,
        };

        void setPort(int port, DspFloatType value)
        {
            switch(port)
            {
                case PORT_COEFF: filter->set_coefficient(value);
            }
        }
    };

    struct StereoHalfTanhShaper : public ATKFilter
    {
        ATK::HalfTanhShaperFilter<DspFloatType> * filter;

        StereoHalfTanhShaper()
        {
            filter = new ATK::HalfTanhShaperFilter<DspFloatType>(2);
            this->setFilter(filter);
        }
        ~StereoHalfTanhShaper() {
            if(filter) delete filter;
        }
        enum {
            PORT_COEFF,
        };

        void setPort(int port, DspFloatType value)
        {
            switch(port)
            {
                case PORT_COEFF: filter->set_coefficient(value);
            }
        }
    };
    struct SD1Overdrive : public ATKMonoFilter
    {
        ATK::SD1OverdriveFilter<DspFloatType> * filter;

        SD1Overdrive()
        {
            filter = new ATK::SD1OverdriveFilter<DspFloatType>;
            this->setFilter(filter);
        }
        ~SD1Overdrive() {
            if(filter) delete filter;
        }
        enum
        {
            PORT_DRIVE
        };
        void setPort(int port, DspFloatType value) {
            switch(port) {
                case PORT_DRIVE: filter->set_drive(value);
            }
        }
    };
    struct TS9Overdrive : public ATKMonoFilter
    {
        ATK::TS9OverdriveFilter<DspFloatType> * filter;
        TS9Overdrive() : ATKMonoFilter()
        {
            filter = new ATK::TS9OverdriveFilter<DspFloatType>;
            this->setFilter(filter);
        }
        ~TS9Overdrive() {
            if(filter) delete filter;
        }
        enum
        {
            PORT_DRIVE
        };
        void setPort(int port, DspFloatType value) {
            switch(port) {
                case PORT_DRIVE: filter->set_drive(value);
            }
        }
    };
    struct MonoTanhShaper : public ATKMonoFilter
    {
        ATK::TanhShaperFilter<DspFloatType> * filter;
        MonoTanhShaper()
        {
            filter = new ATK::TanhShaperFilter<DspFloatType>(1);
            this->setFilter(filter);
        }
        ~MonoTanhShaper() {
            if(filter) delete filter;
        }
        enum {
            PORT_COEFF,
        };

        void setPort(int port, DspFloatType value)
        {
            switch(port)
            {
                case PORT_COEFF: filter->set_coefficient(value);
            }
        }
    };
    
    struct StereoTanhShaper : public ATKFilter
    {
        ATK::TanhShaperFilter<DspFloatType> * filter;
        StereoTanhShaper()
        {
            filter = new ATK::TanhShaperFilter<DspFloatType>(2);
            this->setFilter(filter);
        }
        ~StereoTanhShaper() {
            if(filter) delete filter;
        }
        enum {
            PORT_COEFF,
        };

        void setPort(int port, DspFloatType value)
        {
            switch(port)
            {
                case PORT_COEFF: filter->set_coefficient(value);
            }
        }
    };
}