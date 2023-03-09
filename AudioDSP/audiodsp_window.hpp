#pragma once

namespace AudioDSP
{
    template<typename T>
    struct Window : public sample_vector<T>
    {   
        using sample_vector<T>::operator [];
        using sample_vector<T>::size;
        using sample_vector<T>::data;
        using sample_vector<T>::resize;

        enum WindowType {
            RECTANGLE,
            HANNING,
            HAMMING,
            BLACKMAN,
            BLACKMANHARRIS,
            GAUSSIAN,
            WELCH,
            PARZEN,
            TUKEY,
        } type = HANNING;

        Window() = default;
        Window(WindowType type, size_t i) { resize(i); make_window(i,type); }
        virtual ~Window() = default;

        void make_window(size_t n, WindowType type)
        {
            this->type = type;
            this->resize(n);
            switch(type)
            {
                case RECTANGLE: rectangle(n); break;
                case HANNING: hanning(n); break;
                case HAMMING: hamming(n); break;
                case BLACKMAN: blackman(n); break;
                case BLACKMANHARRIS: blackmanharris(n); break;
                case GAUSSIAN: gaussian(n); break;
                case WELCH: welch(n); break;
                case PARZEN: parzen(n); break;
                case TUKEY: throw std::runtime_error("Can't init tukey window with make_window");
            }            
        }
        void rectangle(size_t n) {            
            std::fill(this->begin(),this->end(),(T)1.0);
        }
        void hamming(size_t n) {            
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)
            {
                (*this)[i] = 0.54 - (0.46 * std::cos(2*M_PI*(double)i/(double)n));
            }        
        }
        void hanning(size_t n)
        {            
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)
            {
                (*this)[i] = 0.5*(1 - std::cos(2*M_PI*(double)i/(double)n));
            }        
        }
        void blackman(size_t n)
        {                    
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)                    
                (*this)[i] = 0.42 - (0.5* std::cos(2*M_PI*i/(n)) + (0.08*std::cos(4*M_PI*i/n)));        
        }
        void blackmanharris(size_t n)
        {            
            #pragma omp simd        
            for(size_t i = 0; i < this->size(); i++)            
            {   
                double ci = (double) i / (double) n;
                (*this)[i] = 0.35875 
                        - 0.48829*std::cos(2*M_PI*(ci))
                        + 0.14128*std::cos(4.0*M_PI*(ci)) 
                        - 0.01168*std::cos(6.0*M_PI*(ci));
            }
        }
        void gaussian(size_t n)
        {            
            T a,b,c=0.5;
            #pragma omp simd        
            for(size_t i = 0; i < this->size(); i++)
            {
                a = ((double)i - c*(this->size()-1)/(std::sqrt(c)*this->size()-1));
                b = -c * std::sqrt(a);
                (*this)[i] = std::exp(b);
            }    
        }
        void welch(size_t n)
        {            
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)
                (*this)[i] = 1.0 - std::sqrt((2.0*(double)i-(double)this->size()-1)/((double)this->size()));        
        }
        void parzen(size_t n)
        {            
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)
                (*this)[i] = 1.0 - std::abs((2.0*(double)i-this->size()-1)/(this->size()));        
        }
        void tukey(size_t num_samples, T alpha)
        {            
            T value = (-1*(num_samples/2)) + 1;
            double n2 = (double)num_samples / 2.0;
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)
            {    
                if(value >= 0 && value <= (alpha * (n2))) 
                    (*this)[i] = 1.0; 
                else if(value <= 0 && (value >= (-1*alpha*(n2)))) 
                    (*this)[i] = 1.0;
                else 
                    (*this)[i] = 0.5 * (1 + std::cos(M_PI *(((2.0*value)/(alpha*(double)num_samples))-1)))        ;
                value = value + 1;
            }     
        }

        #ifdef SWIG
        %extend {
            sample_vector<T> __mul__(const sample_vector<T> & v)
            {
                assert($self->size() == v.size());
                sample_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++)
                    r[i] = (*$self)[i] * v[i];
                return r;
            }
        }
        #endif
    };
}