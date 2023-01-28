#pragma once

#include <valarray>
#include <numeric>
#include <algorithm>
#include <iostream>


namespace qmdsp
{

    class ChromaVector : public std::valarray<double>
    {
    public:
        ChromaVector(size_t uSize = 12) : std::valarray<double>() {
            resize(uSize, 0.0f);
        }
            
        virtual ~ChromaVector() {};
            
        void printDebug() {
            for (int i = 0; i < int(size()); i++) {
                std::cout <<  (*this)[i] << ";";
            }
            std::cout << std::endl;
        }
            
        void normalizeL1() {
            // normalize the chroma vector (L1 norm)
            double dSum = 0.0;
            
            for (size_t i = 0; i < 12; (dSum += std::abs((*this)[i++]))) ;
            for (size_t i = 0; i < 12; dSum > 0.0000001?((*this)[i] /= dSum):(*this)[i]=0.0, i++) ;
        }

        void clear() {
            for (size_t i = 0; i < 12; ++i) (*this)[i] = 0.0;
        }
    };

    class TCSVector : public std::valarray<double>
    {
    public:
        TCSVector() : std::valarray<double>() {
            resize(6, 0.0f);
        }
            
        virtual ~TCSVector() {};

        void printDebug() {
            for (int i = 0; i < int(size()); i++) {
                std::cout <<  (*this)[i] << ";";
            }
            std::cout << std::endl;
        }
            
        double magnitude() const {
            double dMag = 0.0;
                    
            for (size_t i = 0; i < 6; i++) {
                dMag += std::pow((*this)[i], 2.0);
            }
                    
            return std::sqrt(dMag);
        }
    };

    class TonalEstimator
    {
    public:
        TonalEstimator();
        virtual ~TonalEstimator();
        TCSVector transform2TCS(const ChromaVector& rVector);
        
    protected:
        std::valarray< std::valarray<double> > m_Basis;
    };


    typedef std::vector<std::pair<long, TCSVector> > vectorlist_t;

    class TCSGram
    {
    public: 
        TCSGram();
        ~TCSGram();
        void getTCSVector(int, TCSVector&) const;
        void addTCSVector(const TCSVector&);
        long getTime(size_t) const;
        long getDuration() const;
        void printDebug();
        int getSize() const { return m_VectorList.size(); }
        void reserve(size_t uSize) { m_VectorList.reserve(uSize); }
        void clear() { m_VectorList.clear(); }
        void setFrameDuration(const double dFrameDurationMS) {
            m_dFrameDurationMS = dFrameDurationMS;
        }
        void setNumBins(const unsigned int uNumBins) {
            m_uNumBins = uNumBins;
        }
        void normalize();
    protected:
        vectorlist_t m_VectorList;
        unsigned int m_uNumBins;
        double m_dFrameDurationMS;
    };    

    typedef std::valarray<double> ChangeDistance;

    struct ChangeDFConfig
    {
        int smoothingWidth;
    };

    class ChangeDetectionFunction
    {
    public:
        ChangeDetectionFunction(ChangeDFConfig);
        ~ChangeDetectionFunction();
        ChangeDistance process(const TCSGram& rTCSGram);
    private:
        void setFilterWidth(const int iWidth);
            
    private:
        std::valarray<double> m_vaGaussian;
        double m_dFilterSigma;
        int m_iFilterWidth;
    };    


    ChangeDetectionFunction::ChangeDetectionFunction(ChangeDFConfig config) :
        m_dFilterSigma(0.0), m_iFilterWidth(0)
    {
        setFilterWidth(config.smoothingWidth);
    }

    ChangeDetectionFunction::~ChangeDetectionFunction()
    {
    }

    void ChangeDetectionFunction::setFilterWidth(const int iWidth)
    {
        m_iFilterWidth = iWidth*2+1;
            
        // it is assumed that the gaussian is 0 outside of +/- FWHM
        // => filter width = 2*FWHM = 2*2.3548*sigma
        m_dFilterSigma = double(m_iFilterWidth) / double(2*2.3548);
        m_vaGaussian.resize(m_iFilterWidth);
            
        double dScale = 1.0 / (m_dFilterSigma*sqrt(2*M_PI));
            
        for (int x = -(m_iFilterWidth-1)/2; x <= (m_iFilterWidth-1)/2; x++) {
            double w = dScale * std::exp ( -(x*x)/(2*m_dFilterSigma*m_dFilterSigma) );
            m_vaGaussian[x + (m_iFilterWidth-1)/2] = w;
        }
            
    #ifdef DEBUG_CHANGE_DETECTION_FUNCTION
        std::cerr << "Filter sigma: " << m_dFilterSigma << std::endl;
        std::cerr << "Filter width: " << m_iFilterWidth << std::endl;
    #endif
    }


    ChangeDistance ChangeDetectionFunction::process(const TCSGram& rTCSGram)
    {
        ChangeDistance retVal;
        retVal.resize(rTCSGram.getSize(), 0.0);
            
        TCSGram smoothedTCSGram;

        for (int iPosition = 0; iPosition < rTCSGram.getSize(); iPosition++) {
            
            int iSkipLower = 0;
            
            int iLowerPos = iPosition - (m_iFilterWidth-1)/2;
            int iUpperPos = iPosition + (m_iFilterWidth-1)/2;
            
            if (iLowerPos < 0) {
                iSkipLower = -iLowerPos;
                iLowerPos = 0;
            }
            
            if (iUpperPos >= rTCSGram.getSize()) {
                int iMaxIndex = rTCSGram.getSize() - 1;
                iUpperPos = iMaxIndex;
            }
            
            TCSVector smoothedVector;

            // for every bin of the vector, calculate the smoothed value
            for (int iPC = 0; iPC < 6; iPC++) {       

                size_t j = 0;
                double dSmoothedValue = 0.0;
                TCSVector rCV;
                    
                for (int i = iLowerPos; i <= iUpperPos; i++) {
                    rTCSGram.getTCSVector(i, rCV);
                    dSmoothedValue += m_vaGaussian[iSkipLower + j++] * rCV[iPC];
                }

                smoothedVector[iPC] = dSmoothedValue;
            }
                    
            smoothedTCSGram.addTCSVector(smoothedVector);
        }

        for (int iPosition = 0; iPosition < rTCSGram.getSize(); iPosition++) {
            
            /*
            TODO: calculate a confidence measure for the current estimation
            if the current estimate is not confident enough, look further into the future/the past
            e.g., High frequency content, zero crossing rate, spectral flatness
            */
                    
            TCSVector nextTCS;
            TCSVector previousTCS;
                    
            int iWindow = 1;

            // while (previousTCS.magnitude() < 0.1 && (iPosition-iWindow) > 0)
            {
                smoothedTCSGram.getTCSVector(iPosition-iWindow, previousTCS);
                // std::cout << previousTCS.magnitude() << std::endl;
                iWindow++;
            }
                    
            iWindow = 1;
                    
            // while (nextTCS.magnitude() < 0.1 && (iPosition+iWindow) < (rTCSGram.getSize()-1) )
            {
                smoothedTCSGram.getTCSVector(iPosition+iWindow, nextTCS);
                iWindow++;
            }

            double distance = 0.0;
            // Euclidean distance
            for (size_t j = 0; j < 6; j++) {
                distance += std::pow(nextTCS[j] - previousTCS[j], 2.0);
            }
            
            retVal[iPosition] = std::pow(distance, 0.5);
        }

        return retVal;
    }    


    TCSGram::TCSGram() :
        m_uNumBins(6)
    {
    }

    TCSGram::~TCSGram()
    {
    }


    void TCSGram::getTCSVector(int iPosition, TCSVector& rTCSVector) const
    {
        if (iPosition < 0) {
            rTCSVector = TCSVector();
        } else if (iPosition >= int(m_VectorList.size())) {
            rTCSVector = TCSVector();
        } else {
            rTCSVector = m_VectorList[iPosition].second;
        }
    }

    long TCSGram::getTime(size_t uPosition) const
    {
        return m_VectorList[uPosition].first;
    }


    void TCSGram::addTCSVector(const TCSVector& rTCSVector)
    {
        size_t uSize = m_VectorList.size();
        long lMilliSeconds = static_cast<long>(uSize*m_dFrameDurationMS);
        std::pair<long, TCSVector> p; 
        p.first = lMilliSeconds;
        p.second = rTCSVector;
            
        m_VectorList.push_back(p);
    }

    long TCSGram::getDuration() const
    {
        size_t uSize = m_VectorList.size();
        return static_cast<long>(uSize*m_dFrameDurationMS);
    }

    void TCSGram::printDebug()
    {
        vectorlist_t::iterator vectorIterator = m_VectorList.begin();
            
        while (vectorIterator != m_VectorList.end()) {
            vectorIterator->second.printDebug();
            vectorIterator++;
        }
    }

    TonalEstimator::TonalEstimator()
    {
        m_Basis.resize(6);

        int i = 0;
            
            
        // circle of fifths
        m_Basis[i].resize(12);
        for (int iP = 0; iP < 12; iP++) {
            m_Basis[i][iP] = std::sin( (7.0 / 6.0) * iP * M_PI);
        }
            
        i++;

        m_Basis[i].resize(12);
        for (int iP = 0; iP < 12; iP++) {
            m_Basis[i][iP] = std::cos( (7.0 / 6.0) * iP * M_PI);
        }
            
        i++;
            
            
        // circle of major thirds
        m_Basis[i].resize(12);
        for (int iP = 0; iP < 12; iP++) {
            m_Basis[i][iP] = 0.6 * std::sin( (2.0 / 3.0) * iP * M_PI);
        }
            
        i++;

        m_Basis[i].resize(12);
        for (int iP = 0; iP < 12; iP++) {
            m_Basis[i][iP] = 0.6 * std::cos( (2.0 / 3.0) * iP * M_PI);
        }

        i++;


        // circle of minor thirds
        m_Basis[i].resize(12);
        for (int iP = 0; iP < 12; iP++) {
            m_Basis[i][iP] = 1.1 * std::sin( (3.0 / 2.0) * iP * M_PI);
        }
            
        i++;

        m_Basis[i].resize(12);
        for (int iP = 0; iP < 12; iP++) {
            m_Basis[i][iP] = 1.1 * std::cos( (3.0 / 2.0) * iP * M_PI);
        }

    }

    TonalEstimator::~TonalEstimator()
    {
    }

    TCSVector TonalEstimator::transform2TCS(const ChromaVector& rVector)
    {
        TCSVector vaRetVal;
        vaRetVal.resize(6, 0.0);
                    
        for (int i = 0; i < 6; i++) {
            for (int iP = 0; iP < 12; iP++) {
                vaRetVal[i] += m_Basis[i][iP] * rVector[iP];
            }
        }
            
        return vaRetVal;
    }
}