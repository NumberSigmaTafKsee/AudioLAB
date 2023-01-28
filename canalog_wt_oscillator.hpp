#pragma once
#include "canalog_oscillator.hpp"

#define WT_LENGTH 512
#define NUM_TABLES 9

class CWTOscillator : public COscillator
{
public:
    CWTOscillator(void);
    ~CWTOscillator(void);

protected:	
    // oscillator
    DspFloatType m_dReadIndex;
    DspFloatType m_dWT_inc;
        
    // our tables
    DspFloatType m_dSineTable[WT_LENGTH];
    DspFloatType* m_pSawTables[NUM_TABLES];
    DspFloatType* m_pTriangleTables[NUM_TABLES];

    // for storing current table
    DspFloatType* m_pCurrentTable;
    int m_nCurrentTableIndex; //0 - 9

    // correction factor table sum-of-sawtooth
    DspFloatType m_dSquareCorrFactor[NUM_TABLES];

    // find the table with the proper number of harmonics for our pitch
    int getTableIndex();
    void selectTable();

    // create an destory tables
    void createWaveTables();
    void destroyWaveTables();

    // do the selected wavetable
    DspFloatType doWaveTable(DspFloatType& dReadIndex, DspFloatType dWT_inc);

    // for square wave 
    DspFloatType doSquareWave();

public:
    // --- while loop is for phase modulation
    inline void checkWrapIndex(DspFloatType& dIndex)
    {
        while(dIndex < 0.0)
            dIndex += WT_LENGTH;

        while(dIndex >= WT_LENGTH)
            dIndex -= WT_LENGTH;
    }

    // typical overrides
    virtual void reset();
    virtual void startOscillator();
    virtual void stopOscillator();
    
    // render a sample
    // for LFO:		pAuxOutput = QuadPhaseOutput
    //     Pitched: pAuxOutput = Right channel (return value is left Channel
    virtual DspFloatType doOscillate(DspFloatType* pAuxOutput = NULL);

    // wave table specific
    virtual void setSampleRate(DspFloatType dFs);
    virtual void update();
};


CWTOscillator::CWTOscillator(void)
{
    // --- clear out arrays
    memset(m_pSawTables, 0, NUM_TABLES*sizeof(DspFloatType*));
    memset(m_pTriangleTables, 0, NUM_TABLES*sizeof(DspFloatType*));

    // --- init variables
    m_dReadIndex = 0.0;
    m_dWT_inc = 0.0;
    m_nCurrentTableIndex = 0;

    // --- setup correction factors (empirical)
    m_dSquareCorrFactor[0] = 0.5;
    m_dSquareCorrFactor[1] = 0.5;
    m_dSquareCorrFactor[2] = 0.5;
    m_dSquareCorrFactor[3] = 0.49;
    m_dSquareCorrFactor[4] = 0.48;
    m_dSquareCorrFactor[5] = 0.468;
    m_dSquareCorrFactor[6] = 0.43;
    m_dSquareCorrFactor[7] = 0.34;
    m_dSquareCorrFactor[8] = 0.25;

    // --- default to SINE
    m_pCurrentTable = &m_dSineTable[0];
}

CWTOscillator::~CWTOscillator(void)
{
    destroyWaveTables();
}

void CWTOscillator::reset()
{
    COscillator::reset();

    // --- back to top of buffer
    m_dReadIndex = 0.0;
}

void CWTOscillator::startOscillator()
{
    reset();
    m_bNoteOn = true;
}

void CWTOscillator::stopOscillator()
{
    m_bNoteOn = false;
}

void CWTOscillator::update()
{
    COscillator::update();

    // --- calculate the inc value
    m_dWT_inc = WT_LENGTH*m_dInc;

    // --- select the table
    selectTable();
}

void CWTOscillator::setSampleRate(DspFloatType dFs)
{
    // Test for change
    bool bNewSR = m_dSampleRate != dFs ? true : false; // Base class first
    COscillator::setSampleRate(dFs); // Recreate the tables only if sample rate has changed

    if (bNewSR)
    {
        // --- then recreate
        destroyWaveTables();
        createWaveTables();   
        return;
    }   
    
    createWaveTables();
}


void CWTOscillator::createWaveTables()
{
    // create the tables
    //
    // SINE: only need one table
    for(int i = 0; i < WT_LENGTH; i++)
    {
        // sample the sinusoid, WT_LENGTH points
        // sin(wnT) = sin(2pi*i/WT_LENGTH)
        m_dSineTable[i] = sin(((DspFloatType)i/WT_LENGTH)*(2*M_PI));
    }

    // SAW, TRIANGLE: need 10 tables
    DspFloatType dSeedFreq = 27.5; // Note A0, bottom of piano
    for(int j = 0; j < NUM_TABLES; j++)
    {
        DspFloatType* pSawTable = new DspFloatType[WT_LENGTH];
        memset(pSawTable, 0, WT_LENGTH*sizeof(DspFloatType));

        DspFloatType* pTriTable = new DspFloatType[WT_LENGTH];
        memset(pTriTable, 0, WT_LENGTH*sizeof(DspFloatType));

        int nHarms = (int)((m_dSampleRate/2.0/dSeedFreq) - 1.0);
        int nHalfHarms = (int)((DspFloatType)nHarms/2.0);

        DspFloatType dMaxSaw = 0;
        DspFloatType dMaxTri = 0;

        for(int i = 0; i < WT_LENGTH; i++)
        {
            // sawtooth: += (-1)^g+1(1/g)sin(wnT)
            for(int g = 1; g <= nHarms; g++)
            {
                // Lanczos Sigma Factor
                DspFloatType x = g*M_PI/nHarms;
                DspFloatType sigma = sin(x)/x;

                // only apply to partials above fundamental
                if(g == 1)
                    sigma = 1.0;

                DspFloatType n = DspFloatType(g);
                pSawTable[i] += std::pow((DspFloatType)-1.0,(DspFloatType)(g+1))*(1.0/n)*sigma*sin(2.0*M_PI*i*n/WT_LENGTH);
            }

            // triangle: += (-1)^g(1/(2g+1+^2)sin(w(2n+1)T)
            // NOTE: the limit is nHalfHarms here because of the way the sum is constructed
            // (look at the (2n+1) components
            for(int g = 0; g <= nHalfHarms; g++)
            {
                DspFloatType n = DspFloatType(g);
                pTriTable[i] += std::pow((DspFloatType)-1.0, (DspFloatType)n)*(1.0/std::pow((DspFloatType)(2*n + 1),(DspFloatType)2.0))*sin(2.0*M_PI*(2.0*n + 1)*i/WT_LENGTH);
            }

            // store the max values
            if(i == 0)
            {
                dMaxSaw = pSawTable[i];
                dMaxTri = pTriTable[i];
            }
            else
            {
                // test and store
                if(pSawTable[i] > dMaxSaw)
                    dMaxSaw = pSawTable[i];

                if(pTriTable[i] > dMaxTri)
                    dMaxTri = pTriTable[i];
            }
        }
        // normalize
        for(int i = 0; i < WT_LENGTH; i++)
        {
            // normalize it
            pSawTable[i] /= dMaxSaw;
            pTriTable[i] /= dMaxTri;
        }

        // store
        m_pSawTables[j] = pSawTable;
        m_pTriangleTables[j] = pTriTable;

        dSeedFreq *= 2.0;
    }
}

void CWTOscillator::destroyWaveTables()
{
    for(int i = 0; i < NUM_TABLES; i++)
    {
        DspFloatType* p = m_pSawTables[i];
        if(p)
        {
            delete [] p;
            m_pSawTables[i] = 0;
        }

        p = m_pTriangleTables[i];
        if(p)
        {
            delete [] p;
            m_pTriangleTables[i] = 0;
        }
    }
}

// get table index based on current m_dFo
int CWTOscillator::getTableIndex()
{
    if(m_uWaveform == SINE)
        return -1;

    DspFloatType dSeedFreq = 27.5; // Note A0, bottom of piano
    for(int j = 0; j < NUM_TABLES; j++)
    {
        if(m_dFo <= dSeedFreq)
        {
            return j;
        }

        dSeedFreq *= 2.0;
    }

    return -1;
}

void CWTOscillator::selectTable()
{
    m_nCurrentTableIndex = getTableIndex();

    // if the frequency is high enough, the sine table will be returned
    // even for non-sinusoidal waves; anything about 10548 Hz is one
    // harmonic only (sine)
    if(m_nCurrentTableIndex < 0)
    {
        m_pCurrentTable = &m_dSineTable[0];
        return;
    }

    // choose table
    if(m_uWaveform == SAW1 || m_uWaveform == SAW2 || m_uWaveform == SAW3 || m_uWaveform == SQUARE)
        m_pCurrentTable = m_pSawTables[m_nCurrentTableIndex];
    else if(m_uWaveform == TRI)
        m_pCurrentTable = m_pTriangleTables[m_nCurrentTableIndex];
}


DspFloatType CWTOscillator::doWaveTable(DspFloatType& dReadIndex, DspFloatType dWT_inc)
{
    DspFloatType dOut = 0;

    // apply phase modulation, if any
    DspFloatType dModReadIndex = dReadIndex + m_dPhaseMod*WT_LENGTH;

    // check for multi-wrapping on new read index
    checkWrapIndex(dModReadIndex);

    // get INT part
    int nReadIndex = abs((int)dModReadIndex);

    // get FRAC part
    DspFloatType fFrac = dModReadIndex - nReadIndex;

    // setup second index for interpolation; wrap the buffer if needed
    int nReadIndexNext = nReadIndex + 1 > WT_LENGTH-1 ? 0 :  nReadIndex + 1;

    // interpolate the output
    dOut = dLinTerp(0, 1, m_pCurrentTable[nReadIndex], m_pCurrentTable[nReadIndexNext], fFrac);

    // add the increment for next time
    dReadIndex += dWT_inc;

    // check for wrap
    checkWrapIndex(dReadIndex);

    return dOut;
}

DspFloatType CWTOscillator::doSquareWave()
{
    DspFloatType dPW = m_dPulseWidth/100.0;
    DspFloatType dPWIndex = m_dReadIndex + dPW*WT_LENGTH;

    // --- render first sawtooth using dReadIndex
    DspFloatType dSaw1 = doWaveTable(m_dReadIndex, m_dWT_inc);

    // --- find the phase shifted output
    if(m_dWT_inc >= 0)
    {
        if(dPWIndex >= WT_LENGTH)
            dPWIndex = dPWIndex - WT_LENGTH;
    }
    else
    {
        if(dPWIndex < 0)
            dPWIndex = WT_LENGTH + dPWIndex;
    }

    // --- render second sawtooth using dPWIndex (shifted)
    DspFloatType dSaw2 = doWaveTable(dPWIndex, m_dWT_inc);

    // --- find the correction factor from the table
    DspFloatType dSqAmp = m_dSquareCorrFactor[m_nCurrentTableIndex];

    // --- then subtract
    DspFloatType dOut = dSqAmp*dSaw1 -  dSqAmp*dSaw2;

    // --- calculate the DC correction factor
    DspFloatType dCorr = 1.0/dPW;
    if(dPW < 0.5)
        dCorr = 1.0/(1.0-dPW);

    // --- apply correction
    dOut *= dCorr;

    return dOut;
}

DspFloatType CWTOscillator::doOscillate(DspFloatType* pAuxOutput)
{
    if(!m_bNoteOn)
    {
        if(pAuxOutput)
            *pAuxOutput = 0.0;

        return 0.0;
    }

    // if square, it has its own routine
    if(m_uWaveform == SQUARE && m_nCurrentTableIndex >= 0)
    {
        DspFloatType dOut = doSquareWave();
        if(pAuxOutput)
            *pAuxOutput = dOut;

        return dOut;
    }

    // --- get output
    DspFloatType dOutSample = doWaveTable(m_dReadIndex, m_dWT_inc);

    // mono oscillator
    if(pAuxOutput)
        *pAuxOutput = dOutSample*m_dAmplitude*m_dAmpMod;

    return dOutSample*m_dAmplitude*m_dAmpMod;
}
