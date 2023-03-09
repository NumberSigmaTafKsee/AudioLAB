#pragma once
#include "COscillator.h"

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

