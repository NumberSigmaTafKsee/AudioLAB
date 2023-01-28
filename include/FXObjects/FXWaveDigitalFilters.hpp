#pragma once

// ------------------------------------------------------------------ //
// --- WDF LIBRARY -------------------------------------------------- //
// ------------------------------------------------------------------ //

/**
\class IComponentAdaptor
\ingroup Interfaces
\brief
Use this interface for objects in the WDF Ladder Filter library; see book for more information.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class IComponentAdaptor
{
public:
	/** initialize with source resistor R1 */
	virtual void initialize(DspFloatType _R1) {}

	/** initialize all downstream adaptors in the chain */
	virtual void initializeAdaptorChain() {}

	/** set input value into component port  */
	virtual void setInput(DspFloatType _in) {}

	/** get output value from component port  */
	virtual DspFloatType getOutput() { return 0.0; }

	// --- for adaptors
	/** ADAPTOR: set input port 1  */
	virtual void setInput1(DspFloatType _in1) = 0;

	/** ADAPTOR: set input port 2  */
	virtual void setInput2(DspFloatType _in2) = 0;

	/** ADAPTOR: set input port 3 */
	virtual void setInput3(DspFloatType _in3) = 0;

	/** ADAPTOR: get output port 1 value */
	virtual DspFloatType getOutput1() = 0;

	/** ADAPTOR: get output port 2 value */
	virtual DspFloatType getOutput2() = 0;

	/** ADAPTOR: get output port 3 value */
	virtual DspFloatType getOutput3() = 0;

	/** reset the object with new sample rate */
	virtual void reset(DspFloatType _sampleRate) {}

	/** get the commponent resistance from the attached object at Port3 */
	virtual DspFloatType getComponentResistance() { return 0.0; }

	/** get the commponent conductance from the attached object at Port3 */
	virtual DspFloatType getComponentConductance() { return 0.0; }

	/** update the commponent resistance at Port3 */
	virtual void updateComponentResistance() {}

	/** set an individual component value (may be R, L, or C */
	virtual void setComponentValue(DspFloatType _componentValue) { }

	/** set LC combined values */
	virtual void setComponentValue_LC(DspFloatType componentValue_L, DspFloatType componentValue_C) { }

	/** set RL combined values */
	virtual void setComponentValue_RL(DspFloatType componentValue_R, DspFloatType componentValue_L) { }

	/** set RC combined values */
	virtual void setComponentValue_RC(DspFloatType componentValue_R, DspFloatType componentValue_C) { }

	/** get a component value */
	virtual DspFloatType getComponentValue() { return 0.0; }
};

// ------------------------------------------------------------------ //
// --- WDF COMPONENTS & COMMBO COMPONENTS --------------------------- //
// ------------------------------------------------------------------ //
/**
\class WdfResistor
\ingroup WDF-Objects
\brief
The WdfResistor object implements the reflection coefficient and signal flow through
a WDF simulated resistor.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfResistor : public IComponentAdaptor
{
public:
	WdfResistor(DspFloatType _componentValue) { componentValue = _componentValue; }
	WdfResistor() { }
	virtual ~WdfResistor() {}

	/** set sample rate and update component */
	void setSampleRate(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		updateComponentResistance();
	}

	/** get component's value as a resistance */
	virtual DspFloatType getComponentResistance() { return componentResistance; }

	/** get component's value as a conducatance (or admittance) */
	virtual DspFloatType getComponentConductance() { return 1.0 / componentResistance; }

	/** get the component value */
	virtual DspFloatType getComponentValue() { return componentValue; }

	/** set the component value */
	virtual void setComponentValue(DspFloatType _componentValue)
	{
		componentValue = _componentValue;
		updateComponentResistance();
	}

	/** change the resistance of component */
	virtual void updateComponentResistance() { componentResistance = componentValue; }

	/** reset the component; clear registers */
	virtual void reset(DspFloatType _sampleRate) { setSampleRate(_sampleRate);  zRegister = 0.0; }

	/** set input value into component; NOTE: resistor is dead-end energy sink so this function does nothing */
	virtual void setInput(DspFloatType in) {}

	/** get output value; NOTE: a WDF resistor produces no reflected output */
	virtual DspFloatType getOutput() { return 0.0; }

	/** get output1 value; only one resistor output (not used) */
	virtual DspFloatType getOutput1() { return  getOutput(); }

	/** get output2 value; only one resistor output (not used) */
	virtual DspFloatType getOutput2() { return  getOutput(); }

	/** get output3 value; only one resistor output (not used) */
	virtual DspFloatType getOutput3() { return  getOutput(); }

	/** set input1 value; not used for components */
	virtual void setInput1(DspFloatType _in1) {}

	/** set input2 value; not used for components */
	virtual void setInput2(DspFloatType _in2) {}

	/** set input3 value; not used for components */
	virtual void setInput3(DspFloatType _in3) {}

protected:
	DspFloatType zRegister = 0.0;			///< storage register (not used with resistor)
	DspFloatType componentValue = 0.0;	///< component value in electronic form (ohm, farad, henry)
	DspFloatType componentResistance = 0.0;///< simulated resistance
	DspFloatType sampleRate = 0.0;		///< sample rate
};


/**
\class WdfCapacitor
\ingroup WDF-Objects
\brief
The WdfCapacitor object implements the reflection coefficient and signal flow through
a WDF simulated capacitor.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfCapacitor : public IComponentAdaptor
{
public:
	WdfCapacitor(DspFloatType _componentValue) { componentValue = _componentValue; }
	WdfCapacitor() { }
	virtual ~WdfCapacitor() {}

	/** set sample rate and update component */
	void setSampleRate(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		updateComponentResistance();
	}

	/** get component's value as a resistance */
	virtual DspFloatType getComponentResistance() { return componentResistance; }

	/** get component's value as a conducatance (or admittance) */
	virtual DspFloatType getComponentConductance() { return 1.0 / componentResistance; }

	/** get the component value */
	virtual DspFloatType getComponentValue() { return componentValue; }

	/** set the component value */
	virtual void setComponentValue(DspFloatType _componentValue)
	{
		componentValue = _componentValue;
		updateComponentResistance();
	}

	/** change the resistance of component */
	virtual void updateComponentResistance()
	{
		componentResistance = 1.0 / (2.0*componentValue*sampleRate);
	}

	/** reset the component; clear registers */
	virtual void reset(DspFloatType _sampleRate) { setSampleRate(_sampleRate); zRegister = 0.0; }

	/** set input value into component; NOTE: capacitor sets value into register*/
	virtual void setInput(DspFloatType in) { zRegister = in; }

	/** get output value; NOTE: capacitor produces reflected output */
	virtual DspFloatType getOutput() { return zRegister; }	// z^-1

	/** get output1 value; only one capacitor output (not used) */
	virtual DspFloatType getOutput1() { return  getOutput(); }

	/** get output2 value; only one capacitor output (not used) */
	virtual DspFloatType getOutput2() { return  getOutput(); }

	/** get output3 value; only one capacitor output (not used) */
	virtual DspFloatType getOutput3() { return  getOutput(); }

	/** set input1 value; not used for components */
	virtual void setInput1(DspFloatType _in1) {}

	/** set input2 value; not used for components */
	virtual void setInput2(DspFloatType _in2) {}

	/** set input3 value; not used for components */
	virtual void setInput3(DspFloatType _in3) {}

protected:
	DspFloatType zRegister = 0.0;			///< storage register (not used with resistor)
	DspFloatType componentValue = 0.0;	///< component value in electronic form (ohm, farad, henry)
	DspFloatType componentResistance = 0.0;///< simulated resistance
	DspFloatType sampleRate = 0.0;		///< sample rate
};

/**
\class WdfInductor
\ingroup WDF-Objects
\brief
The WdfInductor object implements the reflection coefficient and signal flow through
a WDF simulated inductor.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfInductor : public IComponentAdaptor
{
public:
	WdfInductor(DspFloatType _componentValue) { componentValue = _componentValue; }
	WdfInductor() { }
	virtual ~WdfInductor() {}

	/** set sample rate and update component */
	void setSampleRate(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		updateComponentResistance();
	}

	/** get component's value as a resistance */
	virtual DspFloatType getComponentResistance() { return componentResistance; }

	/** get component's value as a conducatance (or admittance) */
	virtual DspFloatType getComponentConductance() { return 1.0 / componentResistance; }

	/** get the component value */
	virtual DspFloatType getComponentValue() { return componentValue; }

	/** set the component value */
	virtual void setComponentValue(DspFloatType _componentValue)
	{
		componentValue = _componentValue;
		updateComponentResistance();
	}

	/** change the resistance of component R(L) = 2Lfs */
	virtual void updateComponentResistance(){ componentResistance = 2.0*componentValue*sampleRate;}

	/** reset the component; clear registers */
	virtual void reset(DspFloatType _sampleRate) { setSampleRate(_sampleRate); zRegister = 0.0; }

	/** set input value into component; NOTE: inductor sets value into storage register */
	virtual void setInput(DspFloatType in) { zRegister = in; }

	/** get output value; NOTE: a WDF inductor produces reflected output that is inverted */
	virtual DspFloatType getOutput() { return -zRegister; } // -z^-1

	/** get output1 value; only one resistor output (not used) */
	virtual DspFloatType getOutput1() { return  getOutput(); }

	/** get output2 value; only one resistor output (not used) */
	virtual DspFloatType getOutput2() { return  getOutput(); }

	/** get output3 value; only one resistor output (not used) */
	virtual DspFloatType getOutput3() { return  getOutput(); }

	/** set input1 value; not used for components */
	virtual void setInput1(DspFloatType _in1) {}

	/** set input2 value; not used for components */
	virtual void setInput2(DspFloatType _in2) {}

	/** set input3 value; not used for components */
	virtual void setInput3(DspFloatType _in3) {}

protected:
	DspFloatType zRegister = 0.0;			///< storage register (not used with resistor)
	DspFloatType componentValue = 0.0;	///< component value in electronic form (ohm, farad, henry)
	DspFloatType componentResistance = 0.0;///< simulated resistance
	DspFloatType sampleRate = 0.0;		///< sample rate
};


/**
\class WdfSeriesLC
\ingroup WDF-Objects
\brief
The WdfSeriesLC object implements the reflection coefficient and signal flow through
a WDF simulated series LC pair.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
SEE: p143 "Design of Wave Digital Filters" Psenicka, Ugalde, Romero M.
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfSeriesLC : public IComponentAdaptor
{
public:
	WdfSeriesLC() {}
	WdfSeriesLC(DspFloatType _componentValue_L, DspFloatType _componentValue_C)
	{
		componentValue_L = _componentValue_L;
		componentValue_C = _componentValue_C;
	}
	virtual ~WdfSeriesLC() {}

	/** set sample rate and update component */
	void setSampleRate(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		updateComponentResistance();
	}

	/** get component's value as a resistance */
	virtual DspFloatType getComponentResistance() { return componentResistance; }

	/** get component's value as a conducatance (or admittance) */
	virtual DspFloatType getComponentConductance() { return 1.0 / componentResistance; }

	/** change the resistance of component; see FX book for details */
	virtual void updateComponentResistance()
	{
		RL = 2.0*componentValue_L*sampleRate;
		RC = 1.0 / (2.0*componentValue_C*sampleRate);
		componentResistance = RL + (1.0 / RC);
	}

	/** set both LC components at once */
	virtual void setComponentValue_LC(DspFloatType _componentValue_L, DspFloatType _componentValue_C)
	{
		componentValue_L = _componentValue_L;
		componentValue_C = _componentValue_C;
		updateComponentResistance();
	}

	/** set L component */
	virtual void setComponentValue_L(DspFloatType _componentValue_L)
	{
		componentValue_L = _componentValue_L;
		updateComponentResistance();
	}

	/** set C component */
	virtual void setComponentValue_C(DspFloatType _componentValue_C)
	{
		componentValue_C = _componentValue_C;
		updateComponentResistance();
	}

	/** get L component value */
	virtual DspFloatType getComponentValue_L() { return componentValue_L; }

	/** get C component value */
	virtual DspFloatType getComponentValue_C() { return componentValue_C; }

	/** reset the component; clear registers */
	virtual void reset(DspFloatType _sampleRate) { setSampleRate(_sampleRate); zRegister_L = 0.0; zRegister_C = 0.0; }

	/** set input value into component; NOTE: K is calculated here */
	virtual void setInput(DspFloatType in)
	{
		DspFloatType YC = 1.0 / RC;
		DspFloatType K = (1.0 - RL*YC) / (1.0 + RL*YC);
		DspFloatType N1 = K*(in - zRegister_L);
		zRegister_L = N1 + zRegister_C;
		zRegister_C = in;
	}

	/** get output value; NOTE: utput is located in zReg_L */
	virtual DspFloatType getOutput(){ return zRegister_L; }

	/** get output1 value; only one resistor output (not used) */
	virtual DspFloatType getOutput1() { return  getOutput(); }

	/** get output2 value; only one resistor output (not used) */
	virtual DspFloatType getOutput2() { return  getOutput(); }

	/** get output3 value; only one resistor output (not used) */
	virtual DspFloatType getOutput3() { return  getOutput(); }

	/** set input1 value; not used for components */
	virtual void setInput1(DspFloatType _in1) {}

	/** set input2 value; not used for components */
	virtual void setInput2(DspFloatType _in2) {}

	/** set input3 value; not used for components */
	virtual void setInput3(DspFloatType _in3) {}

protected:
	DspFloatType zRegister_L = 0.0; ///< storage register for L
	DspFloatType zRegister_C = 0.0; ///< storage register for C

	DspFloatType componentValue_L = 0.0; ///< component value L
	DspFloatType componentValue_C = 0.0; ///< component value C

	DspFloatType RL = 0.0; ///< RL value
	DspFloatType RC = 0.0; ///< RC value
	DspFloatType componentResistance = 0.0; ///< equivalent resistance of pair of components
	DspFloatType sampleRate = 0.0; ///< sample rate
};

/**
\class WdfParallelLC
\ingroup WDF-Objects
\brief
The WdfParallelLC object implements the reflection coefficient and signal flow through
a WDF simulated parallel LC pair.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
SEE: p143 "Design of Wave Digital Filters" Psenicka, Ugalde, Romero M.
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfParallelLC : public IComponentAdaptor
{
public:
	WdfParallelLC() {}
	WdfParallelLC(DspFloatType _componentValue_L, DspFloatType _componentValue_C)
	{
		componentValue_L = _componentValue_L;
		componentValue_C = _componentValue_C;
	}
	virtual ~WdfParallelLC() {}

	/** set sample rate and update component */
	void setSampleRate(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		updateComponentResistance();
	}

	/** get component's value as a resistance */
	virtual DspFloatType getComponentResistance() { return componentResistance; }

	/** get component's value as a conducatance (or admittance) */
	virtual DspFloatType getComponentConductance() { return 1.0 / componentResistance; }

	/** change the resistance of component; see FX book for details */
	virtual void updateComponentResistance()
	{
		RL = 2.0*componentValue_L*sampleRate;
		RC = 1.0 / (2.0*componentValue_C*sampleRate);
		componentResistance = (RC + 1.0 / RL);
	}

	/** set both LC components at once */
	virtual void setComponentValue_LC(DspFloatType _componentValue_L, DspFloatType _componentValue_C)
	{
		componentValue_L = _componentValue_L;
		componentValue_C = _componentValue_C;
		updateComponentResistance();
	}

	/** set L component */
	virtual void setComponentValue_L(DspFloatType _componentValue_L)
	{
		componentValue_L = _componentValue_L;
		updateComponentResistance();
	}

	/** set C component */
	virtual void setComponentValue_C(DspFloatType _componentValue_C)
	{
		componentValue_C = _componentValue_C;
		updateComponentResistance();
	}

	/** get L component value */
	virtual DspFloatType getComponentValue_L() { return componentValue_L; }

	/** get C component value */
	virtual DspFloatType getComponentValue_C() { return componentValue_C; }

	/** reset the component; clear registers */
	virtual void reset(DspFloatType _sampleRate) { setSampleRate(_sampleRate); zRegister_L = 0.0; zRegister_C = 0.0; }

	/** set input value into component; NOTE: K is calculated here */
	virtual void setInput(DspFloatType in)
	{
		DspFloatType YL = 1.0 / RL;
		DspFloatType K = (YL*RC - 1.0) / (YL*RC + 1.0);
		DspFloatType N1 = K*(in - zRegister_L);
		zRegister_L = N1 + zRegister_C;
		zRegister_C = in;
	}

	/** get output value; NOTE: output is located in -zReg_L */
	virtual DspFloatType getOutput(){ return -zRegister_L; }

	/** get output1 value; only one resistor output (not used) */
	virtual DspFloatType getOutput1() { return  getOutput(); }

	/** get output2 value; only one resistor output (not used) */
	virtual DspFloatType getOutput2() { return  getOutput(); }

	/** get output3 value; only one resistor output (not used) */
	virtual DspFloatType getOutput3() { return  getOutput(); }

	/** set input1 value; not used for components */
	virtual void setInput1(DspFloatType _in1) {}

	/** set input2 value; not used for components */
	virtual void setInput2(DspFloatType _in2) {}

	/** set input3 value; not used for components */
	virtual void setInput3(DspFloatType _in3) {}

protected:
	DspFloatType zRegister_L = 0.0; ///< storage register for L
	DspFloatType zRegister_C = 0.0; ///< storage register for C

	DspFloatType componentValue_L = 0.0; ///< component value L
	DspFloatType componentValue_C = 0.0; ///< component value C

	DspFloatType RL = 0.0; ///< RL value
	DspFloatType RC = 0.0; ///< RC value
	DspFloatType componentResistance = 0.0; ///< equivalent resistance of pair of components
	DspFloatType sampleRate = 0.0; ///< sample rate
};


/**
\class WdfSeriesRL
\ingroup WDF-Objects
\brief
The WdfSeriesRL object implements the reflection coefficient and signal flow through
a WDF simulated series RL pair.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
SEE: p143 "Design of Wave Digital Filters" Psenicka, Ugalde, Romero M.
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfSeriesRL : public IComponentAdaptor
{
public:
	WdfSeriesRL() {}
	WdfSeriesRL(DspFloatType _componentValue_R, DspFloatType _componentValue_L)
	{
		componentValue_L = _componentValue_L;
		componentValue_R = _componentValue_R;
	}
	virtual ~WdfSeriesRL() {}

	/** set sample rate and update component */
	void setSampleRate(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		updateComponentResistance();
	}

	/** get component's value as a resistance */
	virtual DspFloatType getComponentResistance() { return componentResistance; }

	/** get component's value as a conducatance (or admittance) */
	virtual DspFloatType getComponentConductance() { return 1.0 / componentResistance; }

	/** change the resistance of component; see FX book for details */
	virtual void updateComponentResistance()
	{
		RR = componentValue_R;
		RL = 2.0*componentValue_L*sampleRate;
		componentResistance = RR + RL;
		K = RR / componentResistance;
	}

	/** set both RL components at once */
	virtual void setComponentValue_RL(DspFloatType _componentValue_R, DspFloatType _componentValue_L)
	{
		componentValue_L = _componentValue_L;
		componentValue_R = _componentValue_R;
		updateComponentResistance();
	}

	/** set L component */
	virtual void setComponentValue_L(DspFloatType _componentValue_L)
	{
		componentValue_L = _componentValue_L;
		updateComponentResistance();
	}

	/** set R component */
	virtual void setComponentValue_R(DspFloatType _componentValue_R)
	{
		componentValue_R = _componentValue_R;
		updateComponentResistance();
	}

	/** get L component value */
	virtual DspFloatType getComponentValue_L() { return componentValue_L; }

	/** get R component value */
	virtual DspFloatType getComponentValue_R() { return componentValue_R; }

	/** reset the component; clear registers */
	virtual void reset(DspFloatType _sampleRate) { setSampleRate(_sampleRate); zRegister_L = 0.0; zRegister_C = 0.0; }

	/** set input value into component */
	virtual void setInput(DspFloatType in){ zRegister_L = in; }

	/** get output value; NOTE: see FX book for details */
	virtual DspFloatType getOutput()
	{
		DspFloatType NL = -zRegister_L;
		DspFloatType out = NL*(1.0 - K) - K*zRegister_C;
		zRegister_C = out;

		return out;
	}

	/** get output1 value; only one resistor output (not used) */
	virtual DspFloatType getOutput1() { return  getOutput(); }

	/** get output2 value; only one resistor output (not used) */
	virtual DspFloatType getOutput2() { return  getOutput(); }

	/** get output3 value; only one resistor output (not used) */
	virtual DspFloatType getOutput3() { return  getOutput(); }

	/** set input1 value; not used for components */
	virtual void setInput1(DspFloatType _in1) {}

	/** set input2 value; not used for components */
	virtual void setInput2(DspFloatType _in2) {}

	/** set input3 value; not used for components */
	virtual void setInput3(DspFloatType _in3) {}

protected:
	DspFloatType zRegister_L = 0.0; ///< storage register for L
	DspFloatType zRegister_C = 0.0;///< storage register for C (not used)
	DspFloatType K = 0.0;

	DspFloatType componentValue_L = 0.0;///< component value L
	DspFloatType componentValue_R = 0.0;///< component value R

	DspFloatType RL = 0.0; ///< RL value
	DspFloatType RC = 0.0; ///< RC value
	DspFloatType RR = 0.0; ///< RR value

	DspFloatType componentResistance = 0.0; ///< equivalent resistance of pair of componen
	DspFloatType sampleRate = 0.0; ///< sample rate
};

/**
\class WdfParallelRL
\ingroup WDF-Objects
\brief
The WdfParallelRL object implements the reflection coefficient and signal flow through
a WDF simulated parallel RL pair.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfParallelRL : public IComponentAdaptor
{
public:
	WdfParallelRL() {}
	WdfParallelRL(DspFloatType _componentValue_R, DspFloatType _componentValue_L)
	{
		componentValue_L = _componentValue_L;
		componentValue_R = _componentValue_R;
	}
	virtual ~WdfParallelRL() {}

	/** set sample rate and update component */
	void setSampleRate(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		updateComponentResistance();
	}

	/** get component's value as a resistance */
	virtual DspFloatType getComponentResistance() { return componentResistance; }

	/** get component's value as a conducatance (or admittance) */
	virtual DspFloatType getComponentConductance() { return 1.0 / componentResistance; }

	/** change the resistance of component; see FX book for details */
	virtual void updateComponentResistance()
	{
		RR = componentValue_R;
		RL = 2.0*componentValue_L*sampleRate;
		componentResistance = 1.0 / ((1.0 / RR) + (1.0 / RL));
		K = componentResistance / RR;
	}


	/** set both RL components at once */
	virtual void setComponentValue_RL(DspFloatType _componentValue_R, DspFloatType _componentValue_L)
	{
		componentValue_L = _componentValue_L;
		componentValue_R = _componentValue_R;
		updateComponentResistance();
	}

	/** set L component */
	virtual void setComponentValue_L(DspFloatType _componentValue_L)
	{
		componentValue_L = _componentValue_L;
		updateComponentResistance();
	}

	/** set R component */
	virtual void setComponentValue_R(DspFloatType _componentValue_R)
	{
		componentValue_R = _componentValue_R;
		updateComponentResistance();
	}

	/** get L component value */
	virtual DspFloatType getComponentValue_L() { return componentValue_L; }

	/** get R component value */
	virtual DspFloatType getComponentValue_R() { return componentValue_R; }

	/** reset the component; clear registers */
	virtual void reset(DspFloatType _sampleRate) { setSampleRate(_sampleRate); zRegister_L = 0.0; zRegister_C = 0.0; }

	/** set input value into component */
	virtual void setInput(DspFloatType in){ zRegister_L = in; }

	/** get output value; NOTE: see FX book for details */
	virtual DspFloatType getOutput()
	{
		DspFloatType NL = -zRegister_L;
		DspFloatType out = NL*(1.0 - K) + K*zRegister_C;
		zRegister_C = out;
		return out;
	}

	/** get output1 value; only one resistor output (not used) */
	virtual DspFloatType getOutput1() { return  getOutput(); }

	/** get output2 value; only one resistor output (not used) */
	virtual DspFloatType getOutput2() { return  getOutput(); }

	/** get output3 value; only one resistor output (not used) */
	virtual DspFloatType getOutput3() { return  getOutput(); }

	/** set input1 value; not used for components */
	virtual void setInput1(DspFloatType _in1) {}

	/** set input2 value; not used for components */
	virtual void setInput2(DspFloatType _in2) {}

	/** set input3 value; not used for components */
	virtual void setInput3(DspFloatType _in3) {}

protected:
	DspFloatType zRegister_L = 0.0;	///< storage register for L
	DspFloatType zRegister_C = 0.0;	///< storage register for L
	DspFloatType K = 0.0;				///< K value

	DspFloatType componentValue_L = 0.0;	///< component value L
	DspFloatType componentValue_R = 0.0;	///< component value R

	DspFloatType RL = 0.0;	///< RL value
	DspFloatType RC = 0.0;	///< RC value
	DspFloatType RR = 0.0;	///< RR value

	DspFloatType componentResistance = 0.0; ///< equivalent resistance of pair of components
	DspFloatType sampleRate = 0.0; ///< sample rate
};

/**
\class WdfSeriesRC
\ingroup WDF-Objects
\brief
The WdfSeriesRC object implements the reflection coefficient and signal flow through
a WDF simulated series RC pair.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
SEE: p143 "Design of Wave Digital Filters" Psenicka, Ugalde, Romero M.
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfSeriesRC : public IComponentAdaptor
{
public:
	WdfSeriesRC() {}
	WdfSeriesRC(DspFloatType _componentValue_R, DspFloatType _componentValue_C)
	{
		componentValue_C = _componentValue_C;
		componentValue_R = _componentValue_R;
	}
	virtual ~WdfSeriesRC() {}

	/** set sample rate and update component */
	void setSampleRate(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		updateComponentResistance();
	}

	/** get component's value as a resistance */
	virtual DspFloatType getComponentResistance() { return componentResistance; }

	/** get component's value as a conducatance (or admittance) */
	virtual DspFloatType getComponentConductance() { return 1.0 / componentResistance; }

	/** change the resistance of component; see FX book for details */
	virtual void updateComponentResistance()
	{
		RR = componentValue_R;
		RC = 1.0 / (2.0*componentValue_C*sampleRate);
		componentResistance = RR + RC;
		K = RR / componentResistance;
	}

	/** set both RC components at once */
	virtual void setComponentValue_RC(DspFloatType _componentValue_R, DspFloatType _componentValue_C)
	{
		componentValue_R = _componentValue_R;
		componentValue_C = _componentValue_C;
		updateComponentResistance();
	}

	/** set R component */
	virtual void setComponentValue_R(DspFloatType _componentValue_R)
	{
		componentValue_R = _componentValue_R;
		updateComponentResistance();
	}

	/** set C component */
	virtual void setComponentValue_C(DspFloatType _componentValue_C)
	{
		componentValue_C = _componentValue_C;
		updateComponentResistance();
	}

	/** get R component value */
	virtual DspFloatType getComponentValue_R() { return componentValue_R; }

	/** get C component value */
	virtual DspFloatType getComponentValue_C() { return componentValue_C; }

	/** reset the component; clear registers */
	virtual void reset(DspFloatType _sampleRate) { setSampleRate(_sampleRate); zRegister_L = 0.0; zRegister_C = 0.0; }

	/** set input value into component */
	virtual void setInput(DspFloatType in){ zRegister_L = in; }

	/** get output value; NOTE: see FX book for details */
	virtual DspFloatType getOutput()
	{
		DspFloatType NL = zRegister_L;
		DspFloatType out = NL*(1.0 - K) + K*zRegister_C;
		zRegister_C = out;
		return out;
	}

	/** get output1 value; only one resistor output (not used) */
	virtual DspFloatType getOutput1() { return  getOutput(); }

	/** get output2 value; only one resistor output (not used) */
	virtual DspFloatType getOutput2() { return  getOutput(); }

	/** get output3 value; only one resistor output (not used) */
	virtual DspFloatType getOutput3() { return  getOutput(); }

	/** set input1 value; not used for components */
	virtual void setInput1(DspFloatType _in1) {}

	/** set input2 value; not used for components */
	virtual void setInput2(DspFloatType _in2) {}

	/** set input3 value; not used for components */
	virtual void setInput3(DspFloatType _in3) {}

protected:
	DspFloatType zRegister_L = 0.0; ///< storage register for L
	DspFloatType zRegister_C = 0.0; ///< storage register for C
	DspFloatType K = 0.0;

	DspFloatType componentValue_R = 0.0;///< component value R
	DspFloatType componentValue_C = 0.0;///< component value C

	DspFloatType RL = 0.0;	///< RL value
	DspFloatType RC = 0.0;	///< RC value
	DspFloatType RR = 0.0;	///< RR value

	DspFloatType componentResistance = 0.0; ///< equivalent resistance of pair of components
	DspFloatType sampleRate = 0.0; ///< sample rate
};

/**
\class WdfParallelRC
\ingroup WDF-Objects
\brief
The WdfParallelRC object implements the reflection coefficient and signal flow through
a WDF simulated parallal RC pair.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
SEE: p143 "Design of Wave Digital Filters" Psenicka, Ugalde, Romero M.
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfParallelRC : public IComponentAdaptor
{
public:
	WdfParallelRC() {}
	WdfParallelRC(DspFloatType _componentValue_R, DspFloatType _componentValue_C)
	{
		componentValue_C = _componentValue_C;
		componentValue_R = _componentValue_R;
	}
	virtual ~WdfParallelRC() {}

	/** set sample rate and update component */
	void setSampleRate(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		updateComponentResistance();
	}

	/** get component's value as a resistance */
	virtual DspFloatType getComponentResistance() { return componentResistance; }

	/** get component's value as a conducatance (or admittance) */
	virtual DspFloatType getComponentConductance() { return 1.0 / componentResistance; }

	/** change the resistance of component; see FX book for details */
	virtual void updateComponentResistance()
	{
		RR = componentValue_R;
		RC = 1.0 / (2.0*componentValue_C*sampleRate);
		componentResistance = 1.0 / ((1.0 / RR) + (1.0 / RC));
		K = componentResistance / RR;
	}

	/** set both RC components at once */
	virtual void setComponentValue_RC(DspFloatType _componentValue_R, DspFloatType _componentValue_C)
	{
		componentValue_R = _componentValue_R;
		componentValue_C = _componentValue_C;
		updateComponentResistance();
	}

	/** set R component */
	virtual void setComponentValue_R(DspFloatType _componentValue_R)
	{
		componentValue_R = _componentValue_R;
		updateComponentResistance();
	}

	/** set C component */
	virtual void setComponentValue_C(DspFloatType _componentValue_C)
	{
		componentValue_C = _componentValue_C;
		updateComponentResistance();
	}

	/** get R component value */
	virtual DspFloatType getComponentValue_R() { return componentValue_R; }

	/** get C component value */
	virtual DspFloatType getComponentValue_C() { return componentValue_C; }

	/** reset the component; clear registers */
	virtual void reset(DspFloatType _sampleRate) { setSampleRate(_sampleRate); zRegister_L = 0.0; zRegister_C = 0.0; }

	/** set input value into component; */
	virtual void setInput(DspFloatType in){ zRegister_L = in; }

	/** get output value; NOTE: output is located in zRegister_C */
	virtual DspFloatType getOutput()
	{
		DspFloatType NL = zRegister_L;
		DspFloatType out = NL*(1.0 - K) - K*zRegister_C;
		zRegister_C = out;
		return out;
	}

	/** get output1 value; only one resistor output (not used) */
	virtual DspFloatType getOutput1() { return  getOutput(); }

	/** get output2 value; only one resistor output (not used) */
	virtual DspFloatType getOutput2() { return  getOutput(); }

	/** get output3 value; only one resistor output (not used) */
	virtual DspFloatType getOutput3() { return  getOutput(); }

	/** set input1 value; not used for components */
	virtual void setInput1(DspFloatType _in1) {}

	/** set input2 value; not used for components */
	virtual void setInput2(DspFloatType _in2) {}

	/** set input3 value; not used for components */
	virtual void setInput3(DspFloatType _in3) {}

protected:
	DspFloatType zRegister_L = 0.0; ///< storage register for L
	DspFloatType zRegister_C = 0.0; ///< storage register for C
	DspFloatType K = 0.0;

	DspFloatType componentValue_C = 0.0;	///< component value C
	DspFloatType componentValue_R = 0.0;	///< component value R

	DspFloatType RL = 0.0; ///< RL value
	DspFloatType RC = 0.0; ///< RC value
	DspFloatType RR = 0.0; ///< RR value

	DspFloatType componentResistance = 0.0; ///< equivalent resistance of pair of components
	DspFloatType sampleRate = 0.0; ///< sample rate
};


// ------------------------------------------------------------------ //
// --- WDF ADAPTORS ------------------------------------------------- //
// ------------------------------------------------------------------ //

/**
\enum wdfComponent
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set the wdfComponent type
- enum class wdfComponent { R, L, C, seriesLC, parallelLC, seriesRL, parallelRL, seriesRC, parallelRC };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum class wdfComponent { R, L, C, seriesLC, parallelLC, seriesRL, parallelRL, seriesRC, parallelRC };

/**
\struct WdfComponentInfo
\ingroup WDF-Objects
\brief
Custom structure to hold component information.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct WdfComponentInfo
{
	WdfComponentInfo() { }

	WdfComponentInfo(wdfComponent _componentType, DspFloatType value1 = 0.0, DspFloatType value2 = 0.0)
	{
		componentType = _componentType;
		if (componentType == wdfComponent::R)
			R = value1;
		else if (componentType == wdfComponent::L)
			L = value1;
		else if (componentType == wdfComponent::C)
			C = value1;
		else if (componentType == wdfComponent::seriesLC || componentType == wdfComponent::parallelLC)
		{
			L = value1;
			C = value2;
		}
		else if (componentType == wdfComponent::seriesRL || componentType == wdfComponent::parallelRL)
		{
			R = value1;
			L = value2;
		}
		else if (componentType == wdfComponent::seriesRC || componentType == wdfComponent::parallelRC)
		{
			R = value1;
			C = value2;
		}
	}

	DspFloatType R = 0.0; ///< value of R component
	DspFloatType L = 0.0;	///< value of L component
	DspFloatType C = 0.0;	///< value of C component
	wdfComponent componentType = wdfComponent::R; ///< selected component type
};


/**
\class WdfAdaptorBase
\ingroup WDF-Objects
\brief
The WdfAdaptorBase object acts as the base class for all WDF Adaptors; the static members allow
for simplified connection of components. See the FX book for more details.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfAdaptorBase : public IComponentAdaptor
{
public:
	WdfAdaptorBase() {}
	virtual ~WdfAdaptorBase() {}

	/** set the termainal (load) resistance for terminating adaptors */
	void setTerminalResistance(DspFloatType _terminalResistance) { terminalResistance = _terminalResistance; }

	/** set the termainal (load) resistance as open circuit for terminating adaptors */
	void setOpenTerminalResistance(bool _openTerminalResistance = true)
	{
		// --- flag overrides value
		openTerminalResistance = _openTerminalResistance;
		terminalResistance = 1.0e+34; // avoid /0.0
	}

	/** set the input (source) resistance for an input adaptor */
	void setSourceResistance(DspFloatType _sourceResistance) { sourceResistance = _sourceResistance; }

	/** set the component or connected adaptor at port 1; functions is generic and allows extending the functionality of the WDF Library */
	void setPort1_CompAdaptor(IComponentAdaptor* _port1CompAdaptor) { port1CompAdaptor = _port1CompAdaptor; }

	/** set the component or connected adaptor at port 2; functions is generic and allows extending the functionality of the WDF Library */
	void setPort2_CompAdaptor(IComponentAdaptor* _port2CompAdaptor) { port2CompAdaptor = _port2CompAdaptor; }

	/** set the component or connected adaptor at port 3; functions is generic and allows extending the functionality of the WDF Library */
	void setPort3_CompAdaptor(IComponentAdaptor* _port3CompAdaptor) { port3CompAdaptor = _port3CompAdaptor; }

	/** reset the connected component */
	virtual void reset(DspFloatType _sampleRate)
	{
		if (wdfComponent)
			wdfComponent->reset(_sampleRate);
	}

	/** creates a new WDF component and connects it to Port 3 */
	void setComponent(wdfComponent componentType, DspFloatType value1 = 0.0, DspFloatType value2 = 0.0)
	{
		// --- decode and set
		if (componentType == wdfComponent::R)
		{
			wdfComponent = new WdfResistor;
			wdfComponent->setComponentValue(value1);
			port3CompAdaptor = wdfComponent;
		}
		else if (componentType == wdfComponent::L)
		{
			wdfComponent = new WdfInductor;
			wdfComponent->setComponentValue(value1);
			port3CompAdaptor = wdfComponent;
		}
		else if (componentType == wdfComponent::C)
		{
			wdfComponent = new WdfCapacitor;
			wdfComponent->setComponentValue(value1);
			port3CompAdaptor = wdfComponent;
		}
		else if (componentType == wdfComponent::seriesLC)
		{
			wdfComponent = new WdfSeriesLC;
			wdfComponent->setComponentValue_LC(value1, value2);
			port3CompAdaptor = wdfComponent;
		}
		else if (componentType == wdfComponent::parallelLC)
		{
			wdfComponent = new WdfParallelLC;
			wdfComponent->setComponentValue_LC(value1, value2);
			port3CompAdaptor = wdfComponent;
		}
		else if (componentType == wdfComponent::seriesRL)
		{
			wdfComponent = new WdfSeriesRL;
			wdfComponent->setComponentValue_RL(value1, value2);
			port3CompAdaptor = wdfComponent;
		}
		else if (componentType == wdfComponent::parallelRL)
		{
			wdfComponent = new WdfParallelRL;
			wdfComponent->setComponentValue_RL(value1, value2);
			port3CompAdaptor = wdfComponent;
		}
		else if (componentType == wdfComponent::seriesRC)
		{
			wdfComponent = new WdfSeriesRC;
			wdfComponent->setComponentValue_RC(value1, value2);
			port3CompAdaptor = wdfComponent;
		}
		else if (componentType == wdfComponent::parallelRC)
		{
			wdfComponent = new WdfParallelRC;
			wdfComponent->setComponentValue_RC(value1, value2);
			port3CompAdaptor = wdfComponent;
		}
	}

	/** connect two adapters together upstreamAdaptor --> downstreamAdaptor */
	static void connectAdaptors(WdfAdaptorBase* upstreamAdaptor, WdfAdaptorBase* downstreamAdaptor)
	{
		upstreamAdaptor->setPort2_CompAdaptor(downstreamAdaptor);
		downstreamAdaptor->setPort1_CompAdaptor(upstreamAdaptor);
	}

	/** initialize the chain of adaptors from upstreamAdaptor --> downstreamAdaptor */
	virtual void initializeAdaptorChain()
	{
		initialize(sourceResistance);
	}

	/** set value of single-component adaptor */
	virtual void setComponentValue(DspFloatType _componentValue)
	{
		if (wdfComponent)
			wdfComponent->setComponentValue(_componentValue);
	}

	/** set LC value of mjulti-component adaptor */
	virtual void setComponentValue_LC(DspFloatType componentValue_L, DspFloatType componentValue_C)
	{
		if (wdfComponent)
			wdfComponent->setComponentValue_LC(componentValue_L, componentValue_C);
	}

	/** set RL value of mjulti-component adaptor */
	virtual void setComponentValue_RL(DspFloatType componentValue_R, DspFloatType componentValue_L)
	{
		if (wdfComponent)
			wdfComponent->setComponentValue_RL(componentValue_R, componentValue_L);
	}

	/** set RC value of mjulti-component adaptor */
	virtual void setComponentValue_RC(DspFloatType componentValue_R, DspFloatType componentValue_C)
	{
		if (wdfComponent)
			wdfComponent->setComponentValue_RC(componentValue_R, componentValue_C);
	}

	/** get adaptor connected at port 1: for extended functionality; not used in WDF ladder filter library */
	IComponentAdaptor* getPort1_CompAdaptor() { return port1CompAdaptor; }

	/** get adaptor connected at port 2: for extended functionality; not used in WDF ladder filter library */
	IComponentAdaptor* getPort2_CompAdaptor() { return port2CompAdaptor; }

	/** get adaptor connected at port 3: for extended functionality; not used in WDF ladder filter library */
	IComponentAdaptor* getPort3_CompAdaptor() { return port3CompAdaptor; }

protected:
	// --- can in theory connect any port to a component OR adaptor;
	//     though this library is setup with a convention R3 = component
	IComponentAdaptor* port1CompAdaptor = nullptr;	///< componant or adaptor connected to port 1
	IComponentAdaptor* port2CompAdaptor = nullptr;	///< componant or adaptor connected to port 2
	IComponentAdaptor* port3CompAdaptor = nullptr;	///< componant or adaptor connected to port 3
	IComponentAdaptor* wdfComponent = nullptr;		///< WDF componant connected to port 3 (default operation)

	// --- These hold the input (R1), component (R3) and output (R2) resistances
	DspFloatType R1 = 0.0; ///< input port resistance
	DspFloatType R2 = 0.0; ///< output port resistance
	DspFloatType R3 = 0.0; ///< component resistance

	// --- these are input variables that are stored;
	//     not used in this implementation but may be required for extended versions
	DspFloatType in1 = 0.0;	///< stored port 1 input;  not used in this implementation but may be required for extended versions
	DspFloatType in2 = 0.0;	///< stored port 2 input;  not used in this implementation but may be required for extended versions
	DspFloatType in3 = 0.0;	///< stored port 3 input;  not used in this implementation but may be required for extended versions

	// --- these are output variables that are stored;
	//     currently out2 is the only one used as it is y(n) for this library
	//     out1 and out2 are stored; not used in this implementation but may be required for extended versions
	DspFloatType out1 = 0.0;	///< stored port 1 output; not used in this implementation but may be required for extended versions
	DspFloatType out2 = 0.0;	///< stored port 2 output; it is y(n) for this library
	DspFloatType out3 = 0.0;	///< stored port 3 output; not used in this implementation but may be required for extended versions

	// --- terminal impedance
	DspFloatType terminalResistance = 600.0; ///< value of terminal (load) resistance
	bool openTerminalResistance = false; ///< flag for open circuit load

	// --- source impedance, OK for this to be set to 0.0 for Rs = 0
	DspFloatType sourceResistance = 600.0; ///< source impedance; OK for this to be set to 0.0 for Rs = 0
};

/**
\class WdfSeriesAdaptor
\ingroup WDF-Objects
\brief
The WdfSeriesAdaptor object implements the series reflection-free (non-terminated) adaptor
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfSeriesAdaptor : public WdfAdaptorBase
{
public:
	WdfSeriesAdaptor() {}
	virtual ~WdfSeriesAdaptor() {}

	/** get the resistance at port 2; R2 = R1 + component (series)*/
	virtual DspFloatType getR2()
	{
		DspFloatType componentResistance = 0.0;
		if (getPort3_CompAdaptor())
			componentResistance = getPort3_CompAdaptor()->getComponentResistance();

		R2 = R1 + componentResistance;
		return R2;
	}

	/** initialize adaptor with input resistance */
	virtual void initialize(DspFloatType _R1)
	{
		// --- R1 is source resistance for this adaptor
		R1 = _R1;

		DspFloatType componentResistance = 0.0;
		if (getPort3_CompAdaptor())
			componentResistance = getPort3_CompAdaptor()->getComponentResistance();

		// --- calculate B coeff
		B = R1 / (R1 + componentResistance);

		// --- init downstream adaptor
		if (getPort2_CompAdaptor())
			getPort2_CompAdaptor()->initialize(getR2());

		// --- not used in this implementation but saving for extended use
		R3 = componentResistance;
	}

	/** push audio input sample into incident wave input*/
	virtual void setInput1(DspFloatType _in1)
	{
		// --- save
		in1 = _in1;

		// --- read component value
		N2 = 0.0;
		if (getPort3_CompAdaptor())
			N2 = getPort3_CompAdaptor()->getOutput();

		// --- form output
		out2 = -(in1 + N2);

		// --- deliver downstream
		if (getPort2_CompAdaptor())
			getPort2_CompAdaptor()->setInput1(out2);
	}

	/** push audio input sample into reflected wave input */
	virtual void setInput2(DspFloatType _in2)
	{
		// --- save
		in2 = _in2;

		// --- calc N1
		N1 = -(in1 - B*(in1 + N2 + in2) + in2);

		// --- calc out1
		out1 = in1 - B*(N2 + in2);

		// --- deliver upstream
		if (getPort1_CompAdaptor())
			getPort1_CompAdaptor()->setInput2(out1);

		// --- set component state
		if (getPort3_CompAdaptor())
			getPort3_CompAdaptor()->setInput(N1);
	}

	/** set input 3 always connects to component */
	virtual void setInput3(DspFloatType _in3){ }

	/** get OUT1 = reflected output pin on Port 1 */
	virtual DspFloatType getOutput1() { return out1; }

	/** get OUT2 = incident (normal) output pin on Port 2 */
	virtual DspFloatType getOutput2() { return out2; }

	/** get OUT3 always connects to component */
	virtual DspFloatType getOutput3() { return out3; }

private:
	DspFloatType N1 = 0.0;	///< node 1 value, internal use only
	DspFloatType N2 = 0.0;	///< node 2 value, internal use only
	DspFloatType B = 0.0;		///< B coefficient value
};

/**
\class WdfSeriesTerminatedAdaptor
\ingroup WDF-Objects
\brief
The WdfSeriesTerminatedAdaptor object implements the series terminated (non-reflection-free) adaptor
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
// --- Series terminated adaptor
class WdfSeriesTerminatedAdaptor : public WdfAdaptorBase
{
public:
	WdfSeriesTerminatedAdaptor() {}
	virtual ~WdfSeriesTerminatedAdaptor() {}

	/** get the resistance at port 2; R2 = R1 + component (series)*/
	virtual DspFloatType getR2()
	{
		DspFloatType componentResistance = 0.0;
		if (getPort3_CompAdaptor())
			componentResistance = getPort3_CompAdaptor()->getComponentResistance();

		R2 = R1 + componentResistance;
		return R2;
	}

	/** initialize adaptor with input resistance */
	virtual void initialize(DspFloatType _R1)
	{
		// --- source impedance
		R1 = _R1;

		DspFloatType componentResistance = 0.0;
		if (getPort3_CompAdaptor())
			componentResistance = getPort3_CompAdaptor()->getComponentResistance();

		B1 = (2.0*R1) / (R1 + componentResistance + terminalResistance);
		B3 = (2.0*terminalResistance) / (R1 + componentResistance + terminalResistance);

		// --- init downstream
		if (getPort2_CompAdaptor())
			getPort2_CompAdaptor()->initialize(getR2());

		// --- not used in this implementation but saving for extended use
		R3 = componentResistance;
	}

	/** push audio input sample into incident wave input*/
	virtual void setInput1(DspFloatType _in1)
	{
		// --- save
		in1 = _in1;

		N2 = 0.0;
		if (getPort3_CompAdaptor())
			N2 = getPort3_CompAdaptor()->getOutput();

		DspFloatType N3 = in1 + N2;

		// --- calc out2 y(n)
		out2 = -B3*N3;

		// --- form output1
		out1 = in1 - B1*N3;

		// --- form N1
		N1 = -(out1 + out2 + N3);

		// --- deliver upstream to input2
		if (getPort1_CompAdaptor())
			getPort1_CompAdaptor()->setInput2(out1);

		// --- set component state
		if (getPort3_CompAdaptor())
			getPort3_CompAdaptor()->setInput(N1);
	}

	/** push audio input sample into reflected wave input
	    for terminated adaptor, this is dead end, just store it */
	virtual void setInput2(DspFloatType _in2) { in2 = _in2;}

	/** set input 3 always connects to component */
	virtual void setInput3(DspFloatType _in3) { in3 = _in3;}

	/** get OUT1 = reflected output pin on Port 1 */
	virtual DspFloatType getOutput1() { return out1; }

	/** get OUT2 = incident (normal) output pin on Port 2 */
	virtual DspFloatType getOutput2() { return out2; }

	/** get OUT3 always connects to component */
	virtual DspFloatType getOutput3() { return out3; }

private:
	DspFloatType N1 = 0.0;	///< node 1 value, internal use only
	DspFloatType N2 = 0.0;	///< node 2 value, internal use only
	DspFloatType B1 = 0.0;	///< B1 coefficient value
	DspFloatType B3 = 0.0;	///< B3 coefficient value
};

/**
\class WdfParallelAdaptor
\ingroup WDF-Objects
\brief
The WdfParallelAdaptor object implements the parallel reflection-free (non-terminated) adaptor
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfParallelAdaptor : public WdfAdaptorBase
{
public:
	WdfParallelAdaptor() {}
	virtual ~WdfParallelAdaptor() {}

	/** get the resistance at port 2;  R2 = 1.0/(sum of admittances) */
	virtual DspFloatType getR2()
	{
		DspFloatType componentConductance = 0.0;
		if (getPort3_CompAdaptor())
			componentConductance = getPort3_CompAdaptor()->getComponentConductance();

		// --- 1 / (sum of admittances)
		R2 = 1.0 / ((1.0 / R1) + componentConductance);
		return R2;
	}

	/** initialize adaptor with input resistance */
	virtual void initialize(DspFloatType _R1)
	{
		// --- save R1
		R1 = _R1;

		DspFloatType G1 = 1.0 / R1;
		DspFloatType componentConductance = 0.0;
		if (getPort3_CompAdaptor())
			componentConductance = getPort3_CompAdaptor()->getComponentConductance();

		// --- calculate B coeff
		A = G1 / (G1 + componentConductance);

		// --- now, do we init our downstream??
		if (getPort2_CompAdaptor())
			getPort2_CompAdaptor()->initialize(getR2());

		// --- not used in this implementation but saving for extended use
		R3 = 1.0/ componentConductance;
	}

	/** push audio input sample into incident wave input*/
	virtual void setInput1(DspFloatType _in1)
	{
		// --- save
		in1 = _in1;

		// --- read component
		N2 = 0.0;
		if (getPort3_CompAdaptor())
			N2 = getPort3_CompAdaptor()->getOutput();

		// --- form output
		out2 = N2 - A*(-in1 + N2);

		// --- deliver downstream
		if (getPort2_CompAdaptor())
			getPort2_CompAdaptor()->setInput1(out2);
	}

	/** push audio input sample into reflected wave input*/
	virtual void setInput2(DspFloatType _in2)
	{
		// --- save
		in2 = _in2;

		// --- calc N1
		N1 = in2 - A*(-in1 + N2);

		// --- calc out1
		out1 = -in1 + N2 + N1;

		// --- deliver upstream
		if (getPort1_CompAdaptor())
			getPort1_CompAdaptor()->setInput2(out1);

		// --- set component state
		if (getPort3_CompAdaptor())
			getPort3_CompAdaptor()->setInput(N1);
	}

	/** set input 3 always connects to component */
	virtual void setInput3(DspFloatType _in3) { }

	/** get OUT1 = reflected output pin on Port 1 */
	virtual DspFloatType getOutput1() { return out1; }

	/** get OUT2 = incident (normal) output pin on Port 2 */
	virtual DspFloatType getOutput2() { return out2; }

	/** get OUT3 always connects to component */
	virtual DspFloatType getOutput3() { return out3; }

private:
	DspFloatType N1 = 0.0;	///< node 1 value, internal use only
	DspFloatType N2 = 0.0;	///< node 2 value, internal use only
	DspFloatType A = 0.0;		///< A coefficient value
};


/**
\class WdfParallelTerminatedAdaptor
\ingroup WDF-Objects
\brief
The WdfParallelTerminatedAdaptor object implements the parallel  terminated (non-reflection-free) adaptor
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WdfParallelTerminatedAdaptor : public WdfAdaptorBase
{
public:
	WdfParallelTerminatedAdaptor() {}
	virtual ~WdfParallelTerminatedAdaptor() {}

	/** get the resistance at port 2;  R2 = 1.0/(sum of admittances) */
	virtual DspFloatType getR2()
	{
		DspFloatType componentConductance = 0.0;
		if (getPort3_CompAdaptor())
			componentConductance = getPort3_CompAdaptor()->getComponentConductance();

		// --- 1 / (sum of admittances)
		R2 = 1.0 / ((1.0 / R1) + componentConductance);
		return R2;
	}

	/** initialize adaptor with input resistance */
	virtual void initialize(DspFloatType _R1)
	{
		// --- save R1
		R1 = _R1;

		DspFloatType G1 = 1.0 / R1;
		if (terminalResistance <= 0.0)
			terminalResistance = 1e-15;

		DspFloatType G2 = 1.0 / terminalResistance;
		DspFloatType componentConductance = 0.0;
		if (getPort3_CompAdaptor())
			componentConductance = getPort3_CompAdaptor()->getComponentConductance();

		A1 = 2.0*G1 / (G1 + componentConductance + G2);
		A3 = openTerminalResistance ? 0.0 : 2.0*G2 / (G1 + componentConductance + G2);

		// --- init downstream
		if (getPort2_CompAdaptor())
			getPort2_CompAdaptor()->initialize(getR2());

		// --- not used in this implementation but saving for extended use
		R3 = 1.0 / componentConductance;
	}

	/** push audio input sample into incident wave input*/
	virtual void setInput1(DspFloatType _in1)
	{
		// --- save
		in1 = _in1;

		N2 = 0.0;
		if (getPort3_CompAdaptor())
			N2 = getPort3_CompAdaptor()->getOutput();

		// --- form N1
		N1 = -A1*(-in1 + N2) + N2 - A3*N2;

		// --- form output1
		out1 = -in1 + N2 + N1;

		// --- deliver upstream to input2
		if (getPort1_CompAdaptor())
			getPort1_CompAdaptor()->setInput2(out1);

		// --- calc out2 y(n)
		out2 = N2 + N1;

		// --- set component state
		if (getPort3_CompAdaptor())
			getPort3_CompAdaptor()->setInput(N1);
	}

	/** push audio input sample into reflected wave input; this is a dead end for terminated adaptorsthis is a dead end for terminated adaptors  */
	virtual void setInput2(DspFloatType _in2){ in2 = _in2;}

	/** set input 3 always connects to component */
	virtual void setInput3(DspFloatType _in3) { }

	/** get OUT1 = reflected output pin on Port 1 */
	virtual DspFloatType getOutput1() { return out1; }

	/** get OUT2 = incident (normal) output pin on Port 2 */
	virtual DspFloatType getOutput2() { return out2; }

	/** get OUT3 always connects to component */
	virtual DspFloatType getOutput3() { return out3; }

private:
	DspFloatType N1 = 0.0;	///< node 1 value, internal use only
	DspFloatType N2 = 0.0;	///< node 2 value, internal use only
	DspFloatType A1 = 0.0;	///< A1 coefficient value
	DspFloatType A3 = 0.0;	///< A3 coefficient value
};

// ------------------------------------------------------------------------------ //
// --- WDF Ladder Filter Design  Examples --------------------------------------- //
// ------------------------------------------------------------------------------ //
//
// --- 3rd order Butterworth LPF designed with Elsie www.TonneSoftware.comm
//
/*
	3rd Order Inductor-Leading LPF
	Rs = Rload = 600 ohms
	Series(L1) -> Parallel(C1) -> Series(L2)
	--L1-- | --L2--
		   C1
		   |
	fc = 1kHz
		L1 = 95.49e-3;
		C1 = 0.5305e-6;
		L2 = 95.49e-3;
	fc = 10kHz
		L1 = 9.549e-3;
		C1 = 0.05305e-6;
		L2 = 9.549e-3;
*/

/**
\class WDFButterLPF3
\ingroup WDF-Objects
\brief
The WDFButterLPF3 object implements a 3rd order Butterworth ladder filter.
NOTE: designed with Elsie www.TonneSoftware.comm
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- none - this object is hard-wired.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WDFButterLPF3 : public IAudioSignalProcessor
{
public:
	WDFButterLPF3(void) { createWDF(); }	/* C-TOR */
	~WDFButterLPF3(void) {}	/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		// --- rest WDF components (flush state registers)
		seriesAdaptor_L1.reset(_sampleRate);
		parallelAdaptor_C1.reset(_sampleRate);
		seriesTerminatedAdaptor_L2.reset(_sampleRate);

		// --- intialize the chain of adapters
		seriesAdaptor_L1.initializeAdaptorChain();
		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the WDF ladder filter to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- push audio sample into series L1
		seriesAdaptor_L1.setInput1(xn);

		// --- output is at terminated L2's output2
		return seriesTerminatedAdaptor_L2.getOutput2();
	}

	/** create the WDF structure for this object - may be called more than once */
	void createWDF()
	{
		// --- actual component values fc = 1kHz
		DspFloatType L1_value = 95.49e-3;		// 95.5 mH
		DspFloatType C1_value = 0.5305e-6;	// 0.53 uF
		DspFloatType L2_value = 95.49e-3;		// 95.5 mH

										// --- set adapter components
		seriesAdaptor_L1.setComponent(wdfComponent::L, L1_value);
		parallelAdaptor_C1.setComponent(wdfComponent::C, C1_value);
		seriesTerminatedAdaptor_L2.setComponent(wdfComponent::L, L2_value);

		// --- connect adapters
		WdfAdaptorBase::connectAdaptors(&seriesAdaptor_L1, &parallelAdaptor_C1);
		WdfAdaptorBase::connectAdaptors(&parallelAdaptor_C1, &seriesTerminatedAdaptor_L2);

		// --- set source resistance
		seriesAdaptor_L1.setSourceResistance(600.0); // --- Rs = 600

		// --- set terminal resistance
		seriesTerminatedAdaptor_L2.setTerminalResistance(600.0); // --- Rload = 600
	}

protected:
	// --- three adapters
	WdfSeriesAdaptor seriesAdaptor_L1;			///< adaptor for L1
	WdfParallelAdaptor parallelAdaptor_C1;		///< adaptor for C1
	WdfSeriesTerminatedAdaptor seriesTerminatedAdaptor_L2;	///< adaptor for L2
};

/**
\class WDFTunableButterLPF3
\ingroup WDF-Objects
\brief
The WDFTunableButterLPF3 object implements a tunable 3rd order Butterworth ladder filter.
NOTE: designed with Elsie www.TonneSoftware.comm
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- setUsePostWarping(bool b) to enable/disable warping (see book)
- setFilterFc(DspFloatType fc_Hz) to set the tunable fc value
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WDFTunableButterLPF3 : public IAudioSignalProcessor
{
public:
	WDFTunableButterLPF3(void) { createWDF(); }	/* C-TOR */
	~WDFTunableButterLPF3(void) {}	/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		// --- rest WDF components (flush state registers)
		seriesAdaptor_L1.reset(_sampleRate);
		parallelAdaptor_C1.reset(_sampleRate);
		seriesTerminatedAdaptor_L2.reset(_sampleRate);

		// --- intialize the chain of adapters
		seriesAdaptor_L1.initializeAdaptorChain();
		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the WDF ladder filter to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- push audio sample into series L1
		seriesAdaptor_L1.setInput1(xn);

		// --- output is at terminated L2's output2
		return seriesTerminatedAdaptor_L2.getOutput2();
	}

	/** create the filter structure; may be called more than once */
	void createWDF()
	{
		// --- create components, init to noramlized values fc = 1Hz
		seriesAdaptor_L1.setComponent(wdfComponent::L, L1_norm);
		parallelAdaptor_C1.setComponent(wdfComponent::C, C1_norm);
		seriesTerminatedAdaptor_L2.setComponent(wdfComponent::L, L2_norm);

		// --- connect adapters
		WdfAdaptorBase::connectAdaptors(&seriesAdaptor_L1, &parallelAdaptor_C1);
		WdfAdaptorBase::connectAdaptors(&parallelAdaptor_C1, &seriesTerminatedAdaptor_L2);

		// --- set source resistance
		seriesAdaptor_L1.setSourceResistance(600.0); // --- Rs = 600

		// --- set terminal resistance
		seriesTerminatedAdaptor_L2.setTerminalResistance(600.0); // --- Rload = 600
	}

	/** parameter setter for warping */
	void setUsePostWarping(bool b) { useFrequencyWarping = b; }

	/** parameter setter for fc */
	void setFilterFc(DspFloatType fc_Hz)
	{
		if (useFrequencyWarping)
		{
			DspFloatType arg = (kPi*fc_Hz) / sampleRate;
			fc_Hz = fc_Hz*(tan(arg) / arg);
		}

		seriesAdaptor_L1.setComponentValue(L1_norm / fc_Hz);
		parallelAdaptor_C1.setComponentValue(C1_norm / fc_Hz);
		seriesTerminatedAdaptor_L2.setComponentValue(L2_norm / fc_Hz);
	}

protected:
	// --- three adapters
	WdfSeriesAdaptor seriesAdaptor_L1;		///< adaptor for L1
	WdfParallelAdaptor parallelAdaptor_C1;	///< adaptor for C1
	WdfSeriesTerminatedAdaptor seriesTerminatedAdaptor_L2;	///< adaptor for L2

	DspFloatType L1_norm = 95.493;		// 95.5 mH
	DspFloatType C1_norm = 530.516e-6;	// 0.53 uF
	DspFloatType L2_norm = 95.493;		// 95.5 mH

	bool useFrequencyWarping = false;	///< flag for freq warping
	DspFloatType sampleRate = 1.0;			///< stored sample rate
};

/**
\class WDFBesselBSF3
\ingroup WDF-Objects
\brief
The WDFBesselBSF3 object implements a 3rd order Bessel BSF
NOTE: designed with Elsie www.TonneSoftware.comm
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- none - object is hardwired
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WDFBesselBSF3 : public IAudioSignalProcessor
{
public:
	WDFBesselBSF3(void) { createWDF(); }	/* C-TOR */
	~WDFBesselBSF3(void) {}	/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		// --- rest WDF components (flush state registers)
		seriesAdaptor_L1C1.reset(_sampleRate);
		parallelAdaptor_L2C2.reset(_sampleRate);
		seriesTerminatedAdaptor_L3C3.reset(_sampleRate);

		// --- intialize the chain of adapters
		seriesAdaptor_L1C1.initializeAdaptorChain();

		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the WDF ladder filter to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- push audio sample into series L1
		seriesAdaptor_L1C1.setInput1(xn);

		// --- output is at terminated L2's output2
		return seriesTerminatedAdaptor_L3C3.getOutput2();
	}

	/** create the WDF structure; may be called more than once*/
	void createWDF()
	{
		// --- set component values
		// --- fo = 5kHz
		//     BW = 2kHz or Q = 2.5
		seriesAdaptor_L1C1.setComponent(wdfComponent::parallelLC, 16.8327e-3, 0.060193e-6);	/* L, C */
		parallelAdaptor_L2C2.setComponent(wdfComponent::seriesLC, 49.1978e-3, 0.02059e-6);	/* L, C */
		seriesTerminatedAdaptor_L3C3.setComponent(wdfComponent::parallelLC, 2.57755e-3, 0.393092e-6);	/* L, C */

		// --- connect adapters
		WdfAdaptorBase::connectAdaptors(&seriesAdaptor_L1C1, &parallelAdaptor_L2C2);
		WdfAdaptorBase::connectAdaptors(&parallelAdaptor_L2C2, &seriesTerminatedAdaptor_L3C3);

		// --- set source resistance
		seriesAdaptor_L1C1.setSourceResistance(600.0); // Ro = 600

		// --- set terminal resistance
		seriesTerminatedAdaptor_L3C3.setTerminalResistance(600.0);
	}

protected:
	// --- three adapters
	WdfSeriesAdaptor seriesAdaptor_L1C1;		///< adaptor for L1 and C1
	WdfParallelAdaptor parallelAdaptor_L2C2;	///< adaptor for L2 and C2
	WdfSeriesTerminatedAdaptor seriesTerminatedAdaptor_L3C3;	///< adaptor for L3 and C3
};


/**
\class WDFConstKBPF6
\ingroup WDF-Objects
\brief
The WDFConstKBPF6 object implements a 6th order constant K BPF
NOTE: designed with Elsie www.TonneSoftware.comm
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- none - object is hardwired
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WDFConstKBPF6 : public IAudioSignalProcessor
{
public:
	WDFConstKBPF6(void) { createWDF(); }	/* C-TOR */
	~WDFConstKBPF6(void) {}	/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		// --- rest WDF components (flush state registers)
		seriesAdaptor_L1C1.reset(_sampleRate);
		parallelAdaptor_L2C2.reset(_sampleRate);

		seriesAdaptor_L3C3.reset(_sampleRate);
		parallelAdaptor_L4C4.reset(_sampleRate);

		seriesAdaptor_L5C5.reset(_sampleRate);
		parallelTerminatedAdaptor_L6C6.reset(_sampleRate);

		// --- intialize the chain of adapters
		seriesAdaptor_L1C1.initializeAdaptorChain();
		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the WDF ladder filter to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- push audio sample into series L1
		seriesAdaptor_L1C1.setInput1(xn);

		// --- output is at terminated L6C6 output2
		DspFloatType output = parallelTerminatedAdaptor_L6C6.getOutput2();

		return output;
	}

	/** create the WDF structure */
	void createWDF()
	{
		// --- fo = 5kHz
		//     BW = 2kHz or Q = 2.5
		seriesAdaptor_L1C1.setComponent(wdfComponent::seriesLC, 47.7465e-3, 0.02122e-6);
		parallelAdaptor_L2C2.setComponent(wdfComponent::parallelLC, 3.81972e-3, 0.265258e-6);

		seriesAdaptor_L3C3.setComponent(wdfComponent::seriesLC, 95.493e-3, 0.01061e-6);
		parallelAdaptor_L4C4.setComponent(wdfComponent::parallelLC, 3.81972e-3, 0.265258e-6);

		seriesAdaptor_L5C5.setComponent(wdfComponent::seriesLC, 95.493e-3, 0.01061e-6);
		parallelTerminatedAdaptor_L6C6.setComponent(wdfComponent::parallelLC, 7.63944e-3, 0.132629e-6);

		// --- connect adapters
		WdfAdaptorBase::connectAdaptors(&seriesAdaptor_L1C1, &parallelAdaptor_L2C2);
		WdfAdaptorBase::connectAdaptors(&parallelAdaptor_L2C2, &seriesAdaptor_L3C3);
		WdfAdaptorBase::connectAdaptors(&seriesAdaptor_L3C3, &parallelAdaptor_L4C4);
		WdfAdaptorBase::connectAdaptors(&parallelAdaptor_L4C4, &seriesAdaptor_L5C5);
		WdfAdaptorBase::connectAdaptors(&seriesAdaptor_L5C5, &parallelTerminatedAdaptor_L6C6);

		// --- set source resistance
		seriesAdaptor_L1C1.setSourceResistance(600.0); // Ro = 600

		// --- set terminal resistance
		parallelTerminatedAdaptor_L6C6.setTerminalResistance(600.0);
	}

protected:
	// --- six adapters
	WdfSeriesAdaptor seriesAdaptor_L1C1;		///< adaptor for L1 and C1
	WdfParallelAdaptor parallelAdaptor_L2C2;	///< adaptor for L2 and C2

	WdfSeriesAdaptor seriesAdaptor_L3C3;		///< adaptor for L3 and C3
	WdfParallelAdaptor parallelAdaptor_L4C4;	///< adaptor for L4 and C4

	WdfSeriesAdaptor seriesAdaptor_L5C5;		///< adaptor for L5 and C5
	WdfParallelTerminatedAdaptor parallelTerminatedAdaptor_L6C6;///< adaptor for L6 and C6
};


/**
\struct WDFParameters
\ingroup WDF-Objects
\brief
Custom parameter structure for the WDF filter examples.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct WDFParameters
{
	WDFParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	WDFParameters& operator=(const WDFParameters& params)
	{
		if (this == &params)
			return *this;

		fc = params.fc;
		Q = params.Q;
		boostCut_dB = params.boostCut_dB;
		frequencyWarping = params.frequencyWarping;
		return *this;
	}

	// --- individual parameters
	DspFloatType fc = 100.0;				///< filter fc
	DspFloatType Q = 0.707;				///< filter Q
	DspFloatType boostCut_dB = 0.0;		///< filter boost or cut in dB
	bool frequencyWarping = true;	///< enable frequency warping
};

/**
\class WDFIdealRLCLPF
\ingroup WDF-Objects
\brief
The WDFIdealRLCLPF object implements an ideal RLC LPF using the WDF library.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use WDFParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WDFIdealRLCLPF : public IAudioSignalProcessor
{
public:
	WDFIdealRLCLPF(void) { createWDF(); }	/* C-TOR */
	~WDFIdealRLCLPF(void) {}	/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;

		// --- rest WDF components (flush state registers)
		seriesAdaptor_RL.reset(_sampleRate);
		parallelTerminatedAdaptor_C.reset(_sampleRate);

		// --- intialize the chain of adapters
		seriesAdaptor_RL.initializeAdaptorChain();
		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the WDF Ideal RLC filter to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- push audio sample into series L1
		seriesAdaptor_RL.setInput1(xn);

		// --- output is at terminated L2's output2
		//     note compensation scaling by -6dB = 0.5
		//     because of WDF assumption about Rs and Rload
		return 0.5*parallelTerminatedAdaptor_C.getOutput2();
	}

	/** create the WDF structure; may be called more than once */
	void createWDF()
	{
		// --- create components, init to noramlized values fc =
		//	   initial values for fc = 1kHz Q = 0.707
		//     Holding C Constant at 1e-6
		//			   L = 2.533e-2
		//			   R = 2.251131 e2
		seriesAdaptor_RL.setComponent(wdfComponent::seriesRL, 2.251131e2, 2.533e-2);
		parallelTerminatedAdaptor_C.setComponent(wdfComponent::C, 1.0e-6);

		// --- connect adapters
		WdfAdaptorBase::connectAdaptors(&seriesAdaptor_RL, &parallelTerminatedAdaptor_C);

		// --- set source resistance
		seriesAdaptor_RL.setSourceResistance(0.0); // --- Rs = 600

		// --- set open ckt termination
		parallelTerminatedAdaptor_C.setOpenTerminalResistance(true);
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return WDFParameters custom data structure
	*/
	WDFParameters getParameters() { return wdfParameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param WDFParameters custom data structure
	*/
	void setParameters(const WDFParameters& _wdfParameters)
	{
		if (_wdfParameters.fc != wdfParameters.fc ||
			_wdfParameters.Q != wdfParameters.Q ||
			_wdfParameters.boostCut_dB != wdfParameters.boostCut_dB ||
			_wdfParameters.frequencyWarping != wdfParameters.frequencyWarping)
		{
			wdfParameters = _wdfParameters;
			DspFloatType fc_Hz = wdfParameters.fc;

			if (wdfParameters.frequencyWarping)
			{
				DspFloatType arg = (kPi*fc_Hz) / sampleRate;
				fc_Hz = fc_Hz*(tan(arg) / arg);
			}

			DspFloatType inductorValue = 1.0 / (1.0e-6 * pow((2.0*kPi*fc_Hz), 2.0));
			DspFloatType resistorValue = (1.0 / wdfParameters.Q)*(pow(inductorValue / 1.0e-6, 0.5));

			seriesAdaptor_RL.setComponentValue_RL(resistorValue, inductorValue);
			seriesAdaptor_RL.initializeAdaptorChain();
		}
	}

protected:
	WDFParameters wdfParameters;	///< object parameters

	// --- adapters
	WdfSeriesAdaptor				seriesAdaptor_RL;				///< adaptor for series RL
	WdfParallelTerminatedAdaptor	parallelTerminatedAdaptor_C;	///< adaptopr for parallel C

	DspFloatType sampleRate = 1.0;

};

/**
\class WDFIdealRLCHPF
\ingroup WDF-Objects
\brief
The WDFIdealRLCHPF object implements an ideal RLC HPF using the WDF library.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use WDFParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WDFIdealRLCHPF : public IAudioSignalProcessor
{
public:
	WDFIdealRLCHPF(void) { createWDF(); }	/* C-TOR */
	~WDFIdealRLCHPF(void) {}	/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		// --- rest WDF components (flush state registers)
		seriesAdaptor_RC.reset(_sampleRate);
		parallelTerminatedAdaptor_L.reset(_sampleRate);

		// --- intialize the chain of adapters
		seriesAdaptor_RC.initializeAdaptorChain();
		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the WDF Ideal RLC filter to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- push audio sample into series L1
		seriesAdaptor_RC.setInput1(xn);

		// --- output is at terminated L2's output2
		//     note compensation scaling by -6dB = 0.5
		//     because of WDF assumption about Rs and Rload
		return 0.5*parallelTerminatedAdaptor_L.getOutput2();
	}

	/** create WDF structure; may be called more than once */
	void createWDF()
	{
		// --- create components, init to noramlized values fc =
		//	   initial values for fc = 1kHz Q = 0.707
		//     Holding C Constant at 1e-6
		//			   L = 2.533e-2
		//			   R = 2.251131 e2
		seriesAdaptor_RC.setComponent(wdfComponent::seriesRC, 2.251131e2, 1.0e-6);
		parallelTerminatedAdaptor_L.setComponent(wdfComponent::L, 2.533e-2);

		// --- connect adapters
		WdfAdaptorBase::connectAdaptors(&seriesAdaptor_RC, &parallelTerminatedAdaptor_L);

		// --- set source resistance
		seriesAdaptor_RC.setSourceResistance(0.0); // --- Rs = 600

		// --- set open ckt termination
		parallelTerminatedAdaptor_L.setOpenTerminalResistance(true);
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return WDFParameters custom data structure
	*/
	WDFParameters getParameters() { return wdfParameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param WDFParameters custom data structure
	*/
	void setParameters(const WDFParameters& _wdfParameters)
	{
		if (_wdfParameters.fc != wdfParameters.fc ||
			_wdfParameters.Q != wdfParameters.Q ||
			_wdfParameters.boostCut_dB != wdfParameters.boostCut_dB ||
			_wdfParameters.frequencyWarping != wdfParameters.frequencyWarping)
		{
			wdfParameters = _wdfParameters;
			DspFloatType fc_Hz = wdfParameters.fc;

			if (wdfParameters.frequencyWarping)
			{
				DspFloatType arg = (kPi*fc_Hz) / sampleRate;
				fc_Hz = fc_Hz*(tan(arg) / arg);
			}

			DspFloatType inductorValue = 1.0 / (1.0e-6 * pow((2.0*kPi*fc_Hz), 2.0));
			DspFloatType resistorValue = (1.0 / wdfParameters.Q)*(pow(inductorValue / 1.0e-6, 0.5));

			seriesAdaptor_RC.setComponentValue_RC(resistorValue, 1.0e-6);
			parallelTerminatedAdaptor_L.setComponentValue(inductorValue);
			seriesAdaptor_RC.initializeAdaptorChain();
		}
	}


protected:
	WDFParameters wdfParameters;	///< object parameters

	// --- three
	WdfSeriesAdaptor				seriesAdaptor_RC;				///< adaptor for RC
	WdfParallelTerminatedAdaptor	parallelTerminatedAdaptor_L;	///< adaptor for L

	DspFloatType sampleRate = 1.0;	///< sample rate storage
};

/**
\class WDFIdealRLCBPF
\ingroup WDF-Objects
\brief
The WDFIdealRLCBPF object implements an ideal RLC BPF using the WDF library.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use WDFParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WDFIdealRLCBPF : public IAudioSignalProcessor
{
public:
	WDFIdealRLCBPF(void) { createWDF(); }	/* C-TOR */
	~WDFIdealRLCBPF(void) {}	/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		// --- rest WDF components (flush state registers)
		seriesAdaptor_LC.reset(_sampleRate);
		parallelTerminatedAdaptor_R.reset(_sampleRate);

		// --- intialize the chain of adapters
		seriesAdaptor_LC.initializeAdaptorChain();
		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the WDF Ideal RLC filter to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- push audio sample into series L1
		seriesAdaptor_LC.setInput1(xn);

		// --- output is at terminated L2's output2
		//     note compensation scaling by -6dB = 0.5
		//     because of WDF assumption about Rs and Rload
		return 0.5*parallelTerminatedAdaptor_R.getOutput2();
	}

	/** create the WDF structure*/
	void createWDF()
	{
		// --- create components, init to noramlized values fc =
		//	   initial values for fc = 1kHz Q = 0.707
		//     Holding C Constant at 1e-6
		//			   L = 2.533e-2
		//			   R = 2.251131 e2
		seriesAdaptor_LC.setComponent(wdfComponent::seriesLC, 2.533e-2, 1.0e-6);
		parallelTerminatedAdaptor_R.setComponent(wdfComponent::R, 2.251131e2);

		// --- connect adapters
		WdfAdaptorBase::connectAdaptors(&seriesAdaptor_LC, &parallelTerminatedAdaptor_R);

		// --- set source resistance
		seriesAdaptor_LC.setSourceResistance(0.0); // --- Rs = 600

		// --- set open ckt termination
		parallelTerminatedAdaptor_R.setOpenTerminalResistance(true);
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return WDFParameters custom data structure
	*/
	WDFParameters getParameters() { return wdfParameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param WDFParameters custom data structure
	*/
	void setParameters(const WDFParameters& _wdfParameters)
	{
		if (_wdfParameters.fc != wdfParameters.fc ||
			_wdfParameters.Q != wdfParameters.Q ||
			_wdfParameters.boostCut_dB != wdfParameters.boostCut_dB ||
			_wdfParameters.frequencyWarping != wdfParameters.frequencyWarping)
		{
			wdfParameters = _wdfParameters;
			DspFloatType fc_Hz = wdfParameters.fc;

			if (wdfParameters.frequencyWarping)
			{
				DspFloatType arg = (kPi*fc_Hz) / sampleRate;
				fc_Hz = fc_Hz*(tan(arg) / arg);
			}

			DspFloatType inductorValue = 1.0 / (1.0e-6 * pow((2.0*kPi*fc_Hz), 2.0));
			DspFloatType resistorValue = (1.0 / wdfParameters.Q)*(pow(inductorValue / 1.0e-6, 0.5));

			seriesAdaptor_LC.setComponentValue_LC(inductorValue, 1.0e-6);
			parallelTerminatedAdaptor_R.setComponentValue(resistorValue);
			seriesAdaptor_LC.initializeAdaptorChain();
		}
	}

protected:
	WDFParameters wdfParameters;	///< object parameters

	// --- adapters
	WdfSeriesAdaptor				seriesAdaptor_LC; ///< adaptor for LC
	WdfParallelTerminatedAdaptor	parallelTerminatedAdaptor_R; ///< adaptor for R

	DspFloatType sampleRate = 1.0;
};


/**
\class WDFIdealRLCBSF
\ingroup WDF-Objects
\brief
The WDFIdealRLCBSF object implements an ideal RLC BSF using the WDF library.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use WDFParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class WDFIdealRLCBSF : public IAudioSignalProcessor
{
public:
	WDFIdealRLCBSF(void) { createWDF(); }	/* C-TOR */
	~WDFIdealRLCBSF(void) {}	/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		// --- rest WDF components (flush state registers)
		seriesAdaptor_R.reset(_sampleRate);
		parallelTerminatedAdaptor_LC.reset(_sampleRate);

		// --- intialize the chain of adapters
		seriesAdaptor_R.initializeAdaptorChain();
		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the WDF Ideal RLC filter to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- push audio sample into series L1
		seriesAdaptor_R.setInput1(xn);

		// --- output is at terminated L2's output2
		//     note compensation scaling by -6dB = 0.5
		//     because of WDF assumption about Rs and Rload
		return 0.5*parallelTerminatedAdaptor_LC.getOutput2();
	}

	/** create WDF structure */
	void createWDF()
	{
		// --- create components, init to noramlized values fc =
		//	   initial values for fc = 1kHz Q = 0.707
		//     Holding C Constant at 1e-6
		//			   L = 2.533e-2
		//			   R = 2.251131 e2
		seriesAdaptor_R.setComponent(wdfComponent::R, 2.533e-2);
		parallelTerminatedAdaptor_LC.setComponent(wdfComponent::seriesLC, 2.533e-2, 1.0e-6);

		// --- connect adapters
		WdfAdaptorBase::connectAdaptors(&seriesAdaptor_R, &parallelTerminatedAdaptor_LC);

		// --- set source resistance
		seriesAdaptor_R.setSourceResistance(0.0); // --- Rs = 600

		// --- set open ckt termination
		parallelTerminatedAdaptor_LC.setOpenTerminalResistance(true);
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return WDFParameters custom data structure
	*/
	WDFParameters getParameters() { return wdfParameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param WDFParameters custom data structure
	*/
	void setParameters(const WDFParameters& _wdfParameters)
	{
		if (_wdfParameters.fc != wdfParameters.fc ||
			_wdfParameters.Q != wdfParameters.Q ||
			_wdfParameters.boostCut_dB != wdfParameters.boostCut_dB ||
			_wdfParameters.frequencyWarping != wdfParameters.frequencyWarping)
		{
			wdfParameters = _wdfParameters;
			DspFloatType fc_Hz = wdfParameters.fc;

			if (wdfParameters.frequencyWarping)
			{
				DspFloatType arg = (kPi*fc_Hz) / sampleRate;
				fc_Hz = fc_Hz*(tan(arg) / arg);
			}

			DspFloatType inductorValue = 1.0 / (1.0e-6 * pow((2.0*kPi*fc_Hz), 2.0));
			DspFloatType resistorValue = (1.0 / wdfParameters.Q)*(pow(inductorValue / 1.0e-6, 0.5));

			seriesAdaptor_R.setComponentValue(resistorValue);
			parallelTerminatedAdaptor_LC.setComponentValue_LC(inductorValue, 1.0e-6);
			seriesAdaptor_R.initializeAdaptorChain();
		}
	}

protected:
	WDFParameters wdfParameters;	///< object parameters

	// --- adapters
	WdfSeriesAdaptor				seriesAdaptor_R; ///< adaptor for series R
	WdfParallelTerminatedAdaptor	parallelTerminatedAdaptor_LC; ///< adaptor for parallel LC

	DspFloatType sampleRate = 1.0; ///< sample rate storage
};

