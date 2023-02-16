%module amplifiersfuncs
%{
typedef float DspFloatType;
#include "FX/Amplifiers.hpp"

using namespace FX::Distortion;
%}
typedef float DspFloatType;
%include "SoundObject.hpp"
%include "FX/Amplifiers.hpp"

%inline %{
    AmplifierFunction2 AmpUdo1;
    AmplifierFunction1 AmpFold;
    AmplifierFunction1 AmpWrap;
    AmplifierFunction1 AmpSinFold;
    AmplifierFunction2 AmpCheby;
    AmplifierFunction2 AmpChebyPolynomial;
    AmplifierFunction3 AmpClamp;
    AmplifierFunction1 AmpPreamp;
    AmplifierFunction1 AmpPostamp;
    AmplifierFunction3 AmpTanhNormal;
    AmplifierFunction2 AmpSigmoidFunction;
    AmplifierFunction2 AmpSigmoidMinus;
    AmplifierFunction2 AmpBipolarSigmoid;
    AmplifierFunction1 AmpFullWaveRectify;
    AmplifierFunction1 AmpHalfWaveRectify;

    AmplifierFunction1 AmpDistortion1;    
    AmplifierFunction1 AmpAsymmetricSigmoidDistortion;
    AmplifierFunction1 AmpAsymmetricSigmoidDistortion2;
    AmplifierFunction1 AmpHardClip;
    AmplifierFunction1 AmpLogicClip;
    AmplifierFunction1 AmpHardClip2;
    AmplifierFunction1 AmpTanClip;
    AmplifierFunction1 AmpQuintic;
    AmplifierFunction1 AmpCubicBasic;
    AmplifierFunction1 AmpAlgClip;
    AmplifierFunction1 AmpArcClip;
    AmplifierFunction1 AmpSinClip;
    AmplifierFunction1 AmpFuzzCtrTable;
    AmplifierFunction1 AmpAbsf;    
    AmplifierFunction1 AmpCubef;
    AmplifierFunction1 AmpHardClip4;
    AmplifierFunction1 AmpCubicWaveShaper;
    AmplifierFunction1 AmpGloubiBoulga;    
    AmplifierFunction1 AmpGloubiApprox;
    AmplifierFunction1 AmpFuzzEdgeTable;
    AmplifierFunction1 AmpLimitClip;
    AmplifierFunction1 AmpTubeClip;
    AmplifierFunction1 AmpHardClipping;
    AmplifierFunction1 AmpSoftClipping;
    AmplifierFunction1 AmpExponential;
    AmplifierFunction1 AmpFullWaveRectifier;
    AmplifierFunction1 AmpHalfWaveRectifier;
    AmplifierFunction1 AmpArayaAndSuyama;
    AmplifierFunction1 AmpDoidicSymmetric;    
    AmplifierFunction1 AmpDoidicAsymmetric;
    
    AmplifierFunction2 AmpRectify;    
    AmplifierFunction2 AmpHardClip3;
    AmplifierFunction2 AmpSoftCubicClip;
    AmplifierFunction2 AmpSKClip;
    AmplifierFunction2 AmpSoftKnee;
    AmplifierFunction2 AmpSoftClipShape;
    AmplifierFunction2 AmpSinShape;
    AmplifierFunction2 AmpChebyshevShape;
    AmplifierFunction2 AmpArcTangent;
    AmplifierFunction2 AmpSquareLaw;
    AmplifierFunction2 AmpFoldback;
    AmplifierFunction2 AmpWaveshaper1;
    AmplifierFunction2 AmpWaveshaper2;    
    AmplifierFunction2 AmpWaveshaper3;
    
    AmplifierFunction3 AmpSoftCubic;
    AmplifierFunction3 AmpChebyshevRec;

    AmplifierFunction4 AmpFuzzExponential;    
    AmplifierFunction4 AmpSigmoidDistortionFunction;
    AmplifierFunction4 AmpAsymmetricSigmoid;
    AmplifierFunction4 AmpAsymmetricSigmoid2;
    AmplifierFunction4 AmpDistortionFunction;
    AmplifierFunction4 AmpCubicDistortion;
    AmplifierFunction4 AmpASinDistortion;    
    AmplifierFunction4 AmpATanDistortion;
    AmplifierFunction4 AmpASinhDistortion;
    AmplifierFunction4 AmpACoshDistortion;
    AmplifierFunction4 AmpATanhDistortion;
    AmplifierFunction4 AmpExpDistortion;    
    AmplifierFunction4 AmpDcDistortion;
    AmplifierFunction4 AmpBipolarDistortion;
    AmplifierFunction4 AmpQuadraticDistortion;    
    AmplifierFunction4 AmpQuadratic2Distortion;    
    AmplifierFunction4 AmpQuadratic3Distortion;
    AmplifierFunction4 AmpParametricClip;
    AmplifierFunction4 AmpArcTanDistortion;
    AmplifierFunction4 AmpSoftClipper;
    AmplifierFunction4 AmpErrorf;
    AmplifierFunction4 AmpSigmoided;
    AmplifierFunction4 AmpHyperbolicTangent;
    AmplifierFunction4 AmpDiodeClipping;
    AmplifierFunction4 AmpPieceWiseOverdrive;
    AmplifierFunction4 AmpTube;
    AmplifierFunction4 AmpArraya;
    AmplifierFunction4 AmpGallo;    
    AmplifierFunction4 AmpDoubleSoftClipper;
    AmplifierFunction4 AmpCrush;
    AmplifierFunction4 AmpTuboid;    
    AmplifierFunction4 AmpPakarinenYeh;

    /*                    
    AmplifierFunction1 AmpPositiveSignal(positive_signal);
    AmplifierFunction1 AmpNegativeSignal(negative_signal);    
    AmplifierFunction1 AmpModulatedSignal(modulated_signals);
    AmplifierFunction1 AmpCircularModulatedSignal(circular_modulated_signals);
    AmplifierFunction1 AmpPositiveModulatedSignal(positive_modulated_signals);
    AmplifierFunction1 AmpNegativeModulatedSignal(negative_modulated_signals);    
    */
    void initAmps()
    {
        AmpUdo1 = AmplifierFunction2(udo1);
        AmpFold = AmplifierFunction1(Fold);
        AmpWrap = AmplifierFunction1(Wrap);
        AmpSinFold = AmplifierFunction1(SinFold);
        AmpCheby = AmplifierFunction2(cheby);
        AmpChebyPolynomial = AmplifierFunction2(cheby_polynomial);
        AmpClamp = AmplifierFunction3(amp_clamp);
        AmpPreamp = AmplifierFunction1(preamp);
        AmpPostamp = AmplifierFunction1(postamp);
        AmpTanhNormal = AmplifierFunction3(tanh_normal);        
        AmpSigmoidFunction = AmplifierFunction2(sigmoid_function);
        AmpSigmoidMinus = AmplifierFunction2(sigmoid_minus);
        AmpBipolarSigmoid = AmplifierFunction2(bpsigmoid);
        AmpFullWaveRectify = AmplifierFunction1(full_rectify);
        AmpHalfWaveRectify = AmplifierFunction1(half_rectify);
        
        AmpDistortion1 = AmplifierFunction1(distortionFunction);
        AmpAsymmetricSigmoidDistortion = AmplifierFunction1(asymmetricSigmoidDistortionFunction);
        AmpAsymmetricSigmoidDistortion2 = AmplifierFunction1(asymmetricSigmoidDistortionFunction2);
        AmpHardClip = AmplifierFunction1(hardClip);
        AmpLogicClip = AmplifierFunction1(logiclip<DspFloatType>);
        AmpHardClip = AmplifierFunction1(hardclip<DspFloatType>);
        AmpTanClip = AmplifierFunction1(tanclip<DspFloatType>);
        AmpQuintic = AmplifierFunction1(quintic<DspFloatType>);
        AmpCubicBasic = AmplifierFunction1(cubicBasic<DspFloatType>);
        AmpAlgClip = AmplifierFunction1(algClip<DspFloatType>);
        AmpArcClip = AmplifierFunction1(arcClip<DspFloatType>);
        AmpSinClip = AmplifierFunction1(sinclip<DspFloatType>);
        AmpFuzzCtrTable = AmplifierFunction1(FuzzCtrTable);
        AmpAbsf = AmplifierFunction1(absf);
        AmpCubef = AmplifierFunction1(cubef);
        AmpHardClip = AmplifierFunction1(hardClip);
        AmpCubicWaveShaper = AmplifierFunction1(cubicWaveShaper);
        AmpGloubiBoulga = AmplifierFunction1(gloubiBoulga);
        AmpGloubiApprox = AmplifierFunction1(gloubiApprox);
        AmpFuzzEdgeTable = AmplifierFunction1(FuzzEdgeTable);
        AmpLimitClip = AmplifierFunction1(limitclip<DspFloatType>);
        AmpTubeClip = AmplifierFunction1(tubeclip);
        AmpHardClipping = AmplifierFunction1(hardClipping);
        AmpSoftClipping = AmplifierFunction1(softClipping);
        AmpExponential = AmplifierFunction1(exponential);
        AmpFullWaveRectifier = AmplifierFunction1(fullWaveRectifier);
        AmpHalfWaveRectifier = AmplifierFunction1(halfWaveRectifier);
        AmpArayaAndSuyama = AmplifierFunction1(ArayaAndSuyama);
        AmpDoidicSymmetric = AmplifierFunction1(doidicSymmetric);
        AmpDoidicAsymmetric = AmplifierFunction1(doidicAssymetric);

        AmpRectify = AmplifierFunction2(Rectify);
        AmpHardClip3 = AmplifierFunction2(HardClip);
        AmpSoftCubicClip = AmplifierFunction2(SoftCubicClip);   
        AmpSKClip = AmplifierFunction2(SKClip);
        AmpSoftKnee = AmplifierFunction2(SoftKnee);
        AmpSoftClipShape = AmplifierFunction2(softClipShape);
        AmpSinShape = AmplifierFunction2(sineShape);    
        AmpChebyshevShape = AmplifierFunction2(chebyshevShape);
        AmpArcTangent = AmplifierFunction2(arctangent);    
        AmpSquareLaw = AmplifierFunction2(squareLaw);
        AmpFoldback = AmplifierFunction2(foldback);
        AmpWaveshaper1 = AmplifierFunction2(waveShaper1);
        AmpWaveshaper2 = AmplifierFunction2(waveShaper2);
        AmpWaveshaper3 = AmplifierFunction2(waveShaper3);

        AmpSoftCubic = AmplifierFunction3(SoftCubic);
        AmpChebyshevRec = AmplifierFunction3(chebyshevRec);

        AmpASinDistortion = AmplifierFunction4(asin_distortion);
        AmpATanDistortion = AmplifierFunction4(atan_distortion);
        AmpASinhDistortion = AmplifierFunction4(asinh_distortion);
        AmpACoshDistortion = AmplifierFunction4(acosh_distortion);
        AmpATanhDistortion = AmplifierFunction4(atanh_distortion);
        AmpExpDistortion = AmplifierFunction4(exp_distortion);
        AmpDcDistortion = AmplifierFunction4(dc_distortion);
        AmpBipolarDistortion = AmplifierFunction4(bipolar_distortion);
        AmpQuadraticDistortion = AmplifierFunction4(quadratic_distortion);
        AmpQuadratic2Distortion = AmplifierFunction4(quadratic2_distortion);
        AmpQuadratic3Distortion = AmplifierFunction4(quadratic3_distortion);
        AmpParametricClip = AmplifierFunction4(parametric_clip);
        AmpArcTanDistortion = AmplifierFunction4(arcTanDistortion);
        AmpSoftClipper = AmplifierFunction4(softClipper);
        AmpErrorf = AmplifierFunction4(errorf);
        AmpSigmoided = AmplifierFunction4(sigmoided);
        AmpHyperbolicTangent = AmplifierFunction4(hyperbolicTangent);
        AmpDiodeClipping = AmplifierFunction4(diodeClipping);    
        AmpPieceWiseOverdrive = AmplifierFunction4(pieceWiseOverdrive);    
        AmpTube = AmplifierFunction4(tube);    
        AmpArraya = AmplifierFunction4(arraya);
        AmpGallo = AmplifierFunction4(gallo);
        AmpDoubleSoftClipper = AmplifierFunction4(doubleSoftClipper);
        AmpCrush = AmplifierFunction4(crush);
        AmpTuboid = AmplifierFunction4(tuboid);
        AmpPakarinenYeh = AmplifierFunction4(pakarinen_Yeh);
    }      
%}