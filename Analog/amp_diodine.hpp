#pragma once
// https://github.com/pwmagro/diodine    

namespace Diodine
{
    class RingBuffer
    {
    public:
        RingBuffer() = default;
        ~RingBuffer() = default;

        void prepare(juce::dsp::ProcessSpec& spec);
        void writeSamples(juce::dsp::AudioBlock<float>& buffer);
        typedef struct { float max; float min; } maxmin_t;
        maxmin_t readSamples();
        size_t getBufferSize() { return buffer.size(); }
        float getSample(int i);

    private:
        std::vector<float> buffer;
        std::atomic<size_t> writeAtomic{ 0 };
        size_t read{ 0 };
    };

    class Diode {
    public:
        Diode();
        ~Diode();
        void init(double sampleRate, int samplesPerBlock);
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out);
        static float waveshape(float x);
        RingBuffer::maxmin_t readRingBuffer();
        RingBuffer* getRingBuffer() { return &ringBuffer; }
        typedef struct { float left; float right; } rrStatus_t;
        rrStatus_t getRrStatus() { return rrStatus; };

    private:
        typedef struct {
            bool diode1 = true;
            bool diode2 = false;
            float vf = 0;
            float vb = 0;
            float trr = 0;
            float skew = 0.2;
            float charge = 0;
            float gain = 0;
            float sat = 0;
            float mix = 1;
        } DiodeProperties_t;

        static float wsAsym(float x, DiodeProperties_t& diodeProperties);
        static float wsSym(float x, DiodeProperties_t& diodeProperties);
        void setTrr(float trr);
        void recover(juce::AudioBuffer<float>& buffer);

        DiodeProperties_t diodeProperties;
        
        float samplesPerMs;
        int recoverScanner[2] = { 0, 0 };
        float lastSamples[2] = { 0, 0 };

        //juce::dsp::LinkwitzRileyFilter<float> dcOffset;
        RingBuffer ringBuffer;

        std::vector<float> rr;
        float lastTrr;
        float lastSkew;
        rrStatus_t rrStatus = { 0, 0 };
    };

    /*
    ==============================================================================
        RingBuffer.cpp
        Created: 6 Sep 2022 3:33:30pm
        Author:  thesp
    ==============================================================================
    */

    void RingBuffer::prepare(juce::dsp::ProcessSpec& spec)
    {
        buffer.resize(spec.sampleRate);
        writeAtomic.store(0, std::memory_order_relaxed);
        read = 0;
    }
    void RingBuffer::writeSamples(size_t n, DspFloatType ** newBuffer)
    {
        // Only works with 2 channels
        assert(n == 2);
        size_t write = writeAtomic.load(std::memory_order_relaxed);

        const size_t numSamples = n;
        const size_t limit = numSamples + write;
        auto channel1 = newBuffer[0];
        auto channel2 = newBuffer[1];

        // if samples fit within buffer
        if (limit < buffer.size())
        {
            size_t i = 0;
            for (size_t s = write; s < limit; ++s)
            {
                buffer[s] = (channel1[i] + channel2[i]) * 0.5f;
                ++i;
            }

            writeAtomic.store(limit % buffer.size(), std::memory_order_release);

            return;
        }
        
        // if samples go past buffer
        size_t i = 0;
        for (size_t s = write; s < buffer.size(); ++s)
        {
            buffer[s] = (channel1[i] + channel2[i]) * 0.5f;
            ++i;
        }

        const size_t newLimit = limit - buffer.size();

        for (size_t s = 0; s < newLimit; ++s)
        {
            buffer[s] = (channel1[i] + channel2[i]) * 0.5f;
            ++i;
        }

        writeAtomic.store(newLimit % buffer.size(), std::memory_order_release);
        return;
    }
    RingBuffer::maxmin_t RingBuffer::readSamples()
    {
        size_t write = writeAtomic.load(std::memory_order_acquire);

        if (read == write) return { 0.f, 0.f };

        // If read is before write
        if (read < write)
        {
            float loudest = 0.f;
            float quietest = 0.f;

            for (read; read < write; ++read) {
                loudest = std::max(buffer[read], loudest);
                quietest = std::min(buffer[read], quietest);
            }

            return { loudest, quietest };
        }

        // If read is NOT before write
        float loudest = 0.f;
        float quietest = 0.f;

        for (read; read < buffer.size(); ++read) {
            loudest = std::max(buffer[read], loudest);
            quietest = std::min(buffer[read], quietest);
        }

        for (read = 0; read < write; ++read) {
            loudest = std::max(buffer[read], loudest);
            quietest = std::min(buffer[read], quietest);
        }

        return { loudest, quietest };
    }

    Diode::Diode()
    {
        diodeProperties = {};
        samplesPerMs = 0;
    }

    Diode::~Diode() {

    }

    void Diode::init(double sampleRate, int samplesPerBlock) {
        samplesPerMs = sampleRate / 1000.f;

        // TODO switch process() block to use contexts + implement filter
        //dcOffset.setCutoffFrequency(15.f);
        //dcOffset.setType(juce::dsp::LinkwitzRileyFilterType::highpass);

        ringBuffer.prepare(sampleRate,samplesPerBlock);
        rr.reserve(samplesPerMs * 20);
    }

    void Diode::process(juce::AudioBuffer<float>& buffer, juce::AudioProcessorValueTreeState& apvts) {
        ringBuffer.writeSamples(juce::dsp::AudioBlock<float>(buffer));
        auto ch = buffer.getNumChannels();
        if (ch != 2) {
            throw std::exception("Only two-channel audio is supported.");
        }

        auto sm = buffer.getNumSamples();

        diodeProperties.diode1 = apvts.getRawParameterValue(DIODE_1_ID)->load();
        diodeProperties.diode2 = apvts.getRawParameterValue(DIODE_2_ID)->load();
        diodeProperties.vf = apvts.getRawParameterValue(VF_ID)->load();
        diodeProperties.vb = apvts.getRawParameterValue(VB_ID)->load();
        diodeProperties.trr = apvts.getRawParameterValue(TRR_ID)->load();
        diodeProperties.skew = apvts.getRawParameterValue(TRR_SKEW_ID)->load();
        diodeProperties.charge = apvts.getRawParameterValue(TRR_MAG_ID)->load();
        diodeProperties.gain = apvts.getRawParameterValue(GAIN_ID)->load();
        diodeProperties.sat = apvts.getRawParameterValue(SAT_ID)->load() * 10;
        diodeProperties.mix = apvts.getRawParameterValue(MIX_ID)->load();

        if (diodeProperties.trr != lastTrr || diodeProperties.skew != lastSkew) {
            setTrr(diodeProperties.trr);
        }

        if (diodeProperties.diode1 && diodeProperties.diode2) {
            for (int c = 0; c < ch; c++) {
                auto channelWr = buffer.getWritePointer(c);
                auto channelRd = buffer.getReadPointer(c);

                for (int s = 0; s < sm; s++) {
                    channelWr[s] = ((1 - diodeProperties.mix) * channelRd[s]) + (diodeProperties.mix * wsSym(channelRd[s], diodeProperties));
                }
            }
        }
        else if (diodeProperties.diode1) {
            for (int c = 0; c < ch; c++) {
                auto channelWr = buffer.getWritePointer(c);
                auto channelRd = buffer.getReadPointer(c);

                for (int s = 0; s < sm; s++) {
                    channelWr[s] = ((1 - diodeProperties.mix) * channelRd[s]) + (diodeProperties.mix * wsAsym(channelRd[s], diodeProperties));
                }
            }
        }
        else if (diodeProperties.diode2) {
            for (int c = 0; c < ch; c++) {
                auto channelWr = buffer.getWritePointer(c);
                auto channelRd = buffer.getReadPointer(c);

                for (int s = 0; s < sm; s++) {
                    channelWr[s] = (((1 - diodeProperties.mix) * channelRd[s]) - (diodeProperties.mix * wsAsym(-channelRd[s], diodeProperties)));
                }
            }
        }
        else {
            for (int c = 0; c < ch; c++) {
                auto channelWr = buffer.getWritePointer(c);
                auto channelRd = buffer.getReadPointer(c);

                for (int s = 0; s < sm; s++) {
                    channelWr[s] = (((1 - diodeProperties.mix) * channelRd[s]));
                }
            }
        }
        lastSamples[0] = buffer.getReadPointer(0)[buffer.getNumSamples() - 1];
        lastSamples[1] = buffer.getReadPointer(1)[buffer.getNumSamples() - 1];

        recover(buffer);

    }

    float Diode::wsAsym(float x, DiodeProperties_t& diodeProperties) {
        // Saturate
        x = std::tanh((diodeProperties.sat + 0.1) * x) / (std::tanh(diodeProperties.sat + 0.1));

        // Gain
        float s = x * (1 + diodeProperties.gain);
        s = diodeProperties.diode1 ? s : -s;

        // Rectify
        float o1 = 0;
        if (s > diodeProperties.vf)
            o1 = s - diodeProperties.vf;
        float o2;
        if (s > diodeProperties.vb)
            o2 = o1;
        else
            o2 = s - diodeProperties.vb * 0.9;

        // Normalize
        float o3 = o2 / (diodeProperties.gain + 1);

        // Clip
        x = std::max(o3, -1.f);
        x = std::min(x, 1.f);



        return x;

    }

    float Diode::wsSym(float x, DiodeProperties_t& diodeProperties) {
        // Saturate
        x = std::tanh((diodeProperties.sat + 0.1) * x) / (std::tanh(diodeProperties.sat + 0.1));

        int sign = x > 0 ? 1 : -1;

        // Gain
        float s = x * (1 + diodeProperties.gain);

        // Rectify
        s = abs(s);
        float o1 = 0;
        if (s > diodeProperties.vf)
            o1 = s - diodeProperties.vf;
        else
            s = 0;

        // Normalize
        float o2 = o1 / (diodeProperties.gain + 1);

        // Clip
        x = std::min(o2, 1.f);


        return x * sign;
    }

    void Diode::recover(juce::AudioBuffer<float>& buffer) {
        float c = diodeProperties.charge;
        juce::AudioBuffer<float> dry;
        dry.makeCopyOf(buffer);

        for (int n = 0; n < 2; n++) {

            auto chr = dry.getReadPointer(n);
            auto ch = buffer.getWritePointer(n);

            if (lastSamples[n] * diodeProperties.gain> (diodeProperties.vf) && chr[0] * diodeProperties.gain < diodeProperties.vf) {
                recoverScanner[n] = 0;
                n == 0 ? rrStatus.left = -rr[0] : rrStatus.right = -rr[0];
            }
            for (int s = 1; s < buffer.getNumSamples(); s++) {
                if (chr[s - 1] * diodeProperties.gain > (diodeProperties.vf) && chr[s] * diodeProperties.gain < diodeProperties.vf) {
                    recoverScanner[n] = 0;
                    n == 0 ? rrStatus.left = -rr[0] : rrStatus.right = -rr[0];
                }
                if (recoverScanner[n] < juce::roundFloatToInt(diodeProperties.trr * samplesPerMs - 1)) {
                    ch[s] = std::min(1.f, ch[s] + c * rr[recoverScanner[n]] * diodeProperties.vf / diodeProperties.gain);
                    if ((abs(rr[recoverScanner[n]])) > (n == 0 ? rrStatus.left : rrStatus.right))
                        n == 0 ? rrStatus.left = -rr[recoverScanner[n]] : rrStatus.right = -rr[recoverScanner[n]];

                    recoverScanner[n]++;
                }
            }
            rrStatus.left *= 0.97;
            rrStatus.right *= 0.97;
        }
    }

    void Diode::setTrr(float trr) {
        float charge = 1.f;        // not an actual charge value but yknow
        lastTrr = trr;

        int trrInSamples = juce::roundFloatToInt(trr * samplesPerMs);
        if (trrInSamples == 0) {
            rr.push_back(0.f);
            return;
        }

        rr.clear();
        float x;
        float tension = diodeProperties.skew;
        for (float i = 0; i < trrInSamples * tension; i++) {
            x = i / (float)trrInSamples;
            
            rr.push_back((x / pow(tension, 2)) * (x - 2 * tension));
        }

        for (float i = trrInSamples * tension; i <= trrInSamples; i++) {
            x = i / (float)trrInSamples;
            rr.push_back(-pow((x - 1) * (x + (1 - 2 * tension)) / ((tension - 1) * (1 - tension)), 2));
        }
    }

    float Diode::waveshape(float x, juce::AudioProcessorValueTreeState& apvts) {
        DiodeProperties_t diodeProperties;
        diodeProperties.diode1 = apvts.getRawParameterValue(DIODE_1_ID)->load();
        diodeProperties.diode2 = apvts.getRawParameterValue(DIODE_2_ID)->load();
        diodeProperties.vf = apvts.getRawParameterValue(VF_ID)->load();
        diodeProperties.vb = apvts.getRawParameterValue(VB_ID)->load();
        diodeProperties.trr = apvts.getRawParameterValue(TRR_ID)->load();
        diodeProperties.gain = apvts.getRawParameterValue(GAIN_ID)->load();
        diodeProperties.sat = apvts.getRawParameterValue(SAT_ID)->load() * 10;
        diodeProperties.mix = apvts.getRawParameterValue(MIX_ID)->load();

        if (diodeProperties.diode1 && diodeProperties.diode2) {
            return wsSym(x, diodeProperties);
        }
        if (diodeProperties.diode1) {
            return wsAsym(x, diodeProperties);
        }
        if (diodeProperties.diode2) {
            return -wsAsym(x, diodeProperties);
        }
        else return 0;
    }

    RingBuffer::maxmin_t Diode::readRingBuffer() {
        return ringBuffer.readSamples();
    }
}