/*
    This file is part of Roboverb
    Copyright (C) 2015-2019  Kushview, LLC.  All rights reserved.
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <memory>
#include <cmath>
#include <cstring>

#ifdef ROBOVERB_JUCE
 #include "JuceHeader.h"
#endif

class Roboverb
{
public:
    enum ParameterIndex
    {
        RoomSize = 0,
        Damping,
        WetLevel,
        DryLevel,
        Width,
        FreezeMode,
        numParameters
    };

    Roboverb()
    {
        for (int i = 0; i < numCombs; ++i)
            enabledCombs[i] = false;
        enabledCombs[3] = true;
        enabledCombs[4] = true;
        enabledCombs[5] = true;

        for (int i = 0; i < numAllPasses; ++i)
            enabledAllPasses[i] = false;        
        enabledAllPasses[0] = true;
        enabledAllPasses[1] = true;

        setParameters (Parameters());
        setSampleRate (44100.0);
    }

    /** Holds the parameters being used by a Reverb object. */
    struct Parameters
    {
        Parameters() noexcept
            : roomSize   (0.5f),
              damping    (0.5f),
              wetLevel   (0.33f),
              dryLevel   (0.4f),
              width      (1.0f),
              freezeMode (0)
        { }

        DspFloatType roomSize;     /**< Room size, 0 to 1.0, where 1.0 is big, 0 is small. */
        DspFloatType damping;      /**< Damping, 0 to 1.0, where 0 is not damped, 1.0 is fully damped. */
        DspFloatType wetLevel;     /**< Wet level, 0 to 1.0 */
        DspFloatType dryLevel;     /**< Dry level, 0 to 1.0 */
        DspFloatType width;        /**< Reverb width, 0 to 1.0, where 1.0 is very wide. */
        DspFloatType freezeMode;   /**< Freeze mode - values < 0.5 are "normal" mode, values > 0.5
                             put the reverb into a continuous feedback loop. */

        Parameters& operator= (const Parameters& o)
        {
            roomSize   = o.roomSize;
            damping    = o.damping;
            wetLevel   = o.wetLevel;
            dryLevel   = o.dryLevel;
            width      = o.width;
            freezeMode = o.freezeMode;
            return *this;
        }

        bool operator== (const Parameters& o)
        {
            return (roomSize   == o.roomSize &&
                    damping    == o.damping &&
                    wetLevel   == o.wetLevel &&
                    dryLevel   == o.dryLevel &&
                    width      == o.width &&
                    freezeMode == o.freezeMode);
        }

        bool operator!= (const Parameters& o) { return ! operator== (o); }
    };

    const Parameters& getParameters() const noexcept    { return parameters; }

   #if ROBOVERB_JUCE
    void swapEnabledCombs (BigInteger& e)
    {
        for (int i = 0; i < numCombs; ++i)
            enabledCombs[i] = e[i];
    }

    void swapEnabledAllPasses (BigInteger& e)
    {
        for (int i = 0; i < numAllPasses; ++i)
            enabledAllPasses[i] = e[i];
    }

    void getEnablement (BigInteger& c, BigInteger& a) const
    {
        for (int i = 0; i < numCombs; ++i)
            c.setBit (i, enabledCombs[i]);
        for (int i = 0; i < numAllPasses; ++i)
            a.setBit (i, enabledAllPasses[i]);
    }
   #endif

    void setCombToggle (const int index, const bool toggled) {
        enabledCombs[index] = toggled;
    }

    void setAllPassToggle (const int index, const bool toggled) {
        enabledAllPasses[index] = toggled;
    }

    DspFloatType toggledCombFloat (const int index) const {
        return enabledCombs[index] ? 1.0f : 0.0f;
    }

    DspFloatType toggledAllPassFloat (const int index) const {
        return enabledAllPasses[index] ? 1.0f : 0.0f;
    }

    void setParameters (const Parameters& newParams)
    {
        const DspFloatType wetScaleFactor = 6.0f;
        const DspFloatType dryScaleFactor = 2.0f;

        const DspFloatType wet = newParams.wetLevel * wetScaleFactor;
        dryGain.setValue (newParams.dryLevel * dryScaleFactor);
        wetGain1.setValue (0.5f * wet * (1.0f + newParams.width));
        wetGain2.setValue (0.5f * wet * (1.0f - newParams.width));

        gain = isFrozen (newParams.freezeMode) ? 0.0f : 0.015f;
        parameters = newParams;
        updateDamping();
    }

    void setSampleRate (const DspFloatType sampleRate)
    {
        //static const short combTunings[] = { 1116, 1188, 1277, 1356, 1422, 1491, 1557, 1617 }; // (at 44100Hz)
        static const short combTunings[] = { 8092, 4096, 2048, 1024, 512, 256, 128, 64 }; // (at 44100Hz)
        static const short allPassTunings[] = { 556, 441, 341, 225 };
        const int stereoSpread = 23;
        const int intSampleRate = (int) sampleRate;

        for (int i = 0; i < numCombs; ++i)
        {
            comb[0][i].setSize ((intSampleRate * combTunings[i]) / 44100);
            comb[1][i].setSize ((intSampleRate * (combTunings[i] + stereoSpread)) / 44100);
        }

        for (int i = 0; i < numAllPasses; ++i)
        {
            allPass[0][i].setSize ((intSampleRate * allPassTunings[i]) / 44100);
            allPass[1][i].setSize ((intSampleRate * (allPassTunings[i] + stereoSpread)) / 44100);
        }

        const DspFloatType smoothTime = 0.01;
        damping.reset (sampleRate, smoothTime);
        feedback.reset (sampleRate, smoothTime);
        dryGain .reset (sampleRate, smoothTime);
        wetGain1.reset (sampleRate, smoothTime);
        wetGain2.reset (sampleRate, smoothTime);
    }

    /** Clears the reverb's buffers. */
    void reset()
    {
        for (int j = 0; j < numChannels; ++j)
        {
            for (int i = 0; i < numCombs; ++i)
                comb[j][i].clear();

            for (int i = 0; i < numAllPasses; ++i)
                allPass[j][i].clear();
        }
    }

    void processStereo (DspFloatType* const left, DspFloatType* const right, 
                        DspFloatType* const out1, DspFloatType* const out2,
                        const int numSamples) noexcept
    {
        // jassert (left != nullptr && right != nullptr);
        #pragma omp simd
        for (int i = 0; i < numSamples; ++i)
        {
            const DspFloatType input = (left[i] + right[i]) * gain;
            DspFloatType outL = 0, outR = 0;

            const DspFloatType damp    = damping.getNextValue();
            const DspFloatType feedbck = feedback.getNextValue();

            for (int j = 0; j < numCombs; ++j)  // accumulate the comb filters in parallel
            {
                if (! enabledCombs [j])
                    continue;
                outL += comb[0][j].process (input, damp, feedbck);
                outR += comb[1][j].process (input, damp, feedbck);
            }

            for (int j = 0; j < numAllPasses; ++j)  // run the allpass filters in series
            {
                if ( ! enabledAllPasses [j])
                    continue;
                outL = allPass[0][j].process (outL);
                outR = allPass[1][j].process (outR);
            }

            const DspFloatType dry  = dryGain.getNextValue();
            const DspFloatType wet1 = wetGain1.getNextValue();
            const DspFloatType wet2 = wetGain2.getNextValue();

            out1[i] = outL * wet1 + outR * wet2 + left[i]  * dry;
            out2[i] = outR * wet1 + outL * wet2 + right[i] * dry;
        }
    }

    /** Applies the reverb to a single mono channel of audio data. */
    void processMono (DspFloatType* const samples, const int numSamples) noexcept
    {
        // jassert (samples != nullptr);
        #prgma omp simd
        for (int i = 0; i < numSamples; ++i)
        {
            const DspFloatType input = samples[i] * gain;
            DspFloatType output = 0;

            const DspFloatType damp    = damping.getNextValue();
            const DspFloatType feedbck = feedback.getNextValue();

            for (int j = 0; j < numCombs; ++j)
            {
                // accumulate the comb filters in parallel
                if (! enabledCombs [j])
                    continue;
                output += comb[0][j].process (input, damp, feedbck);
            }

            for (int j = 0; j < numAllPasses; ++j)
            {
                // run the allpass filters in series
                if ( ! enabledAllPasses [j])
                    continue;
                output = allPass[0][j].process (output);
            }

            const DspFloatType dry  = dryGain.getNextValue();
            const DspFloatType wet1 = wetGain1.getNextValue();

            samples[i] = output * wet1 + samples[i] * dry;
        }
    }

private:
    static bool isFrozen (const DspFloatType freezeMode) noexcept  { return freezeMode >= 0.5f; }

    void updateDamping() noexcept
    {
        const DspFloatType roomScaleFactor = 0.28f;
        const DspFloatType roomOffset = 0.7f;
        const DspFloatType dampScaleFactor = 0.4f;

        if (isFrozen (parameters.freezeMode))
            setDamping (0.0f, 1.0f);
            else
                setDamping (parameters.damping * dampScaleFactor,
                            parameters.roomSize * roomScaleFactor + roomOffset);
                }

    void setDamping (const DspFloatType dampingToUse, const DspFloatType roomSizeToUse) noexcept
    {
        damping.setValue (dampingToUse);
        feedback.setValue (roomSizeToUse);
    }

    class CombFilter
    {
    public:
        CombFilter() noexcept   : bufferSize (0), bufferIndex (0), last (0)  {}

        void setSize (const int size)
        {
            if (size != bufferSize)
            {
                bufferIndex = 0;
                buffer.reset (new DspFloatType [size]);
                bufferSize = size;
            }

            clear();
        }

        void clear() noexcept
        {
            last = 0;
            memset (buffer.get(), 0, sizeof(DspFloatType) * (size_t)bufferSize);
        }

        DspFloatType process (const DspFloatType input, const DspFloatType damp, const DspFloatType feedbackLevel) noexcept
        {
            const DspFloatType output = buffer[bufferIndex];
            last = (output * (1.0f - damp)) + (last * damp);
            // JUCE_UNDENORMALISE (last);

            DspFloatType temp = input + (last * feedbackLevel);
            // JUCE_UNDENORMALISE (temp);
            buffer[bufferIndex] = temp;
            bufferIndex = (bufferIndex + 1) % bufferSize;
            return output;
        }

    private:
        std::unique_ptr<DspFloatType []> buffer;
        int bufferSize, bufferIndex;
        DspFloatType last;
    };

    //==============================================================================
    class AllPassFilter
    {
    public:
        AllPassFilter() noexcept  : bufferSize (0), bufferIndex (0) {}

        void setSize (const int size)
        {
            if (size != bufferSize)
            {
                bufferIndex = 0;
                buffer.reset (new DspFloatType [size]);
                bufferSize = size;
            }

            clear();
        }

        void clear() noexcept
        {
            memset (buffer.get(), 0, sizeof(DspFloatType) * (size_t)bufferSize);
        }

        DspFloatType process (const DspFloatType input) noexcept
        {
            const DspFloatType bufferedValue = buffer [bufferIndex];
            DspFloatType temp = input + (bufferedValue * 0.5f);
            // JUCE_UNDENORMALISE (temp);
            buffer [bufferIndex] = temp;
            bufferIndex = (bufferIndex + 1) % bufferSize;
            return bufferedValue - input;
        }

    private:
        std::unique_ptr<DspFloatType []> buffer;
        int bufferSize, bufferIndex;
    };

    class LinearSmoothedValue
    {
    public:
        LinearSmoothedValue() noexcept
            : currentValue (0), target (0), step (0), countdown (0), stepsToTarget (0)
        {}

        void reset (DspFloatType sampleRate, DspFloatType fadeLengthSeconds) noexcept
        {
            // jassert (sampleRate > 0 && fadeLengthSeconds >= 0);
            stepsToTarget = (int) std::floor (fadeLengthSeconds * sampleRate);
            currentValue = target;
            countdown = 0;
        }

        void setValue (DspFloatType newValue) noexcept
        {
            if (target != newValue)
            {
                target = newValue;
                countdown = stepsToTarget;

                if (countdown <= 0)
                    currentValue = target;
                else
                    step = (target - currentValue) / (DspFloatType) countdown;
            }
        }

        DspFloatType getNextValue() noexcept
        {
            if (countdown <= 0)
                return target;

            --countdown;
            currentValue += step;
            return currentValue;
        }

    private:
        DspFloatType currentValue, target, step;
        int countdown, stepsToTarget;
    };

    //==============================================================================
    enum { numCombs = 8, numAllPasses = 4, numChannels = 2 };
    bool enabledCombs [numCombs];
    bool enabledAllPasses [numAllPasses];

    Parameters parameters;
    DspFloatType gain;

    CombFilter comb [numChannels][numCombs];
    AllPassFilter allPass [numChannels][numAllPasses];

    LinearSmoothedValue damping, feedback, dryGain, wetGain1, wetGain2;
};