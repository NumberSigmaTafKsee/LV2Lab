/*
==============================================================================
EnvelopeDetector.h
Created: 16 Dec 2017 9:49:56pm
Author:  Toby
==============================================================================
*/
#pragma once
// turn off non-critical warnings
#pragma warning(disable : 4244)//DspFloatType to DspFloatType
#pragma warning(disable : 4996)//strncpy
#pragma warning(disable : 4305)//DspFloatType DspFloatType truncation

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>

// For WIN vs MacOS
// XCode requires these be defined for compatibility
#if defined _WINDOWS || defined _WINDLL
#include <windows.h>
#else // MacOS
typedef unsigned int        UINT;
typedef unsigned long       DWORD;
typedef unsigned char		UCHAR;
typedef unsigned char       BYTE;
#endif

#include <cmath>



// constants for dealing with overflow or underflow
#define FLT_EPSILON_PLUS      1.192092896e-07         /* smallest such that 1.0+FLT_EPSILON != 1.0 */
#define FLT_EPSILON_MINUS    -1.192092896e-07         /* smallest such that 1.0-FLT_EPSILON != 1.0 */
#define FLT_MIN_PLUS          1.175494351e-38         /* std::min positive value */
#define FLT_MIN_MINUS        -1.175494351e-38         /* std::min negative value */

const DspFloatType DIGITAL_TC = -2.0; // log(1%)
const DspFloatType ANALOG_TC = -0.43533393574791066201247090699309; // (log(36.7%)
const DspFloatType METER_UPDATE_INTERVAL_MSEC = 15.0;
const DspFloatType METER_MIN_DB = -60.0;

namespace FX::Filters
{
	class CEnvelopeDetector
	{
	public:
		CEnvelopeDetector(void);
		~CEnvelopeDetector(void);

	public:

		// Call the Init Function to initialize and setup all at once; this can be called as many times
		// as you want
		void init(DspFloatType samplerate, DspFloatType attack_in_ms, DspFloatType release_in_ms, bool bAnalogTC, UINT uDetect, bool bLogDetector);

		// these functions allow you to change modes and attack/release one at a time during
		// realtime operation
		void setTCModeAnalog(bool bAnalogTC); // {m_bAnalogTC = bAnalogTC;}

											// THEN do these after init
		void setAttackTime(DspFloatType attack_in_ms);
		void setReleaseTime(DspFloatType release_in_ms);

		// Use these "codes"
		// DETECT PEAK   = 0
		// DETECT MS	 = 1
		// DETECT RMS	 = 2
		//
		void setDetectMode(UINT uDetect) { m_uDetectMode = uDetect; }

		void setSampleRate(DspFloatType f) { m_fSampleRate = f; }

		void setLogDetect(bool b) { m_bLogDetector = b; }

		// call this to detect; it returns the peak ms or rms value at that instant
		DspFloatType detect(DspFloatType fInput);

		// call this from your prepareForPlay() function each time to reset the detector
		void prepareForPlay();

	protected:
		int  m_nSample;
		DspFloatType m_fAttackTime;
		DspFloatType m_fReleaseTime;
		DspFloatType m_fAttackTime_mSec;
		DspFloatType m_fReleaseTime_mSec;
		DspFloatType m_fSampleRate;
		DspFloatType m_fEnvelope;
		UINT  m_uDetectMode;
		bool  m_bAnalogTC;
		bool  m_bLogDetector;
	};

	inline DspFloatType lagrpol(DspFloatType* x, DspFloatType* y, int n, DspFloatType xbar)
	{
		int i, j;
		DspFloatType fx = 0.0;
		DspFloatType l = 1.0;
		for (i = 0; i<n; i++)
		{
			l = 1.0;
			for (j = 0; j<n; j++)
			{
				if (j != i)
					l *= (xbar - x[j]) / (x[i] - x[j]);
			}
			fx += l * y[i];
		}
		return (fx);
	}


	CEnvelopeDetector::CEnvelopeDetector(void)
	{
		m_fAttackTime_mSec = 0.0;
		m_fReleaseTime_mSec = 0.0;
		m_fAttackTime = 0.0;
		m_fReleaseTime = 0.0;
		m_fSampleRate = 44100;
		m_fEnvelope = 0.0;
		m_uDetectMode = 0;
		m_nSample = 0;
		m_bAnalogTC = false;
		m_bLogDetector = false;
	}

	CEnvelopeDetector::~CEnvelopeDetector(void)
	{
	}

	void CEnvelopeDetector::prepareForPlay()
	{
		m_fEnvelope = 0.0;
		m_nSample = 0;
	}

	void CEnvelopeDetector::init(DspFloatType samplerate, DspFloatType attack_in_ms, DspFloatType release_in_ms, bool bAnalogTC, UINT uDetect, bool bLogDetector)
	{
		m_fEnvelope = 0.0;
		m_fSampleRate = samplerate;
		m_bAnalogTC = bAnalogTC;
		m_fAttackTime_mSec = attack_in_ms;
		m_fReleaseTime_mSec = release_in_ms;
		m_uDetectMode = uDetect;
		m_bLogDetector = bLogDetector;

		// set themm_uDetectMode = uDetect;
		setAttackTime(attack_in_ms);
		setReleaseTime(release_in_ms);
	}

	void CEnvelopeDetector::setAttackTime(DspFloatType attack_in_ms)
	{
		m_fAttackTime_mSec = attack_in_ms;

		if (m_bAnalogTC)
			m_fAttackTime = exp(ANALOG_TC / (attack_in_ms * m_fSampleRate * 0.001));
		else
			m_fAttackTime = exp(DIGITAL_TC / (attack_in_ms * m_fSampleRate * 0.001));
	}

	void CEnvelopeDetector::setReleaseTime(DspFloatType release_in_ms)
	{
		m_fReleaseTime_mSec = release_in_ms;

		if (m_bAnalogTC)
			m_fReleaseTime = exp(ANALOG_TC / (release_in_ms * m_fSampleRate * 0.001));
		else
			m_fReleaseTime = exp(DIGITAL_TC / (release_in_ms * m_fSampleRate * 0.001));
	}

	void CEnvelopeDetector::setTCModeAnalog(bool bAnalogTC)
	{
		m_bAnalogTC = bAnalogTC;
		setAttackTime(m_fAttackTime_mSec);
		setReleaseTime(m_fReleaseTime_mSec);
	}


	DspFloatType CEnvelopeDetector::detect(DspFloatType fInput)
	{
		switch (m_uDetectMode)
		{
		case 0:
			fInput = fabs(fInput);
			break;
		case 1:
			fInput = fabs(fInput) * fabs(fInput);
			break;
		case 2:
			fInput = pow((DspFloatType)fabs(fInput) * (DspFloatType)fabs(fInput), (DspFloatType)0.5);
			break;
		default:
			fInput = (DspFloatType)fabs(fInput);
			break;
		}

		DspFloatType fOld = m_fEnvelope;
		if (fInput> m_fEnvelope)
			m_fEnvelope = m_fAttackTime * (m_fEnvelope - fInput) + fInput;
		else
			m_fEnvelope = m_fReleaseTime * (m_fEnvelope - fInput) + fInput;

		if (m_fEnvelope > 0.0 && m_fEnvelope < FLT_MIN_PLUS) m_fEnvelope = 0;
		if (m_fEnvelope < 0.0 && m_fEnvelope > FLT_MIN_MINUS) m_fEnvelope = 0;

		// bound them; can happen when using pre-detector gains of more than 1.0
		m_fEnvelope = std::min((DspFloatType)m_fEnvelope, (DspFloatType)1.0);
		m_fEnvelope = std::max((DspFloatType)m_fEnvelope, (DspFloatType)0.0);

		//16-bit scaling!
		if (m_bLogDetector)
		{
			if (m_fEnvelope <= 0)
				return -96.0; // 16 bit noise floor

			return 20 * log10(m_fEnvelope);
		}

		return m_fEnvelope;
	}
}