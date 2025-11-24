/// @file

#include <iostream>
#include <math.h>
#include "pitch_analyzer.h"

const float PI_M = 3.14169265f;

using namespace std;

/// Name space of UPC
namespace upc {
  void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {

    for (unsigned int l = 0; l < r.size(); ++l) {
  		/**  \DONE Compute the autocorrelation r[l]
       \f[
       r[l] = \sum_{i=l}^{N-1} x[i]*x[i-l]
       \f]
      */
      r[l] = 0.0f;
      for (unsigned int i=l; i<x.size(); ++i) {
        r[l] += x[i]*x[i-l];
      }
    }

    if (r[0] == 0.0F) //to avoid log() and divide zero 
      r[0] = 1e-10; 
      
  }

  void PitchAnalyzer::set_window(Window win_type) {
    if (frameLen == 0)
      return;

    window.resize(frameLen);

    switch (win_type) {
    case HAMMING:
      window.assign(frameLen, 0);
      for(size_t i = 0; i<frameLen; i++) {
        window[i] = 0.54f - 0.46f*cos((2*PI_M*i)/(frameLen-1));
      }
      
      
    case RECT:
      window.assign(frameLen, 1);
    default:
      window.assign(frameLen, 1);
    }
  }

  void PitchAnalyzer::set_f0_range(float min_F0, float max_F0) {
    npitch_min = (unsigned int) samplingFreq/max_F0;
    if (npitch_min < 2)
      npitch_min = 2;  // samplingFreq/2

    npitch_max = 1 + (unsigned int) samplingFreq/min_F0;

    //frameLen should include at least 2*T0
    if (npitch_max > frameLen/2)
      npitch_max = frameLen/2;
  }


  bool PitchAnalyzer::unvoiced(float pot, float r1norm, float rmaxnorm, float zcrnorm) {
    /// \DONE Automatic Noise Floor Tracking + Hysteresis Logic

    // --- Constants ---
    // How much louder than the noise floor must the signal be? (e.g., 6dB or 10dB)
    const float kSignalToNoiseThreshold = 6.0f; 
    
    // Hysteresis thresholds (Strict to Start, Loose to Keep)
    const float kStart_RMax = 0.6f; 
    const float kStart_ZCR  = 0.2f;
    const float kKeep_RMax  = 0.4f; 
    const float kKeep_ZCR   = 0.4f;

    // --- 1. Automatic Noise Floor Tracking (Adaptive) ---
    
    // Initialize noiseFloor if it looks uninitialized (e.g., 0 or very low depending on your scale)
    // Assuming 'pot' is in dB (negative values), -100 is usually a safe start.
    if (this->noiseFloordB == 0.0f) this->noiseFloordB = pot;

    if (pot < this->noiseFloordB) {
        // Case A: Current signal is QUIETER than our tracked noise.
        // Action: Fast Decay. Adapt quickly to this new, lower noise level.
        // Formula: 90% old, 10% new
        this->noiseFloordB = 0.90f * this->noiseFloordB + 0.10f * pot;
    } else {
        // Case B: Current signal is LOUDER than noise (likely speech).
        // Action: Slow Attack. Rise very slowly to account for drifting background noise 
        // without capturing the speech energy as noise.
        // Formula: 99.9% old, 0.1% new
        this->noiseFloordB = 0.999f * this->noiseFloordB + 0.001f * pot;
    }

    // --- 2. The "Gate" Check ---
    
    // Calculate actual Signal-to-Noise Ratio
    float snr = pot - this->noiseFloordB;

    // If SNR is too low, it's automatically unvoiced/silence
    if (snr < kSignalToNoiseThreshold) {
        this->prevState = true; // True = Unvoiced/Silence
        return true;
    }

    // --- 3. Voicing Classification with Hysteresis ---
    
    // Determine if we are currently considering the signal voiced (inverted logic of prevState)
    bool isVoiced = !this->prevState;

    if (isVoiced) {
        // We are currently VOICED. Only stop if signal degrades significantly.
        // Condition to drop to Unvoiced: Low Correlation OR High Frequency Noise
        if (rmaxnorm < kKeep_RMax || zcrnorm > kKeep_ZCR) {
            isVoiced = false;
        }
    } else {
        // We are currently UNVOICED. Only start if signal is very clean.
        // Condition to start Voiced: High Correlation AND Low Frequency Noise
        if (rmaxnorm >= kStart_RMax && zcrnorm <= kStart_ZCR) {
            isVoiced = true;
        }
    }

    // --- 4. State Updates ---
    this->prevZcr = zcrnorm;
    this->prevR1 = r1norm;
    this->prevRMax = rmaxnorm;
    this->prevPot = pot;
    
    // Store current state (Remember: function returns true for unvoiced)
    this->prevState = !isVoiced;

    return !isVoiced; 
  }

  float PitchAnalyzer::compute_zcr(vector<float> &x) {
    if (x.size() != frameLen)
      return -1.0F;

    unsigned int abszcr = 0;
    
    for(size_t i=1; i<frameLen; i++) {
      if(x[i]*x[i-1]<=0) {
        abszcr++;
      }
    }

    return (float)((float)abszcr/(float)frameLen); 
  }

  float PitchAnalyzer::compute_pitch(vector<float> & x) {
    if (x.size() != frameLen)
      return -1.0F;

    //Window input frame
    for (unsigned int i=0; i<x.size(); ++i)
      x[i] *= window[i];

    vector<float> r(npitch_max);

    //Compute correlation
    autocorrelation(x, r);

    unsigned int lag = 0;
    float rMax = 0.0f;

    for(size_t i = npitch_min; i<npitch_max; i++) {
      if(r[i] > rMax) {
        rMax = r[i];
        lag = i;
      }
    }

    

    float pot = 10 * log10(r[0]);

 
#if 0
    if (r[0] > 0.0F)
      cout << pot << '\t' << r[1]/r[0] << '\t' << r[lag]/r[0] << endl;
#endif
    
    if (unvoiced(pot, r[1]/r[0], r[lag]/r[0], compute_zcr(x)))
      return 0;
    else
      return (float) samplingFreq/(float) lag;
  }
}
