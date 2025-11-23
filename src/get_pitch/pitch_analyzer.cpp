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
        window[i] = 0.54f - 0.46*cos((2*PI_M*i)/(frameLen-1));
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
    /// \DONE Implement a rule to decide whether the sound is voiced or not.
    /// * You can use the standard features (pot, r1norm, rmaxnorm),
    ///   or compute and use other ones.

    

    if(r1norm >= 0.5 && rmaxnorm >= 0.4 && pot>=noiseFloordB && zcrnorm<=0.25f) {
      
      return false;
    }
    this->noiseFloordB = 1.2*pot;
    return true;
    
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
