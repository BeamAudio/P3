/// @file

#include <iostream>
#include <math.h>
#include "pitch_analyzer.h"

const float PI_M = 3.14169265f;

using namespace std;

/// Name space of UPC
namespace upc {
  void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {

    //Preprocessing Normalitzation
    vector<float> x_norm=x;
    
   
    for (unsigned int l = 0; l < r.size(); ++l) {
  	
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
    for (size_t i = 0; i < frameLen; i++) {
      window[i] = 0.54f - 0.46f * cos((2 * PI_M * i) / (frameLen - 1));
    }
    break;

  case RECT:
    window.assign(frameLen, 1.0f);
    break;

  default:
    window.assign(frameLen, 1.0f);
    break;
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

    // ========= 1) Normalización de autocorrelación =========
    // rmaxnorm = r[lag] / r[0]  → mide periodicidad (independiente del volumen)
    if (rmaxnorm < 0.0f) rmaxnorm = 0.0f;
    if (rmaxnorm > 1.0f) rmaxnorm = 1.0f;

    // ========= 2) Umbrales =========
    const float SNR_START  = 18.0f;
    const float SNR_KEEP   = -10.0f;

    const float RMAX_START = 0.42f;   // entrar voiced
    const float RMAX_KEEP  = 0.34f;   // mantener voiced

    const float ZCR_START  = 0.14f;
    const float ZCR_KEEP   = 0.55f;


    // ========= 3) Limitar potencia =========
    if (pot < -100.0f) pot = -100.0f;
    if (pot >   20.0f) pot =  20.0f;

    // ========= 4) Inicializar ruido =========
    if (noiseFloordB < -1e8f)
        noiseFloordB = pot;

    float snr = pot - noiseFloordB;

    bool isVoiced = !prevState;

    // ========= 5) Histéresis =========
    if (isVoiced) {
        // salir de voiced
        if (snr < SNR_KEEP || rmaxnorm < RMAX_KEEP || zcrnorm > ZCR_KEEP)
            isVoiced = false;
    } else {
        // entrar en voiced
        if (snr >= SNR_START && rmaxnorm >= RMAX_START && zcrnorm <= ZCR_START)
            isVoiced = true;
    }

    bool isUnvoiced = !isVoiced;

    // ========= 6) Actualizar ruido SOLO en unvoiced =========
    if (isUnvoiced) {
        noiseFloordB = (pot < noiseFloordB)
            ? 0.90f  * noiseFloordB + 0.10f  * pot
            : 0.999f * noiseFloordB + 0.001f * pot;
    }

    prevState = isUnvoiced;
    return isUnvoiced;   // true = unvoiced
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




 float PitchAnalyzer::compute_pitch(std::vector<float> &x) {
    if (x.size() != frameLen) return -1.0f;

    // 1) Ventana
    for (unsigned int i = 0; i < x.size(); ++i)
        x[i] *= window[i];

    // 2) Autocorrelación
    std::vector<float> r(npitch_max);
    autocorrelation(x, r);

    // 3) Buscar desde el primer cruce a negativo (para no coger el lóbulo del origen)
    unsigned int start = npitch_min;
    for (unsigned int i = 1; i < npitch_max; ++i) {
        if (r[i] < 0.0f) { start = std::max(start, i); break; }
    }

    // 4) Máximo global en [start, npitch_max)
    unsigned int lag = start;
    float rMax = r[start];
    for (unsigned int i = start; i < npitch_max; ++i) {
        if (r[i] > rMax) { rMax = r[i]; lag = i; }
    }

    // 5) Chequeo armónico simple: si 2*lag también es pico fuerte, quizá el fundamental es 2*lag
    if (2 * lag < npitch_max) {
        if (r[2 * lag] >= 0.90f * r[lag]) { // 0.90 ajustable
            lag = 2 * lag;
        }
    }

    // 6) Potencia + decisión voiced/unvoiced
    float pot = 10.0f * log10(r[0]);

    if (unvoiced(pot, r[1] / r[0], r[lag] / r[0], compute_zcr(x)))
        return 0.0f;

    return (float)samplingFreq / (float)lag;
}




  
    

   

    
      


 

}