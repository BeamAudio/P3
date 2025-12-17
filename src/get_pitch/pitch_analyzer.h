/// @file

#ifndef PITCH_ANALYZER_H
#define PITCH_ANALYZER_H

#include <vector>
#include <algorithm>

namespace upc {

  const float MIN_F0 = 50.0F;   ///< Minimum value of pitch in Hertz
  const float MAX_F0 = 500.0F;  ///< Maximum value of pitch in Hertz

  ///
  /// PitchAnalyzer: class that computes the pitch (in Hz) from a signal frame.
  /// No pre-processing or post-processing has been included
  ///
  class PitchAnalyzer {
  public:
    /// Window type
    enum Window {
      RECT,    ///< Rectangular window
      HAMMING  ///< Hamming window
    };

    void set_window(Window type); ///< pre-compute window

  private:
    std::vector<float> window; ///< precomputed window

    unsigned int frameLen;     ///< length of frame (in samples)
    unsigned int samplingFreq; ///< sampling rate (in samples per second)
    unsigned int npitch_min;   ///< minimum value of pitch period, in samples
    unsigned int npitch_max;   ///< maximum value of pitch period, in samples

    float noiseFloordB;
    bool  prevState;   ///< true = unvoiced, false = voiced
    float prevZcr, prevR1, prevRMax, prevPot;

    // opcional: suavizados (por si los usas)
    float smoothedPot   = -1e9f;
    float smoothedRmax  = 0.0f;
    float smoothedZcr   = 0.0f;

    // Umbral para central clipping (fracción del máximo)
    const float centralClippingThreshold = 0.7f;

    ///
    /// Computes correlation from lag=0 to r.size()
    ///
    void autocorrelation(const std::vector<float> &x, std::vector<float> &r) const;

    ///
    /// Returns the pitch (in Hz) of input frame x
    ///
    float compute_pitch(std::vector<float> & x);

    ///
    /// Returns normalized frame ZCR
    ///
    float compute_zcr(std::vector<float> &x);
	
    ///
    /// Returns true if the frame is unvoiced
    ///
    bool unvoiced(float pot, float r1norm, float rmaxnorm, float zcrnorm);

  public:
    PitchAnalyzer(unsigned int fLen,             ///< Frame length in samples
                  unsigned int sFreq,            ///< Sampling rate in Hertz
                  Window w = PitchAnalyzer::HAMMING, ///< Window type
                  float min_F0 = MIN_F0,         ///< Min pitch (Hz)
                  float max_F0 = MAX_F0)         ///< Max pitch (Hz)
    {
      frameLen     = fLen;
      samplingFreq = sFreq;

      set_f0_range(min_F0, max_F0);
      set_window(w);

      // Inicialización estado V/UV
      noiseFloordB = -1e9f;
      prevState    = true;   // empezamos como unvoiced
      prevZcr      = 0.0f;
      prevR1       = 0.0f;
      prevRMax     = 0.0f;
      prevPot      = -1e9f;
    }

    ///
    /// Operator (): computes the pitch for the given vector x
    ///
    float operator()(const std::vector<float> & _x) {
      if (_x.size() != frameLen)
        return -1.0F;

      std::vector<float> x(_x); //local copy of input frame
      return compute_pitch(x);
    }

    ///
    /// Operator (): computes the pitch for the given "C" vector (float *).
    /// N is the size of the vector pointed by pt.
    ///
    float operator()(const float * pt, unsigned int N) {
      if (N != frameLen)
        return -1.0F;

      std::vector<float> x(N); //local copy of input frame, size N
      std::copy(pt, pt+N, x.begin()); ///copy input values into local vector x
      return compute_pitch(x);
    }

    ///
    /// Operator (): computes the pitch for the given vector, expressed
    /// by the begin and end iterators
    ///
    float operator()(std::vector<float>::const_iterator begin,
                     std::vector<float>::const_iterator end) {

      if (static_cast<unsigned int>(end - begin) != frameLen)
        return -1.0F;

      std::vector<float> x(end - begin); //local copy of input frame, size N
      std::copy(begin, end, x.begin());  //copy input values into local vector x
      return compute_pitch(x);
    }
    
    ///
    /// Sets pitch range: takes min_F0 and max_F0 in Hz,
    /// sets npitch_min and npitch_max in samples
    ///
    void set_f0_range(float min_F0, float max_F0);

   
    
  };

} // namespace upc

#endif // PITCH_ANALYZER_H
