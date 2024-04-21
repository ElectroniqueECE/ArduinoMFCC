/*
	MFCC library
	Copyright (C) 2023 Foued DERRAZ

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

  updated from https://github.com/FouedDrz/arduinoMFCC by Raphaël Jeantet

*/
#include<Arduino.h>
#include<arduinoFFT.h>
#include "arduinoMFCC.h"
#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif
#if !defined(__SAM3X8E__)
#error "Ce code est conçu pour être utilisé uniquement avec Arduino Due. Veuillez utiliser un Arduino Due pour le projet Neurospeech."
#endif

arduinoMFCC::arduinoMFCC(uint8_t mfcc_size,uint8_t dct_mfcc_size, uint16_t frame_size, float samplerate) {
    _mfcc_size = mfcc_size;
    _frame_size = frame_size;
    _dct_mfcc_size = dct_mfcc_size;
    _samplerate = samplerate;

    _frame = (float*)malloc(_frame_size * sizeof(float));

    _hamming_window = (float*)malloc(_frame_size * sizeof(float));

    _mel_filter_bank  = (float **)malloc(_mfcc_size * sizeof(float *));
      for (uint8_t i = 0; i < _mfcc_size; i++) {
        _mel_filter_bank[i] = (float *)malloc((_frame_size/2) * sizeof(float));
    }

    _dct_matrix  = (float **)malloc(_dct_mfcc_size * sizeof(float *));
        for (uint8_t i = 0; i < _dct_mfcc_size; i++) {
             _dct_matrix[i] = (float *)malloc(_mfcc_size * sizeof(float));//
             }      

    _dct_mfcc_coeffs = (float*)malloc(_dct_mfcc_size * sizeof(float)); 
    _mfcc_coeffs = (float*)malloc(_mfcc_size* sizeof(float)); 	
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

arduinoMFCC::~arduinoMFCC() {
    free(_frame);
    free(_hamming_window);
    for (int i = 0; i < _mfcc_size; i++) {
          free(_mel_filter_bank[i]);
          }
    free(_mel_filter_bank);
    //
    for (int i = 0; i < _dct_mfcc_size; i++) {
          free(_dct_matrix[i]);
          }
    free(_dct_matrix);
    //
    free(_dct_mfcc_coeffs);
    free(_mfcc_coeffs);

}
/////////////////////////////////////////////////////////////////////////////
void arduinoMFCC::pre_emphasis()
{
// pre-emphasis
    for (uint16_t j = 1; j < _frame_size; j++)
          _frame[j] = _frame[j] - 0.95 * _frame[j-1];
  
}       
/////////////////////////////////////////////////////////////////////////
void arduinoMFCC::compute(float* frame, float* mfcc_coeffs) {
    for(int i = 0; i < _frame_size; i++) {
        _frame[i] = frame[i];
    }
    
    pre_emphasis();
    apply_hamming_window();
    apply_fft();
    apply_mel_filter_bank();
    for(int i = 0; i < _mfcc_size; i++) {
        mfcc_coeffs[i] = _mfcc_coeffs[i];
    }
}

void arduinoMFCC::computeWithDCT(float* frame, float* mfcc_coeffs) {
    for(int i = 0; i < _frame_size; i++) {
        _frame[i] = frame[i];
    }
    
    pre_emphasis();
    apply_hamming_window();
    apply_fft();
    apply_mel_filter_bank();
    apply_dct();

    for(int i = 0; i < _dct_mfcc_size; i++) {
        mfcc_coeffs[i] = _dct_mfcc_coeffs[i];
    }
}
///////////////////////////////////////////////////////////////////////////
// Fonction publique pour créer la fenêtre de Hamming
void arduinoMFCC::create_hamming_window() {
    // ... code pour créer la fenêtre de Hamming ...
	// Fonction publique pour créer la fenêtre de Hamming
    for (uint16_t  i = 0; i < _frame_size; i++) {
        _hamming_window[i] = 0.54 - 0.46 * cos(2 * PI * i / (_frame_size - 1));
    }
}


///////////////////////////////////////////////////////////////////////////
// Fonction publique pour créer la fenêtre de Hanning
void arduinoMFCC::create_hanning_window() {
    // ... code pour créer la fenêtre de Hanning ...
	// Fonction publique pour créer la fenêtre de Hanning
    for (uint16_t  i = 0; i < _frame_size; i++) {
        _hanning_window[i] = 0.5 - 0.5 * cos(2 * PI * i / (_frame_size - 1));
    }
}
///////////////////////////////////////////////////////////////////////////
// Fonction publique pour créer les filtres de Mel
void arduinoMFCC::create_mel_filter_bank()
 {
    float f_low          = 300.;  
    float f_high         = _samplerate; // Nyquist Frequency
    float mel_low_freq   = 2595. * log10f(1 + (f_low / 2.) / 700.);
    float mel_high_freq  = 2595. * log10f(1 + (f_high / 2.) / 700.);
    float* mel_f         = (float*)malloc((_mfcc_size + 2) * sizeof(float));
    float* hzPoints      = (float*)malloc((_mfcc_size + 2) * sizeof(float)); // Corresponding Hz scale points
  // Calculate Mel and Hz scale points
    float mel_freq_delta = (mel_high_freq - mel_low_freq) / (_mfcc_size + 1);
  for (uint8_t i = 0; i < _mfcc_size + 2; i++) {
    mel_f[i] = mel_low_freq + i * mel_freq_delta;
    hzPoints[i] = 700.0 * (powf(10, mel_f[i] / 2595.0) - 1);
  }
  // Create the filter bank
  for (uint16_t i = 0; i < _mfcc_size; i++) {
    for (uint16_t j = 0; j < _frame_size/2; j++) {
      float freq = (float)j * (_samplerate / 2) / (_frame_size / 2);
      if (freq  < hzPoints[i])
        _mel_filter_bank[i][j] = 0;
      else if ( freq >= hzPoints[i] && freq < hzPoints[i+1])
        _mel_filter_bank[i][j] = (freq - hzPoints[i]) / (hzPoints[i+1] - hzPoints[i]);
      else if ( freq >= hzPoints[i+1] && freq <= hzPoints[i+2])
        _mel_filter_bank[i][j] = (hzPoints[i+2] - freq) / (hzPoints[i+2] - hzPoints[i+1]);
      else
        _mel_filter_bank[i][j] = 0;
    
    }
    
  }
  free(mel_f);
  free(hzPoints);
}
//////////////////////////////////////////////////////////////////////////////////
// Fonction publique pour créer la matrice de transformée de cosinus discrète (DCT)
// Function is OK
void arduinoMFCC::create_dct_matrix() {
    // ... code pour créer la matrice DCT ...
	float sqrt_2_over_n = sqrt(2.0 / _mfcc_size);
	for (uint8_t  i = 0; i < _dct_mfcc_size; i++) {
    for (uint8_t  j = 0; j < _mfcc_size; j++) {

        _dct_matrix[i][j] = cos((PI * i * (j + 0.5)) / _dct_mfcc_size);
}
  }
}
/////////////////////////////////////////////////////////////////////////////
// Fonction publique pour appliquer la fenêtre de Hamming au signal audio
void arduinoMFCC::apply_hamming_window() {
    for (uint16_t  n = 0; n < _frame_size; n++) {
        _frame[n] = _frame[n] * _hamming_window[n];
    }
}
/////////////////////////////////////////////////////////////////////////////
// Fonction publique pour appliquer la fenêtre de Hamming au signal audio
void arduinoMFCC::apply_hanning_window() {
    for (uint16_t  n = 0; n < _frame_size; n++) {
        _frame[n] = _frame[n] * _hanning_window[n];
    }
}
/////////////////////////////////////////////////////////////////////////////
// Fonction publique pour appliquer les filtres de Mel au signal audio
void arduinoMFCC::apply_mel_filter_bank() {
    for (uint16_t  i = 0; i < _mfcc_size; i++) {
        float output = 0.0f;
        for (uint16_t  j = 0; j < _frame_size/2; j++) {
            // Serial.print(_mel_filter_bank[i][j]);
            // Serial.print(" ");
            output += _frame[j] * _mel_filter_bank[i][j];
        }
        _mfcc_coeffs[i] = log10f(output);
        // Serial.println();
    }
}
/////////////////////////////////////////////////////////////////////////////

void arduinoMFCC::apply_fft() {
  arduinoFFT myFFTframe;
  double *_vframe= (double*)malloc(_frame_size * sizeof(double));
   double *_rframe= (double*)malloc(_frame_size * sizeof(double));
  for(uint16_t i=0;i<_frame_size;i++){
    _rframe[i]=_frame[i];
    _vframe[i]=0.0;
  }
  myFFTframe= arduinoFFT(_rframe, _vframe,  _frame_size, (double)_samplerate);
  myFFTframe.DCRemoval();
  myFFTframe.Compute(FFT_FORWARD);
  myFFTframe.ComplexToMagnitude();
  for(uint16_t i=0;i<_frame_size;i++)  _frame[i]=(float)_rframe[i];

  free(_rframe);
  free(_vframe);
}

/////////////////////////////////////////////////////////////////////////////
// Fonction publique pour appliquer la transformée de cosinus discrète (DCT) au signal audio
void arduinoMFCC::apply_dct() {
      for (uint16_t  i = 0; i < _dct_mfcc_size; i++) {
        _dct_mfcc_coeffs[i] = 0.0;
        for (uint8_t  j = 0; j < _mfcc_size; j++) {
            _dct_mfcc_coeffs[i] += _mfcc_coeffs[j] *_dct_matrix[i][j];
        }
    }
}


