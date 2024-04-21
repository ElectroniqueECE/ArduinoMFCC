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
#ifndef arduinoMFCC_H
#define arduinoMFCC_H

#include <Arduino.h>
#include <math.h>

class arduinoMFCC {
    public:
        // Constructeur
        arduinoMFCC();
        arduinoMFCC(uint8_t mfcc_size,uint8_t dct_mfcc_size, uint16_t frame_size, float samplerate);
        ~arduinoMFCC();

        // Fonctions publiques
        void compute(float* frame, float* mfcc_coeffs);
        void computeWithDCT(float* frame, float* mfcc_coeffs);

        //
        void pre_emphasis();
        void apply_mel_filter_bank();    

        void create_mel_filter_bank();

        void apply_fft();
        void apply_dct();

        void create_hamming_window();
        void create_hanning_window();
        
        void apply_hamming_window();
        void apply_hanning_window();

        void create_dct_matrix();
        // 
        uint16_t _samplerate;
        float* _frame;

        float** _mel_filter_bank;
        uint16_t  _frame_size;

        float* _mfcc_coeffs;
        float *_dct_mfcc_coeffs;

        uint8_t _mfcc_size;
        uint8_t _dct_mfcc_size;

        float* _hamming_window;
        float* _hanning_window;
        

        float** _dct_matrix;

};

#endif
