/*
 * Copyright (c) 2008-2011 Zhang Ming (M. Zhang), zmjerry@163.com
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 2 or any later version.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details. A copy of the GNU General Public License is available at:
 * http://www.fsf.org/licensing/licenses
 */


/*****************************************************************************
 *                                  fir.h
 *
 * Finite Impulse Response Digital Filter.
 *
 * This class is designed for designing the FIR filter using window function
 * method. The unit of frequency and gain parameters are "Hz" and "dB",
 * respectively.
 *
 * The valid filter types are:
 *      lowpass     highpass    bandpass    bandstop
 * and the valid windows are:
 *      Rectangle   Bartlett    Blackman    Hanning     Hamming
 *      Gauss       Kaiser
 *
 * The length of filter( filter's order plus one ) are multiples of 4, which
 * is the least number(L=4n) satisfying the design specifies.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef FIR_H
#define FIR_H


#include <iomanip>
#include <window.h>
#include <dfd.h>


namespace splab
{

    using std::setw;
    using std::ios;
    using std::setiosflags;
    using std::setprecision;


    class FIR : public DFD
    {

     public:

        FIR( const string &select, const string &win );
        FIR( const string &select, const string &win, double a );
        ~FIR();

        void    design();
        void    dispInfo() const;
        Vector<double> getCoefs() const;

    private:

        void    orderEst();
        void    idealCoef();
        void    calcCoef();
        double  frqeResp( double freq );
        void    calcGain();
        bool    isSatisfy();

        string  windType;           // window type

        Vector<double>  wind,       // window function
                        coefs,      // coefficients
                        edgeGain;   // gain at edge frquency
        double  alpha;              // window parameter

    };
    // class FIR


    #include <fir-impl.h>

}
// namespace splab


#endif
// FIR_H
