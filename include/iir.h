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
 *                                  iir.h
 *
 * Infinite Impulse Response.
 *
 * This class is designed for designing the IIR filter using bilinear
 * transform method. The unit of frequency and gain parameters are "Hz" and
 * "dB", respectively.
 *
 * The valid filter types are:
 *      lowpass         highpass    bandpass            bandstop
 * and the valid approximation method are:
 *      Butterworth     Chebyshev   Inverse Chebyshev   Elliptic
 *
 * The filter's coefficients consist of quadratic terms( the linear term can
 * be viewed as: a + b*z^-1 + 0*z^-2 ). The numerator and denominator
 * coefficients are stored in "aCoefs" and "bCoefs", respectively. The kth
 * quadratic term is as follow:
 *
 *             aCoefs[3k] + aCoefs[3k+1]*z^-1 + aCoefs[3k+2]*z^-2
 *           -------------------------------------------------------
 *             bCoefs[3k] + bCoefs[3k+1]*z^-1 + bCoefs[3k+2]*z^-2
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef IIR_H
#define IIR_H


#include <iomanip>
#include <cstdio>
#include <dfd.h>
#include <advmath.h>


namespace splab
{

    using std::ios;
    using std::setw;
    using std::setiosflags;
    using std::setprecision;

    class IIR : public DFD
    {

     public:

        IIR( const string &select, const string &appMethod );
        ~IIR();

        void    design();
        void    dispInfo() const;
        Vector<double> getNumCoefs() const;
        Vector<double> getDenCoefs() const;

    private:

        void    freqWrap();
        bool    orderEst();
        void    normCoefs();
        void    unnormCoefs();
        void    bilinTran();
        void    freqUnwrap();

        void    butterCoefs();
        void    chebyCoefs();
        void    invChebyCoefs();
        void    elliptCoefs();

        void    unnormLP( double freq );
        void    unnormHP( double freq );
        void    unnormBP( double BW, double Wo );
        void    unnormBS( double BW, double Wo );

        double  frqeResp( double freq );
        void    calcGain();

        string  approx;             // approximation method
        double  gain;               // overall gain

        Vector<double>  aCoefs,     // numerator
                        bCoefs,     // denominator
                        edgeGain;   // gain at edge frquency

    };
    // class IIR


    #include <iir-impl.h>

}
// namespace splab


#endif
// IIR_H
