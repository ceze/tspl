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
 *                               fir-impl.h
 *
 * Implementation for FIR class.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
FIR::FIR( const string &select, const string &win )
        : DFD( select ), windType(win)
{
    bool cond = ( filtType=="lowpass" || filtType=="highpass" ||
                  filtType=="bandpass" || filtType=="bandstop" );
    assert(cond);

    cond = ( windType=="Rectangle" || windType=="Bartlett" ||
             windType=="Blackman" || windType=="Hanning" ||
             windType=="Hamming" || windType=="Gauss" );
    assert(cond);
}

FIR::FIR( const string &select, const string &win, double a )
        : DFD( select ), windType(win), alpha(a)
{
    bool cond = ( filtType=="lowpass" || filtType=="highpass" ||
                  filtType=="bandpass" || filtType=="bandstop" );
    assert(cond);

    cond = ( windType=="Kaiser" || windType=="Gauss" );
    assert(cond);
}

FIR::~FIR()
{
}


/**
 * Design digital FIR filter.
 */
void FIR::design()
{
    // Estimate the length of filter.
    orderEst();
    calcCoef();
    calcGain();

    while( !isSatisfy() )
    {
        order += 4;
        calcCoef();
        calcGain();
    }
}


/**
 * Displays the design result of filter.
 */
void FIR::dispInfo() const
{
    int i, j, k;

    cout << endl;
    cout << "\t\t    Filter selectivity      :  " << filtType << endl;
    cout << "\t\t    Window type             :  " << windType << endl;
    cout << "\t\t    Sampling Frequency (Hz) :  " << fsamp << endl;

    //  gains and edge frequency
    if( filtType=="lowpass" )
    {
        cout << "\t\t    Passband frequency (Hz) :  "
             << wpass1/TWOPI << endl;
        cout << "\t\t    Passband gain      (dB) :  "
             << apass1 << endl;
        cout << "\t\t    Stopband frequency (Hz) :  "
             << wstop1/TWOPI << endl;
        cout << "\t\t    Stopband gain      (dB) :  "
             << astop1 << endl;
    }
    else if( filtType=="highpass" )
    {
        cout << "\t\t    Stopband frequency (Hz) :  "
             << wstop1/TWOPI << endl;
        cout << "\t\t    Stopband gain      (dB) :  "
             << astop1 << endl;
        cout << "\t\t    Passband frequency (Hz) :  "
             << wpass1/TWOPI << endl;
        cout << "\t\t    Passband gain      (dB) :  "
             << apass1 << endl;
    }
    else if( filtType=="bandpass" )
    {
        cout << "\t\t    Lower stopband frequency (Hz) :  "
             << wstop1/TWOPI << endl;
        cout << "\t\t    Lower stopband gain      (dB) :  "
             << astop1 << endl;
        cout << "\t\t    Lower passband frequency (Hz) :  "
             << wpass1/TWOPI << endl;
        cout << "\t\t    Upper passband frequency (Hz) :  "
             << wpass2/TWOPI << endl;
        cout << "\t\t    Passband gain            (dB) :  "
             << apass1 << endl;
        cout << "\t\t    Upper stopband frequency (Hz) :  "
             << wstop2/TWOPI << endl;
        cout << "\t\t    Upper stopband gain      (dB) :  "
             << astop2 << endl;
    }
    else
    {
        cout << "\t\t    Lower passband frequency (Hz) :  "
             << wpass1/TWOPI << endl;
        cout << "\t\t    Lower passband gain      (dB) :  "
             << apass1 << endl;
        cout << "\t\t    Lower stopband frequency (Hz) :  "
             << wstop1/TWOPI << endl;
        cout << "\t\t    Upper stopband frequency (Hz) :  "
             << wstop2/TWOPI << endl;
        cout << "\t\t    Stopband gain            (dB) :  "
             << astop1 << endl;
        cout << "\t\t    Upper passband frequency (Hz) :  "
             << wpass2/TWOPI << endl;
        cout << "\t\t    Upper passband gain      (dB) :  "
             << apass2 << endl;
    }

    // display coefficients
    cout << endl << endl;
    cout << "\t\t\t\t  Filter Coefficients" << endl << endl;
    cout << "   N    [     N + 0            N + 1" ;
    cout << "            N + 2            N + 3     ]" << endl;
    cout << "  ===   ============================" ;
    cout << "========================================" ;
    for( i=0; i<order/4; ++i )
    {
        j = i * 4;
        cout << endl << setw(4) << j << "    ";

        cout << setiosflags(ios::scientific) << setprecision(8);
        for( k=0; k<4; ++k )
            cout << setw(16) << coefs[j+k] << " " ;
    }

    cout << endl << endl << endl
         << "\t ==================== Edge Frequency Response";
    cout << " ====================" << endl;
    cout << setiosflags(ios::fixed);
	if( edgeGain.size() == 2 )
    {
        cout << "\t     Mag(fp) = " << edgeGain[0] << "(dB)";
        cout << "       Mag(fs) = " << edgeGain[1] << "(dB)" << endl;
    }
    else
    {
        cout << "\t     Mag(fp1) = " << edgeGain[0] << "(dB)";
        cout << "       Mag(fp2) = " << edgeGain[1] << "(dB)" << endl;
        cout << "\t     Mag(fs1) = " << edgeGain[2] << "(dB)";
        cout << "       Mag(fs2) = " << edgeGain[3] << "(dB)" << endl;
    }
}


/**
 * Get the FIR filter's coefficients.
 */
inline Vector<double> FIR::getCoefs() const
{
    return coefs;
}


/**
 * Estimates the length of FIR filter.
 */
void FIR::orderEst()
{
    double  deltaFreq,  // delta frequency
            lowerDF,    // delta freq of lower band
            upperDF,    // delta freq of upper band
            errSB,      // stopband error
            errPB,      // passband error
            errMin,     // minimum of sb/pb errors
            errDB;      // minimum error in dB

    // Determine frequency delta
    if( filtType=="lowpass" || filtType=="highpass" )
        deltaFreq = abs(wstop1-wpass1) / fsamp;
    else
    {
        lowerDF = abs(wstop1-wpass1) / fsamp;
        upperDF = abs(wstop2-wpass2) / fsamp;

        if( lowerDF > upperDF )
            deltaFreq = upperDF;
        else
            deltaFreq = lowerDF;
    }

    // Determine stopband and passband errors
    errPB = 1 - pow( 10, 0.05*apass1 );
    errSB = pow( 10, 0.05*astop1 );

    if( errSB < errPB )
        errMin = errSB;
    else
        errMin = errPB;
    errDB = -20 * log10(errMin);

    // Store filter length in pFilt and return beta.
    if( errDB > 21 )
        order = int( ceil( 1 + (errDB-7.95) / (2.285*deltaFreq) ) );
    else
        order = int( ceil( 1 + (5.794/deltaFreq) ) );

    order += 4 - order%4;
}


/**
 * Calculates ideal FIR coefficients.
 */
void FIR::idealCoef()
{
    int     i;
    double  t,
            tau,
            Wc1,
            Wc2;

    // Calculate tau as non-integer if order is even.
    tau =  ( order - 1 ) / 2.0;

    // Adjust cutoff frequency to midway point between stop
    // and pass edge frequency and convert to digital frequency.
    Wc1 = (wstop1 + wpass1) / (2*fsamp);
    Wc2 = (wstop2 + wpass2) / (2*fsamp);

    // Calc coefs based on selectivity of filter.
    if( filtType == "lowpass" )
        for( i=0; i<order; ++i )
        {
            if( i == tau )
                coefs[i] = Wc1/PI;
            else
            {
                t = i - tau;
                coefs[i] = sin(Wc1*t) / (PI*t);
            }
        }
    else if( filtType == "highpass" )
        for( i=0; i<order; ++i )
        {
            if( i == tau )
                coefs[i] = (PI-Wc1) / PI;
            else
            {
                t = i - tau;
                coefs[i] = ( sin(PI*t) - sin(Wc1*t) ) / (PI*t);
            }
        }
    else if( filtType == "bandpass" )
        for( i=0; i<order; ++i )
        {
            if( i == tau )
                coefs[i] = ( Wc2-Wc1 ) / PI;
            else
            {
                t = i - tau;
                coefs[i] = ( sin(Wc2*t) - sin(Wc1*t) ) / (PI*t);
            }
        }
    else
        for( i=0; i<order; ++i )
        {
            if( i == tau )
                coefs[i] = ( PI+Wc1-Wc2 ) / PI;
            else
            {
                t = i - tau;
                coefs[i] = ( sin(PI*t)-sin(Wc2*t)+sin(Wc1*t) ) / (PI*t);
            }
        }
}


/**
 * Calculates the digital FIR coefficients.
 */
void FIR::calcCoef()
{
    coefs.resize(order);
    wind.resize(order);

    //  Calculate the ideal FIR coefficients.
    idealCoef();

    // Get window function.
    if( windType=="Kaiser" || windType=="Gauss" )
        wind = window( windType, order, alpha, 1.0 );
    else
        wind = window( windType, order, 1.0 );

    //  Multiply window and ideal coefficients.
    coefs *= wind;
}


/**
 * Calculate the response at edge freqency.
 */
double FIR::frqeResp( double freq )
{
    double  mag = 1.0,
            rea = 0.0,
            img = 0.0,
            omega = TWOPI*freq / fsamp;

    // Loop through all the coefficients.
    for( int i=0; i<order; ++i )
    {
        double  domega = i * omega;
        rea += coefs[i] * cos(domega);
        img += coefs[i] * sin(domega);
    }

    // Calculate final result and convert to degrees.
    mag = sqrt( rea*rea + img*img );
    if( mag < EPS )
        mag = EPS;
    mag = 20 * log10( mag );

    return mag;
}


/**
 * Calculate the response at edge freqency.
 */
void FIR::calcGain()
{
    // Determine the edge freqency.
    if( filtType=="lowpass" || filtType=="highpass" )
    {
        double  f1 = wpass1 / TWOPI,
                f2 = wstop1 / TWOPI;

        edgeGain.resize(2);
        edgeGain(1) = frqeResp( f1 );
        edgeGain(2) = frqeResp( f2 );
    }
    else
    {
        double  f1 = wpass1 / TWOPI,
                f2 = wpass2 / TWOPI,
                f3 = wstop1 / TWOPI,
                f4 = wstop2 / TWOPI;

        edgeGain.resize(4);
        edgeGain(1) = frqeResp( f1 );
        edgeGain(2) = frqeResp( f2 );
        edgeGain(3) = frqeResp( f3 );
        edgeGain(4) = frqeResp( f4 );
    }
}


/**
 * Design digital FIR filter.
 */
bool FIR::isSatisfy()
{
    if( edgeGain.size() == 2 )
    {
        if( edgeGain(1)<apass1 || edgeGain(2)>astop1 )
            return false;
    }
    else
    {
        if( edgeGain(1)<apass1 || edgeGain(2)<apass2 ||
            edgeGain(3)>astop1 || edgeGain(4)>astop2 )
            return false;
    }

    return true;
}
