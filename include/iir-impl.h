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
 *                               iir-impl.h
 *
 * Implementation for IIR class.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
IIR::IIR( const string &select, const string &app )
        : DFD( select ), approx(app)
{
    bool cond = ( filtType=="lowpass" || filtType=="highpass" ||
                  filtType=="bandpass" || filtType=="bandstop" );
    assert(cond);

    cond = ( approx=="Butterworth" || approx=="Chebyshev" ||
             approx=="InvChebyshev" || approx=="Elliptic" );
    assert(cond);
}

IIR::~IIR()
{
}


/**
 * Calculates the digital FIR coefficients.
 */
void IIR::design()
{
    freqWrap();
    if( orderEst() )
    {
        normCoefs();
        unnormCoefs();
        bilinTran();
        freqUnwrap();
        calcGain();
    }
    else
        cerr << "The filter's order exceeds exceed the range!" << endl;
}


/**
 * Displays the design result of filter.
 */
void IIR::dispInfo() const
{
    int i, j;

    cout << endl;
    cout << "\t\t    Filter selectivity      :  " << filtType << endl;
    cout << "\t\t    Approximation method    :  " << approx << endl;
    cout << "\t\t    Filter order            :  " << order << endl;
    cout << "\t\t    Overall gain            :  " << gain << endl;
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
    cout << endl << endl << "           Numerator Coefficients    ";
    cout << "            Denominator Coefficients" << endl << endl ;
    cout << "  N   [  1       z^-1           z^-2  ]";
    cout  << "     [  1       z^-1           z^-2  ]" << endl;
    cout << " ===  =================================";
    cout  << "     =================================" << endl ;
    for( i=0; i<(order+1)/2; ++i )
    {
        j = i * 3;
        cout << setw(3) << i+1 << " ";

        cout << setiosflags(ios::fixed );
        cout << setw(6) << setprecision(1) << aCoefs[j+0];
        cout << setw(14) << setprecision(8) << aCoefs[j+1];
        cout << setw(14) << setprecision(8) << aCoefs[j+2];
        cout << setw(10) << setprecision(1) << bCoefs[j+0];
        cout << setw(14) << setprecision(8) << bCoefs[j+1];
        cout << setw(14) << setprecision(8) << bCoefs[j+2] << endl ;
    }

    // display edge frequency response
    cout << endl << endl
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
inline Vector<double> IIR::getNumCoefs() const
{
    return aCoefs;
}

inline Vector<double> IIR::getDenCoefs() const
{
    return bCoefs;
}


/**
 * Pre-warps analog frequency for IIR filtrs.
 */
void IIR::freqWrap( )
{
    // Pre-warps critical frequency based on selectivity.
    if( filtType=="lowpass" || filtType=="highpass" )
    {
        wpass1 = 2*fsamp * tan( wpass1/(2*fsamp) );
        wstop1 = 2*fsamp * tan( wstop1/(2*fsamp) );
    }
    else
    {
        wpass1 = 2*fsamp * tan( wpass1/(2*fsamp) );
        wpass2 = 2*fsamp * tan( wpass2/(2*fsamp) );
        wstop1 = 2*fsamp * tan( wstop1/(2*fsamp) );
        wstop2 = 2*fsamp * tan( wstop2/(2*fsamp) );
    }
}


/**
 * Calculate the order of digital IIR filter.
 */
bool IIR::orderEst()
{
    double  kernel, ratio,
            n,
            rt, kn,
            CEIrt, CEIrtq,
            CEIkn, CEIknq,
            wp1, wp2, ws1, ws2;

    // Simplify the variables.
    wp1 = wpass1;
    wp2 = wpass2;
    ws1 = wstop1;
    ws2 = wstop2;

    // Calculate frequency ratio based on filter selectivity.
    if( filtType == "lowpass" )
        ratio = ws1 / wp1;
    else if( filtType == "highpass" )
        ratio = wp1 / ws1;
    else if( filtType == "bandpass" )
    {
        if( ws1 > (wp1*wp2)/ws2 )
        {
            ws2 = (wp1*wp2) / ws1;
            wstop2 = ws2;
        }
        else
        {
            ws1 = (wp1*wp2) / ws2;
            wstop1 = ws1;
        }
        ratio = (ws2-ws1) / (wp2-wp1);
    }
    else
    {
        if( wp1 > (ws1*ws2)/wp2 )
        {
            wp2 = (ws1*ws2) / wp1;
            wpass2 = wp2;
        }
        else
        {
            wp1 = (ws1*ws2) / wp2;
            wpass1 = wp1;
        }
        ratio = (wp2-wp1) / (ws2-ws1);
    }

    // Calculate kernel which is used in all calculations.
    kernel = ( pow(10.0,-0.1*astop1) - 1 )
           / ( pow(10.0,-0.1*apass1) - 1 );

    // Calculate order and round up to integer value.
    if( approx == "Butterworth" )
        n = log10(kernel) / (2*log10(ratio));
    else if( approx == "Chebyshev" || approx == "InvChebyshev" )
        n = acosh(sqrt(kernel)) / acosh(ratio);
    else
    {
        rt = 1 / ratio;
        if( rt > 0.9999 )
            return false;
        kn = 1 / sqrt(kernel);
        if( kn < 2e-8 )
            return false;

        CEIrt   = ellipticIntegral(rt);
        CEIrtq  = ellipticIntegral(sqrt(1-rt*rt));
        CEIkn   = ellipticIntegral(kn);
        CEIknq  = ellipticIntegral(sqrt(1-kn*kn));
        n = ( CEIrt*CEIknq ) / ( CEIrtq * CEIkn );
    }

    if( n > 200 )
        return false;

    order = int( std::ceil(n) );
    return true;
}


/**
 * Calculate normal filter coefficients.
 */
void IIR::normCoefs()
{
    // There are 3 coefficients for each quadratic.
    // First order factors are considered as quadratics.
    int coefsNum = 3 * ( (order+1) / 2 );
    aCoefs.resize(coefsNum);
    bCoefs.resize(coefsNum);

    // Calculate coefs based on approximation.
    if( approx == "Butterworth" )
        butterCoefs();
    else if( approx == "Chebyshev" )
        chebyCoefs();
    else if( approx == "InvChebyshev" )
        invChebyCoefs();
    else
        elliptCoefs();
}


/**
 * Converts normal lowpass coefficients to unnormalized LP/HP/BP/BS.
 */
void IIR::unnormCoefs()
{
    double  freq, BW, Wo;

    // Calculate frequency, Wo and BW based on approx method.
    if( approx == "InvChebyshev" )
    {
        freq = wstop1;
        Wo = sqrt( wstop1 * wstop2 );
        BW = wstop2 - wstop1;
    }
    else
    {
        freq = wpass1;
        Wo = sqrt( wpass1 * wpass2 );
        BW = wpass2 - wpass1;
    }

    // Call unnormal function based on selectivity.
    if( filtType == "lowpass" )
        unnormLP( freq );
    else if( filtType == "highpass" )
        unnormHP( freq );
    else if( filtType == "bandpass" )
        unnormBP( BW, Wo );
    else
        unnormBS( BW, Wo );
}


/**
 * Use bilinear transform to convert transfer function
 * from s-domain to z-domain.
 */
void IIR::bilinTran()
{
    int     i, j, start, numQuads;
    double  f2, f4,
            N0, N1, N2,
            D0, D1, D2;

    f2 = 2 * fsamp;
    f4 = f2 * f2;
    numQuads = ( order + 1 ) / 2;

    // Handle first order factor if present start
    // indicates if one quad factor handled.
    start = 0;
    if( order % 2 )
    {
        N0 = aCoefs[2] + aCoefs[1] * f2;
        N1 = aCoefs[2] - aCoefs[1] * f2;
        D0 = bCoefs[2] + bCoefs[1] * f2;
        D1 = bCoefs[2] - bCoefs[1] * f2;

        aCoefs[0] = 1.0;
        aCoefs[1] = N1 / N0;
        aCoefs[2] = 0.0;
        bCoefs[0] = 1.0;
        bCoefs[1] = D1 / D0;
        bCoefs[2] = 0.0;
        gain *= (N0 / D0);
        start = 1;
    }

    // Handle quadratic factors. j indexes quad factors
    // since each quad factor has three coefs.
    for( i=start; i<numQuads; ++i )
    {
        j = 3 * i;
        N0 = aCoefs[j]*f4 + aCoefs[j+1]*f2 + aCoefs[j+2];
        N1 = 2 * ( aCoefs[j+2] - aCoefs[j]*f4 );
        N2 = aCoefs[j]*f4 - aCoefs[j+1]*f2 + aCoefs[j+2];
        D0 = bCoefs[j]*f4 + bCoefs[j+1]*f2 + bCoefs[j+2];
        D1 = 2 * ( bCoefs[j+2] - bCoefs[j]*f4 );
        D2 = bCoefs[j]*f4 - bCoefs[j+1]*f2 + bCoefs[j+2];

        aCoefs[j] = 1.0;
        aCoefs[j+1] = N1 / N0;
        aCoefs[j+2] = N2 / N0;
        bCoefs[j] = 1.0;
        bCoefs[j+1] = D1 / D0;
        bCoefs[j+2] = D2 / D0;
        gain *= (N0 / D0);
    }
}


/**
 * Unwarps analog frequency for IIR filtrs.
 */
void IIR::freqUnwrap()
{
    if( filtType=="lowpass" || filtType=="highpass" )
    {
        wpass1 = 2*fsamp * atan( wpass1 / (2*fsamp) );
        wstop1 = 2*fsamp * atan( wstop1 / (2*fsamp) );
    }
    else
    {
        wpass1 = 2*fsamp * atan( wpass1 / (2*fsamp) );
        wpass2 = 2*fsamp * atan( wpass2 / (2*fsamp) );
        wstop1 = 2*fsamp * atan( wstop1 / (2*fsamp) );
        wstop2 = 2*fsamp * atan( wstop2 / (2*fsamp) );
    }
}


/**
 * Calculates normal Butterworth filter coefficients.
 */
void IIR::butterCoefs()
{
    int     m, a, b;
    double  R, epsilon, theta, sigma, omega;

    // Make calculations of necessary constants.
    epsilon = sqrt( pow(10.0,-0.1*apass1) - 1.0 );
    R = pow( epsilon, -1.0/order );

    // Initialize gain to 1.0, start indices at 0.
    gain = 1.0;
    a = 0;
    b = 0;

    // Handle odd order if necessary.
    if( order % 2 )
    {
        aCoefs[a++] = 0.0;
        aCoefs[a++] = 0.0;
        aCoefs[a++] = R;
        bCoefs[b++] = 0.0;
        bCoefs[b++] = 1.0;
        bCoefs[b++] = R;
    }

    // Handle all quadratic terms.
    for( m=0; m<order/2; ++m )
    {
        // Calculate angle first, then real and imag position.
        theta = PI*( 2*m + order+1 ) / ( 2*order );
        sigma = R * cos(theta);
        omega = R * sin(theta);

        // Set the quadratic coefficients.
        aCoefs[a++] = 0.0;
        aCoefs[a++] = 0.0;
        aCoefs[a++] = sigma*sigma + omega*omega;
        bCoefs[b++] = 1.0;
        bCoefs[b++] = -2*sigma;
        bCoefs[b++] = sigma*sigma + omega*omega;
    }
}


/**
 * Calculaets normal Chebyshev coefficients.
 */
void IIR::chebyCoefs()
{
    int     m, a, b;
    double  D, epsilon, phi, sigma, omega;

    // Make calculations of necessary constants.
    epsilon = sqrt( pow(10.0,-0.1*apass1) - 1.0 );
    D = asinh(1/epsilon) / order;

    // Handle odd order if necessary.
    a = 0;
    b = 0;

    if( order % 2 )
    {
        gain = 1.0;
        aCoefs[a++] = 0.0;
        aCoefs[a++] = 0.0;
        aCoefs[a++] = sinh(D);
        bCoefs[b++] = 0.0;
        bCoefs[b++] = 1.0;
        bCoefs[b++] = sinh(D);
    }
    else
        gain = pow( 10.0, 0.05*apass1 );

    // Handle all quadratic terms.
    for( m=0; m<order/2; m++ )
    {
        // Calculate angle first, then real and imag position.
        phi = PI * ( 2*m+1 ) / ( 2*order );
        sigma = -1 * sinh(D) * sin(phi);
        omega = cosh(D) * cos(phi);

        // Set the quadratic coefs.
        aCoefs[a++] = 0.0;
        aCoefs[a++] = 0.0;
        aCoefs[a++] = sigma*sigma + omega*omega;
        bCoefs[b++] = 1.0;
        bCoefs[b++] = -2 * sigma;
        bCoefs[b++] = sigma*sigma + omega*omega;
    }
}


/**
 * Calculates normal inverse Chebyshev coefficients.
 */
void IIR::invChebyCoefs()
{
    int     m, a, b;
    double  D, epsilon, mag2, phi, zero, sigma, omega;

    // Make calculations of necessary constants.
    epsilon = 1 / sqrt( pow(10.0,-0.1*astop1) - 1.0 );
    D = asinh( 1/epsilon ) / order;

    // Initialize gain to 1.0, start indices at 0.
    gain = 1.0;
    a = 0;
    b = 0;

    // Handle odd order if necessary.
    if( order % 2 )
    {
        aCoefs[a++] = 0.0;
        aCoefs[a++] = 0.0;
        aCoefs[a++] = 1 / sinh(D);
        bCoefs[b++] = 0.0;
        bCoefs[b++] = 1.0;
        bCoefs[b++] = 1 / sinh(D);
    }

    // Handle all quadratic terms.
    for( m=0; m<order/2; m++ )
    {
        // Calculate angle first, then real and imag position.
        phi = PI * ( 2*m+1) / ( 2*order );
        sigma = -1 * sinh(D) * sin(phi);
        omega = cosh(D) * cos(phi);

        // Calculate inverse of pole locations and zero location.
        mag2 = omega*omega + sigma*sigma;
        omega = -1 * omega / mag2;
        sigma = sigma / mag2;
        zero = 1. / cos(phi);

        // Set the quadratic coefficients.
        mag2 = omega*omega + sigma*sigma;
        aCoefs[a++] = 1.0;
        aCoefs[a++] = 0.0;
        aCoefs[a++] = zero * zero;
        bCoefs[b++] = 1.0;
        bCoefs[b++] = -2 * sigma;
        bCoefs[b++] = mag2;
        gain *= mag2 / (zero*zero);
    }
}


/**
 * Calculates normal elliptic coefficients.
 */
void IIR::elliptCoefs()
{
    int     m, a, b, odd;
    double  kernel, ratio, epsilon, kn, rt, sigma, omega, zero,
            CEIrt, CEIkn, vo, fm,
            SP, CP, DP, SN, CN, DN;

    // Make calculations of necessary constants.
    epsilon = sqrt( pow(10.0,-0.1*apass1) - 1.0 );

    // Calculate frequency ratio based on filter selectivity.
    if( filtType == "lowpass" )
        ratio = wstop1 / wpass1;
    else if( filtType == "highpass" )
        ratio = wpass1 / wstop1;
    else if( filtType == "bandpass" )
        ratio = (wstop2-wstop1) / (wpass2-wpass1);
    else
        ratio = (wpass2-wpass1) / (wstop2-wstop1);

    // Calculate values used in future calculations.
    kernel = ( pow(10,-0.1*astop1) - 1 ) / ( pow(10,-0.1*apass1) - 1 );
    rt = 1 / ratio;
    kn = 1 / sqrt(kernel);

    // Calculate ellip integrals, vo and ellipFun.
    CEIrt = ellipticIntegral(rt);
    CEIkn = ellipticIntegral(kn);
    vo = ( CEIrt / (CEIkn*order) ) * arcsc( 1/epsilon, kn );
    ellipticFun( vo, sqrt(1-rt*rt), &SP, &CP, &DP );

    // Handle odd order if necessary.
    a = 0;
    b = 0;
    odd = order % 2;
    if( odd )
    {
        gain = 1.0;
        aCoefs[a++] = 0.0;
        aCoefs[a++] = 0.0;
        aCoefs[a++] = SP * CP / (1-SP*SP);
        bCoefs[b++] = 0.0;
        bCoefs[b++] = 1.0;
        bCoefs[b++] = SP * CP / (1-SP*SP);
    }
    else
        gain = pow( 10.0, 0.05*apass1 );

    // Handle all quadratic terms.
    for( m=0; m<order/2; ++m )
    {
        // Make intermediate calculations.
        fm = CEIrt * ( 2*m+1+odd ) / order;
        ellipticFun( fm, rt, &SN, &CN, &DN );

        // Calculate real and imag coordinates of poles.
        sigma = -1*CN*DN*SP*CP / ( 1 - DN*DN*SP*SP );
        omega = SN*DP / ( 1 - DN*DN*SP*SP );
        zero = 1 / (rt*SN);

        // Set the quadratic coefficients.
        aCoefs[a++] = 1.0;
        aCoefs[a++] = 0.0;
        aCoefs[a++] = zero*zero;
        bCoefs[b++] = 1.0;
        bCoefs[b++] = -2*sigma;
        bCoefs[b++] = sigma*sigma + omega*omega;
        gain *= (sigma*sigma + omega*omega) / (zero*zero);
    }
}


/**
 * Converts normal lowpass coefs to unnormal lowpass
 * coefficients at a specific frequency.
 */
void IIR::unnormLP( double freq )
{
    int qd, cf, qdStart;

    // Handle first order, if odd; set qdStart.
    if( order % 2 )
    {
        aCoefs[2] *= freq;
        bCoefs[2] *= freq;
        qdStart = 1;
    }
    else
        qdStart = 0;

    // Handle quadratic factors, qd indexes through
    // quadratic factors, cf converts to coef number.
    for( qd=qdStart; qd<(order+1)/2; ++qd )
    {
        cf = qd * 3;
        aCoefs[cf+1] *= freq;
        aCoefs[cf+2] *= (freq*freq);
        bCoefs[cf+1] *= freq;
        bCoefs[cf+2] *= (freq*freq);
    }
}


/**
 * Converts normal lowpass coefs to unnormal highpass
 * coefficients at a specific frequency.
 */
void IIR::unnormHP( double freq )
{
    int qd, cf, qdStart;

    // Handle first order, if odd; set qdStart
    if( order % 2 )
    {
        gain *= ( aCoefs[2] / bCoefs[2] );
        aCoefs[2] = freq*aCoefs[1] / aCoefs[2];
        aCoefs[1] = 1.0;
        bCoefs[2] = freq / bCoefs[2];
        qdStart = 1;
    }
    else
        qdStart = 0;

    // Handle quadratic factors, qd indexes through
    // quadratic factors, cf converts to coef number.
    for( qd=qdStart; qd<(order+1)/2; ++qd )
    {
        cf = qd * 3;
        gain *= ( aCoefs[cf+2] / bCoefs[cf+2] );
        aCoefs[cf+1] *= ( freq / aCoefs[cf+2] );
        aCoefs[cf+2] = freq * freq *aCoefs[cf]
                            / aCoefs[cf+2];
        aCoefs[cf] = 1.0;
        bCoefs[cf+1] *= ( freq / bCoefs[cf+2] );
        bCoefs[cf+2] = freq * freq * bCoefs[cf]
                            / bCoefs[cf+2];
        bCoefs[cf] = 1.0;
    }
}


/**
 * Converts normal lowpass coefficients to unnormal
 * bandpass coefficients at a specific frequency.
 */
void IIR::unnormBP( double BW, double Wo )
{
    int     qd, ocf, ncf, qdStart,
            numCoefs, orgQuads, orgOrder;
    complex<double> A, B, C, D, E;
    Vector<double>  orgNum,
                    orgDen,
                    newNum,
                    newDen;

    // Store orig number of quadratics and order,
    // new order will be twice original.
    orgOrder = order;
    orgQuads = (orgOrder + 1)/2;
    order = orgOrder * 2;
    orgNum = aCoefs;
    orgDen = bCoefs;

    numCoefs = 3 * orgOrder;
    newNum.resize(numCoefs);
    newDen.resize(numCoefs);

    // If orgOrder odd, convert first order factor to
    // quadratic, qdStart indic start pt for loop.
    if( orgOrder % 2 )
    {
        newNum[0] = orgNum[1];
        newNum[1] = BW * orgNum[2];
        newNum[2] = orgNum[1] * Wo * Wo;
        newDen[0] = orgDen[1];
        newDen[1] = BW * orgDen[2];
        newDen[2] = orgDen[1] * Wo * Wo;
        qdStart = 1;
    }
    else
        qdStart = 0;

    // Each orig quad term will be converted to two new
    // quadratics via complex quadratic factoring.
    for( qd=qdStart; qd<orgQuads; ++qd )
    {
        ocf = qd * 3;
        ncf = qd * 6 - qdStart * 3;

        // For numers which DON'T have s^2 or s terms.
        if( orgNum[ocf] == 0.0 )
        {
            newNum[ncf] = 0.0;
            newNum[ncf+1] = sqrt(orgNum[ocf+2]) * BW;
            newNum[ncf+2] = 0.0;
            newNum[ncf+3] = 0.0;
            newNum[ncf+4] = sqrt(orgNum[ocf+2]) * BW;
            newNum[ncf+5] = 0.0;
        }

        // For numers which DO have s^2 and s terms.
        else
        {
            A = complex<double>( orgNum[ocf], 0.0 );
            B = complex<double>( orgNum[ocf+1], 0.0 );
            C = complex<double>( orgNum[ocf+2], 0.0 );
            quadradicRoot( A, B, C, D, E );

            A = complex<double>( 1.0, 0.0 );
            B = - D * complex<double>( BW, 0.0 );
            C = complex<double>( Wo*Wo, 0.0 );
            quadradicRoot( A, B, C, D, E );

            newNum[ncf] = 1.0;
            newNum[ncf+1] = -2.0 * real(D);
            newNum[ncf+2] = norm(D);
            newNum[ncf+3] = 1.0;
            newNum[ncf+4] = -2.0 * real(E);
            newNum[ncf+5] = norm(E);
        }

        A = complex<double>( orgDen[ocf], 0.0 );
        B = complex<double>( orgDen[ocf+1], 0.0 );
        C = complex<double>( orgDen[ocf+2], 0.0 );
        quadradicRoot( A, B, C, D, E );

        // Make required substitutions, factor again.
        A = complex<double>( 1.0, 0.0 );
        B = - D * complex<double>( BW, 0.0 );
        C = complex<double>( Wo*Wo, 0.0 );
        quadradicRoot( A, B, C, D, E );

        // Make required substitutions, factor again.
        newDen[ncf] = 1.0;
        newDen[ncf+1] = -2.0 * real(D);
        newDen[ncf+2] = norm(D);
        newDen[ncf+3] = 1.0;
        newDen[ncf+4] = -2.0 * real(E);
        newDen[ncf+5] = norm(E);
    }

    aCoefs = newNum;
    bCoefs = newDen;
}


/**
 * Converts normal lowpass coefficients to unnormal
 * bandstop coefficients at a specific frequency.
 */
void IIR::unnormBS( double BW, double Wo )
{
    int     qd, ocf, ncf, qdStart,
            numCoefs, orgQuads, orgOrder;
    complex<double> A, B, C, D, E;
    Vector<double>  orgNum,
                    orgDen,
                    newNum,
                    newDen;

    // Store orig number of quadratics and order,
    // new order will be twice original.
    orgOrder = order;
    orgQuads = (orgOrder + 1)/2;
    order = orgOrder * 2;
    orgNum = aCoefs;
    orgDen = bCoefs;

    numCoefs = 3 * orgOrder;
    newNum.resize(numCoefs);
    newDen.resize(numCoefs);

    // If orgOrder odd, convert first order factor to
    // quadratic, qdStart indic start pt for loop.
    if( orgOrder % 2 )
    {
        gain *= (orgNum[2] / orgDen[2]);
        newNum[0] = 1.0;
        newNum[1] = BW * orgNum[1] / orgNum[2];
        newNum[2] = Wo * Wo;
        newDen[0] = 1.0;
        newDen[1] = BW * orgDen[1] / orgDen[2];
        newDen[2] = Wo * Wo;
        qdStart = 1;
    }
    else
        qdStart = 0;

    // Each orig quad term will be converted to two new
    // quadratics via complex quadratic factoring.
    for( qd=qdStart; qd<orgQuads; qd++ )
    {
        ocf = qd * 3;
        ncf = qd * 6 - qdStart * 3;
        gain *= (orgNum[ocf+2] / orgDen[ocf+2]);

        // For numers which DON'T have s^2 or s terms.
        if(orgNum[ocf] == 0)
        {
            newNum[ncf] = 1.0;
            newNum[ncf+1] = 0.0;
            newNum[ncf+2] = Wo * Wo;
            newNum[ncf+3] = 1.0;
            newNum[ncf+4] = 0.0;
            newNum[ncf+5] = Wo * Wo;
        }

        // For numers which DO have s^2 and s terms.
        else
        {
            A = complex<double>( orgNum[ocf], 0.0 );
            B = complex<double>( orgNum[ocf+1], 0.0 );
            C = complex<double>( orgNum[ocf+2], 0.0 );
            quadradicRoot( A, B, C, D, E );

            A = complex<double>( 1.0, 0.0 );
            B = - 1.0/D * complex<double>( BW, 0.0 );
            C = complex<double>( Wo*Wo, 0.0 );
            quadradicRoot( A, B, C, D, E );

            newNum[ncf] = 1.0;
            newNum[ncf+1] = -2.0 * real(D);
            newNum[ncf+2] = norm(D);
            newNum[ncf+3] = 1.0;
            newNum[ncf+4] = -2.0 * real(E);
            newNum[ncf+5] = norm(E);
        }

        A = complex<double>( orgDen[ocf], 0.0 );
        B = complex<double>( orgDen[ocf+1], 0.0 );
        C = complex<double>( orgDen[ocf+2], 0.0 );
        quadradicRoot( A, B, C, D, E );

        // Make required substitutions, factor again.
        A = complex<double>( 1.0, 0.0 );
        B = - 1.0/D * complex<double>( BW, 0.0 );
        C = complex<double>( Wo*Wo, 0.0 );
        quadradicRoot( A, B, C, D, E );

        // Make required substitutions, factor again.
        newDen[ncf] = 1.0;
        newDen[ncf+1] = -2.0 * real(D);
        newDen[ncf+2] = norm(D);
        newDen[ncf+3] = 1.0;
        newDen[ncf+4] = -2.0 * real(E);
        newDen[ncf+5] = norm(E);
    }

    aCoefs = newNum;
    bCoefs = newDen;
}


/**
 * Calculates response for IIR filters.
 */
double IIR::frqeResp( double freq )
{
    int     c, q;
    double  omega, omega2,
            rea, img, tmp,
            mag;

    // Initialize magnatude.
    mag = gain;

    // Pre-calculate omega and omega squared.
    omega = TWOPI * freq / fsamp;
    omega2 = 2 * omega;

    // Loop through coefficients for each quadratic.
    for( q=0; q<(order+1)/2; q++ )
    {
        c = 3*q;

        rea = aCoefs[c] + aCoefs[c+1]*cos(omega) + aCoefs[c+2]*cos(omega2);
        img = -aCoefs[c+1]*sin(omega) - aCoefs[c+2]*sin(omega2);
        tmp = sqrt( rea*rea + img*img );
        mag *= tmp;

        rea = bCoefs[c] + bCoefs[c+1]*cos(omega) + bCoefs[c+2]*cos(omega2);
        img = -bCoefs[c+1]*sin(omega) - bCoefs[c+2]*sin(omega2);
        tmp = sqrt(rea*rea + img*img);
        mag /= tmp;
    }

    // Convert magnitude response to dB.
    if( mag < EPS )
        mag = EPS;
    mag = 20 * log10(mag);

    return mag;
}


/**
 * Calculate the response at edge freqency.
 */
void IIR::calcGain()
{

    // Determine the edge frequency.
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
