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
 *                                 fftw-impl.h
 *
 * Implementation for C++ interface for FFTW.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Real to complex DFT of 1D signal. If "xn" has N points, "Xk"
 * should has N/2+1 points.
 */
inline void fftw( Vector<double> &xn, Vector< complex<double> > &Xk )
{
    if( Xk.size() < xn.size()/2+1 )
        Xk.resize( xn.size()/2+1 );

    fftw_plan r2cP;
    r2cP = fftw_plan_dft_r2c_1d( xn.dim(), xn.begin(),
           reinterpret_cast<fftw_complex*>(Xk.begin()),
           FFTW_ESTIMATE );
    fftw_execute( r2cP );
    fftw_destroy_plan( r2cP );
}

inline void fftw( Vector<float> &xn, Vector< complex<float> > &Xk )
{
    if( Xk.size() < xn.size()/2+1 )
        Xk.resize( xn.size()/2+1 );

    fftwf_plan r2cP;
    r2cP = fftwf_plan_dft_r2c_1d( xn.dim(), xn.begin(),
           reinterpret_cast<fftwf_complex*>(Xk.begin()),
           FFTW_ESTIMATE );
    fftwf_execute( r2cP );
    fftwf_destroy_plan( r2cP );
}

inline void fftw( Vector<long double> &xn,
                  Vector< complex<long double> > &Xk )
{
    if( Xk.size() < xn.size()/2+1 )
        Xk.resize( xn.size()/2+1 );

    fftwl_plan r2cP;
    r2cP = fftwl_plan_dft_r2c_1d( xn.dim(), xn.begin(),
           reinterpret_cast<fftwl_complex*>(Xk.begin()),
           FFTW_ESTIMATE );
    fftwl_execute( r2cP );
    fftwl_destroy_plan( r2cP );
}


/**
 * Complex to complex DFT of 1D signal.
 */
inline void fftw( Vector< complex<double> > &xn,
                  Vector< complex<double> > &Xk )
{
    if( Xk.size() < xn.size() )
        Xk.resize( xn.size() );

    fftw_plan c2cP;
    c2cP = fftw_plan_dft_1d( xn.dim(),
           reinterpret_cast<fftw_complex*>(xn.begin()),
           reinterpret_cast<fftw_complex*>(Xk.begin()),
           FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute( c2cP );
    fftw_destroy_plan( c2cP );
}

inline void fftw( Vector< complex<float> > &xn,
                  Vector< complex<float> > &Xk )
{
    if( Xk.size() < xn.size() )
        Xk.resize( xn.size() );

    fftwf_plan c2cP;
    c2cP = fftwf_plan_dft_1d( xn.dim(),
           reinterpret_cast<fftwf_complex*>(xn.begin()),
           reinterpret_cast<fftwf_complex*>(Xk.begin()),
           FFTW_FORWARD, FFTW_ESTIMATE );
    fftwf_execute( c2cP );
    fftwf_destroy_plan( c2cP );
}

inline void fftw( Vector< complex<long double> > &xn,
                  Vector< complex<long double> > &Xk )
{
    if( Xk.size() < xn.size() )
        Xk.resize( xn.size() );

    fftwl_plan c2cP;
    c2cP = fftwl_plan_dft_1d( xn.dim(),
           reinterpret_cast<fftwl_complex*>(xn.begin()),
           reinterpret_cast<fftwl_complex*>(Xk.begin()),
           FFTW_FORWARD, FFTW_ESTIMATE );
    fftwl_execute( c2cP );
    fftwl_destroy_plan( c2cP );
}


/**
 * Complex to real IDFT of 1D. If "xn" has N points, "Xk" should
 * has N/2+1 points.
 */
inline void ifftw( Vector< complex<double> > &Xk, Vector<double> &xn )
{
    assert( Xk.size() == xn.size()/2+1 );

    fftw_plan c2rP;
    c2rP = fftw_plan_dft_c2r_1d( xn.dim(),
           reinterpret_cast<fftw_complex*>(Xk.begin()),
           xn.begin(), FFTW_ESTIMATE );
    fftw_execute( c2rP );
    fftw_destroy_plan( c2rP );

    xn /= double( xn.dim() );
}

inline void ifftw( Vector< complex<float> > &Xk, Vector<float> &xn )
{
    assert( Xk.size() == xn.size()/2+1 );

    fftwf_plan c2rP;
    c2rP = fftwf_plan_dft_c2r_1d( xn.dim(),
           reinterpret_cast<fftwf_complex*>(Xk.begin()),
           xn.begin(), FFTW_ESTIMATE );
    fftwf_execute( c2rP );
    fftwf_destroy_plan( c2rP );

    xn /= float( xn.dim() );
}

inline void ifftw( Vector< complex<long double> > &Xk,
                   Vector<long double> &xn )
{
    assert( Xk.size() == xn.size()/2+1 );

    fftwl_plan c2rP;
    c2rP = fftwl_plan_dft_c2r_1d( xn.dim(),
           reinterpret_cast<fftwl_complex*>(Xk.begin()),
           xn.begin(), FFTW_ESTIMATE );
    fftwl_execute( c2rP );
    fftwl_destroy_plan( c2rP );

    xn /= (long double)( xn.dim() );
}


/**
 * Complex to complex IDFT of 1D.
 */
inline void ifftw( Vector< complex<double> > &Xk,
                   Vector< complex<double> > &xn )
{
    if( xn.size() < Xk.size() )
        xn.resize( Xk.size() );

    fftw_plan c2cP;
    c2cP = fftw_plan_dft_1d( xn.dim(),
           reinterpret_cast<fftw_complex*>(Xk.begin()),
           reinterpret_cast<fftw_complex*>(xn.begin()),
           FFTW_BACKWARD, FFTW_ESTIMATE );
    fftw_execute( c2cP );
    fftw_destroy_plan( c2cP );

    xn /= complex<double>( xn.dim(), 0.0 );
}

inline void ifftw( Vector< complex<float> > &Xk,
                   Vector< complex<float> > &xn )
{
    if( xn.size() < Xk.size() )
        xn.resize( Xk.size() );

    fftwf_plan c2cP;
    c2cP = fftwf_plan_dft_1d( xn.dim(),
           reinterpret_cast<fftwf_complex*>(Xk.begin()),
           reinterpret_cast<fftwf_complex*>(xn.begin()),
           FFTW_BACKWARD, FFTW_ESTIMATE );
    fftwf_execute( c2cP );
    fftwf_destroy_plan( c2cP );

    xn /= complex<float>( xn.dim(), 0.0 );
}

inline void ifftw( Vector< complex<long double> > &Xk,
                   Vector< complex<long double> > &xn )
{
    if( xn.size() < Xk.size() )
        xn.resize( Xk.size() );

    fftwl_plan c2cP;
    c2cP = fftwl_plan_dft_1d( xn.dim(),
           reinterpret_cast<fftwl_complex*>(Xk.begin()),
           reinterpret_cast<fftwl_complex*>(xn.begin()),
           FFTW_BACKWARD, FFTW_ESTIMATE );
    fftwl_execute( c2cP );
    fftwl_destroy_plan( c2cP );

    xn /= complex<long double>( xn.dim(), 0.0 );
}
