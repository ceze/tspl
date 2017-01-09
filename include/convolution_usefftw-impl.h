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
 *                          convolution_usefftw-impl.h
 *
 * Implementation for linear convolution by using FFTW.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Fast convolution by FFT.
 */
template<typename Type>
Vector<Type> fastConvFFTW( const Vector<Type> &xn, const Vector<Type> &yn )
{
    int M = xn.dim(),
        N = yn.dim(),
        L = M + N - 1;

    Vector<Type> zn(L),
                 xnPadded = wextend( xn, N-1, "right", "zpd" ),
                 ynPadded = wextend( yn, M-1, "right", "zpd" );

    Vector< complex<Type> > Xk( L/2+1 ),
                            Yk( L/2+1 ),
                            Zk( L/2+1 );
    fftw( xnPadded, Xk );
    fftw( ynPadded, Yk );
    Zk = Xk * Yk;
    ifftw( Zk, zn );

    return zn;
}


//    void fastConv( const Vector<double> &xn, const Vector<double> &yn,
//                   Vector<double> &zn )
//    {
//        int M = xn.dim(),
//            N = yn.dim(),
//            L = M + N - 1;
//
//        zn.resize( L );
//        Vector<double>  xnPadded( L ),
//                        ynPadded( L );
//
//        for( int i=0; i<M; ++i )
//            xnPadded[i] = xn[i];
//        for( int i=0; i<N; ++i )
//            ynPadded[i] = yn[i];
//
//        Vector< complex<double> >   Xk( L/2+1 ),
//                                    Yk( L/2+1 ),
//                                    Zk( L/2+1 );
//        fftw_plan r2cP,
//                  c2rP;
//        r2cP = fftw_plan_dft_r2c_1d( L, xnPadded.begin(),
//               reinterpret_cast<fftw_complex*>(Xk.begin()),
//               FFTW_ESTIMATE );
//        fftw_execute( r2cP );
//        r2cP = fftw_plan_dft_r2c_1d( L, ynPadded.begin(),
//               reinterpret_cast<fftw_complex*>(Yk.begin()),
//               FFTW_ESTIMATE );
//        fftw_execute( r2cP );
//
//        Zk = Xk * Yk;
//
//        c2rP = fftw_plan_dft_c2r_1d( L,
//               reinterpret_cast<fftw_complex*>(Zk.begin()),
//               zn.begin(), FFTW_ESTIMATE );
//        fftw_execute( c2rP );
//        zn /= double(L);
//
//        fftw_destroy_plan( r2cP );
//        fftw_destroy_plan( c2rP );
//    }
