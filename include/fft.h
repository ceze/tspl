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
 *                                    fft.h
 *
 * Some convenient interface for FFT algorithm. If signal length "N" is power
 * of 2, then calls FFTR2 class to compute, else calls FFTPF class.
 * Forward:     Xk = fftr2c(xn);    Xk = fftc2c(xn);
 * Inverse:     xn = fftc2r(Xk);    xn = fftc2c(Xk);
 *
 * These routines don't need FFTW lib, but less efficiency than FFTW.
 *
 * Zhang Ming, 2010-09, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef FFT_H
#define FFT_H


#include <fftmr.h>
#include <fftpf.h>


namespace splab
{

    template<typename Type>
    Vector< complex<Type> > fft( const Vector<Type>& );
    template<typename Type>
    Vector< complex<Type> > fft( const Vector< complex<Type> >& );
    template<typename Type>
    Vector< complex<Type> > ifft( const Vector< complex<Type> >& );

    template<typename Type>
    Vector< complex<Type> > fftr2c( const Vector<Type>& );
    template<typename Type>
    Vector< complex<Type> > fftc2c( const Vector< complex<Type> >& );

    template<typename Type>
    Vector<Type> ifftc2r( const Vector< complex<Type> >& );
    template<typename Type>
    Vector< complex<Type> > ifftc2c( const Vector< complex<Type> >& );


    #include <fft-impl.h>

}
// namespace splab


#endif
// FFT_H
