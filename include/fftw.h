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
 *                                    fftw.h
 *
 * A simple C++ interface for FFTW. This file only provides the one dimension
 * DFT and IDFT for data type as "Vector<Type>" in SPLab, where "Type" can
 * be "double", "complex<double>, "float", "complex<float>, "long double",
 * "complex<long double>.
 *
 * If you need execute FFT many times in your program, it's not a good idear
 * for using this interface because of losing efficiency.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef FFTW_H
#define FFTW_H


#include <complex>
#include <fftw3.h>
#include <vector.h>


namespace splab
{

    void fftw( Vector<double>&,      Vector< complex<double> >& );
    void fftw( Vector<float>&,       Vector< complex<float> >& );
    void fftw( Vector<long double>&, Vector< complex<long double> >& );

    void fftw( Vector< complex<double> >&, Vector< complex<double> >& );
    void fftw( Vector< complex<float> >&,  Vector< complex<float> >& );
    void fftw( Vector< complex<long double> >&,
               Vector< complex<long double> >& );

    void ifftw( Vector< complex<double> >&, Vector<double>& );
    void ifftw( Vector< complex<float> >&,  Vector<float>& );
    void ifftw( Vector< complex<long double> >&, Vector<long double>& );

    void ifftw( Vector< complex<double> >&, Vector< complex<double> >& );
    void ifftw( Vector< complex<float> >&,  Vector< complex<float> >& );
    void ifftw( Vector< complex<long double> >&,
                Vector< complex<long double> >& );


    #include <fftw-impl.h>

}
// namespace splab


#endif
// FFTW_H
