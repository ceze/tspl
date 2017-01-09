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
 *                                 cwt_usefftw.h
 *
 * Continuous Wavelet Transform by Using FFTW.
 *
 * Class for continuous wavelet transform, which is designed for computing
 * the continuous wavelet transform and it's inverse transform of 1D signals.
 * For now this class only supports "Mexican hat"(real) and "Morlet"(complex)
 * wavlet.
 *
 * The mother wavelet doubles and scale parameters are specified by users. The
 * inverse transform can not achive perfect reconstruction, but with a
 * sufficient precision in practice. Of course you can improve the accurate
 * by extend the range of scale parameter.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef CWT_USEFFTW_H
#define CWT_USEFFTW_H


#include <cstdlib>
#include <string>
#include <fftw.h>
#include <matrix.h>


namespace splab
{

    template <typename Type>
	class CWTFFTW
	{

	public:

		CWTFFTW( const string &name );
		~CWTFFTW();

        void setScales( Type fs, Type fmin, Type fmax, Type dj=Type(0.25) );

        Matrix<Type> cwtR( Vector<Type> &signal );
		Vector<Type> icwtR( const Matrix<Type> &coefs );
		Matrix< complex<Type> > cwtC( Vector<Type> &signal );
		Vector<Type> icwtC( const Matrix< complex<Type> > &coefs );

	private:

        string waveType;
        Type delta;
        Vector<Type> scales;
        Matrix<Type> table;

        void setTable( int N );
		Type constDelta();

	};
	// class CWTFFTW

    #include <cwt_usefftw-impl.h>

}
// namespace splab


#endif
// CWT_USEFFTW_H
