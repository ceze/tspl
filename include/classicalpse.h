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
 *                               classicalpse.h
 *
 * Classical Power Spectrum Estimation Mothods.
 *
 * The goal of spectral density estimation is to estimate the spectral density
 * of a random signal from a sequence of time samples of the signal. The
 * purpose of estimating the spectral density is to detect any periodicities
 * in the data, by observing peaks at the frequencies corresponding to these
 * periodicities.
 *
 * This file provides 5 usually used Classical-Specturm-Estimation methods:
 *    correlogram method,    periodogram method,
 *    smoothed periodogram method (Barteltt method and Welch method),
 *    and Blackman-Tukey method
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef CLASSICALPSE_H
#define CLASSICALPSE_H


#include <window.h>
#include <utilities.h>
#include <fft.h>
#include <correlation.h>


namespace splab
{

    template<typename Type> Vector<Type> correlogramPSE( const Vector<Type>&,
                                                         int );

    template<typename Type> Vector<Type> periodogramPSE( const Vector<Type>&,
                                                         const Vector<Type>&,
                                                         int );
    template<typename Type> Vector<Type> bartlettPSE( const Vector<Type>&,
                                                      int, int );
    template<typename Type> Vector<Type> welchPSE( const Vector<Type>&,
                                                  const Vector<Type>&,
                                                  int, int );

    template<typename Type> Vector<Type> btPSE( const Vector<Type>&,
                                               const Vector<Type>&, int );

    #include <classicalpse-impl.h>

}
// namespace splab


#endif
// CLASSICALSE_H
