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
 *                               dfd-impl.h
 *
 * Implementation for DFD class.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
DFD::DFD( const string &select ) : filtType(select)
{
}

DFD::~DFD()
{
}


/**
 * Get a number which falls between min and max.
 */
double DFD::getValue( double x, double min, double max )
{

    if( (x > min) && (x < max) )
        return x;
    else
    {
        cerr << "The parameter is out of range!" << endl;
        return 0;
    }
}


/**
 * Set filter's parameters
 */
void DFD::setParams( double fs, double f1, double a1, double f2, double a2 )
{
    fsamp = getValue( fs, FREQ_MIN, FREQ_MAX );
    double maxFreq = fsamp / 2;

    if( filtType == "lowpass" )
    {
        wpass1 = getValue( f1, FREQ_MIN, maxFreq );
        apass1 = getValue( a1, GAIN_TRAN, GAIN_PASS );
        wstop1 = getValue( f2, wpass1, maxFreq );
        astop1 = getValue( a2, GAIN_STOP, GAIN_TRAN );
    }
    else if( filtType == "highpass" )
    {
        wstop1 = getValue( f1, FREQ_MIN, maxFreq );
        astop1 = getValue( a1, GAIN_STOP, GAIN_TRAN );
        wpass1 = getValue( f2, wstop1, maxFreq );
        apass1 = getValue( a2, GAIN_TRAN, GAIN_PASS );
    }
    else
        cerr << "Parameters setting has failed!" << endl;

    // Convert all edge frequencies to rad/sec
    wpass1 *= TWOPI;
    wstop1 *= TWOPI;

}

void DFD::setParams( double fs, double f1, double a1, double f2,
                     double f3, double a2, double f4, double a3 )
{
    fsamp = getValue( fs, FREQ_MIN, FREQ_MAX );
    double maxFreq = fsamp / 2;

    if( filtType == "bandpass" )
    {
        wstop1 = getValue( f1, FREQ_MIN, maxFreq );
        astop1 = getValue( a1, GAIN_STOP, GAIN_TRAN );
        wpass1 = getValue( f2, wstop1, maxFreq );
        wpass2 = getValue( f3, wpass1, maxFreq );
        apass1 = getValue( a2, GAIN_TRAN, GAIN_PASS );
        wstop2 = getValue( f4, wpass2, maxFreq );
        astop2 = getValue( a3, GAIN_STOP, GAIN_TRAN );
        apass2 = apass1;
        if( astop1 < astop2 )
            astop2 = astop1;
        else
            astop1 = astop2;
    }
    else if( filtType == "bandstop" )
    {
        wpass1 = getValue( f1, FREQ_MIN, maxFreq );
        apass1 = getValue( a1, GAIN_TRAN, GAIN_PASS );
        wstop1 = getValue( f2, wpass1, maxFreq );
        wstop2 = getValue( f3, wstop1, maxFreq );
        astop1 = getValue( a2, GAIN_STOP, GAIN_TRAN );
        wpass2 = getValue( f4, wstop2, maxFreq );
        apass2 = getValue( a3, GAIN_TRAN, GAIN_PASS );
        astop2 = astop1;
        if( apass1 < apass2 )
            apass1 = apass2;
        else
            apass2 = apass1;
    }
    else
        cerr << "Parameters setting has failed!" << endl;

    // Convert all edge frequencies to rad/sec
    wpass1 *= TWOPI;
    wpass2 *= TWOPI;
    wstop1 *= TWOPI;
    wstop2 *= TWOPI;
}
