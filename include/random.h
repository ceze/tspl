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
 *                                    random.h
 *
 * Random Number Generator.
 *
 * The class Random can generate an "long int" random number. The default
 * initial state is set to 1, of cause it can be set to any integer through
 * routine "seed(int)".
 *
 * Based on the "Random" class, we can generate some usually used distributed
 * pseudorandom number or pseudorandom sequence. These distributions include:
 *      Uniform     Normal      Exponential     Rayleigh        Poisson
 *      Bernoulli
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef RANDOM_H
#define RANDOM_H


#include <vector.h>


namespace splab
{

    class Random
    {

    public:

        explicit Random( long int initValue=1 );
        ~Random();

        void seed( long int value );
        long int random();
        long int getM() const;

    private:

        static const long int A = 48271;
        static const long int M = 2147483647;
        static const long int Q = 44488;
        static const long int R = 3399;

        long int state;

    };


    template<typename Type> Type randu( int, const Type&, const Type& );
    template<typename Type> Type randn( int, const Type&, const Type& );
    template<typename Type> Type rande( int, const Type& );
    template<typename Type> Type randr( int, const Type& );
    template<typename Type> int randp( int, const Type& );
    template<typename Type> int randb( int, const Type& );

    template<typename Type> Vector<Type> randu( int, const Type&,
                                                const Type&, int );
    template<typename Type> Vector<Type> randn( int, const Type&,
                                                const Type&, int );
    template<typename Type> Vector<Type> rande( int, const Type&, int );
    template<typename Type> Vector<Type> randr( int, const Type&, int );
    template<typename Type> Vector<int>  randp( int, const Type&, int );
    template<typename Type> Vector<int>  randb( int, const Type&, int );


    #include <random-impl.h>

}
// namespace splab


#endif
// RANDOM_H
