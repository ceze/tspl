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
 *                                    sort.h
 *
 * Some sort algorithms.
 *
 * This file includes several usually used sorting algorithm, such as: bubble
 * sorting, selection sorting, insertion sorting, quick sorting, merging
 * sorting, and heap sorting.
 *
 * Zhang Ming, 2010-07, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef SORT_H
#define SORT_H


#include <vector.h>


using namespace std;


namespace splab
{

    template<typename Type> void bubbleSort( Vector<Type>&, int, int );
    template<typename Type> void selectSort( Vector<Type>&, int, int );
    template<typename Type> void insertSort( Vector<Type>&, int, int );
    template<typename Type> void quickSort( Vector<Type>&, int, int );
    template<typename Type> void mergSort( Vector<Type>&, int, int );
    template<typename Type> void heapSort( Vector<Type>&, int, int );

    template<typename Type> const Type& median3( Vector<Type>&, int, int );
    template<typename Type> void merg( Vector<Type>&, int, int, int, int );
    template<typename Type> void filterDown( Vector<Type>&, int, int );


    #include <sort-impl.h>

}
// namespace splab


#endif
// SORT_H
