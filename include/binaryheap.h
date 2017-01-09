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
 *                                    binaryheap.h
 *
 * Minimum binary heap implemented by C++ template class.
 *
 * This class provides the following operatins of a minimum binary heap:
 *      build a heap
 *      insert an element into the heap
 *      find the minimum element in the heap
 *      delete the minimum element in the heap
 *
 * The defualt initial size of the heap is set to 20. If the elements number
 * exceed initial size, then it will be extended by a factor of 2.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef BINARYHEAP_H
#define BINARYHEAP_H


#include <iostream>
#include <cstdlib>
#include <constants.h>


using namespace std;


namespace splab
{

    template <typename Type>
    class BinaryHeap
    {

    public:

        explicit BinaryHeap( int maxSize = INITSIZE );
        BinaryHeap( Type *array, int length );
        BinaryHeap( const BinaryHeap<Type> &rhs );
        ~BinaryHeap();
        BinaryHeap<Type>& operator=( const BinaryHeap<Type> &rhs );

        bool isEmpty() const;
        void makeEmpty();
        int size() const;

        void insert( const Type &x );
        void findMin( Type &x );
        void deleteMin();
        void deleteMin( Type &minItem );

    private:

        Type *elements;

        int currentSize;
        int capacity;

        void filterDown( int hole );
        void filterUp( int hole );

        void handleOverflow();
        void handleUnderflow();

    };
    // class BinaryHeap


    #include <binaryheap-impl.h>

}
// namespace splab


#endif
// BINARYHEAP_H
