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
 *                               queue.h
 *
 * A simple queue implemented by C++ template class .
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef QUEUE_H
#define QUEUE_H


#include <iostream>
#include <cstdlib>


using namespace std;


namespace splab
{

    /**
     * Queue Node
     */
    template <typename Type>
    struct LinkNode
    {
        Type data;
        LinkNode<Type> *next;

        LinkNode() : next(NULL)
        { }
        LinkNode( const Type &x, LinkNode<Type> *p=NULL ) : data(x), next(p)
        { }
    };


    /**
     * Queue
     */
    template <typename Type>
    class Queue
    {

    public:

        Queue();
        ~Queue();

//        Queue( const Queue<Type> &rhs);
//        Queue<Type>& operator=( const Queue<Type> &rhs);

        bool isEmpty() const;
        void makeEmpty();

        void enqueue( const Type &x );
        void dequeue();
        void dequeue( Type &x );
        void getFront( Type &x );

    private:

        LinkNode<Type> *front,
                       *rear;

        void handleUnderflow();

    };


    #include <queue-impl.h>

}
// namespace splab


#endif
// QUEUE_H
