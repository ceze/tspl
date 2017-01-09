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
 *                                    stack.h
 *
 * A C++ template class implemented stack.
 *
 * The defualt initial size of the stack is set to 20. If the elements number
 * exceed initial size, then it will be extended by a factor of 2.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef STACK_H
#define STACK_H


#include <iostream>
#include <cstdlib>
#include <constants.h>


using namespace std;


namespace splab
{

    template <typename Type>
    class Stack
    {

    public:

        explicit Stack( int maxSize = INITSIZE );
        Stack( const Stack<Type> &rhs );
        ~Stack();
        Stack<Type>& operator=( const Stack<Type> &rhs );

        bool isEmpty() const;
        void makeEmpty();

        void push( const Type &x );
        void pop();
        void pop( Type &x );
        void getTop( Type &x );

    private:

        int top;
        int capacity;

        Type *elements;

        void handleOverflow();
        void handleUnderflow();

    };
    // class Stack


    #include <stack-impl.h>

}
// namespace splab


#endif
// STACK_H
