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
 *                               stack-impl.h
 *
 * Implementation for Stack class.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template <typename Type>
Stack<Type>::Stack( int maxSize ) : top (-1), capacity(maxSize)
{
    elements = new Type[maxSize];
    if( !elements )
    {
        cerr << "Out of memory!" << endl;
        exit(1);
    }
};

template <typename Type>
Stack<Type>::Stack( const Stack<Type> &rhs )
{
    top = rhs.top;
    capacity = rhs.capacity;

    elements = new Type[capacity];
    for( int i=0; i<=top; ++i )
        elements[i] = rhs.elements[i];
}

template <typename Type>
Stack<Type>::~Stack()
{
    top = -1;
    capacity = INITSIZE;
    delete []elements;
}


/**
 * Overload copy assignment operation.
 */
template <typename Type>
Stack<Type>& Stack<Type>::operator=( const Stack<Type> &rhs )
{
    top = rhs.top;
    capacity = rhs.capacity;

    delete []elements;
    elements = new Type[capacity];
    for( int i=0; i<=top; ++i )
        elements[i] = rhs.elements[i];

    return *this;
}


/**
 * If the stack is empty, return true.
 */
template <typename Type>
inline bool Stack<Type>::isEmpty() const
{
    return ( top == -1 );
}


/**
 * Make the stack empty.
 */
template <typename Type>
inline void Stack<Type>::makeEmpty()
{
    top = -1;
}


/**
 * Push an element into the stack.
 */
template <typename Type>
inline void Stack<Type>::push( const Type &x )
{
    elements[++top] = x;
    if( top == capacity-1 )
        handleOverflow();
};


/**
 * Pop an element out of stack.
 */
template <typename Type>
inline void Stack<Type>::pop()
{
    if( !isEmpty() )
        top--;
    else
        handleUnderflow();
};

template <typename Type>
inline void Stack<Type>::pop( Type &x )
{
    if( !isEmpty() )
        x = elements[top--];
    else
        handleUnderflow();
};


/**
 * Get the top element of the stack.
 */
template <typename Type>
inline void Stack<Type>::getTop( Type &x )
{
    if( !isEmpty() )
        x = elements[top];
    else
        handleUnderflow();
};


/**
 * If the capability of the stack exceeds the initial size, make it double.
 */
template <typename Type>
void Stack<Type>::handleOverflow()
{
    capacity = EXTFACTOR * capacity ;

    Type *newArray = new Type[capacity];
    if( newArray == NULL )
    {
        cerr << "Out of memory!" << endl;
        exit(1);
    }

    for( int i=0; i<=top; ++i )
        newArray[i] = elements[i];

    delete []elements;
    elements = newArray;
};


/**
 * Handle the error of get element from an empty stack.
 */
template <typename Type>
inline void Stack<Type>::handleUnderflow()
{
    cerr << "The stack is empty!" << endl << endl;
    exit( 1 );
};
