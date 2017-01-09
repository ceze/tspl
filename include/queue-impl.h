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
 *                               queue-impl.h
 *
 * Implementation for Queue class.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template <typename Type>
Queue<Type>::Queue() : front(NULL), rear(NULL)
{
}

template <typename Type>
Queue<Type>::~Queue()
{
    makeEmpty();
}


/**
 * If the queue is empty, return true.
 */
template <typename Type>
inline bool Queue<Type>::isEmpty() const
{
    return front == NULL;
}


/**
 * Make the queue empty.
 */
template <typename Type>
inline void Queue<Type>::makeEmpty()
{
    LinkNode<Type> *p;
	while( front != NULL )
	{
		p = front;
		front = front->next;
		delete p;
	}
}


/**
 * Enter an element into the queue.
 */
template <typename Type>
void Queue<Type>::enqueue( const Type &x )
{
    if( front == NULL )
	{
		front = rear = new LinkNode<Type>( x );
        if( !front )
        {
            cerr << "Out of memory!" << endl;
            exit( 1 );
        }
    }
	else
	{
		rear->next = new LinkNode<Type>( x );
		if( !rear )
		{
            cerr << "Out of memory!" << endl;
            exit( 1 );
        }
		rear = rear->next;
	}
}


/**
 * Pop an element from the queue.
 */
template <typename Type>
void Queue<Type>::dequeue()
{
    if( !isEmpty() )
    {
        LinkNode<Type> *p = front;
        front = front->next;
        delete p;
    }
    else
        handleUnderflow();
}

template <typename Type>
void Queue<Type>::dequeue( Type &x )
{
    if( !isEmpty() )
    {
        LinkNode<Type> *p = front;
        x = front->data;
        front = front->next;
        delete p;
    }
    else
        handleUnderflow();
}


/**
 * Get the front element of the queue.
 */
template <typename Type>
inline void Queue<Type>::getFront( Type &x )
{
    if( !isEmpty() )
        x = front->data;
    else
        handleUnderflow();
}


/**
 * Handle the error of get element from an empty queue.
 */
template <typename Type>
inline void Queue<Type>::handleUnderflow()
{
    cerr << "The queue is empty!" << endl << endl;
    exit( 1 );
}
