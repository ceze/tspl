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
 *                                    hashtable.h
 *
 * Hash table implemented by C++ template class.
 *
 * This class provides "search", "insert" and "remove" operations by using a
 * Hash Table. We use the quadratic probing method to prevent collision.
 * To ensure a new element can always be inserted, the table at least half
 * empty.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef HASHTABLE_H
#define HASHTABLE_H


#include <iostream>
#include <cstdlib>
#include <string>


using namespace std;


namespace splab
{

    template <typename Object, typename Key>
    class HashTable
    {

    public:

        explicit HashTable( int size=7 );
        HashTable( const HashTable<Object, Key> &rhs );
        ~HashTable();
        HashTable<Object, Key>& operator=( const HashTable<Object, Key> &rhs );

        void makeEmpty();

        bool search( const Key k, Object &x ) const;
        bool insert( const Object &x );
        bool remove( const Key k, Object &x );

        enum EntryFlag { ACTIVE, EMPTY, DELETED };

    private:

        struct HashEntry
        {
            Object      element;
            EntryFlag   info;

            HashEntry( const Object &x=Object(), EntryFlag i=EMPTY )
              : element(x), info(i)
            { }
        };

        HashEntry *array;

        int currentSize;
        int tableSize;

        bool isActive( int currentPos ) const;
        int findPos( const Key k ) const;
        void rehash( );

        int myhash( const Key k ) const;

    };


    static bool isPrime( int n );
    static int  nextPrime( int n );


    #include <hashtable-impl.h>

}
// namespace splab


#endif
// HASHTABLE_H
