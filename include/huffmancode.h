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
 *                               huffmancode.h
 *
 * Huffman coding algorighm implemented by C++ template.
 *
 * This class is designed for Huffman codeing and decoding algorighm. The
 * priority queue is used in building Huffman tree. For comparing pointers
 * of Huffman Tree Node by their pointing contents, I use the smart pointer
 * replacing nomal pointer. And for saving memory, the coded word is stored
 * bit-by-bit in an unsigned char array.
 *
 * Zhang Ming, 20010-07, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef HUFFMANCODE_H
#define HUFFMANCODE_H


#include <iostream>
#include <binaryheap.h>


namespace splab
{

    /**
     * Object node to be coded
     */
    template <typename Object, typename Weight>
    struct CodeObject
    {
        Object data;
        Weight cost;
    };
    // struct CodeObject


    /**
     * Huffman codeword node
     */
    const int CODESIZE = 2;
    template <typename Object>
    struct HCNode
    {
        Object data;
        unsigned char bits[CODESIZE];
        unsigned int length;
    };
    // struct HCNode


    /**
     * Huffman tree node
     */
    template <typename Object, typename Weight>
    struct HTNode
    {
        Object data;
        Weight cost;
        HTNode<Object,Weight> *left,
                              *right;

        HTNode() : data(Object()), cost(Weight(0)), left(NULL), right(NULL)
        { }
        HTNode( const Object &x, const Weight w,
                HTNode<Object,Weight> *lp=NULL,
                HTNode<Object,Weight> *rp=NULL )
        {   data = x; cost = w; left = lp; right = rp;    }

        operator Weight() const
        {   return cost;    }

    };
    // struct HTNode


    /**
     * Smart pointer for Huffman tree bode
     */
    template <typename Object, typename Weight>
    class HTNSmartPtr
    {

    public:

        HTNSmartPtr( HTNode<Object,Weight> *p = NULL )
        {   ptr = p;    }

        HTNode<Object,Weight> operator*()
        {   return *ptr;    }
        HTNode<Object,Weight>* operator->()
        {   return ptr;    }

        operator HTNode<Object,Weight> *()
        {   return ptr;    }
        bool operator<( HTNSmartPtr<Object,Weight> p )
        {   return ptr->cost < p->cost;    }
        bool operator>( HTNSmartPtr<Object,Weight> p )
        {   return ptr->cost > p->cost;    }

    private:

        HTNode<Object,Weight> *ptr;

    };
    // class HTNSmartPtr


    /**
     * Huffman tree
     */
    template <typename Object, typename Weight>
    class HuffmanTree
    {

    public:

        HuffmanTree();
        ~HuffmanTree();

//        HuffmanTree( const HuffmanTree &rhs );
//        HuffmanTree & operator=( const HuffmanTree &rhs );

        void code( CodeObject<Object,Weight> *codeArray, int length );
        bool decode( unsigned char bits[CODESIZE], unsigned int length,
                     Object &decodeword );
        void printCode( unsigned char bits[CODESIZE], unsigned int length );
        void printCodeTable();

    private:

        HTNSmartPtr<Object,Weight> root;
        int arraySize;
        HCNode<Object> *codeTable;

        void createHuffmanTree( CodeObject<Object,Weight> *codeArray );
        void createCodeTable();
        void createCodeTableRecursive( HTNSmartPtr<Object,Weight> ht,
                                       unsigned char *code,
                                       int pos, int &index );

        void setBit( unsigned char *bits, unsigned int pos,
                     unsigned int state );
        unsigned int getBit( unsigned char *codeword, unsigned int pos );

        void destroy( HTNode<Object,Weight> *&r );

    };
    // class HuffmanTree


    #include <huffmancode-impl.h>

}
// namespace splab


#endif
// HUFFMANCODE_H
