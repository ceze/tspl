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
 *                             huffmancode-impl.h
 *
 * Implementation for Huffman code.
 *
 * Zhang Ming, 2010-07, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template <typename Object, typename Weight>
HuffmanTree<Object,Weight>::HuffmanTree() : root(NULL)
{
    arraySize = 0;
    codeTable = NULL;
}

template <typename Object, typename Weight>
HuffmanTree<Object,Weight>::~HuffmanTree()
{
    HTNode<Object,Weight> *r = root;
    destroy(r);
    root = NULL;
    delete []codeTable;
}


/**
 * Huffman codeing.
 */
template <typename Object, typename Weight>
void HuffmanTree<Object,Weight>::code( CodeObject<Object,Weight> *codeArray,
                                       int length )
{
    arraySize = length;

    if( codeTable != NULL )
        delete []codeTable;
    codeTable = new HCNode<Object>[arraySize];

    createHuffmanTree( codeArray );
    createCodeTable();
}


/**
 * If decodeing successful return ture, else return false;
 */
template <typename Object, typename Weight>
bool HuffmanTree<Object,Weight>::decode( unsigned char codeword[CODESIZE],
                                         unsigned int length,
                                         Object &decodeword )
{
    unsigned int pos = 0;
    HTNSmartPtr<Object,Weight> p = root;
    for( pos=0; pos<length; ++pos )
        if( p != NULL )
        {
            if( getBit( codeword, pos ) == 0 )
                p = p->left;
            else
                p = p->right;
        }
        else
            break;

    if( pos == length && p != NULL && p->left == NULL && p->right == NULL  )
    {
        decodeword = p->data;
        return true;
    }

    else
        return false;
}


/**
 * Print the codedword.
 */
template <typename Object, typename Weight>
inline void
HuffmanTree<Object,Weight>::printCode( unsigned char codewrod[CODESIZE],
                                       unsigned int length )
{
    for( unsigned int j=0; j<length; ++j )
        cout << getBit( codewrod, j );
}


/**
 * Print the coded table.
 */
template <typename Object, typename Weight>
void HuffmanTree<Object,Weight>::printCodeTable()
{
	cout << "Object\tCode\tSize" << endl;
	for( int i=0; i<arraySize; ++i )
    {
        cout << codeTable[i].data  << "\t";
        printCode( codeTable[i].bits, codeTable[i].length );
        cout << "\t" << codeTable[i].length << endl;
    }
	cout << endl;
}


/**
 * Create the Huffman tree.
 */
template <typename Object, typename Weight>
void HuffmanTree<Object,Weight>::createHuffmanTree(
                                 CodeObject<Object,Weight> *codeArray )
{
    BinaryHeap<HTNSmartPtr<Object,Weight> > heap( arraySize );
    HTNSmartPtr<Object,Weight>  tree = NULL,
                                subtreeL = NULL,
                                subtreeR = NULL;

    for( int i=0; i<arraySize; ++i )
    {
        tree = new HTNode<Object,Weight>( codeArray[i].data,
                                          codeArray[i].cost,
                                          NULL, NULL );
        heap.insert( tree );
    }

    while( heap.size() > 1 )
    {
        heap.deleteMin( subtreeL );
        heap.deleteMin( subtreeR );
        tree = new HTNode<Object,Weight>( Object(),
                                          subtreeL->cost+subtreeR->cost,
                                          subtreeL, subtreeR );
        heap.insert( tree );
    }
    heap.deleteMin( root );
}


/**
 * Create the coded character table.
 */
template <typename Object, typename Weight>
void HuffmanTree<Object,Weight>::createCodeTable()
{
    for( int i=0; i<arraySize; ++i )
    {
        codeTable[i].data = Object();
        codeTable[i].length = 0;
    }

    unsigned char code[CODESIZE];
    int index = 0;
    createCodeTableRecursive( root, code, 0, index );
}


/**
 * Create the coded character table.
 */
template <typename Object, typename Weight>
void HuffmanTree<Object,Weight>::createCodeTableRecursive(
                                 HTNSmartPtr<Object,Weight> ht,
                                 unsigned char *code,
                                 int pos, int &index )
{
    if( ht->left )
    {
        setBit( code, pos, 0 );
        createCodeTableRecursive( ht->left, code, pos+1, index );
    }

    if( ht->right )
    {
        setBit( code, pos, 1 );
        createCodeTableRecursive( ht->right, code, pos+1, index );
    }

    if( ht->left==NULL && ht->right==NULL )
    {
        codeTable[index].data = ht->data;
        for( int i=0; i<CODESIZE; ++i )
            codeTable[index].bits[i] = code[i];
        codeTable[index].length = pos;
        index++;
    }
}


/**
 * Set the bit at position pos in the array bits to the value state.
 */
template <typename Object, typename Weight>
void HuffmanTree<Object,Weight>::setBit( unsigned char *bits,
                                         unsigned int pos,
                                         unsigned int state )
{
	unsigned char mask = 0x80;
	for( unsigned int i=0; i<(pos % 8); ++i )
		mask = mask >> 1;

	if( state == 1 )
		bits[pos/8] = bits[pos/8] | mask;
	else if( state == 0 )
		bits[pos/8] = bits[pos/8] & (~mask);
    else
        cerr << endl << "The bit to be set should be '1' or '0'!" << endl;

	return;
}


/**
 * Get the state of the bit at position pos in the array bits.
 */
template <typename Object, typename Weight>
unsigned int HuffmanTree<Object,Weight>::getBit( unsigned char *bits,
                                                 unsigned int pos )
{
	unsigned char mask = 0x80;
	for( unsigned int i=0; i<(pos%8); ++i )
		mask = mask >> 1;

	return ( ((mask & bits[(int)(pos/8)]) == mask) ? 1 : 0 );
}


/**
 * Destroy the tree.
 */
template <typename Object, typename Weight>
void HuffmanTree<Object,Weight>::destroy( HTNode<Object,Weight> *&r )
{
    if( r != NULL )
    {
        destroy( r->left );
        destroy( r->right );
        delete r;
    }
    r = NULL;
}
