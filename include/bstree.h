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
 *                                 bstree.h
 *
 * Binary Search Tree
 *
 * This class provides "traversal(preorder,inorder or postorder), "search",
 * "insert" and "remove" operations for Binary Search Tree.
 *
 * All of the algorithms are implemented by the nonrecursion method.
 *
 * Zhang Ming, 20010-07, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef BSTREE_H
#define BSTREE_H


#include <iostream>
#include <stack.h>


using namespace std;


namespace splab
{

    /**
     * Node in Binary Search Tree
     */
    template <typename Object>
    struct BSNode
    {
        Object data;

        BSNode<Object>  *left,
                        *right;

        BSNode() : left(NULL), right(NULL)
        { }

        BSNode( const Object &x, BSNode<Object> *lp=NULL,
                BSNode<Object> *rp=NULL )
        {   data = x; left = lp; right = rp;    }

        bool operator<( BSNode<Object> &x )
        {   return data < x.data;    }

        bool operator>( BSNode<Object> &x )
        {   return data > x.data;    }

        bool operator==( BSNode<Object> &x )
        {	return data == x.data;    }

    };
    // struct BSNode


    template <typename Object, typename Key>
    class BSTree
    {

    public:

        BSTree();
        ~BSTree();

//        BSTree( const BSTree &rhs );
//        BSTree & operator=( const BSTree &rhs );

        bool isEmpty();
        void makeEmpty();

        void preTraverse();
        void inTraverse();
        void postTraverse();

        bool insert( const Object &x );
        bool remove( Key k );

        BSNode<Object>* search( Key k );
        Object minItem();
        Object maxItem();

    private:

        BSNode<Object> *root;

        void preTravNonrecursion( BSNode<Object> *t );
        void inTravNonrecursion( BSNode<Object> *t );
        void postTravNonrecursion( BSNode<Object> *t );

        void destroy( BSNode<Object> *&r );
        bool insert( BSNode<Object> *&r, const Object &x );
        bool remove( BSNode<Object> *&r, Key k );

        BSNode<Object>* search( BSNode<Object> *r, Key k );
        BSNode<Object>* minItem( BSNode<Object> *r );
        BSNode<Object>* maxItem( BSNode<Object> *r );
        void visit( BSNode<Object> *p );

    };
    // class BSTree


    #include <bstree-impl.h>

}
// namespace splab


#endif
// BSTREE_H
