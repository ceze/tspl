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
 *                                 avltree.h
 *
 * AVL Tree implemented by C++ class template.
 *
 * This class provides "traversal(preorder,inorder or postorder), "search",
 * "insert" and "remove" operations by using an AVL Tree (Balanced Binary
 * Search Tree). The rotation technology is used to keep the tree balanced.
 * Of cause, you can get the height of the tree by the subroutine "height".
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef AVLTREE_H
#define AVLTREE_H


#include <cstdlib>
#include <stack.h>


namespace splab
{

    /**
     * Node in AVL Tree
     */
    template <typename Object, typename Key>
    struct AVLNode
    {
        Object data;
        int balance;

        AVLNode<Object, Key> *left,
                             *right;

        AVLNode() : left(NULL), right(NULL), balance(0)
        { }

        AVLNode( const Object &x, AVLNode<Object, Key> *lp=NULL,
                 AVLNode<Object, Key> *rp=NULL ) : balance(0)
        {   data = x; left = lp; right = rp;    }

        bool operator<( AVLNode<Object, Key> &t )
        {   return data < t.data;    }

        bool operator>( AVLNode<Object, Key> &t )
        {   return data > t.data;    }

        bool operator==( AVLNode<Object, Key> &t )
        {	return data == t.data;    }

    };
    // class AVLNode


    /**
     * AVL Tree
     */
    template <typename Object, typename Key>
    class AVLTree
    {

    public:

        AVLTree();
        ~AVLTree();

    //	AVLTree( const AVLTree<Object, Type> &rhs );
    //	AVLTree<Object, Type>& operator=( const AVLTree<Object, Type> &rhs );

        bool isEmpty() const;
        void makeEmpty();

        void print( const string &mode );

        int height() const;
        AVLNode<Object, Key>* search( Key k );
        bool insert( Object &x );
        bool remove( Key k, Object &x );

    private:

        AVLNode<Object, Key> *root;

        void preTraversal( AVLNode<Object, Key> *t );
        void inTraversal( AVLNode<Object, Key> *t );
        void postTraversal( AVLNode<Object, Key> *t );

        void makeEmpty( AVLNode<Object, Key> * &t );

        int height( AVLNode<Object, Key> *t ) const;
        AVLNode<Object, Key>* search( Key k, AVLNode<Object, Key> *t );
        bool insert( AVLNode<Object, Key> * &ptr, Object &x );
        bool remove( AVLNode<Object, Key> * &ptr, Key k, Object &x );

        void rotateL( AVLNode<Object, Key> * &ptr );
        void rotateR( AVLNode<Object, Key> * &ptr );
        void rotateLR( AVLNode<Object, Key> * &ptr );
        void rotateRL( AVLNode<Object, Key> * &ptr );

        void handleUnderflow();

    };
    // class AVLTree


    #include <avltree-impl.h>

}
// namespace splab


#endif
// AVLTREE_H
