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
 *                                 bstree-impl.h
 *
 * Implementation for binary search rree.
 *
 * Zhang Ming, 2010-07, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template <typename Object, typename Key>
BSTree<Object, Key>::BSTree() : root(NULL)
{ }

template <typename Object, typename Key>
BSTree<Object, Key>::~BSTree()
{
    destroy(root);
}


/**
 * If the tree is empty, return true.
 */
template <typename Object, typename Key>
inline bool BSTree<Object, Key>::isEmpty()
{
    return ( root == NULL );
}


/**
 * Make the tree empty.
 */
template <typename Object, typename Key>
inline void BSTree<Object, Key>::makeEmpty()
{
    destroy( root );
}


/**
 * preorder traversal
 */
template <typename Object, typename Key>
inline void BSTree<Object, Key>::preTraverse()
{
    preTravNonrecursion( root );
}


/**
 * inorder traversal
 */
template <typename Object, typename Key>
inline void BSTree<Object, Key>::inTraverse()
{
    inTravNonrecursion( root );
}


/**
 * postorder traversal
 */
template <typename Object, typename Key>
inline void BSTree<Object, Key>::postTraverse()
{
    postTravNonrecursion( root );
}


/**
 * Insert a new element x in the tree.
 * If "x" isn't existent, return "true", else return "false".
 */
template <typename Object, typename Key>
inline bool BSTree<Object, Key>::insert( const Object &x )
{
    return insert( root, x );
}


/**
 * Delete the element x, which key is "k", in the tree.
 * If such key is existent, copy the deleted element to "x",
 * and return "true"; else return "false".
 */
template <typename Object, typename Key>
inline bool BSTree<Object, Key>::remove( Key k )
{
    return remove( root, k );
}


/**
 * Finding an element with key value of "k" in the tree.
 */
template <typename Object, typename Key>
inline BSNode<Object>* BSTree<Object, Key>::search( Key k )
{
    return search( root, k );
}


/**
 * Return smallest item in the tree.
 */
template <typename Object, typename Key>
Object BSTree<Object, Key>::minItem()
{
    return minItem(root)->data;
}


/**
 * Return largest item in the tree.
 */
template <typename Object, typename Key>
Object BSTree<Object, Key>::maxItem()
{
    return maxItem(root)->data;
}


/**
 * Destroy the tree.
 */
template <typename Object, typename Key>
void BSTree<Object, Key>::destroy( BSNode<Object> *&r )
{
    if( r != NULL )
    {
        destroy( r->left );
        destroy( r->right );
        delete r;
    }
    r = NULL;
}


/**
 * Preorder traversa the tree by nonrecursion method.
 */
template <typename Object, typename Key>
void BSTree<Object, Key>::preTravNonrecursion( BSNode<Object> *t )
{
    Stack< BSNode<Object>* > s;

    while( (t != NULL) || !s.isEmpty() )
    {
        if( t != NULL )
        {
            visit( t );
            s.push( t );
            t = t->left;
        }
        else
        {
            s.getTop( t );
            s.pop();
            t = t->right;
        }
    }
}


/**
 * Inorder traversa the tree by nonrecursion method.
 */
template <typename Object, typename Key>
void BSTree<Object, Key>::inTravNonrecursion( BSNode<Object> *t )
{
    Stack< BSNode<Object>* > s;

    while( (t != NULL) || !s.isEmpty() )
    {
        if( t != NULL )
        {
            s.push( t );
            t = t->left;
        }
        else
        {
            s.getTop( t );
            visit(t);
            s.pop();
            t = t->right;
        }
    }
}


/**
 * Postorder traversa the tree by nonrecursion method.
 */
template <typename Object, typename Key>
void BSTree<Object, Key>::postTravNonrecursion( BSNode<Object> *t )
{
    Stack< BSNode<Object>* > s;
    BSNode<Object>  *pre = NULL,        //last traversed node
                    *top = NULL;        //current traversed node

    while( (t != NULL) || !s.isEmpty() )
    {
        if( t != NULL )
        {
            s.push( t );
            t = t->left;
        }
        else
        {
            s.getTop( top );
            if( top->right != NULL && top->right != pre )
                t = top->right;
            else
            {
                visit( top );
                pre = top;
                s.pop();
            }
        }
    }
}


/**
 * Insert elemetn "x" into the tree "r".
 */
template <typename Object, typename Key>
bool BSTree<Object, Key>::insert( BSNode<Object> *&r, const Object &x )
{
    BSNode<Object>  *cur = r,
                    *par = NULL;

    while( cur != NULL )
    {
        if( cur->data == x )
            break;
        else
        {
            par = cur;
            if( cur->data < x )
                cur = cur->right;
            else if( cur->data > x )
                cur = cur->left;
        }
    }

    if( cur != NULL )
        return false;
    else
    {
        BSNode<Object> *tmp = new BSNode<Object>( x, NULL, NULL );
        if( par == NULL )
            r = tmp;
        else
        {
            if( par->data < x )
                par->right = tmp;
            else
                par->left = tmp;
        }
        return true;
    }
}


/**
 * Remove elemetn with key "k" from the tree "r".
 */
template <typename Object, typename Key>
bool BSTree<Object, Key>::remove( BSNode<Object> *&r, Key k )
{
    if( r == NULL )
        return false;

    if( k < r->data.key )
        return remove( r->left, k );
    else if( r->data.key < k )
        return remove( r->right, k );
    else if( r->left != NULL && r->right != NULL )
    {
        r->data = minItem( r->right )->data;
        return remove( r->right, r->data.key );
    }
    else
    {
        BSNode<Object> *oldNode = r;
        r = ( r->left != NULL ) ? r->left : r->right;
        delete oldNode;
        return true;
    }
}


/**
 * Remove elemetn with key "k" from the tree "r".
 */
template <typename Object, typename Key>
BSNode<Object>* BSTree<Object, Key>::search( BSNode<Object> *r, Key k )
{
    while( r != NULL )
        if( k < r->data.key )
            r = r->left;
        else if( k > r->data.key )
            r = r->right;
        else
            return r;

    return NULL;
}


/**
 * Find the smallest element in the tree "r".
 */
template <typename Object, typename Key>
BSNode<Object>* BSTree<Object, Key>::minItem( BSNode<Object> *r )
{
    if( r == NULL )
        return NULL;

    while( r->left != NULL )
        r = r->left;
    return r;
}


/**
 * Find the maximum element in the tree "r".
 */
template <typename Object, typename Key>
BSNode<Object>* BSTree<Object, Key>::maxItem( BSNode<Object> *r )
{
    if( r == NULL )
        return NULL;

    while( r->right != NULL )
        r = r->right;
    return r;
}


/**
 * Visit the item pointed by "p".
 */
template <typename Object, typename Key>
inline void BSTree<Object, Key>::visit( BSNode<Object> *p )
{
    cout << p->data;
}
