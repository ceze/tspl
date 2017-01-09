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
 *                               avltree-impl.h
 *
 * Implementation for AVLTree class.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template <typename Object, typename Key>
AVLTree<Object, Key>::AVLTree() : root(NULL)
{ }

template <typename Object, typename Key>
AVLTree<Object, Key>::~AVLTree()
{
    delete root;
}


/**
 * If the tree is empty, return true.
 */
template <typename Object, typename Key>
inline bool AVLTree<Object, Key>::isEmpty() const
{
    return ( root == NULL );
}


/**
 * Make the tree empty.
 */
template <typename Object, typename Key>
inline void AVLTree<Object, Key>::makeEmpty()
{
    makeEmpty( root );
}


/**
 * Traverse the tree, the order is specified by parameter
 * "mode", which can be "preorder", "inorder" or "postorder".
 */
template <typename Object, typename Key>
void AVLTree<Object, Key>::print( const string &mode )
{
    if( isEmpty() )
        cout << "The tree is empty!" << endl;
    else
    {
        if( mode == "preorder" )
            preTraversal( root );
        else if( mode == "inorder" )
            inTraversal( root );
        else if( mode == "postorder" )
            postTraversal( root );
        else
            cout << "Invalid travesal method!" << endl;
    }
}


/**
 * Return the high of the tree.
 */
template <typename Object, typename Key>
inline int AVLTree<Object, Key>::height() const
{
    return height( root );
}


/**
 * Searching a node with key "k" and then return a pointer
 * pointed to it.
 */
template <typename Object, typename Key>
inline AVLNode<Object, Key>* AVLTree<Object, Key>::search( Key k )
{
    if( !isEmpty() )
        return search( k, root );
    else
    {
        handleUnderflow();
        return NULL;
    }
}


/**
 * Inserting an element into the tree.
 */
template <typename Object, typename Key>
inline bool AVLTree<Object, Key>::insert( Object &x )
{
    return insert( root, x );
}


/**
 * Removing an element with key "k" and then assign it to "x".
 */
template <typename Object, typename Key>
inline bool AVLTree<Object, Key>::remove( Key k, Object &x )
{
    if( !isEmpty() )
        return remove( root, k, x );
    else
    {
        handleUnderflow();
        return false;
    }
}


/**
 * Make the tree empty.
 */
template <typename Object, typename Key>
inline void AVLTree<Object, Key>::makeEmpty( AVLNode<Object, Key> * &t )
{
    if( t != NULL )
    {
        makeEmpty( t->left );
        makeEmpty( t->right );
        delete t;
    }
    t = NULL;
}


/**
 * preorder traversal
 */
template <typename Object, typename Key>
inline void AVLTree<Object, Key>::preTraversal( AVLNode<Object, Key> *t )
{
	if( t != NULL )
	{
        cout << "  " << t->balance << "\t" << t->data;
        preTraversal( t->left );
		preTraversal( t->right );
	}
};


/**
 * inorder traversal
 */
template <typename Object, typename Key>
inline void AVLTree<Object, Key>::inTraversal( AVLNode<Object, Key> *t )
{
	if( t != NULL )
	{
		inTraversal( t->left );
		cout << "  " << t->data;
		inTraversal( t->right );
	}
};


/**
 * postorder traversal
 */
template <typename Object, typename Key>
inline void AVLTree<Object, Key>::postTraversal( AVLNode<Object, Key> *t )
{
	if( t != NULL )
	{
        postTraversal( t->left );
		postTraversal( t->right );
		cout << "  " << t->data;
	}
};


/**
 * Compute the high of the subtree with root of "t".
 */
template <typename Object, typename Key>
inline int AVLTree<Object, Key>::height( AVLNode<Object, Key> *t ) const
{
	if( !t )
		return -1;

	int heightL = height( t->left );
	int heightR = height( t->right );

	return ( heightL < heightR ) ? heightR+1 : heightL+1;
}


/**
 * Finding an element with key value of "k" in the
 * subtree with root of "t".
 */
template <typename Object, typename Key>
AVLNode<Object, Key>* AVLTree<Object, Key>::search( Key k,
                      AVLNode<Object, Key> *t )
{
	if( t == NULL )
		return NULL;

	if( k == t->data.key )
		return t;
    else if( k < t->data.key )
        return search( k, t->left );
    else
        return search( k, t->right );
}


/**
 * Insert a new element x in the tree "t".
 * If "x" isn't existent, return "true", else return "false".
 */
template <typename Object, typename Key>
bool AVLTree<Object, Key>::insert( AVLNode<Object, Key> * &t,
                                   Object &x )
{
    int flag;       // flag for sigle or double rotation

	AVLNode<Object, Key> *child = t,
                         *parent = NULL,
                         *grandpa = NULL;

	Stack< AVLNode<Object, Key> * > st;

    // find the inserting position
	while( child != NULL )
	{
		if( x == child->data )
            return false;

		parent = child;
		st.push( parent );

		if( x < parent->data )
            child = parent->left;
		else
            child = parent->right;
	}

	child = new AVLNode<Object, Key>( x );
	if( child == NULL )
	{
	    cerr << "Out of memory!" << endl;
	    exit( 1 );
    }

    // The tree is empty, the new element is the root of the tree.
	if( parent == NULL )
	{
	    t = child;
	    return true;
    }

	if( x < parent->data )
        parent->left = child;
	else
        parent->right = child;

    // equalization
	while( !st.isEmpty() )
	{
		st.pop( parent );

        // Modify the balance factor of node "parent".
		if( child == parent->left )
            parent->balance--;
		else
            parent->balance++;

        // have equalized
		if( parent->balance == 0 )
		{
            break;
		}

        // Current node is balance, but need backtracking upward.
		else if( parent->balance == 1 || parent->balance == -1 )
		{
			child = parent;
		}

        // Current node is unbalance, doing equalization.
		else
		{
			flag = ( parent->balance < 0 ) ? -1 : 1;

            // Current node and its parent have the same balance
            // factor, so making sigle rotation.
			if( child->balance == flag )
			{
				if( flag == -1 )
                    rotateR( parent );
				else
                    rotateL( parent );
			}

			// Current node and its parent have the opposite balance
            // factor, so making double rotation.
			else
			{
				if( flag == -1 )
                    rotateLR( parent );
				else
                    rotateRL( parent );
			}

			// Over the equalization.
			break;
		}
	}

    // The equalized node is the root of the tree.
	if( st.isEmpty() )
	{
        t = parent;
	}

	// The equalized node is at the middle of the tree.
	else
	{
		st.getTop( grandpa );
		if( grandpa->data > parent->data )
            grandpa->left = parent;
		else
            grandpa->right = parent;
	}

	return true;
};


/**
 * Delete the element x, which key is "k", in the tree "t".
 * If such key is existent, copy the deleted element to "x",
 * and return "true"; else return "false".
 */
template <typename Object, typename Key>
bool AVLTree<Object, Key>::remove( AVLNode<Object, Key> * &t,
                                   Key k, Object &x )
{
    int flag = 0;       // 1 : right subtree higher than left subtre
                        // -1 : left subtree higher than right subtre

    int leftRight = 0;  // 1 : right subtree; -1 : left subtree

	AVLNode<Object, Key> *current = t,
                         *parent = NULL,
                         *grandpa = NULL,
                         *child = NULL;

	Stack< AVLNode<Object, Key> * > st;

    // Finding the deleted position.
	while( current != NULL )
	{
		if( k == current->data.key )
		{
            x = current->data;
            break;
		}

		parent = current;
		st.push( parent );

		if( k < parent->data.key )
            current = parent->left;
		else
            current = parent->right;
	}

    // No such node in the tree, deletion false.
	if( current == NULL )
        return false;

    // The deleted node has tow children.
	if( ( current->left != NULL ) && ( current->right != NULL) )
	{
		parent = current;
		st.push( parent );

        // Find the immediate predecessor of "current" pointing
        // node in its left subtree.
		child = current->left;
		while( child->right != NULL )
		{
		    parent = child;
		    st.push( parent );
		    child = parent->right;
        }

        // Copy "child" to "current", and the deleted node have
        // translated into the new node "child".
		current->data = child->data;
		current = child;
	}

    // The deleted node has no more than one child.
	if( current->left != NULL )
        child = current->left;
	else
        child = current->right;

    // The deleted node is the root.
	if( parent == NULL )
	{
        t = child;
	}

	// The deleted node is leaf or has only on child.
	else
	{
	    // relinking
		if( parent->left == current )
            parent->left = child;
		else
            parent->right = child;

        // equalization
		while( !st.isEmpty() )
		{
			st.pop( parent );

            // Modify the balance factor of "parent" node.
			if( parent->right == child )
                parent->balance--;
			else
                parent->balance++;

            // The stack is empty, so do not need link upward
            // after rotation.
			if( st.isEmpty() )
			{
			    leftRight = 0;
			}

			// The stack isn't empty, so need link upward after rotation.
			else
			{
			    st.getTop( grandpa );
				leftRight = ( grandpa->left == parent ) ? -1 : 1;
			}

            // The balance of the node pointed by "parent" is 0,
            // so do not need modifying after deletion.
			if( parent->balance == 1 || parent->balance == -1 )
                break;

            // The node pointed by "parent" is unbalance.
			if( parent->balance != 0 )
			{
                // "child" pointing the higher subtree.
				if( parent->balance < 0)
				{
				    flag = -1;
				    child = parent->left;
                }
				else
				{
				    flag = 1;
				    child = parent->right;
                }

                // The left and right subtree of the higher subtree have
                // the same height.
				if( child->balance == 0 )
				{
				    // The higher subtree is the left subtree of "parent".
					if( flag == -1 )
					{
					    rotateR( parent );
					    parent->balance = 1;
                        parent->right->balance = -1;
                    }

                    // The higher subtree is the right subtree of "parent".
					else
					{
					    rotateL( parent );
					    parent->balance = -1;
                        parent->left->balance = 1;
                    }
					break;
				}

                // The higher subtree have the same balance factor
                // of its parent, so need sigle rotation.
				if( child->balance == flag )
				{
					if( flag == -1 )
                        rotateR( parent );
					else
                        rotateL( parent );
				}

				// The higher subtree have the opposite balance factor
                // of its parent, so need double rotation.
				else
				{
					if( flag == -1 )
                        rotateLR( parent );
					else
                        rotateRL( parent );
				}

                // Linking upward after rotatin.
				if( leftRight == -1 )
                    grandpa->left = parent;
				else if( leftRight == 1 )
                    grandpa->right = parent;
			}

            // backtracking upward
			child = parent;
		}

        // Have modified to the root.
		if( st.isEmpty() )
            t = parent;
	}

	delete current;
	return true;
};


/**
 * Making left sigle rotation fo node pointed by "t",
 * and the root of the rotated tree is still pointed by "t".
 */
template <typename Object, typename Key>
void AVLTree<Object, Key>::rotateL( AVLNode<Object, Key> * &t )
{
	AVLNode<Object, Key> *leftChild = t;

	t = leftChild->right;
	leftChild->right = t->left;
	t->left = leftChild;

	t->balance = leftChild->balance = 0;
};


/**
 * Making right sigle rotation fo node pointed by "t",
 * and the root of the rotated tree is still pointed by "t".
 */
template <typename Object, typename Key>
void AVLTree<Object, Key>::rotateR( AVLNode<Object, Key> * &t )
{
	AVLNode<Object, Key> *rightChild = t;

	t = rightChild->left;
	rightChild->left = t->right;
	t->right = rightChild;

	t->balance = rightChild->balance = 0;
};


/**
 * Making left-right double rotation fo node pointed by "t",
 * and the root of the rotated tree is still pointed by "t".
 */
template <typename Object, typename Key>
void AVLTree<Object, Key>::rotateLR( AVLNode<Object, Key> * &t)
{
	AVLNode<Object, Key> *rightChild = t;
	AVLNode<Object, Key> *leftChild = rightChild->left;
	t = leftChild->right;

	leftChild->right = t->left;
	t->left = leftChild;
	if( t->balance <= 0 )
        leftChild->balance = 0;
	else
        leftChild->balance = -1;

	rightChild->left = t->right;
	t->right = rightChild;
	if( t->balance == -1 )
        rightChild->balance = 1;
	else
        rightChild->balance = 0;

    t->balance = 0;
};


/**
 * Making right-left double rotation fo node pointed by "t",
 * and the root of the rotated tree is still pointed by "t".
 */
template <typename Object, typename Key>
void AVLTree<Object, Key>::rotateRL( AVLNode<Object, Key> * &t )
{
	AVLNode<Object, Key> *leftChild = t;
	AVLNode<Object, Key> *rightChild = leftChild->right;
	t = rightChild->left;

	rightChild->left = t->right;
	t->right = rightChild;
	if( t->balance >= 0 )
        rightChild->balance = 0;
	else
        rightChild->balance = 1;

	leftChild->right = t->left;
	t->left = leftChild;
	if( t->balance == 1 )
        leftChild->balance = -1;
	else
        leftChild->balance = 0;

	t->balance = 0;
};


/**
 * Handle the error of get element from an empty tree.
 */
template <typename Object, typename Key>
inline void AVLTree<Object, Key>::handleUnderflow()
{
    cerr << "The tree is empty!" << endl << endl;
    exit( 1 );
}
