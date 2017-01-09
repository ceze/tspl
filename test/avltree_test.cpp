/*****************************************************************************
 *                               avltree_test.cpp
 *
 * AVL tree class testing.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <student.h>
#include <avltree.h>


using namespace std;
using namespace splab;


int main()
{
    int x[16] = { 3, 2, 1, 4, 5, 6, 7, 16, 15, 14, 13, 12, 11, 10, 8, 9 };
    int y[16] = { 9, 4, 12, 6, 5, 2, 3, 15, 14, 7, 8, 1, 1, 3, 20, 12 };
    int z[4] = { 10, 13, 5, 1 };

    Student stu;
    AVLNode<Student, int> *pNode;
    AVLTree<Student, int> stuTree;

    for( int i=0; i<16; ++i )
    {
        stu.key = x[i];
        stuTree.insert( stu );
    }

    cout << "Preorder Travesal: " << endl;
    stuTree.print( "preorder" );
    cout << endl << endl;

    for( int i=0; i<16; ++i )
    {
        if( stuTree.remove( y[i], stu ) )
        {
            cout << "The removed item is:  "<< stu;
            cout << "Preorder Travesal: " << endl;
            stuTree.print( "preorder" );
        }
        else
        {
            cout << "No such item (key=" << y[i] << ") in the tree!";
        }
        cout << endl << endl;

    }

    cout << endl;
    for( int i=0; i<4; ++i)
    {
        pNode = stuTree.search( z[i] );
        if( pNode )
            cout << "Have finding the element: " << pNode->data;
        else
            cout << "No such item (key=" << z[i] << ") in the tree!";
        cout << endl;
    }
    cout << endl;

    stuTree.makeEmpty();
    pNode = stuTree.search( 10 );

	return 0;
}
