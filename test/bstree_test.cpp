/*****************************************************************************
 *                               bstree_test.cpp
 *
 * Binary search tree class testing.
 *
 * Zhang Ming, 2010-07, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <student.h>
#include <bstree.h>


using namespace std;
using namespace splab;


int main()
{
    const int N = 8,
              M = 4;

    int x[N] = { 3, 1, 2, 7, 5, 8, 4, 6 };
    int y[N] = { 3, 2, 5, 7, 7, 9, 5, 0 };
    int z[M] = { 3, 2, 1, 7 };
    Student stu;
    BSNode<Student> *pNode = NULL;
    BSTree<Student, int> stuTree;

    for( int i=0; i<N; ++i )
    {
        stu.key = x[i];
        if( stuTree.insert( stu ) )
            cout << "Insert " << x[i] << " success." << endl;
        else
            cout << "Insert failing." << endl;
    }
    cout << endl;

    cout << "Preorder Travesal: " << endl;
    stuTree.preTraverse();
    cout << endl << "Inorder Travesal: " << endl;
    stuTree.inTraverse();
    cout << endl  << "Postorder Travesal: " << endl;
    stuTree.postTraverse();
    cout << endl;

    for( int i=0; i<N; ++i )
    {
        if( stuTree.remove( y[i] ) )
        {
            cout << "Remove " << y[i] <<" sucess. Preorder Travesal: "<< endl;
//            stuTree.preTraverse();
            stuTree.inTraverse();
//            stuTree.postTraverse();
        }
        else
            cout << "No such item (key=" << y[i] << ") in the tree!" << endl;
    }

    cout << endl;
    for( int i=0; i<M; ++i)
    {
        pNode = stuTree.search( z[i] );
        if( pNode )
            cout << "Finding the element: " << pNode->data;
        else
            cout << "No such item (key=" << z[i] << ") in the tree!" << endl;
    }
    cout << endl;

    stuTree.makeEmpty();
    if( stuTree.isEmpty() )
        cout << "The tree is empty." << endl;
    else
        cerr << "Making tree empty failing!" << endl;

    cout << endl;
    return 0;
}
