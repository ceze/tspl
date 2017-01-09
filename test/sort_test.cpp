/*****************************************************************************
 *                                sort_test.cpp
 *
 * Sorting algorithm testing.
 *
 * Zhang Ming, 2010-07, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <random.h>
#include <sort.h>


using namespace std;
using namespace splab;


const int SIZE = 10;


template <typename Type>
void printVector( const Vector<Type> &a )
{
    Vector<int>::const_iterator itr = a.begin();
    while( itr != a.end() )
        cout << *itr++ << "\t";

    cout << endl;
}


int main()
{
    Vector<int> a( SIZE );

    cout << "Unsorted Numbers : " << endl;
    a = randu( 127, 1, 10*SIZE, SIZE );
    printVector( a );
    cout << "Bubble Sorted Numbers : " << endl;
    bubbleSort( a, 0, a.size()-1 );
    printVector( a );
    cout << endl;

    cout << "Unsorted Numbers : " << endl;
    a = randu( 13579, 1, 10*SIZE, SIZE );
    printVector( a );
    cout << "Select Sorted Numbers : " << endl;
    selectSort( a, 0, a.size()-1 );
    printVector( a );
    cout << endl;

    cout << "Unsorted Numbers : " << endl;
    a = randu( 37, 1, 10*SIZE, SIZE );
    printVector( a );
    cout << "Insert Sorted Numbers : " << endl;
    insertSort( a, 0, a.size()-1 );
    printVector( a );
    cout << endl;

    cout << "Unsorted Numbers : " << endl;
    for( int i=0; i<a.size(); ++i )
        a[i] = randu( 127, 1, 10*SIZE );
    printVector( a );
    cout << "Quick Sorted Numbers : " << endl;
    quickSort( a, 0, a.size()-1 );
    printVector( a );
    cout << endl;

    cout << "Unsorted Numbers : " << endl;
    for( int i=0; i<a.size(); ++i )
        a[i] = randu( 127, 1, 10*SIZE );
    printVector( a );
    cout << "Merg Sorted Numbers : " << endl;
    mergSort( a, 0, a.size()-1 );
    printVector( a );
    cout << endl;

    cout << "Unsorted Numbers : " << endl;
    for( int i=0; i<a.size(); ++i )
        a[i] = randu( 127, 1, 10*SIZE );
    printVector( a );
    cout << "Heap Sorted Numbers : " << endl;
    heapSort( a, 0, a.size()-1 );
    printVector( a );
    cout << endl;

    return 0;
}
