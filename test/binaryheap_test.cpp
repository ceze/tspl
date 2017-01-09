/*****************************************************************************
 *                              binaryheap_test.cpp
 *
 * Minimum binary heap class testing.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <random.h>
#include <binaryheap.h>


using namespace std;
using namespace splab;


int main()
{
    int x, a[100];
    BinaryHeap<int> h;

    for( int i=1; i<=25; ++i )
        h.insert( randu( 1, 1, 100 ) );

    BinaryHeap<int> h1(h);
    while( !h1.isEmpty() )
    {
        h1.findMin( x );
        cout << "  " << x << "\t";
        h1.deleteMin();
    }
    cout << endl;

    cout << endl << endl;
    for( int i=0; i<40; ++i )
    {
        a[i] = randu( 1, 1, 100 );
        cout << "  " << a[i] << "\t";
    }
    cout << endl;

    BinaryHeap<int> hh( a, 40 ), hh1;
    hh1 = hh;
    while( !hh1.isEmpty() )
    {
        hh1.deleteMin( x );
        cout << "  " << x << "\t";
    }
    cout << endl;

    for( int i=1; i<=5; ++i )
        h1.insert( randu( 1, 1, 100 ) );
    h1.makeEmpty();
    h1.deleteMin();

    cout << endl << endl;
    return 0;
}
