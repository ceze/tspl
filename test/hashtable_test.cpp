/*****************************************************************************
 *                               hashtable_test.cpp
 *
 * Hash table class testing.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <string>
#include <student.h>
#include <hashtable.h>


using namespace std;
using namespace splab;


int const N = 10;


int main()
{
    int x[N] = { 3, 2, 1, 4, 5, 6, 7, 10, 9, 8 };
    int y[N] = { 3, 2, 1, 4, 5, 16, 17, 20, 19, 18 };
    Student stu;
    HashTable<Student, int> ht;

    for( int i=0; i<N; ++i )
    {
        stu.key = x[i];
        if( ht.insert( stu ) )
            cout << "Inserting success: " << stu;
        else
            cout << "Inserting failure." << endl;
    }
    cout << endl;

    HashTable<Student, int> ht1(ht), ht2;
    for( int i=0; i<N; ++i )
        if( ht1.remove( y[i], stu ) )
            cout << "Removing success: " << stu;
        else
            cout << "Removing failure." << endl;
    cout << endl;

    ht2 = ht1;
    for( int i=0; i<N; ++i )
        if( ht2.search( x[i], stu ) )
            cout << "Searching success: " << stu;
        else
            cout << "Searching failure." << endl;
    cout << endl;

    return 0;
}
