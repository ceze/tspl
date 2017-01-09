/*****************************************************************************
 *                               queue_test.cpp
 *
 * Queue class testing.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <string>
#include <queue.h>


using namespace std;
using namespace splab;


struct Student
{
    int     stuNum;
    string  stuName;
    string  department;

    Student( int number = 0, const string &name = "Tom & Jerry",
             const string &dpt = "Information" )
    : stuNum(number), stuName(name), department(dpt)
    { }

};


int main()
{
    Student stu;
    Queue<Student> q;

    q.enqueue( Student(11, "Zhang Ming", "Information") );
    q.enqueue( Student(16, "Hu Zhaojun") );
    q.enqueue( Student(11, "Zhang Yuyang", "AutoControl") );
    q.enqueue( Student() );


    q.getFront( stu );
    cout << "  " << stu.stuNum << "\t" << stu.stuName << "\t"
         << stu.department << endl;

    cout << endl;
    q.dequeue();
    while( !q.isEmpty() )
    {
        q.dequeue( stu );
        cout << "  " << stu.stuNum << "\t" << stu.stuName << "\t"
             << stu.department << endl;
    }
    cout << endl;

    q.getFront( stu );
    cout << "  " << stu.stuNum << "\t" << stu.stuName << "\t"
         << stu.department;

    cout << endl << endl;
    return 0;
}
