/*****************************************************************************
 *                               stack_test.cpp
 *
 * Stack class testing.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <string>
#include <utility>
#include <stack.h>


using namespace std;
using namespace splab;


int main()
{
    pair<int, string > stu;

    Stack< pair<int, string> > s;

    s.push( make_pair(1, "XiaoMing") );
    s.push( make_pair(2, "XiaoHong") );
    s.push( make_pair(3, "XiaoLing") );

    s.getTop( stu );
    cout << stu.first << "\t" << stu.second << endl;

    Stack< pair<int, string> > s1(s), s2;
    s2 = s1;

    cout << endl;
    while( !s2.isEmpty() )
    {
        s2.pop( stu );
        cout << stu.first << "\t" << stu.second <<endl;
    }
    cout << endl;

    s2.getTop( stu );
    cout << stu.first << "\t" << stu.second <<endl;

    cout << endl;
    return 0;
}
