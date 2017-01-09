/*****************************************************************************
 *                              huffmancode_test.cpp
 *
 * Huffman Code class testing.
 *
 * Zhang Ming, 2010-07, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <string>
#include <huffmancode.h>


using namespace std;
using namespace splab;


const int CHARNUM = 128;
template <typename Object, typename Weight>
void creatCodeInfo( string &str, CodeObject<Object,Weight> *&codeArr, int &length );
template <typename Object, typename Weight>
void printCodeInfo( CodeObject<Object,Weight> *&codeArr, int length );


int main()
{
    int length = 0;
    CodeObject<char, int> *codeArr = NULL;
    string str = "AAAABBCD";
    creatCodeInfo( str, codeArr, length );
    printCodeInfo( codeArr, length );

    HuffmanTree<char, int> ht;
    ht.code( codeArr, length );
    ht.printCodeTable();

    char decodedResult;
    unsigned char codeword[CODESIZE];
    cout << "Coding\t" << "Decodeding" << endl;
    codeword[0] = 0xC0;
    if( ht.decode( codeword, 3, decodedResult ) )
    {
        ht.printCode( codeword, 3 );
        cout << "\t" << decodedResult << endl;
    }
    else
        cerr << "invalid codeword! " << endl;
    codeword[0] = 0xC0;
    if( ht.decode( codeword, 4, decodedResult ) )
    {
        ht.printCode( codeword, 3 );
        cout << "\t" << decodedResult << endl;
    }
    else
    {
        ht.printCode( codeword, 4 );
        cerr << "\t" << "invalid codeword! " << endl;
    }
    codeword[0] = 0x40;
    if( ht.decode( codeword, 3, decodedResult ) )
    {
        ht.printCode( codeword, 3 );
        cout << "\t" << decodedResult << endl;
    }
    else
    {
        ht.printCode( codeword, 3 );
        cerr << "\t" << "invalid codeword! " << endl;
    }

    cout << endl;
    return 0;
}


/**
 * Generate the array for counting character frequency.
 */
template <typename Object, typename Weight>
void creatCodeInfo( string &str, CodeObject<Object,Weight> *&codeArr, int &length )
{
    char charFreq[CHARNUM];
    int index[CHARNUM];
    for( int i=0; i<CHARNUM; ++i )
    {
        charFreq[i] = 0;
        index[i] = -1;
    }

    length = 0;
    for( unsigned int i=0; i<str.size(); ++i )
        charFreq[int(str[i])] += 1;

    for( int i=0; i<CHARNUM; ++i )
        if( charFreq[i] )
            index[length++] = i;


    codeArr = new CodeObject<Object,Weight>[length];
    for( int i=0; i<length; ++i )
    {
        codeArr[i].data = index[i];
        codeArr[i].cost = charFreq[index[i]];
    }
}


/**
 * Print the characters and there frequency.
 */
template <typename Object, typename Weight>
void printCodeInfo( CodeObject<Object,Weight> *&codeArr, int length )
{
    cout << "Object\tFrequency" << endl;
    for( int i=0; i<length; ++i )
        cout << codeArr[i].data << "\t" << codeArr[i].cost << endl;
    cout << endl;
}
