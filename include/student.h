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
 *                                    student.h
 *
 * A student class including comparison operators.
 *
 * Zhang Ming, 2009-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef STUDENT_H
#define STUDENT_H


#include <iostream>
#include <string>


using namespace std;


namespace splab
{

    class Student
    {

    public:

        Student() : key(0), firstName("Ming"), lastName("Zhang")
        { }

        Student( int number, const string &name1="Ming",
                             const string &name2="Zhang" )
        {
            key = number;
            firstName = name1;
            lastName = name2;
        }

        ~Student()
        { }

        inline void operator=( const Student &stu )
        {
            key = stu.key;
            firstName = stu.firstName;
            lastName = stu.lastName;
        }

        inline bool operator<( const Student &stu )
        {   return key < stu.key;   }

        inline bool operator>( const Student &stu )
        {   return key > stu.key;   }

        inline bool operator==( const Student &stu )
        {   return key == stu.key;   }

        friend istream& operator>>( istream &in, Student &stu )
        {
            in >> stu.key >> stu.lastName >> stu.firstName;
            return in;
        }

        friend ostream& operator<<( ostream &out, Student &stu )
        {
            out << stu.key << "\t"
                << stu.lastName << " " << stu.firstName << endl;
            return out;
        }

        int key;

    private:

        string firstName;
        string lastName;

    };
    // class Student

}
// namespace splab


#endif
// STUDENT_H
