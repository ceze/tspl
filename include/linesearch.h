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
 *                               linesearch.h
 *
 * Line search algorithm based on quadratic polynomial interpolation.
 *
 * This routine is often used in orther optimation method, so it is designed
 * as a templated base class.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef LINESEARCH_H
#define LINESEARCH_H


#include <vector.h>


namespace splab
{

    template <typename Dtype, typename Ftype>
    class LineSearch
    {

     public:

        LineSearch();
        ~LineSearch();

        Dtype getStep( Ftype &func, Vector<Dtype> &xk, Vector<Dtype> &dk,
                       int maxItr=10 );

        int getFuncNum() const;
        bool isSuccess() const;

    protected:

        int     funcNum;
        bool    success;

    };
    // class LineSearch


    #include <linesearch-impl.h>

}
// namespace splab


#endif
// LINESEARCH_H
