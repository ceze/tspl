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
 *                              matrixallocate-impl.h
 *
 * Implementation matrix allocate and free.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Allocate matrix by specified rows and columns.
 */
template<typename Type>
inline void makeMatrix( Type **&m, int rows, int cols )
{
    m = new Type*[rows];
    for( int i=0; i<rows; ++i )
        m[i] = new Type[cols];
}


/**
 * Delete the matrix with specified rows.
 */
template<typename Type>
inline void deleteMatrix( Type **&m, int rows )
{
    for( int i=0; i<rows; ++i )
        delete []m[i];

    delete []m;
    m = 0;
}
