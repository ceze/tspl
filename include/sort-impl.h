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
 *                                  sort-impl.h
 *
 * Implementation for sort algorithms.
 *
 * Zhang Ming, 2010-07, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Bubble sort algorithm.
 * "a"      ---->   array of Comparable items.
 * "left"   ---->   the left-most index of the subarray.
 * "right"  ---->   the right-most index of the subarray.
 */
 template <typename Type>
 void bubbleSort( Vector<Type> &a, int left, int right )
 {
     bool cond;
     for( int i=left; i<right; ++i )
     {
         cond = false;
         for( int j=right; j>i; --j )
             if( a[j] < a[j-1] )
             {
                 swap( a[j], a[j-1] );
                 cond = true;
             }

        if( !cond  )
            return;
     }
 }


/**
 * Selection sort algorithm.
 * "a"      ---->   array of Comparable items.
 * "left"   ---->   the left-most index of the subarray.
 * "right"  ---->   the right-most index of the subarray.
 */
 template <typename Type>
 void selectSort( Vector<Type> &a, int left, int right )
 {
     Type minPos;
     for( int i=left; i<right; ++i )
     {
         minPos = i;
         for( int j=i+1; j<=right; ++j )
             if( a[j] < a[minPos] )
                 minPos = j;

        if( i != minPos )
            swap( a[i], a[minPos] );
     }
 }


/**
 * Insertion sort algorithm.
 * "a"      ---->   array of Comparable items.
 * "left"   ---->   the left-most index of the subarray.
 * "right"  ---->   the right-most index of the subarray.
 */
template <typename Type>
void insertSort( Vector<Type> &a, int left, int right )
{
    for( int p=left+1; p<=right; p++ )
    {
        Type tmp = a[p];
        int j;

        for( j=p; j>left && tmp<a[j-1]; --j )
            a[j] = a[j-1];
        a[j] = tmp;
    }
}


/**
 * Internal quicksort method that makes recursive calls.
 * Uses median-of-three partitioning and a cutoff of 20.
 * "a"      ---->   array of Comparable items.
 * "left"   ---->   the left-most index of the subarray.
 * "right"  ---->   the right-most index of the subarray.
 */
template <typename Type>
void quickSort( Vector<Type> &a, int left, int right )
{
    if( left+20 <= right )
    {
        Type pivot = median3( a, left, right );

        // begin partitioning
        int i = left, j = right-1;
        for( ; ; )
        {
            while( a[++i] < pivot ) { }
            while( pivot < a[--j] ) { }

            if( i < j )
                swap( a[i], a[j] );
            else
                break;
        }

        // Restore pivot
        swap( a[i], a[right-1] );

        // Sort small elements
        quickSort( a, left, i-1 );

        // Sort large elements
        quickSort( a, i+1, right );
    }
    else
        insertSort( a, left, right );
}


/**
 * Merg sort algorithm (nonrecursion).
 * "a"      ---->   array of Comparable items.
 * "left"   ---->   the left-most index of the subarray.
 * "right"  ---->   the right-most index of the subarray.
 */
template <typename Type>
void mergSort( Vector<Type> &a, int left, int right )
{
    int left1, right1, left2, right2,
        n = right - left + 1,
        size = 1;

    while( size < n )
    {
        left1 = left;

        while( left1+size < n )
        {
            left2 = left1+size;
            right1 = left2-1;
            if( left2+size > n )
                right2 = right;
            else
                right2 = left2+size-1;

            merg( a, left1, right1, left2, right2 );

            left1 = right2+1;
        }

        size *= 2;
    }
}


/**
 * Heap sort algorthm.
 * "a"      ---->   array of Comparable items.
 * "left"   ---->   the left-most index of the subarray.
 * "right"  ---->   the right-most index of the subarray.
 */
template <typename Type>
void heapSort( Vector<Type> &a, int left, int right )
{
    int n = right-left+1;
    Vector<Type> tmp( n );
    for( int i=0; i<n; ++i )
        tmp[i] = a[left+i];

    for( int i=n/2; i>=0; --i )
        filterDown( tmp, i, n );
    for( int j=n-1; j>0; --j )
    {
        swap( tmp[0], tmp[j] );
        filterDown( tmp, 0, j );
    }

    for( int i=0; i<n; ++i )
        a[left+i] = tmp[i];
}


/**
 * Return median of left, center, and right.
 * Order these and hide the pivot.
 */
template <typename Type>
const Type& median3( Vector<Type> &a, int left, int right )
{
    int center = (left+right) / 2;

    if( a[center] < a[left] )
        swap( a[left], a[center] );

    if( a[right] < a[left] )
        swap( a[left], a[right] );

    if( a[right] < a[center] )
        swap( a[center], a[right] );

    swap( a[center], a[right-1] );

    return a[right-1];
}


/**
 * Merg two subsequence to a bigger one.
 * The first subsequence is a[left1] ... a[right1], and
 * The second subsqeuence is a[left2] ... a[right2].
 */
template <typename Type>
void merg( Vector<Type> &a, int left1, int right1, int left2, int right2 )
{
    int k = 0,
        i = left1,
        j = left2,
        n1 = right1-left1+1,
        n2 = right2-left2+1;

    Type *tmp = new Type[n1+n2];

    while( i <= right1 && j <= right2 )
        if( a[i] < a[j] )
            tmp[k++] = a[i++];
        else
            tmp[k++] = a[j++];

    while( i <= right1 )
        tmp[k++] = a[i++];
    while( j <= right2 )
        tmp[k++] = a[j++];

    for( int i=0; i<n1; ++i )
        a[left1++] = tmp[i];
    for( int i=0; i<n2; ++i )
        a[left2++] = tmp[n1+i];

    delete []tmp;
}


/**
 * Percolate down the heap.
 * "i"  ---->  the position from which to percolate down.
 * "n"  ---->  the logical size of the binary heap.
 */
template <typename Type>
void filterDown( Vector<Type> &a, int i, int n )
{
    int child;
    Type tmp;

    for( tmp=a[i]; 2*i+1<n; i=child )
    {
        child = 2*i+1;
        if( child!=n-1 && a[child]<a[child+1] )
            child++;

        if( tmp < a[child] )
            a[i] = a[child];
        else
            break;
    }
    a[i] = tmp;
}
