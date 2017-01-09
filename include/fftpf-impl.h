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
 *                               fftpf-impl.h
 *
 * Implementation for FFTPF class.
 *
 * Zhang Ming, 2010-09, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template<class Type>
FFTPF<Type>::FFTPF()
{
    bAlloc = false;
    mNOldSize = nNewFactorSize = nHalfFactorSize = 0;

    pftwRe = twiddleRe;
    pftwIm = twiddleIm;
	pftgRe = trigRe;
	pftgIm = trigIm;
	pfzRe = zRe;
	pfzIm = zIm;
	pfvRe = vRe;
	pfvIm = vIm;
	pfwRe = wRe;
	pfwIm = wIm;

	tIRT2 = Type(sqrt(2.0)/2);
	tIRT3 = Type(sqrt(3.0)/2);
	tPI   = Type(4*atan(1.0));
	t2PI  = 2*tPI;
	Type tphi = t2PI/5;

	c3_1 = -1.5;
	c5_1 = Type((cos(tphi)-cos(2*tphi))/2);
	c5_2 = Type(sin(tphi));
	c5_3 = Type(sin(2*tphi));
	c5_4 = Type((c5_3+c5_2));
	c5_3 = c5_2-c5_3;
	c5_5 = 1.25;

	tphi = t2PI/7;
	c7_1 = Type(cos(tphi));
	c7_2 = Type(cos(2*tphi));
	c7_3 = Type(cos(3*tphi));
	c7_4 = Type(sin(tphi));
	c7_5 = Type(sin(2*tphi));
	c7_6 = Type(sin(3*tphi));

	tphi = t2PI/9;
	c9_4 = Type(cos(tphi));
	c9_8 = Type(sin(8*tphi));
	c9_3 = Type(cos(2*tphi));
	c9_7 = Type(sin(2*tphi));
	c9_5 = Type(cos(3*tphi));
	c9_2 = Type(cos(4*tphi));
	c9_6 = Type(sin(5*tphi));

	tphi = t2PI/11;
	c11_1 = Type(cos(tphi));
	c11_2 = Type(cos(2*tphi));
	c11_3 = Type(cos(3*tphi));
	c11_4 = Type(cos(4*tphi));
	c11_5 = Type(cos(5*tphi));
	c11_6 = Type(sin(tphi));
	c11_7 = Type(sin(2*tphi));
	c11_8 = Type(sin(3*tphi));
	c11_9 = Type(sin(4*tphi));
	c11_10 = Type(sin(5*tphi));

	tphi = t2PI/13;
	c13_1 = Type(cos(tphi));
	c13_2 = Type(cos(2*tphi));
	c13_3 = Type(cos(3*tphi));
	c13_4 = Type(cos(4*tphi));
	c13_5 = Type(cos(5*tphi));
	c13_6 = Type(cos(6*tphi));
	c13_7 = Type(sin(tphi));
	c13_8 = Type(sin(2*tphi));
	c13_9 = Type(sin(3*tphi));
	c13_10 = Type(sin(4*tphi));
	c13_11 = Type(sin(5*tphi));
	c13_12 = Type(sin(6*tphi));

	tphi = t2PI/16;
	c16_2 = Type(sin(tphi));
	c16_5 = Type(cos(tphi));
	c16_3 = c16_2+c16_5;
	c16_4 = c16_5-c16_2;
}


template<class Type>
FFTPF<Type>::~FFTPF()
{
    releaseMem();
}


template<class Type>
void FFTPF<Type>::releaseMem()
{
	if ( !bAlloc)
        return;

	if ( NULL != pftwRe )
        delete []pftwRe;
	if ( NULL != pftwIm )
        delete []pftwIm;
	if ( NULL != pftgRe )
        delete []pftgRe;
	if ( NULL != pftgIm )
        delete []pftgIm;
	if ( NULL != pfzRe )
        delete []pfzRe;
	if ( NULL != pfzIm )
        delete []pfzIm;
	if ( NULL != pfvRe )
        delete []pfvRe;
	if ( NULL != pfvIm )
        delete []pfvIm;
	if ( NULL != pfwRe )
        delete []pfwRe;
	if ( NULL != pfwIm )
        delete []pfwIm;

	bAlloc = false;
	nNewFactorSize = nHalfFactorSize = 0;
}


/**
 * allocate memory
 */
template<class Type>
inline void FFTPF<Type>::allocateMem()
{
	pftwRe = new Type[nNewFactorSize];
	pftwIm = new Type[nNewFactorSize];
	pftgRe = new Type[nNewFactorSize];
	pftgIm = new Type[nNewFactorSize];
	pfzRe  = new Type[nNewFactorSize];
	pfzIm  = new Type[nNewFactorSize];
	pfvRe  = new Type[nHalfFactorSize];
	pfvIm  = new Type[nHalfFactorSize];
	pfwRe  = new Type[nHalfFactorSize];
	pfwIm  = new Type[nHalfFactorSize];
	bAlloc = true;
}


/**
 * factoring the transformation length.
 */
template<class Type>
void FFTPF<Type>::factorize( int n, int &nFact, int *fact )
{
	int i, j=0, k,
        nRadix = 11, radices[12],
        factors[PRIMECOUNT];

	radices[1] = 2;  radices[2] = 3;   radices[3] = 4;   radices[4] = 5;
	radices[5] = 7;  radices[6] = 8;   radices[7] = 9;   radices[8] = 10;
	radices[9] = 11; radices[10] = 13; radices[11] = 16;

	if( 1 == n )
	{
	    j = 1;
	    factors[1] = 1;
    }

	i = nRadix;
	while( (n>1) && (i>0) )
	{
		if( (n%radices[i]) != 0 )
            --i;
		else
		{
		    n = n/radices[i];
		    j++;
		    factors[j] = radices[i];
        }
	}

	if( 2 == factors[j] )
	{
		i = j-1;
		while( (i>0) && (factors[i]!=8) )
            i--;
		if( i > 0 )
		{
		    factors[j] = 4;
		    factors[i] = 4;
        }
	}

	if( n > 1 )
	{
		for( k=2; k<sqrt(1.0*n)+1; ++k )
            while( 0 == (n%k) )
            {
                n = n/k;
                j++;
                factors[j] = k;
            }
		if( n > 1 )
		{
		    j++;
		    factors[j] = n;
        }
	}

	for( i=1; i<=j; ++i )
        fact[i] = factors[j-i+1];

	nFact = j;
}


/**
 * After N is factored the parameters that control the stages are generated.
 * For each stage we have:
 * mSofarRadix  : the product of the radices so far.
 * mActualRadix : the radix handled in this stage.
 * mRemainRadix : the product of the remaining radices.
 */
template<class Type>
void FFTPF<Type>::primeSetup( int nPoints )
{
	factorize( nPoints, mNFactor, mActualRadix );

	if( mActualRadix[1] > PRIMEFACTOR )
	{
		bool bdone = true;
		if( bAlloc )
		{
			if( mActualRadix[1] <= nNewFactorSize )
                bdone=false;
			else
			{
			    releaseMem();
			    bdone=true;
            }
		}

		if( bdone )
		{
			nNewFactorSize  = mActualRadix[1];
			nHalfFactorSize = (nNewFactorSize+1)/2;
			allocateMem();
		}
	}

    mSofarRadix[1]  = 1;
	mRemainRadix[0] = nPoints;
	mRemainRadix[1] = nPoints/mActualRadix[1];

	for( int i=2; i<=mNFactor; ++i )
	{
		mSofarRadix[i]  = mSofarRadix[i-1] * mActualRadix[i-1];
		mRemainRadix[i] = mRemainRadix[i-1] / mActualRadix[i];
	}

	mNOldSize = nPoints;
}


/**
 * The sequence "yn" is the permuted input sequence "xn" so that the
 * following transformations can be performed in-place,and the final
 * result is the normal order.
 */
template<class Type>
void FFTPF<Type>::permute( const Vector<Type> &xn,
                           Vector< complex<Type> > &yn )
{
	int i, j, k=0, nPoint = xn.size(),
        count[PRIMECOUNT];

	for( i=1; i<=mNFactor; ++i )
        count[i] = 0;

	for( i=0; i<=nPoint-2; ++i )
	{
        yn[i] = xn[k];
		j = 1;
		k += mRemainRadix[j];
		count[1] = count[1]+1;
		while( count[j] >= mActualRadix[j] )
		{
			count[j] = 0;
			k = k - mRemainRadix[j-1] + mRemainRadix[j+1];
			j++;
			count[j] = count[j] + 1;
		}
	}
	yn[nPoint-1] = xn[nPoint-1];
}

template<class Type>
void FFTPF<Type>::permute( const Vector< complex<Type> > &xn,
                           Vector< complex<Type> > &yn,
                           bool bTrans )
{
	int i, j, k=0, nPoint = xn.size(),
        count[PRIMECOUNT];

	for( i=1; i<=mNFactor; ++i )
        count[i] = 0;

	for( i=0; i<=nPoint-2; ++i )
	{
		yn[i] = xn[k];
		j = 1;
		k += mRemainRadix[j];
		count[1] = count[1]+1;
		while( count[j] >= mActualRadix[j] )
		{
			count[j] = 0;
			k = k - mRemainRadix[j-1] + mRemainRadix[j+1];
			j++;
			count[j] = count[j] + 1;
		}
	}
	yn[nPoint-1] = xn[nPoint-1];

	if( !bTrans )
        for( i=0; i< nPoint; ++i )
            yn[i] = conj(yn[i]);
}


/**
 * initialise sine/cosine table.
 */
template<class Type>
inline void FFTPF<Type>::initTrig( int radix )
{
	Type w, xre, xim;

	w = 2*tPI/radix;
	pftgRe[0] = 1;
	pftgIm[0] = 0;
	xre = Type(cos(w));
	xim = -Type(sin(w));
	pftgRe[1] = xre;
	pftgIm[1] = xim;

	for( int i=2; i<radix; ++i )
	{
		pftgRe[i] = xre*pftgRe[i-1] - xim*pftgIm[i-1];
		pftgIm[i] = xim*pftgRe[i-1] + xre*pftgIm[i-1];
	}
}


/**
 * length 2 DFT
 */
template<class Type>
inline void FFTPF<Type>::radix2( Type *aRe, Type *aIm )
{
	t1_re  = aRe[0];
	aRe[0] = t1_re+aRe[1];
	aRe[1] = t1_re-aRe[1];

	t1_im  = aIm[0];
	aIm[0] = t1_im+aIm[1];
	aIm[1] = t1_im-aIm[1];
}


/**
 * length 3 DFT
 */
template<class Type>
inline void FFTPF<Type>::radix3( Type *aRe, Type *aIm )
{
	t1_re  = aRe[1]+aRe[2];
	t2_re  = (aRe[2]-aRe[1])*tIRT3;
	aRe[0] = aRe[0]+t1_re;

	t1_re  = aRe[0]+t1_re*c3_1;
	t1_im  = aIm[1]+aIm[2];
	t2_im  = (aIm[2]-aIm[1])*tIRT3;
	aIm[0] = aIm[0]+t1_im;

	t1_im  = aIm[0]+t1_im*c3_1;
	aRe[1] = t1_re-t2_im;
	aRe[2] = t1_re+t2_im;
	aIm[1] = t1_im+t2_re;
	aIm[2] = t1_im-t2_re;
}


/**
 * length 4 DFT
 */
template<class Type>
inline void FFTPF<Type>::radix4( Type *aRe, Type *aIm )
{
	t1_re  = aRe[0]+aRe[2];
	m2_re  = aRe[0]-aRe[2];
	t2_re  = aRe[1]+aRe[3];
	aRe[0] = t1_re+t2_re;
	aRe[2] = t1_re-t2_re;

	t1_re  = aIm[0]+aIm[2];
	m2_im  = aIm[0]-aIm[2];
	t2_re  = aIm[1]+aIm[3];
	aIm[0] = t1_re+t2_re;
	aIm[2] = t1_re-t2_re;

	t1_re = aRe[1]-aRe[3];
	t2_re = aIm[1]-aIm[3];
	aRe[1] = m2_re+t2_re;
	aRe[3] = m2_re-t2_re;
	aIm[1] = m2_im-t1_re;
	aIm[3] = m2_im+t1_re;
}


/**
 * length 5 DFT
 */
template<class Type>
inline void FFTPF<Type>::radix5( Type *aRe, Type *aIm )
{
	t1_re = aRe[1]+aRe[4];
	t2_re = aRe[2]-aRe[3];
	t3_re = aRe[1]-aRe[4];
	t4_re = aRe[2]+aRe[3];
	t5_re = (t1_re-t4_re)*c5_1;
	t1_re = t1_re+t4_re;

	aRe[0] = aRe[0]+t1_re;

	t1_re = aRe[0]-t1_re*c5_5;
	t4_re = t1_re-t5_re;
	t1_re = t1_re+t5_re;
	t5_re = (t3_re+t2_re)*c5_2;
	t3_re = t5_re-t3_re*c5_4;
	t2_re = t5_re-t2_re*c5_3;
	t1_im = aIm[1]+aIm[4];
    t2_im = aIm[2]-aIm[3];
	t3_im = aIm[1]-aIm[4];
	t4_im = aIm[2]+aIm[3];
	t5_re = (t1_im-t4_im)*c5_1; t1_im = t1_im+t4_im;

	aIm[0] = aIm[0]+t1_im;

	t1_im = aIm[0]-t1_im*c5_5;
	t4_im = t1_im-t5_re;
	t1_im = t1_im+t5_re;
	t5_re = (t3_im+t2_im)*c5_2;
	t3_im = t5_re-t3_im*c5_4;
	t2_im = t5_re-t2_im*c5_3;

	aRe[1] = t1_re+t2_im;
	aIm[1] = t1_im-t2_re;
	aRe[2] = t4_re-t3_im;
	aIm[2] = t4_im+t3_re;
	aRe[3] = t4_re+t3_im;
	aIm[3] = t4_im-t3_re;
	aRe[4] = t1_re-t2_im;
	aIm[4] = t1_im+t2_re;
}


/**
 * length 7 DFT
 */
template<class Type>
inline void FFTPF<Type>::radix7( Type *aRe, Type *aIm )
{
	t1_re = aRe[1]+aRe[6];
	t1_im = aIm[1]+aIm[6];
	t2_re = aRe[2]+aRe[5];
	t2_im = aIm[2]+aIm[5];
	t3_re = aRe[3]+aRe[4];
	t3_im = aIm[3]+aIm[4];
	t4_re = aRe[6]-aRe[1];
	t4_im = aIm[6]-aIm[1];
	t5_re = aRe[5]-aRe[2];
	t5_im = aIm[5]-aIm[2];
	t6_re = aRe[4]-aRe[3];
	t6_im = aIm[4]-aIm[3];
	t7_re = aRe[0]-Type(0.5*t3_re);
	t7_im = aIm[0]-Type(0.5*t3_im);
	t8_re = t1_re-t3_re;
	t8_im = t1_im-t3_im;
	t9_re = t2_re-t3_re;
	t9_im = t2_im-t3_im;

	m1_re = t7_re+c7_1*t8_re+c7_2*t9_re;
	m1_im = t7_im+c7_1*t8_im+c7_2*t9_im;
	m2_re = t7_re+c7_2*t8_re+c7_3*t9_re;
	m2_im = t7_im+c7_2*t8_im+c7_3*t9_im;
	m3_re = t7_re+c7_3*t8_re+c7_1*t9_re;
	m3_im = t7_im+c7_3*t8_im+c7_1*t9_im;
	m4_re = c7_6*t4_re-c7_4*t5_re+c7_5*t6_re;
	m4_im = c7_6*t4_im-c7_4*t5_im+c7_5*t6_im;
	m5_re = c7_5*t4_re-c7_6*t5_re-c7_4*t6_re;
	m5_im = c7_5*t4_im-c7_6*t5_im-c7_4*t6_im;
	m6_re = c7_4*t4_re+c7_5*t5_re+c7_6*t6_re;
	m6_im = c7_4*t4_im+c7_5*t5_im+c7_6*t6_im;

	aRe[0] = aRe[0]+t1_re+t2_re+t3_re;
	aIm[0] = aIm[0]+t1_im+t2_im+t3_im;
	aRe[1] = m1_re-m6_im;
	aIm[1] = m1_im+m6_re;
	aRe[2] = m2_re-m5_im;
	aIm[2] = m2_im+m5_re;
	aRe[3] = m3_re-m4_im;
	aIm[3] = m3_im+m4_re;
	aRe[4] = m3_re+m4_im;
	aIm[4] = m3_im-m4_re;
	aRe[5] = m2_re+m5_im;
	aIm[5] = m2_im-m5_re;
	aRe[6] = m1_re+m6_im;
	aIm[6] = m1_im-m6_re;
}


/**
 * length 8 DFT
 */
template<class Type>
inline void FFTPF<Type>::radix8( Type *aRe, Type *aIm )
{
	t1_re = aRe[0]+aRe[4];
	t2_re = aRe[0]-aRe[4];
	t3_re = aRe[1]-aRe[7];
	t4_re = aRe[1]+aRe[7];
	t7_re = aRe[2]+aRe[6];
	t6_re = aRe[2]-aRe[6];
	t5_re = aRe[3]+aRe[5];
	t8_re = aRe[3]-aRe[5];
	m2_re = t1_re+t7_re;
	m2_im = t1_re-t7_re;
	m1_re = t4_re+t5_re;
	t4_re = (t4_re-t5_re)*tIRT2;

	aRe[0] = m2_re+m1_re;
	aRe[4] = m2_re-m1_re;

	m2_re = t2_re+t4_re;
	m1_re = t2_re-t4_re;
	t1_im = t3_re-t8_re;
	t3_re = (t3_re+t8_re)*tIRT2;
	t2_im = t3_re+t6_re;
	t4_im = t3_re-t6_re;
	t1_re = aIm[0]+aIm[4];
	t2_re = aIm[0]-aIm[4];
	t3_re = aIm[1]-aIm[7];
	t4_re = aIm[1]+aIm[7];
	t5_re = aIm[3]+aIm[5];
	t6_re = aIm[2]-aIm[6];
	t7_re = aIm[2]+aIm[6];
	t8_re = aIm[3]-aIm[5];
	m1_im = t1_re+t7_re;
	t1_re = t1_re-t7_re;
	t7_re = t4_re+t5_re;
	t4_re = (t4_re-t5_re)*tIRT2;

	aIm[0] = m1_im+t7_re;
	aIm[4] = m1_im-t7_re;

	t7_re = t2_re+t4_re;
	t2_re = t2_re-t4_re;
	t4_re = t3_re-t8_re;
	t3_re = (t3_re+t8_re)*tIRT2;
	t5_re = t3_re+t6_re;
	t3_re = t3_re-t6_re;

	aRe[1] = m2_re+t5_re;
	aIm[1] = t7_re-t2_im;
	aRe[2] = m2_im+t4_re;
	aIm[2] = t1_re-t1_im;
	aRe[3] = m1_re+t3_re;
	aIm[3] = t2_re-t4_im;
	aRe[5] = m1_re-t3_re;
	aIm[5] = t2_re+t4_im;
	aRe[6] = m2_im-t4_re;
	aIm[6] = t1_re+t1_im;
	aRe[7] = m2_re-t5_re;
	aIm[7] = t7_re+t2_im;
}


/**
 * length 9 DFT
 */
template<class Type>
inline void FFTPF<Type>::radix9( Type *aRe, Type *aIm )
{
	t1_re = aRe[1]+aRe[8];
	t2_re = aRe[1]-aRe[8];
	t3_re = aRe[2]+aRe[7];
	t4_re = aRe[2]-aRe[7];
	t5_re = aRe[3]+aRe[6];
	ttmp  = aRe[0]+t5_re;
	t7_re = aRe[4]+aRe[5];
	t8_re = aRe[4]-aRe[5];
	t8_im = (aRe[6]-aRe[3])*tIRT3;
	t7_im = aRe[0]+t5_re*c9_5;
	t5_re = t1_re+t3_re+t7_re;

	aRe[0] = ttmp+t5_re;

	t5_im = ttmp+t5_re*c9_5;
	t3_im = (t7_re-t3_re)*c9_2;
	t7_re = (t7_re-t1_re)*c9_3;
	t3_re = (t1_re-t3_re)*c9_4;
	t1_im = t7_im+t3_im+t3_re;
	t3_im = t7_im-t3_im-t7_re;
	t7_im = t7_im+t7_re-t3_re;
	t6_im = (t4_re-t2_re-t8_re)*tIRT3;
	t4_im = (t4_re+t8_re)*c9_6;
	t8_re = (t8_re-t2_re)*c9_7;
	t2_re = (t2_re+t4_re)*c9_8;
	t2_im = t8_im+t4_im+t2_re;
	t4_im = t8_im-t4_im-t8_re;
	t8_im = t8_im+t8_re-t2_re;
	t1_re = aIm[1]+aIm[8];
	t2_re = aIm[1]-aIm[8];
	t3_re = aIm[2]+aIm[7];
	t4_re = aIm[2]-aIm[7];
	t5_re = aIm[3]+aIm[6];
	t6_re = (aIm[6]-aIm[3])*tIRT3;
	t7_re = aIm[4]+aIm[5];
	t8_re = aIm[4]-aIm[5];
	ttmp  = aIm[0]+t5_re;
	t9_im = aIm[0]+t5_re*c9_5;
	t5_re = t1_re+t3_re+t7_re;

	aIm[0] = ttmp+t5_re;

	t5_re = ttmp+t5_re*c9_5;
	ttmp  = (t7_re-t3_re)*c9_2;
	t7_re = (t7_re-t1_re)*c9_3;
	t3_re = (t1_re-t3_re)*c9_4;
	t1_re = t9_im+ttmp+t3_re;
	ttmp  = t9_im-ttmp-t7_re;
	t7_re = t9_im+t7_re-t3_re;
	t9_re = (t4_re-t2_re-t8_re)*tIRT3;
	t3_re = (t4_re+t8_re)*c9_6;
	t8_re = (t8_re-t2_re)*c9_7;
	t4_re = (t2_re+t4_re)*c9_8;
	t2_re = t6_re+t3_re+t4_re;
	t3_re = t6_re-t8_re-t3_re;
	t8_re = t6_re+t8_re-t4_re;

	aRe[1] = t1_im-t2_re;
	aIm[1] = t1_re+t2_im;
	aRe[2] = t3_im+t3_re;
	aIm[2] = ttmp-t4_im;
	aRe[3] = t5_im-t9_re;
	aIm[3] = t5_re+t6_im;
	aRe[4] = t7_im-t8_re;
	aIm[4] = t7_re+t8_im;
	aRe[5] = t7_im+t8_re;
	aIm[5] = t7_re-t8_im;
	aRe[6] = t5_im+t9_re;
	aIm[6] = t5_re-t6_im;
	aRe[7] = t3_im-t3_re;
	aIm[7] = ttmp+t4_im;
	aRe[8] = t1_im+t2_re;
	aIm[8] = t1_re-t2_im;
}


/**
 * length 10 DFT
 */
template<class Type>
inline void FFTPF<Type>::radix10( Type *aRe, Type *aIm )
{
	t1_re = aRe[2]+aRe[8];
	t2_re = aRe[4]-aRe[6];
	t3_re = aRe[2]-aRe[8];
	t4_re = aRe[4]+aRe[6];
	t5_re = (t1_re-t4_re)*c5_1;
	t1_re = t1_re+t4_re;
	m9_re = aRe[0]+t1_re;
	t1_re = m9_re-t1_re*c5_5;
	t4_re = t1_re-t5_re;
	t1_re = t1_re+t5_re;
	t5_re = (t3_re+t2_re)*c5_2;
	t3_re = t5_re-t3_re*c5_4;
	t2_re = t5_re-t2_re*c5_3;
	t1_im = aIm[2]+aIm[8];
	t2_im = aIm[4]-aIm[6];
	t3_im = aIm[2]-aIm[8];
	t4_im = aIm[4]+aIm[6];
	t5_re = (t1_im-t4_im)*c5_1;
	t1_im = t1_im+t4_im;
	m9_im = aIm[0]+t1_im;
	t1_im = m9_im-t1_im*c5_5;
	t4_im = t1_im-t5_re;
	t1_im = t1_im+t5_re;
	t5_re = (t3_im+t2_im)*c5_2;
	t3_im = t5_re-t3_im*c5_4;
	t2_im = t5_re-t2_im*c5_3;
	m1_re = t1_re+t2_im;
	m1_im = t1_im-t2_re;
	m2_re = t4_re-t3_im;
	m2_im = t4_im+t3_re;
	m3_re = t4_re+t3_im;
	m3_im = t4_im-t3_re;
	m4_re = t1_re-t2_im;
	m4_im = t1_im+t2_re;
	t1_re = aRe[7]+aRe[3];
	t2_re = aRe[9]-aRe[1];
	t3_re = aRe[7]-aRe[3];
	t4_re = aRe[9]+aRe[1];
	t5_re = (t1_re-t4_re)*c5_1;
	t1_re = t1_re+t4_re;
	t9_re = aRe[5]+t1_re;
	t1_re = t9_re-t1_re*c5_5;
	t4_re = t1_re-t5_re;
	t1_re = t1_re+t5_re;
	t5_re = (t3_re+t2_re)*c5_2;
	t3_re = t5_re-t3_re*c5_4;
	t2_re = t5_re-t2_re*c5_3;
	t1_im = aIm[7]+aIm[3];
	t2_im = aIm[9]-aIm[1];
	t3_im = aIm[7]-aIm[3];
	t4_im = aIm[9]+aIm[1];
	t5_re = (t1_im-t4_im)*c5_1;
	t1_im = t1_im+t4_im;
	t9_im = aIm[5]+t1_im;
	t1_im = t9_im-t1_im*c5_5;
	t4_im = t1_im-t5_re;
	t1_im = t1_im+t5_re;
	t5_re = (t3_im+t2_im)*c5_2;
	t3_im = t5_re-t3_im*c5_4;
	t2_im = t5_re-t2_im*c5_3;
	m5_re = t1_re+t2_im;
	m5_im = t1_im-t2_re;
	m6_re = t4_re-t3_im;
	m6_im = t4_im+t3_re;
	m7_re = t4_re+t3_im;
	m7_im = t4_im-t3_re;
	m8_re = t1_re-t2_im;
	m8_im = t1_im+t2_re;

	aRe[0] = m9_re+t9_re;
	aIm[0] = m9_im+t9_im;
	aRe[1] = m1_re-m5_re;
	aIm[1] = m1_im-m5_im;
	aRe[2] = m2_re+m6_re;
	aIm[2] = m2_im+m6_im;
	aRe[3] = m3_re-m7_re;
	aIm[3] = m3_im-m7_im;
	aRe[4] = m4_re+m8_re;
	aIm[4] = m4_im+m8_im;
	aRe[5] = m9_re-t9_re;
	aIm[5] = m9_im-t9_im;
	aRe[6] = m1_re+m5_re;
	aIm[6] = m1_im+m5_im;
	aRe[7] = m2_re-m6_re;
	aIm[7] = m2_im-m6_im;
	aRe[8] = m3_re+m7_re;
	aIm[8] = m3_im+m7_im;
	aRe[9] = m4_re-m8_re;
	aIm[9] = m4_im-m8_im;
}


/**
 * length 11 DFT
 */
template<class Type>
inline void FFTPF<Type>::radix11( Type *aRe, Type *aIm )
{
	t1_re = aRe[1]+aRe[10];
	t1_im = aIm[1]+aIm[10];
	t2_re = aRe[2]+aRe[9];
	t2_im = aIm[2]+aIm[9];
	t3_re = aRe[3]+aRe[8];
	t3_im = aIm[3]+aIm[8];
	t4_re = aRe[4]+aRe[7];
	t4_im = aIm[4]+aIm[7];
	t5_re = aRe[5]+aRe[6];
	t5_im = aIm[5]+aIm[6];
	t6_re = aRe[10]-aRe[1];
	t6_im = aIm[10]-aIm[1];
	t7_re = aRe[9]-aRe[2];
	t7_im = aIm[9]-aIm[2];
	t8_re = aRe[8]-aRe[3];
	t8_im = aIm[8]-aIm[3];
	t9_re = aRe[7]-aRe[4];
	t9_im = aIm[7]-aIm[4];
	t10_re = aRe[6]-aRe[5];
	t10_im = aIm[6]-aIm[5];
	t11_re = aRe[0]-Type(0.5*t5_re);
	t11_im = aIm[0]-Type(0.5*t5_im);
	t12_re = t1_re-t5_re;
	t12_im = t1_im-t5_im;
	t13_re = t2_re-t5_re;
	t13_im = t2_im-t5_im;
	t14_re = t3_re-t5_re;
	t14_im = t3_im-t5_im;
	t15_re = t4_re-t5_re;
	t15_im =  t4_im-t5_im;

	m1_re = t11_re+c11_1*t12_re+c11_2*t13_re+c11_3*t14_re+c11_4*t15_re;
	m1_im = t11_im+c11_1*t12_im+c11_2*t13_im+c11_3*t14_im+c11_4*t15_im;
	m2_re = t11_re+c11_2*t12_re+c11_4*t13_re+c11_5*t14_re+c11_3*t15_re;
	m2_im = t11_im+c11_2*t12_im+c11_4*t13_im+c11_5*t14_im+c11_3*t15_im;
	m3_re = t11_re+c11_3*t12_re+c11_5*t13_re+c11_2*t14_re+c11_1*t15_re;
	m3_im = t11_im+c11_3*t12_im+c11_5*t13_im+c11_2*t14_im+c11_1*t15_im;
	m4_re = t11_re+c11_4*t12_re+c11_3*t13_re+c11_1*t14_re+c11_5*t15_re;
	m4_im = t11_im+c11_4*t12_im+c11_3*t13_im+c11_1*t14_im+c11_5*t15_im;
	m5_re = t11_re+c11_5*t12_re+c11_1*t13_re+c11_4*t14_re+c11_2*t15_re;
	m5_im = t11_im+c11_5*t12_im+c11_1*t13_im+c11_4*t14_im+c11_2*t15_im;
	m6_re = c11_10*t6_re-c11_6*t7_re+c11_9*t8_re-c11_7*t9_re+c11_8*t10_re;
	m6_im = c11_10*t6_im-c11_6*t7_im+c11_9*t8_im-c11_7*t9_im+c11_8*t10_im;
	m7_re = c11_9*t6_re-c11_8*t7_re+c11_6*t8_re+c11_10*t9_re-c11_7*t10_re;
	m7_im = c11_9*t6_im-c11_8*t7_im+c11_6*t8_im+c11_10*t9_im-c11_7*t10_im;
	m8_re = c11_8*t6_re-c11_10*t7_re-c11_7*t8_re+c11_6*t9_re+c11_9*t10_re;
	m8_im = c11_8*t6_im-c11_10*t7_im-c11_7*t8_im+c11_6*t9_im+c11_9*t10_im;
	m9_re = c11_7*t6_re+c11_9*t7_re-c11_10*t8_re-c11_8*t9_re-c11_6*t10_re;
	m9_im = c11_7*t6_im+c11_9*t7_im-c11_10*t8_im-c11_8*t9_im-c11_6*t10_im;
	m10_re = c11_6*t6_re+c11_7*t7_re+c11_8*t8_re+c11_9*t9_re+c11_10*t10_re;
	m10_im = c11_6*t6_im+c11_7*t7_im+c11_8*t8_im+c11_9*t9_im+c11_10*t10_im;

	aRe[0] = aRe[0]+t1_re+t2_re+t3_re+t4_re+t5_re;
	aIm[0] = aIm[0]+t1_im+t2_im+t3_im+t4_im+t5_im;
	aRe[1] = m1_re-m10_im;
	aIm[1] = m1_im+m10_re;
	aRe[2] = m2_re-m9_im;
	aIm[2] = m2_im+m9_re;
	aRe[3] = m3_re-m8_im;
	aIm[3] = m3_im+m8_re;
	aRe[4] = m4_re-m7_im;
	aIm[4] = m4_im+m7_re;
	aRe[5] = m5_re-m6_im;
	aIm[5] = m5_im+m6_re;
	aRe[6] = m5_re+m6_im;
	aIm[6] = m5_im-m6_re;
	aRe[7] = m4_re+m7_im;
	aIm[7] = m4_im-m7_re;
	aRe[8] = m3_re+m8_im;
	aIm[8] = m3_im-m8_re;
	aRe[9] = m2_re+m9_im;
	aIm[9] = m2_im-m9_re;
	aRe[10] = m1_re+m10_im;
	aIm[10] = m1_im-m10_re;
}


/**
 * length 13 DFT
 */
template<class Type>
inline void FFTPF<Type>::radix13( Type *aRe, Type *aIm )
{
	t1_re = aRe[1]+aRe[12];
	t1_im = aIm[1]+aIm[12];
	t2_re = aRe[2]+aRe[11];
	t2_im = aIm[2]+aIm[11];
	t3_re = aRe[3]+aRe[10];
	t3_im = aIm[3]+aIm[10];
	t4_re = aRe[4]+aRe[9];
	t4_im = aIm[4]+aIm[9];
	t5_re = aRe[5]+aRe[8];
	t5_im = aIm[5]+aIm[8];
	t6_re = aRe[6]+aRe[7];
	t6_im = aIm[6]+aIm[7];
	t7_re = aRe[12]-aRe[1];
	t7_im = aIm[12]-aIm[1];
	t8_re = aRe[11]-aRe[2];
	t8_im = aIm[11]-aIm[2];
	t9_re = aRe[10]-aRe[3];
	t9_im = aIm[10]-aIm[3];
	t10_re = aRe[9]-aRe[4];
	t10_im = aIm[9]-aIm[4];
	t11_re = aRe[8]-aRe[5];
	t11_im = aIm[8]-aIm[5];
	t12_re = aRe[7]-aRe[6];
	t12_im = aIm[7]-aIm[6];
	t13_re = aRe[0]-Type(0.5*t6_re);
	t13_im = aIm[0]-Type(0.5*t6_im);
	t14_re = t1_re-t6_re;
	t14_im = t1_im-t6_im;
	t15_re = t2_re-t6_re;
	t15_im = t2_im-t6_im;
	t16_re = t3_re-t6_re;
	t16_im = t3_im-t6_im;
	t17_re = t4_re-t6_re;
	t17_im = t4_im-t6_im;
	t18_re = t5_re-t6_re;
	t18_im = t5_im-t6_im;

	m1_re = t13_re + c13_1*t14_re + c13_2*t15_re + c13_3*t16_re +
            c13_4*t17_re + c13_5*t18_re;
	m1_im = t13_im + c13_1*t14_im + c13_2*t15_im + c13_3*t16_im +
            c13_4*t17_im + c13_5*t18_im;
	m2_re = t13_re + c13_2*t14_re + c13_4*t15_re + c13_6*t16_re +
            c13_5*t17_re + c13_3*t18_re;
	m2_im = t13_im + c13_2*t14_im + c13_4*t15_im + c13_6*t16_im +
            c13_5*t17_im + c13_3*t18_im;
	m3_re = t13_re + c13_3*t14_re + c13_6*t15_re + c13_4*t16_re +
            c13_1*t17_re + c13_2*t18_re;
	m3_im = t13_im + c13_3*t14_im + c13_6*t15_im + c13_4*t16_im +
            c13_1*t17_im + c13_2*t18_im;
	m4_re = t13_re + c13_4*t14_re + c13_5*t15_re + c13_1*t16_re +
            c13_3*t17_re + c13_6*t18_re;
	m4_im = t13_im + c13_4*t14_im + c13_5*t15_im + c13_1*t16_im +
            c13_3*t17_im + c13_6*t18_im;
	m5_re = t13_re + c13_5*t14_re + c13_3*t15_re + c13_2*t16_re +
            c13_6*t17_re + c13_1*t18_re;
	m5_im = t13_im + c13_5*t14_im + c13_3*t15_im + c13_2*t16_im +
            c13_6*t17_im + c13_1*t18_im;
	m6_re = t13_re + c13_6*t14_re + c13_1*t15_re + c13_5*t16_re +
            c13_2*t17_re + c13_4*t18_re;
	m6_im = t13_im + c13_6*t14_im + c13_1*t15_im + c13_5*t16_im +
            c13_2*t17_im + c13_4*t18_im;
	m7_re = c13_12*t7_re - c13_7*t8_re + c13_11*t9_re - c13_8*t10_re +
            c13_10*t11_re - c13_9*t12_re;
	m7_im = c13_12*t7_im - c13_7*t8_im + c13_11*t9_im - c13_8*t10_im +
            c13_10*t11_im - c13_9*t12_im;
	m8_re = c13_11*t7_re - c13_9*t8_re + c13_8*t9_re - c13_12*t10_re -
            c13_7*t11_re + c13_10*t12_re;
	m8_im = c13_11*t7_im - c13_9*t8_im + c13_8*t9_im - c13_12*t10_im -
            c13_7*t11_im + c13_10*t12_im;
	m9_re = c13_10*t7_re - c13_11*t8_re - c13_7*t9_re + c13_9*t10_re -
            c13_12*t11_re - c13_8*t12_re;
	m9_im = c13_10*t7_im - c13_11*t8_im - c13_7*t9_im + c13_9*t10_im -
            c13_12*t11_im - c13_8*t12_im;
	m10_re = c13_9*t7_re + c13_12*t8_re - c13_10*t9_re - c13_7*t10_re +
	         c13_8*t11_re + c13_11*t12_re;
	m10_im = c13_9*t7_im + c13_12*t8_im - c13_10*t9_im - c13_7*t10_im +
	         c13_8*t11_im + c13_11*t12_im;
	m11_re = c13_8*t7_re + c13_10*t8_re + c13_12*t9_re - c13_11*t10_re -
	         c13_9*t11_re - c13_7*t12_re;
	m11_im = c13_8*t7_im + c13_10*t8_im + c13_12*t9_im - c13_11*t10_im -
             c13_9*t11_im - c13_7*t12_im;
	m12_re = c13_7*t7_re + c13_8*t8_re + c13_9*t9_re + c13_10*t10_re +
             c13_11*t11_re + c13_12*t12_re;
	m12_im = c13_7*t7_im + c13_8*t8_im + c13_9*t9_im + c13_10*t10_im +
             c13_11*t11_im + c13_12*t12_im;

	aRe[0] = aRe[0]+t1_re+t2_re+t3_re+t4_re+t5_re+t6_re;
	aIm[0] = aIm[0]+t1_im+t2_im+t3_im+t4_im+t5_im+t6_im;
	aRe[1] = m1_re-m12_im;
	aIm[1] = m1_im+m12_re;
	aRe[2] = m2_re-m11_im;
	aIm[2] = m2_im+m11_re;
	aRe[3] = m3_re-m10_im;
	aIm[3] = m3_im+m10_re;
	aRe[4] = m4_re-m9_im;
	aIm[4] = m4_im+m9_re;
	aRe[5] = m5_re-m8_im;
	aIm[5] = m5_im+m8_re;
	aRe[6] = m6_re-m7_im;
	aIm[6] = m6_im+m7_re;
	aRe[7] = m6_re+m7_im;
	aIm[7] = m6_im-m7_re;
	aRe[8] = m5_re+m8_im;
	aIm[8] = m5_im-m8_re;
	aRe[9] = m4_re+m9_im;
	aIm[9] = m4_im-m9_re;
	aRe[10] = m3_re+m10_im;
	aIm[10] = m3_im-m10_re;
	aRe[11] = m2_re+m11_im;
	aIm[11] = m2_im-m11_re;
	aRe[12] = m1_re+m12_im;
	aIm[12] = m1_im-m12_re;
}


/**
 * length 16 DFT
 */
template<class Type>
inline void FFTPF<Type>::radix16( Type *aRe, Type *aIm )
{
	t1_re = aRe[0]+aRe[8];
	t2_re = aRe[0]-aRe[8];
	t3_re = aRe[1]+aRe[9];
	t4_re = aRe[1]-aRe[9];
	t5_re = aRe[2]+aRe[10];
	t6_re = aRe[2]-aRe[10];
	t7_re = aRe[3]+aRe[11];
	t8_re = aRe[3]-aRe[11];
	t9_re = aRe[4]+aRe[12];
	t10_re = aRe[4]-aRe[12];
	t11_re = aRe[5]+aRe[13];
	t12_re = aRe[5]-aRe[13];
	t13_re = aRe[6]+aRe[14];
	t14_re = aRe[6]-aRe[14];
	t15_re = aRe[7]+aRe[15];
	t16_re = aRe[7]-aRe[15];
	t1_im = t1_re+t9_re;
	t2_im = t1_re-t9_re;
	t3_im = t3_re+t11_re;
	t4_im = t3_re-t11_re;
	t5_im = t5_re+t13_re;
	t6_im = t5_re-t13_re;
	t7_im = t7_re+t15_re;
	t8_im = t7_re-t15_re;
	t1_re = t1_im+t5_im;
	t3_re = t1_im-t5_im;
	t5_re = t3_im+t7_im;
	t7_re = t3_im-t7_im;

	aRe[0] = t1_re+t5_re;
	aRe[8] = t1_re-t5_re;

	t1_im = tIRT2*(t4_im+t8_im);
	t5_im = tIRT2*(t4_im-t8_im);
	t9_re = t2_im+t5_im;
	t11_re = t2_im-t5_im;
	t13_re = t6_im+t1_im;
	t15_re = t6_im-t1_im;
	t1_im = t4_re+t16_re;
	t2_im = t4_re-t16_re;
	t3_im = tIRT2*(t6_re+t14_re);
	t4_im = tIRT2*(t6_re-t14_re);
	t5_im = t8_re+t12_re;
	t6_im = t8_re-t12_re;
	t7_im = c16_2*(t2_im-t6_im);
	t2_im = c16_3*t2_im-t7_im;
	t6_im = c16_4*t6_im-t7_im;
	t7_im = t2_re+t4_im;
	t8_im = t2_re-t4_im;
	t2_re = t7_im+t2_im;
	t4_re = t7_im-t2_im;
	t6_re = t8_im+t6_im;
	t8_re = t8_im-t6_im;
	t7_im = c16_5*(t1_im+t5_im);
	t2_im = t7_im-c16_4*t1_im;
	t4_im = t7_im-c16_3*t5_im;
	t6_im = t10_re+t3_im;
	t8_im = t10_re-t3_im;
	t10_re = t6_im+t2_im;
	t12_re = t6_im-t2_im;
	t14_re = t8_im+t4_im;
	t16_re = t8_im-t4_im;
	t1_re = aIm[0]+aIm[8];
	t9_im = aIm[0]-aIm[8];
	t10_im = aIm[1]+aIm[9];
	t11_im = aIm[1]-aIm[9];
	t5_re = aIm[2]+aIm[10];
	t12_im = aIm[2]-aIm[10];
	t13_im = aIm[3]+aIm[11];
	t14_im = aIm[3]-aIm[11];
	t15_im = aIm[4]+aIm[12];
	t16_im = aIm[4]-aIm[12];
	t17_im = aIm[5]+aIm[13];
	t18_im = aIm[5]-aIm[13];
	t19_im = aIm[6]+aIm[14];
	t20_im = aIm[6]-aIm[14];
	t21_im = aIm[7]+aIm[15];
	t22_im = aIm[7]-aIm[15];
	t1_im = t1_re+t15_im;
	t2_im = t1_re-t15_im;
	t3_im = t10_im+t17_im;
	t4_im = t10_im-t17_im;
	t5_im = t5_re+t19_im;
	t6_im = t5_re-t19_im;
	t7_im = t13_im+t21_im;
	t8_im = t13_im-t21_im;
	t1_re = t1_im+t5_im;
	t10_im = t1_im-t5_im;
	t5_re = t3_im+t7_im;
	t13_im = t3_im-t7_im;

	aIm[0] = t1_re+t5_re;
	aIm[8] = t1_re-t5_re;
	aRe[4] = t3_re+t13_im;
	aIm[4] = t10_im-t7_re;
	aRe[12] = t3_re-t13_im;
	aIm[12] = t10_im+t7_re;

	t1_im = tIRT2*(t4_im+t8_im);
	t5_im = tIRT2*(t4_im-t8_im);
	t15_im = t2_im+t5_im;
	t17_im = t2_im-t5_im;
	t19_im = t6_im+t1_im;
	t21_im = t6_im-t1_im;
	t1_im = t11_im+t22_im;
	t2_im = t11_im-t22_im;
	t3_im = tIRT2*(t12_im+t20_im);
	t4_im = tIRT2*(t12_im-t20_im);
	t5_im = t14_im+t18_im;
	t6_im = t14_im-t18_im;
	t7_im = c16_2*(t2_im-t6_im);
	t2_im = c16_3*t2_im-t7_im;
	t6_im = c16_4*t6_im-t7_im;
	t7_im = t9_im+t4_im;
	t8_im = t9_im-t4_im;
	t9_im = t7_im+t2_im;
	t11_im = t7_im-t2_im;
	t12_im = t8_im+t6_im;
	t14_im = t8_im-t6_im;
	t7_im = c16_5*(t1_im+t5_im);
	t2_im = t7_im-c16_4*t1_im;
	t4_im = t7_im-c16_3*t5_im;
	t6_im = t16_im+t3_im;
	t8_im = t16_im-t3_im;
	t16_im = t6_im+t2_im;
	t18_im = t6_im-t2_im;
	t20_im = t8_im+t4_im;
	t22_im = t8_im-t4_im;

	aRe[1] = t2_re+t16_im;
	aIm[1] = t9_im-t10_re;
	aRe[2] = t9_re+t19_im;
	aIm[2] = t15_im-t13_re;
	aRe[3] = t8_re-t22_im;
	aIm[3] = t14_im+t16_re;
	aRe[5] = t6_re+t20_im;
	aIm[5] = t12_im-t14_re;
	aRe[6] = t11_re-t21_im;
	aIm[6] = t17_im+t15_re;
	aRe[7] = t4_re-t18_im;
	aIm[7] = t11_im+t12_re;
	aRe[9] = t4_re+t18_im;
	aIm[9] = t11_im-t12_re;
	aRe[10] = t11_re+t21_im;
	aIm[10] = t17_im-t15_re;
	aRe[11] = t6_re-t20_im;
	aIm[11] = t12_im+t14_re;
	aRe[13] = t8_re+t22_im;
	aIm[13] = t14_im-t16_re;
	aRe[14] = t9_re-t19_im;
	aIm[14] = t15_im+t13_re;
	aRe[15] = t2_re-t16_im;
	aIm[15] = t9_im+t10_re;
}


/**
 * other length DFT
 */
template<class Type>
inline void FFTPF<Type>::radixOther( int radix )
{
	Type    rere, reim, imre, imim;
	int     i, j, k, n, max;

	n = radix;
	max = (n+1)/2;

	for( j=1; j<max; ++j )
	{
		pfvRe[j] = pfzRe[j]+pfzRe[n-j];
		pfvIm[j] = pfzIm[j]-pfzIm[n-j];
		pfwRe[j] = pfzRe[j]-pfzRe[n-j];
		pfwIm[j] = pfzIm[j]+pfzIm[n-j];
	}

	for( j=1; j<max; ++j )
	{
		pfzRe[j] = pfzRe[0];
		pfzIm[j] = pfzIm[0];
		pfzRe[n-j] = pfzRe[0];
		pfzIm[n-j] = pfzIm[0];
		k = j;
		for( i=1; i<max; ++i )
		{
			rere = pftgRe[k]*pfvRe[i];
			imim = pftgIm[k]*pfvIm[i];
			reim = pftgRe[k]*pfwIm[i];
			imre = pftgIm[k]*pfwRe[i];

			pfzRe[n-j] += rere+imim;
			pfzIm[n-j] += reim-imre;
			pfzRe[j] += rere-imim;
			pfzIm[j] += reim+imre;

			k = k+j;
			if( k >= n )
                k = k-n;
		}
	}

	for( j=1; j<max; ++j )
	{
		pfzRe[0] = pfzRe[0]+pfvRe[j];
		pfzIm[0] = pfzIm[0]+pfwIm[j];
	}
}

/**
 * twiddle multiplications and DFT's for one stage.
 */
template<class Type>
void FFTPF<Type>::twiddleFFT( int sofarRadix, int radix, int remainRadix,
                              Vector< complex<Type> > &yn )
{
	Type   cosw, sinw;
	initTrig( radix );

	omega = 2*tPI/(Type)(sofarRadix*radix);
	cosw  = Type(cos(omega));
	sinw  = -Type(sin(omega));
	twRe  = 1;
	twIim = 0;
	dataOffset = 0;
    groupOffset = dataOffset;
    adr = groupOffset;

	for( dataNo=0; dataNo<sofarRadix; ++dataNo )
	{
		if( sofarRadix > 1 )
		{
			pftwRe[0] = 1;
			pftwIm[0] = 0;
			pftwRe[1] = twRe;
			pftwIm[1] = twIim;

			for( twNo=2; twNo<radix; ++twNo )
			{
				pftwRe[twNo] = twRe*pftwRe[twNo-1]-twIim*pftwIm[twNo-1];
				pftwIm[twNo] = twIim*pftwRe[twNo-1]+twRe*pftwIm[twNo-1];
			}
			ttmp = cosw*twRe-sinw*twIim;
			twIim = sinw*twRe+cosw*twIim;
			twRe = ttmp;
		}

		for( groupNo=0; groupNo<remainRadix; ++groupNo )
		{
			if( (sofarRadix>1) && (dataNo>0) )
			{
				pfzRe[0] = yn[adr].real();
				pfzIm[0] = yn[adr].imag();
				blockNo = 1;
				do
				{
					adr += sofarRadix;
					pfzRe[blockNo] = pftwRe[blockNo]*yn[adr].real() -
                                     pftwIm[blockNo]*yn[adr].imag();
					pfzIm[blockNo] = pftwRe[blockNo]*yn[adr].imag() +
                                     pftwIm[blockNo]*yn[adr].real();
					blockNo++;
				}
				while( blockNo < radix );
			}
			else
                for( blockNo=0; blockNo<radix; ++blockNo )
                {
                    pfzRe[blockNo] = yn[adr].real();
                    pfzIm[blockNo] = yn[adr].imag();
                    adr = adr+sofarRadix;
                }

			switch( radix )
			{
                case 2 :
                    radix2( pfzRe, pfzIm );
                    break;
                case 3 :
                    radix3( pfzRe, pfzIm );
                    break;
                case 4 :
                    radix4( pfzRe, pfzIm );
                    break;
                case 5 :
                    radix5( pfzRe, pfzIm );
                    break;
                case 7 :
                    radix7( pfzRe, pfzIm );
                    break;
                case 8 :
                    radix8( pfzRe, pfzIm );
                    break;
                case 9 :
                    radix9( pfzRe, pfzIm );
                    break;
                case 10 :
                    radix10( pfzRe, pfzIm );
                    break;
                case 11 :
                    radix11( pfzRe, pfzIm );
                    break;
                case 13 :
                    radix13( pfzRe, pfzIm );
                    break;
                case 16 :
                    radix16( pfzRe, pfzIm );
                    break;
                default :
                    radixOther( radix );
                    break;
			}

			adr = groupOffset;
			for( blockNo=0; blockNo<radix; ++blockNo )
			{
                yn[adr] = complex<Type>( pfzRe[blockNo], pfzIm[blockNo] );
				adr += sofarRadix;
			}
			groupOffset += sofarRadix*radix;
			adr = groupOffset;
		}

		dataOffset = dataOffset+1;
		groupOffset = dataOffset;
		adr = groupOffset;
	}
}


/**
 * forward transform
 */
template<class Type>
void FFTPF<Type>::fft( const Vector<Type> &xn, Vector< complex<Type> > &Xk )
{
    int n = xn.size();
	if( n != mNOldSize )
        primeSetup( n );

	permute( xn, Xk );

	for( int i=1; i<=mNFactor; ++i )
		twiddleFFT( mSofarRadix[i], mActualRadix[i], mRemainRadix[i], Xk );
}

template<class Type>
void FFTPF<Type>::fft( const Vector< complex<Type> > &xn,
                       Vector< complex<Type> > &Xk )
{
    int n = xn.size();
	if( n != mNOldSize )
        primeSetup( n );

	permute( xn, Xk, true );

	for( int i=1; i<=mNFactor; ++i )
		twiddleFFT( mSofarRadix[i], mActualRadix[i], mRemainRadix[i], Xk );
}


/**
 * inverse transform
 */
template<class Type>
void FFTPF<Type>::ifft( const Vector< complex<Type> > &Xk, Vector<Type> &xn )
{
    int n = xn.size();
	if( n != mNOldSize )
        primeSetup( n );

    Vector< complex<Type> > Yk(n);
	permute( Xk, Yk, false );

	for( int i=1; i<=mNFactor; ++i )
		twiddleFFT( mSofarRadix[i], mActualRadix[i], mRemainRadix[i], Yk );

    for( int i=0; i<n; ++i )
        xn[i] = Yk[i].real() / Type(n);
}

template<class Type>
void FFTPF<Type>::ifft( const Vector< complex<Type> > &Xk,
                        Vector< complex<Type> > &xn )
{
    int n = xn.size();
	if( n != mNOldSize )
        primeSetup( n );

	permute( Xk, xn, false );

	for( int i=1; i<=mNFactor; ++i )
		twiddleFFT( mSofarRadix[i], mActualRadix[i], mRemainRadix[i], xn );

    for( int i=0; i<n; ++i )
        xn[i] = conj(xn[i]) / Type(n);
}
