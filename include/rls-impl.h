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
 *                              rls-impl.h
 *
 * Implementation for RLS Filter.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * The conventional RLS algorighm. The parameter "lambda" is the Forgetting
 * Factor, the smaller "lambda" is, the smaller contribution of previous samples.
 * This makes the filter more sensitive to recent samples, which means more
 * fluctuations in the filter co-efficients. Suggesting range is: [0.8, 1.0].
 * The parametet "delta" is the value to initialeze the inverse of the Auto-
 * Relation Matrix of input signal, which can be chosen as an estimation of
 * the input signal power.
 */
template <typename Type>
Type rls( const Type &xk, const Type &dk, Vector<Type> &wn,
          const Type &lambda, const Type &delta )
{
    assert( Type(0.8) <= lambda );
    assert( lambda <= Type(1.0) );

    int filterLen = wn.size();
    Vector<Type> vP(filterLen);
    Vector<Type> vQ(filterLen);

    static Vector<Type> xn(filterLen);
    static Matrix<Type> invR = eye( filterLen, Type(1.0/delta) );

    // updata input signal
    for( int i=filterLen; i>1; --i )
        xn(i) = xn(i-1);
    xn(1) = xk;

    // priori error
    Type ak = dk - dotProd(wn,xn);

    vQ = invR * xn;
    vP = vQ / (lambda+dotProd(vQ,xn));

    // updata Correlation-Matrix's inverse
    invR = (invR - multTr(vQ,vP)) / lambda;

    // update Weight Vector
    wn += ak * vP;
    //    wn += ak * (invR*xn);

    return dotProd(wn,xn);
}


/**
 * Stabilized Fast Tranversal RLS
 */
template <typename Type>
Type sftrls( const Type &xk, const Type &dk, Vector<Type> &wn,
             const Type &lambda, const Type &epsilon,
             const string &training )
{
    int filterLen = wn.size(),
        L = wn.size()-1;

    assert( Type(1.0-1.0/(2*L+2)) <= lambda );
    assert( lambda <= Type(1.0) );

    static Vector<Type> xn(filterLen), xnPrev(filterLen);

    const Type  k1 = Type(1.5),
                k2 = Type(2.5),
                k3 = Type(1.0);

    // initializing for begin
    Type    e, ep,
            ef, efp,
            eb1, eb2, ebp1, ebp2, ebp31, ebp32, ebp33;

    static Type gamma = 1,
                xiBmin = epsilon,
                xiFminInv = 1/epsilon;

    Vector<Type> phiExt(L+2);

    static Vector<Type> phi(filterLen),
                        wf(filterLen), wb(filterLen);

    // updata input signal
    xnPrev = xn;
    for( int i=1; i<=L; ++i )
        xn[i] = xnPrev[i-1];
    xn[0] = xk;

    if( training == "on" )
    {
        // forward prediction error
        efp = xk - dotProd(wf,xnPrev);
        ef  = gamma * efp;

        phiExt[0] = efp * xiFminInv/lambda;
        for( int i=0; i<filterLen; ++i )
            phiExt[i+1] = phi[i] - phiExt[0]*wf[i];

        // gamma1
        gamma = 1 / ( 1/gamma + phiExt[0]*efp );

        // forward minimum weighted least-squares error
        xiFminInv = xiFminInv/lambda - gamma*phiExt[0]*phiExt[0];

        // forward prediction coefficient vector
        wf += ef * phi;

        // backward prediction errors
        ebp1 = lambda * xiBmin * phiExt[filterLen];
        ebp2 = xnPrev[L] - dotProd(wb,xn);
        ebp31 = (1-k1)*ebp1 + k1*ebp2;
        ebp32 = (1-k2)*ebp1 + k2*ebp2;
        ebp33 = (1-k3)*ebp1 + k3*ebp2;

        // gamma2
        gamma = 1 / ( 1/gamma - phiExt[filterLen]*ebp33 );

        // backward prediction errors
        eb1 = gamma * ebp31;
        eb2 = gamma * ebp32;

        // backward minimum weighted least-squares error
        xiBmin = lambda*xiBmin + eb2*ebp2;

        for( int i=0; i<filterLen; ++i )
            phi[i] = phiExt[i] + phiExt[filterLen]*wb[i];

        // forward prediction coefficient vector
        wb += eb1 * phi;

        // gamma3
        gamma = 1 / ( 1 + dotProd(phi,xn) );

        // Joint-Process Estimation
        ep = dk - dotProd(wn,xn);
        e  = gamma * ep;
        wn += e * phi;
    }

    return dotProd(wn,xn);
}


/**
 * Lattice RLS
 */
template <typename Type>
Type lrls( const Type &xk, const Type &dk, Vector<Type> &vn,
           const Type &lambda, const Type &epsilon,
           const string &training )
{
    assert( Type(0.8) <= lambda );
    assert( lambda <= Type(1.0) );

    int filterLen = vn.size(),
        L = filterLen-1;

    // initializing for begin
    Vector<Type>    gamma(filterLen),
                    eb(filterLen),
                    kb(L), kf(L),
                    xiBmin(filterLen), xiFmin(filterLen);

    static Vector<Type> delta(L), deltaD(filterLen),
                        gammaOld(filterLen,Type(1.0)),
                        ebOld(filterLen),
                        xiBminOld(filterLen,epsilon),
                        xiFminOld(filterLen,epsilon);

    // initializing for Zeor Order
    gamma[0] = 1;
    xiBmin[0]= xk*xk + lambda*xiFminOld[0];
    xiFmin[0] = xiBmin[0];

    Type e = dk;
    Type ef = xk;
    eb[0] = xk;

    for( int j=0; j<L; ++j )
    {
        // auxiliary parameters
        delta[j] = lambda*delta[j] + ebOld[j]*ef/gammaOld[j];
        gamma[j+1] = gamma[j] - eb[j]*eb[j]/xiBmin[j];

        // reflection coefficients
        kb[j] = delta[j] / xiFmin[j];
        kf[j] = delta[j] / xiBminOld[j];

        // prediction errors
        eb[j+1] = ebOld[j] - kb[j]*ef;
        ef -= kf[j]*ebOld[j];

        // minimum least-squares
        xiBmin[j+1] = xiBminOld[j] - delta[j]*kb[j];
        xiFmin[j+1] = xiFmin[j] - delta[j]*kf[j];

        // feedforward filtering
        if( training == "on" )
        {
            deltaD[j] = lambda*deltaD[j] + e*eb[j]/gamma[j];
            vn[j] = deltaD[j] / xiBmin[j];
        }
        e -= vn[j]*eb[j];
    }

    // last order feedforward filtering
    if( training == "on" )
    {
        deltaD[L] = lambda*deltaD[L] + e*eb[L]/gamma[L];
        vn[L] = deltaD[L] / xiBmin[L];
    }
    e -= vn[L]*eb[L];

    // updated parameters
    gammaOld = gamma;
    ebOld = eb;
    xiFminOld = xiFmin;
    xiBminOld = xiBmin;

    return dk-e;
}


/**
 * Error Feedback Lattice RLS
 */
template <typename Type>
Type eflrls( const Type &xk, const Type &dk, Vector<Type> &vn,
             const Type &lambda, const Type &epsilon,
             const string &training )
{
    assert( Type(0.8) <= lambda );
    assert( lambda <= Type(1.0) );

    int filterLen = vn.size(),
        L = filterLen-1;

    // initializing for begin
    Vector<Type>    gamma(filterLen),
                    eb(filterLen),
                    xiBmin(filterLen), xiFmin(filterLen);

    static Vector<Type> delta(L), deltaD(filterLen),
                        gammaOld(filterLen,Type(1.0)),
                        ebOld(filterLen),
                        kb(L), kf(L),
                        xiBminOld2(filterLen,epsilon),
                        xiBminOld(filterLen,epsilon),
                        xiFminOld(filterLen,epsilon);

    // initializing for Zeor Order
    gamma[0] = 1;
    xiBmin[0]= xk*xk + lambda*xiFminOld[0];
    xiFmin[0] = xiBmin[0];

    Type tmp = 0;
    Type e = dk;
    Type ef = xk;
    eb[0] = xk;

    for( int j=0; j<L; ++j )
    {
        // auxiliary parameters
        delta[j] = lambda*delta[j] + ebOld[j]*ef/gammaOld[j];
        gamma[j+1] = gamma[j] - eb[j]*eb[j]/xiBmin[j];

        // reflection coefficients
        tmp = ebOld[j]*ef / gammaOld[j]/lambda;

        kb[j] = gamma[j+1]/gammaOld[j] * ( kb[j] + tmp/xiFminOld[j] );
	    kf[j] = gammaOld[j+1]/gammaOld[j] * ( kf[j] + tmp/xiBminOld2[j] ) ;

        // prediction errors
        eb[j+1] = ebOld[j] - kb[j]*ef;
        ef -= kf[j]*ebOld[j];

        // minimum least-squares
        xiBmin[j+1] = xiBminOld[j] - delta[j]*delta[j]/xiFmin[j];
        xiFmin[j+1] = xiFmin[j] - delta[j]*delta[j]/xiBminOld[j];

        // feedforward filtering
        if( training == "on" )
           vn[j] = gamma[j+1]/gamma[j] *
                   ( vn[j] + e*eb[j]/(lambda*gamma[j]*xiBminOld[j]) );

        e -= vn[j]*eb[j];
    }

    // last order feedforward filtering
    if( training == "on" )
        vn[L] = (gamma[L]-eb[L]*eb[L]/xiBmin[L])/gamma[L] *
                (vn[L]+e*eb[L]/(lambda*gamma[L]*xiBminOld[L]));
    e -= vn[L]*eb[L];

    // updated parameters
    gammaOld = gamma;
    ebOld = eb;
    xiBminOld2 = xiBminOld;
    xiBminOld = xiBmin;
    xiFminOld = xiFmin;

    return dk-e;
}


/**
 * QR-RLS
 */
template <typename Type>
Type qrrls( const Type &xk, const Type &dk, Vector<Type> &wn,
            const Type &lambdaSqrt, const string &training )
{
    int filterLen = wn.size(),
		fL1 = filterLen+1;

    // initializing for begin
	static int k = 1;

	Type	dp, ep,
            gammap,
			c, cosTheta, sinTheta,
			tmp;

    static Type  sx1;

    Vector<Type> xp(filterLen);

    static Vector<Type> xn(filterLen), dn(filterLen), dq2p(filterLen);

	static Matrix<Type> Up(filterLen,filterLen);

	// updata input signal
	for( int i=filterLen; i>1; --i )
            xn(i) = xn(i-1);
        xn(1) = xk;

    if( training == "on" )
    {
        // initializing for 0 to L iterations
        if( k <= filterLen )
        {
            // updata Up
            for( int i=filterLen; i>1; --i )
                for( int j=1; j<=filterLen; ++j )
                    Up(i,j) = lambdaSqrt * Up(i-1,j);
            for( int j=1; j<=filterLen; ++j )
                    Up(1,j) = lambdaSqrt * xn(j);

            dn(k) = dk;
            for( int i=filterLen; i>1; --i )
                dq2p(i) = lambdaSqrt * dq2p(i-1);
            dq2p(1) = lambdaSqrt * dn(k);

            if( k == 1 )
            {
                sx1 = xk;
                if( abs(sx1) > Type(1.0-6) )
                    sx1 = Type(1.0);
            }

            // new Weight Vector
            wn(1) = dn(1) / sx1;
            if( k > 1 )
                for( int i=2; i<=k; ++i )
                {
                    tmp = 0;
                    for( int j=2; j<=i; ++j )
                        tmp += xn(k-j+1)*wn(i-j+1);

                    wn(i) = (-tmp+dn(i)) / sx1;
                }

            ep = dk - dotProd(wn,xn);

            k++;
        }
        else
        {
            xp = xn;
            gammap = 1;
            dp = dk;

            // Givens rotation
            for( int i=1; i<=filterLen; ++i )
            {
                c = sqrt( Up(i,fL1-i)*Up(i,fL1-i) + xp(fL1-i)*xp(fL1-i) );
                cosTheta = Up(i,fL1-i) / c;
                sinTheta = xp(fL1-i) / c;
                gammap *= cosTheta;

                for( int j=1; j<=filterLen; ++j )
                {
                    tmp = xp(j);
                    xp(j) = cosTheta*tmp - sinTheta*Up(i,j);
                    Up(i,j) = sinTheta*tmp + cosTheta*Up(i,j);
                }

                tmp = dp;
                dp = cosTheta*tmp - sinTheta*dq2p(i);
                dq2p(i) = sinTheta*tmp + cosTheta*dq2p(i);
            }

            ep = dp / gammap;

            // new Weight Vector
            wn(1) = dq2p(filterLen) / Up(filterLen,1);
            for( int i=2; i<=filterLen; ++i )
            {
                tmp = 0;
                for( int j=2; j<=i; ++j )
                    tmp += Up(fL1-i,i-j+1) * wn(i-j+1);
                wn(i) = (-tmp+dq2p(fL1-i)) / Up(fL1-i,i);
            }

            // updating internal variables
            Up *= lambdaSqrt;
            dq2p *= lambdaSqrt;
        }

        return dk-ep;
    }
    else
        return dotProd(wn,xn);
}
