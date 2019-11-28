/*************************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Filename:    mtimesx_RealTimesReal.c
 * Programmer:  James Tursa
 * Version:     1.20
 * Date:        February 23, 2010
 * Copyright:   (c) 2009, 2010 by James Tursa, All Rights Reserved
 *
 *  This code uses the BSD License:
 *
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are 
 *  met:
 *
 *     * Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright 
 *       notice, this list of conditions and the following disclaimer in 
 *       the documentation and/or other materials provided with the distribution
 *      
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * mtimesx_RealTimesReal.c is a support file for the mtimesx.c mex routine.
 *
 * Change Log:
 * 2009/Sep/27 --> 1.00, Initial Release
 * 2009/Dec/03 --> 1.01, Fixed scalar * sparse for scalar with inf or NaN
 * 2009/Dec/05 --> 1.02, Fixed bug, added line scalarmultiply = 0;
 * 2009/Dec/08 --> 1.10, Added singleton expansion capability
 * 2009/Dec/10 --> 1.11, Slight code simplification for singleton expansion
 * 2010/Feb/23 --> 1.20, Fixed bug for dgemv and sgemv calls
 *
 ****************************************************************************/

//---------------------------------------------------------------------------------
// Complex type for function return
//---------------------------------------------------------------------------------

struct RealKindComplex {RealKind r; RealKind i;};

//---------------------------------------------------------------------------------
// Function Prototypes
//---------------------------------------------------------------------------------

int AllRealZero(RealKind *x, mwSignedIndex n);
void RealKindEqP1P0TimesRealKindN(RealKind *Cpr, RealKind *Cpi, 
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1P0TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1P0TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1P0TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1P1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1P1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1P1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1P1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1M1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1M1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1M1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1M1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1PxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1PxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1PxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1PxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1P0TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1P0TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1P0TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1P0TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1P1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1P1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1P1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1P1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1M1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1M1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1M1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1M1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1PxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1PxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1PxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1PxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxP1TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxP1TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxP1TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxP1TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxM1TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxM1TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxM1TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxM1TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxP0TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxP0TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxP0TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxP0TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxPxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind ai, RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxPxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxPxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxPxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void xFILLPOS(RealKind *Cpr, mwSignedIndex n);
void xFILLNEG(RealKind *Cpr, mwSignedIndex n);
RealKind xDOT(mwSignedIndex *, RealKind *, mwSignedIndex *, RealKind *, mwSignedIndex *);
void     xGER(mwSignedIndex *, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *,
              RealKind *, mwSignedIndex *, RealKind *, mwSignedIndex *);
void    xGEMV(char *, mwSignedIndex *, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *,
              RealKind *, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *);
void    xGEMM(char *, char *, mwSignedIndex *, mwSignedIndex *, mwSignedIndex *,
              RealKind *, RealKind *, mwSignedIndex *, RealKind *, mwSignedIndex *,
              RealKind *, RealKind *, mwSignedIndex *);
void    xSYRK(char *, char *, mwSignedIndex *, mwSignedIndex *, RealKind *, RealKind *,
              mwSignedIndex *, RealKind *, RealKind *, int*);
void   xSYR2K(char *, char *, int*, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *,
              RealKind *, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *);
struct RealKindComplex RealKindDotProduct(mwSignedIndex, RealKind *, RealKind *, RealKind, 
                                                         RealKind *, RealKind *, RealKind);
void RealKindOuterProduct(mwSignedIndex, mwSignedIndex, RealKind *, RealKind *, char,
                          RealKind *, RealKind *, char, RealKind *, RealKind *);
mxArray *RealScalarTimesReal(mxArray *, char, mwSize, mwSize, mxArray *, char, mwSize, mwSize);

//-------------------------------------------------------------------------------------
// Function for multiplying two MATLAB variables
//-------------------------------------------------------------------------------------

mxArray *RealTimesReal(mxArray *A, char transa, mxArray *B, char transb)
{
    mwSignedIndex inc = 1;
    RealKind Zero = zero;
    RealKind One = one;
    RealKind Minusone = minusone;
    char uplo = 'L';  // Arbitrary choice. Pick the symmetric case as lower triangular

    mwSize m1, n1, m2, n2, Andim, Bndim, Cndim, Ap, Bp, Cp, ip, p, Asize, Bsize, Csize, ndim;
    mwSize Ablock, Bblock;
    mwSize *Adims, *Bdims, *Cdims, *Adimz, *Bdimz, *Cdimz, *Cindx;
    register mwSignedIndex i, j; 
    mwSignedIndex m, n, k, l, lda, ldb, ldc;
    mxArray *C, *result, *rhs[4];
    RealKind *Apr, *Api, *Bpr, *Bpi, *Cpr, *Cpi, *Apr0, *Api0, *Bpr0, *Bpi0;
    RealKind *apr, *api, *bpr, *bpi, *cpr, *cpi;
    RealKind ai, bi, alpha;
    char trans, ptransa, ptransb;
    char transstring[2] = "_";
    struct RealKindComplex z;
    int scalarmultiply;
    
//--------------------------------------------------------------------------------
// Get sizes. Note that in the multi-dimensional case, mxGetN returns the product
// of all of the dimension sizes 2 through end.
//--------------------------------------------------------------------------------
    
    m1 = mxGetM(A);
    n1 = mxGetN(A);
    m2 = mxGetM(B);
    n2 = mxGetN(B);
    
//--------------------------------------------------------------------------------
// Get pointers to the data areas of the operands.
//--------------------------------------------------------------------------------
    
    Apr0 = Apr = mxGetData(A);
    Api0 = Api = mxGetImagData(A);
    Bpr0 = Bpr = mxGetData(B);
    Bpi0 = Bpi = mxGetImagData(B);
    
//--------------------------------------------------------------------------------
// Scalar expansion cases (custom sparse array code for these cases only).
// If there is a inf or NaN in the scalar and the other variable is sparse,
// then don't do the custom code because the sparse zeros will not remain
// zero. So in that case just fall through and call mtimesx_sparse.
// (1 x 1) * (K x N) or (M x K) * (1 x 1)
//--------------------------------------------------------------------------------
    
    scalarmultiply = 0;  // Bug fix added line 2009/12/05
    if( m1 == 1 && n1 == 1 ) {
        scalarmultiply = 1;
        if( mxIsSparse(B) ) {
            if( mxIsInf(*Apr) || mxIsNaN(*Apr) ) {
                scalarmultiply = 0;
            } else if( (Api != NULL) && (mxIsInf(*Api) || mxIsNaN(*Api)) ) {
                scalarmultiply = 0;
            }
        }
    } else if( m2 == 1 && n2 == 1 ) {
        scalarmultiply = 1;
        if( mxIsSparse(A) ) {
            if( mxIsInf(*Bpr) || mxIsNaN(*Bpr) ) {
                scalarmultiply = 0;
            } else if( (Bpi != NULL) && (mxIsInf(*Bpi) || mxIsNaN(*Bpi)) ) {
                scalarmultiply = 0;
            }
        }
    }
    if( scalarmultiply ) {
        return RealScalarTimesReal(A, transa, m1, n1, B, transb, m2, n2);
     }
    
//--------------------------------------------------------------------------------
// Generic sparse matrix or vector operations are not directly supported.
// So just call an m-file to do the work. Won't save any time, but at least
// the function will be supported.
//--------------------------------------------------------------------------------
    
    if( mxIsSparse(A) || mxIsSparse(B) ) {
        rhs[0] = A;
        transstring[0] = transa;
        rhs[1] = mxCreateString(transstring);
        rhs[2] = B;
        transstring[0] = transb;
        rhs[3] = mxCreateString(transstring);
        mexCallMATLAB(1, &result, 4, rhs, "mtimesx_sparse");
        mxDestroyArray(rhs[3]);
        mxDestroyArray(rhs[1]);
        return result;
    }
    
//-------------------------------------------------------------------------------
// Rename array sizes for convenience. Also makes sure that the integer arguments
// to the BLAS routines are of type mwSignedIndex.
//-------------------------------------------------------------------------------
    
    Andim = mxGetNumberOfDimensions(A);
    Adims = mxGetDimensions(A);
    Bndim = mxGetNumberOfDimensions(B);
    Bdims = mxGetDimensions(B);
    m1 = Adims[0];
    n1 = Adims[1];
    m2 = Bdims[0];
    n2 = Bdims[1];
    
    if( transa == 'N' || transa == 'G' ) {
        m = m1;
        k = n1;
    } else {
        m = n1;
        k = m1;
    }
    if( transb == 'N' || transb == 'G' ) {
        l = m2;
        n = n2;
    } else {
        l = n2;
        n = m2;
    }
    lda = m1;
    ldb = m2;
    ldc = m;
    
//-------------------------------------------------------------------------------
// Check for conforming sizes
//-------------------------------------------------------------------------------
    
    if( k != l ) {
        mexErrMsgTxt("Inner matrix dimensions must agree.");
    }
    
    ndim = (Andim <= Bndim) ? Andim : Bndim;
    for( Cp=2; Cp<ndim; Cp++ ) {
        if( Adims[Cp] != Bdims[Cp] && Adims[Cp] != 1 && Bdims[Cp] != 1 ) {
            mexErrMsgTxt("Dimensions 3:end must agree or be 1 in ND case.");
        }
    }
    
//-------------------------------------------------------------------------------
// Construct the dimensions of the result. Also use the p variable to keep track
// of the total number of individual matrix multiples that are involved. The
// first two dimensions are simply the result of a single matrix multiply, with
// accouting for the transa and transb pre-operations. The remaining dimensions
// are copied from A or B, whichever happens to be non-singleton.
//-------------------------------------------------------------------------------
    
    Cndim = (Andim > Bndim) ? Andim : Bndim;
    Cindx = mxMalloc( Cndim * sizeof(*Cindx) );
    Cdims = mxMalloc( Cndim * sizeof(*Cdims) );
    Cdims[0] = m;
    Cdims[1] = n;
    Adimz = mxMalloc( Cndim * sizeof(*Adimz) );
    Adimz[0] = Adims[0];
    Adimz[1] = Adims[1];
    Bdimz = mxMalloc( Cndim * sizeof(*Bdimz) );
    Bdimz[0] = Bdims[0];
    Bdimz[1] = Bdims[1];
    p = 1;
    for( Cp=2; Cp<Cndim; Cp++ ) {
        Adimz[Cp] = (Cp < Andim) ? Adims[Cp] : 1;
        Bdimz[Cp] = (Cp < Bndim) ? Bdims[Cp] : 1;
        Cdims[Cp] = (Adimz[Cp] > Bdimz[Cp]) ? Adimz[Cp] : Bdimz[Cp];
        p *= Cdims[Cp];
    }
    for( Cp=0; Cp<Cndim; Cp++ ) {
        Cindx[Cp] = 0;
    }
    
//------------------------------------------------------------------------------
// Set up conjugate factors
//------------------------------------------------------------------------------
    
    ai = ( transa == 'C' || transa == 'G' ) ? -one : one;
    bi = ( transb == 'C' || transb == 'G' ) ? -one : one;
    
//------------------------------------------------------------------------------
// Create output array
//------------------------------------------------------------------------------
    
    if( mxGetNumberOfElements(A) == 0 || mxGetNumberOfElements(B) == 0 ) {
        result = mxCreateNumericArray(Cndim, Cdims, MatlabReturnType, mxREAL);
        mxFree(Cindx);
        mxFree(Cdims);
        mxFree(Adimz);
        mxFree(Bdimz);
        return result;
    }
    if( mxIsComplex(A) || mxIsComplex(B) ) {
        result = mxCreateNumericArray(Cndim, Cdims, MatlabReturnType, mxCOMPLEX);
    } else {
        result = mxCreateNumericArray(Cndim, Cdims, MatlabReturnType, mxREAL);
    }
    C = result;
    Cpr = mxGetData(C);
    Cpi = mxGetImagData(C);
    
//----------------------------------------------------------------------------
// Outer Loop to process all of the individual matrix multiplies
//----------------------------------------------------------------------------

    Asize = m1 * n1;
    Bsize = m2 * n2;
    Csize = m * n;
    
    for( ip=0; ip<p; ip++ ) {
        ptransa = transa;  // Restore the original transa and transb, because
        ptransb = transb;  // they might have been changed in previous iteration
    
//----------------------------------------------------------------------------
// Vector dot product (1 x K) * (K x 1)
//----------------------------------------------------------------------------
    
    if( m == 1 && n == 1 ) {
        z = RealKindDotProduct(k, Apr, Api, ai, Bpr, Bpi, bi);
        *Cpr = z.r;
        if( mxIsComplex(C) ) {
            *Cpi = z.i;
        }
        
//----------------------------------------------------------------------------
// Vector outer product (M x 1) * (1 x N)
//----------------------------------------------------------------------------
        
    } else if( k == 1 && !matlab ) {
        RealKindOuterProduct(m, n, Apr, Api, transa, Bpr, Bpi, transb, Cpr, Cpi);
        
//----------------------------------------------------------------------------
// Matrix times vector (M x K) * (K x 1)
//----------------------------------------------------------------------------
        
    } else if( n == 1 ) {
        
//----------------------------------------------------------------------------
// If the first matrix is not transposed, use calls to xGEMV. Also use this
// method if running in the 'MATLAB' mode (indicated by matlab variable).
//----------------------------------------------------------------------------
        
        if( transa == 'N' || transa == 'G' || matlab ) {
            if( transa == 'G' ) ptransa = 'N';
            xGEMV(PTRANSA, M1, N1, ONE, Apr, LDA, Bpr, INCX, ZERO, Cpr, INCY);
            if( mxIsComplex(B) ) {
                alpha = bi;
                xGEMV(PTRANSA, M1, N1, ALPHA, Apr, LDA, Bpi, INCX, ZERO, Cpi, INCY);
                if( mxIsComplex(A) ) {                    // (complex matrix) * (complex vector)
                    alpha = -ai * bi;
                    xGEMV(PTRANSA, M1, N1, ALPHA, Api, LDA, Bpi, INCX,  ONE, Cpr, INCY);
                    alpha = ai;
                    xGEMV(PTRANSA, M1, N1, ALPHA, Api, LDA, Bpr, INCX,  ONE, Cpi, INCY);
                } else {                                  // (real matrix) * (complex vector)
                    // already done
                }
            } else {
                if( mxIsComplex(A) ) {                    // (complex matrix) * (real vector)
                    alpha = ai;
                    xGEMV(PTRANSA, M1, N1, ALPHA, Api, LDA, Bpr, INCX, ZERO, Cpi, INCY);
                } else {                                  // (real matrix) * (real vector)
                    // already done
                }
            }

// Alternate method ... doesn't match MATLAB exactly
//
//         if( transa == 'N' || transa == 'G' || matlab ) {
//             if( transa == 'G' ) ptransa = 'N';
//             xGEMV(PTRANSA, M1, N1, ONE, Apr, LDA, Bpr, INCX, ZERO, Cpr, INCY);
//             if( mxIsComplex(A) ) {
//                 alpha = ai;
//                 xGEMV(PTRANSA, M1, N1, ALPHA, Api, LDA, Bpr, INCX, ZERO, Cpi, INCY);
//                 if( mxIsComplex(B) ) {                    // (complex matrix) * (complex vector)
//                     alpha = -ai * bi;
//                     xGEMV(PTRANSA, M1, N1, ALPHA, Api, LDA, Bpi, INCX,  ONE, Cpr, INCY);
//                     alpha = bi;
//                     xGEMV(PTRANSA, M1, N1, ALPHA, Apr, LDA, Bpi, INCX,  ONE, Cpi, INCY);
//                 } else {                                  // (complex matrix) * (real vector)
//                     // already done
//                 }
//             } else {
//                 if( mxIsComplex(B) ) {                    // (real matrix) * (complex vector)
//                     alpha = bi;
//                     xGEMV(PTRANSA, M1, N1, ALPHA, Apr, LDA, Bpi, INCX, ZERO, Cpi, INCY);
//                 } else {                                  // (real matrix) * (real vector)
//                     // already done
//                 }
//             }
            
//-----------------------------------------------------------------------------------------
// Else if the first matrix is transposed, then use calls to xDOT instead (faster) because
// the matrix can be accessed as a series of contiguous column vectors.
//-----------------------------------------------------------------------------------------
            
        } else { // transa == 'T' || transa == 'C'
            apr = Apr;
            api = Api;
            if( mxIsComplex(A) ) {
                for( i=0; i<m; i++ ) {                   // (complex matrix) * (vector)
                    z = RealKindDotProduct(k, apr, api, ai, Bpr, Bpi, bi);
                    Cpr[i] = z.r;
                    Cpi[i] = z.i;
                    apr += k;
                    api += k;
                }
            } else {                                     // (real matrix) * (complex vector)
                if( mxIsComplex(B) ) {
                    for( i=0; i<m; i++ ) {
                        z = RealKindDotProduct(k, apr, Api, ai, Bpr, Bpi, bi);
                        Cpr[i] = z.r;
                        Cpi[i] = z.i;
                        apr += k;
                    }
                } else {                                 // (real matrix) * (real vector)
                    for( i=0; i<m; i++ ) {
                        z = RealKindDotProduct(k, apr, Api, ai, Bpr, Bpi, bi);
                        Cpr[i] = z.r;
                        apr += k;
                    }
                }
            }
        }
        
//----------------------------------------------------------------------------------------
// Vector times matrix (1 x K) * (K x N)
//----------------------------------------------------------------------------------------

    } else if( m == 1 ) {
        
//----------------------------------------------------------------------------------------
// If the second matrix is transposed, then use calls to xGEMV with the arguments reversed.
// Also use this method if running in 'MATLAB' mode (indicated by matlab variable).
//----------------------------------------------------------------------------------------
        
        if( transb == 'C' || transb == 'T' || matlab ) {
            if( transb == 'C' || transb == 'T' ) {
                ptransb = 'N';
            } else {
                ptransb = 'T';
            }
            xGEMV(PTRANSB, M2, N2, ONE, Bpr, LDB, Apr, INCX, ZERO, Cpr, INCY);
            if( mxIsComplex(A) ) {
                alpha = ai;
                xGEMV(PTRANSB, M2, N2, ALPHA, Bpr, LDB, Api, INCX, ZERO, Cpi, INCY);
                if( mxIsComplex(B) ) {                    // (complex matrix) * (complex vector)
                    alpha = -ai * bi;
                    xGEMV(PTRANSB, M2, N2, ALPHA, Bpi, LDB, Api, INCX,  ONE, Cpr, INCY);
                    alpha = bi;
                    xGEMV(PTRANSB, M2, N2, ALPHA, Bpi, LDB, Apr, INCX,  ONE, Cpi, INCY);
                } else {                                  // (complex matrix) * (real vector)
                    // already done
                }
            } else {
                if( mxIsComplex(B) ) {                    // (real matrix) * (complex vector)
                    alpha = bi;
                    xGEMV(PTRANSB, M2, N2, ALPHA, Bpi, LDB, Apr, INCX, ZERO, Cpi, INCY);
                } else {                                  // (real matrix) * (real vector)
                    // already done
                }
            }

// Alternate method ... doesn't match MATLAB exactly
//
//         if( transb == 'C' || transb == 'T' || matlab ) {
//             if( transb == 'C' || transb == 'T' ) {
//                 ptransb = 'N';
//             } else {
//                 ptransb = 'T';
//             }
//             xGEMV(PTRANSB, M2, N2, ONE, Bpr, LDB, Apr, INCX, ZERO, Cpr, INCY);
//             if( mxIsComplex(B) ) {
//                 alpha = bi;
//                 xGEMV(PTRANSB, M2, N2, ALPHA, Bpi, LDB, Apr, INCX, ZERO, Cpi, INCY);
//                 if( mxIsComplex(A) ) {                    // (complex matrix) * (complex vector)
//                     alpha = -ai * bi;
//                     xGEMV(PTRANSB, M2, N2, ALPHA, Bpi, LDB, Api, INCX,  ONE, Cpr, INCY);
//                     alpha = ai;
//                     xGEMV(PTRANSB, M2, N2, ALPHA, Bpr, LDB, Api, INCX,  ONE, Cpi, INCY);
//                 } else {                                  // (real matrix) * (complex vector)
//                     // already done
//                 }
//             } else {
//                 if( mxIsComplex(A) ) {                    // (complex matrix) * (real vector)
//                     alpha = ai;
//                     xGEMV(PTRANSB, M2, N2, ALPHA, Bpr, LDB, Api, INCX, ZERO, Cpi, INCY);
//                 } else {                                  // (real matrix) * (real vector)
//                     // already done
//                 }
//             }
            
//-----------------------------------------------------------------------------------------
// Else if the second matrix is not transposed, then use calls to dot product instead
// (faster) because the matrix can be accessed as a series of contiguous column vectors.
//-----------------------------------------------------------------------------------------
            
        } else {
            bpr = Bpr;
            bpi = Bpi;
            if( mxIsComplex(B) ) {
                for( i=0; i<n; i++ ) {
                    z = RealKindDotProduct(k, Apr, Api, ai, bpr, bpi, bi);
                    Cpr[i] = z.r;
                    Cpi[i] = z.i;
                    bpr += k;
                    bpi += k;
                }
            } else {                                     // (complex vector) * (real matrix)
                if( mxIsComplex(A) ) {
                    for( i=0; i<n; i++ ) {
                        z = RealKindDotProduct(k, Apr, Api, ai, bpr, Bpi, bi);
                        Cpr[i] = z.r;
                        Cpi[i] = z.i;
                        bpr += k;
                    }
                } else {                                 // (real vector) * (real matrix)
                    for( i=0; i<n; i++ ) {
                        z = RealKindDotProduct(k, Apr, Api, ai, bpr, Bpi, bi);
                        Cpr[i] = z.r;
                        bpr += k;
                    }
                }
            }
        }
        
//---------------------------------------------------------------------------------
// Matrix times matrix (M x K) * (K x N) with N small and first matrix transposed.
// Use dot product (faster) because the 1st matrix can be accessed as a series of
// contiguous column vectors. When the column size reaches about 8 then the memory
// access efficiency of the BLAS routines increases and this custom method is no
// longer faster. The number 8 is likely machine / implementation dependent. Only
// use this method if running in the 'SPEED' mode.
//---------------------------------------------------------------------------------
        
    } else if( !matlab && n < 7 && (transa == 'T' || transa == 'C') && (transb == 'N' || transb == 'G') ) {
        bpr = Bpr;
        bpi = Bpi;
        cpr = Cpr;
        cpi = Cpi;
        for( j=0; j<n; j++ ) {
            apr = Apr;
            api = Api;
            for( i=0; i<m; i++ ) {
                z = RealKindDotProduct(k, apr, api, ai, bpr, bpi, bi);
                *cpr++ = z.r;
                if( cpi ) *cpi++ = z.i;
                apr += k;
                if( api ) api += k;
            }
            bpr += k;
            if( bpi ) bpi += k;
        }
        
//---------------------------------------------------------------------------------------------
// Matrix times matrix (M x K) *(K x N)
//---------------------------------------------------------------------------------------------
        
    } else {
        
///--------------------------------------------------------------------------------------------
// If Matrix product is actually the same matrix, use calls to the symmetric routines xSYRK and
// xSYR2K where possible. These only work on the lower or upper triangle part of the matrix, so
// we will have to fill out the other half manually, but even so this will be faster. Some of
// these will not match MATLAB exactly, so only run them in the 'SPEED' mode.
//---------------------------------------------------------------------------------------------
        
        if( A == B && ((transa == 'N' && transb == 'T') || 
                       (transa == 'T' && transb == 'N')) ) {
            xSYRK(UPLO, TRANSA, N, K, ONE, Apr, LDA, ZERO, Cpr, LDC);
            if( mxIsComplex(A) ) {
                xSYRK(UPLO, TRANSA, N, K, MINUSONE, Api, LDA, ONE, Cpr, LDC);
                xSYR2K(UPLO,TRANSA, N, K, ONE, Api, LDA, Apr, LDA, ZERO, Cpi, LDC);
                xFILLPOS(Cpi, n);
            }
            xFILLPOS(Cpr, n);
            
        } else if( A == B && (!matlab || (Api == NULL && Bpi == NULL)) &&
                             ((transa == 'G' && transb == 'C') || 
                              (transa == 'C' && transb == 'G')) ) {
            if( transa == 'G')  ptransa = 'N';
            xSYRK(UPLO, PTRANSA, N, K, ONE, Apr, LDA, ZERO, Cpr, LDC);
            if( mxIsComplex(A) ) {
                xSYRK(UPLO, PTRANSA, N, K, MINUSONE, Api, LDA, ONE, Cpr, LDC);
                xSYR2K(UPLO,PTRANSA, N, K, MINUSONE, Apr, LDA, Api, LDA, ZERO, Cpi, LDC);
                xFILLPOS(Cpi, n);
            }
            xFILLPOS(Cpr, n);
            
        } else if( A == B && ((transa == 'N' && transb == 'C') || 
                              (transa == 'T' && transb == 'G' && (!matlab || (Api == NULL && Bpi == NULL)))) ) {
            if( transb == 'G' ) ptransb = 'N';
            xSYRK(UPLO, TRANSA, N, K, ONE, Apr, LDA, ZERO, Cpr, LDC);
            if( mxIsComplex(A) ) {
                xSYRK(UPLO, TRANSA, N, K, ONE, Api, LDA, ONE, Cpr, LDC);
                xGEMM(TRANSA, PTRANSB, M, N, K, ONE, Api, LDA, Apr, LDA, ZERO, Cpi, LDC);
                xFILLNEG(Cpi, n);
            }
            xFILLPOS(Cpr, n);
            
        } else if( A == B && ((transa == 'C' && transb == 'N') || 
                              (transa == 'G' && transb == 'T' && (!matlab || (Api == NULL && Bpi == NULL)))) ) {
            if( transa == 'G' ) ptransa = 'N';
            xSYRK(UPLO, PTRANSA, N, K, ONE, Apr, LDA, ZERO, Cpr, LDC);
            if( mxIsComplex(A) ) {
                xSYRK(UPLO, PTRANSA, N, K, ONE, Api, LDA, ONE, Cpr, LDC);
                xGEMM(PTRANSA, TRANSB, M, N, K, ONE, Apr, LDA, Api, LDA, ZERO, Cpi, LDC);
                xFILLNEG(Cpi, n);
            }
            xFILLPOS(Cpr, n);
            
//-------------------------------------------------------------------------------------------
// Else this is not a symmetric case, so just call the general matrix multiply routine xGEMM.
//-------------------------------------------------------------------------------------------
            
        } else {
            if( transa == 'G' ) ptransa = 'N';
            if( transb == 'G' ) ptransb = 'N';
            xGEMM(PTRANSA, PTRANSB, M, N, K, ONE, Apr, LDA, Bpr, LDB, ZERO, Cpr, LDC);
            if( mxIsComplex(B) ) {
                alpha = bi;
                xGEMM(PTRANSA, PTRANSB, M, N, K, ALPHA, Apr, LDA, Bpi, LDB, ZERO, Cpi, LDC);
                if( mxIsComplex(A) ) {                    // (complex matrix) * (complex matrix)
                    alpha = -ai * bi;
                    xGEMM(PTRANSA, PTRANSB, M, N, K, ALPHA, Api, LDA, Bpi, LDB,  ONE, Cpr, LDC);
                    alpha = ai;
                    xGEMM(PTRANSA, PTRANSB, M, N, K, ALPHA, Api, LDA, Bpr, LDB,  ONE, Cpi, LDC);
                } else {                                  // (real matrix) * (complex matrix)
                    // already done
                }
            } else {
                if( mxIsComplex(A) ) {                    // (complex matrix) * (real matrix)
                    alpha = ai;
                    xGEMM(PTRANSA, PTRANSB, M, N, K, ALPHA, Api, LDA, Bpr, LDB, ZERO, Cpi, LDC);
                } else {                                  // (real matrix) * (real matrix)
                    // already done
                }
            }
            
// Alternate method ... doesn't match MATLAB exactly
//            
//         } else {
//             if( transa == 'G' ) ptransa = 'N';
//             if( transb == 'G' ) ptransb = 'N';
//             xGEMM(PTRANSA, PTRANSB, M, N, K, ONE, Apr, LDA, Bpr, LDB, ZERO, Cpr, LDC);
//             if( mxIsComplex(A) ) {
//                 alpha = ai;
//                 xGEMM(PTRANSA, PTRANSB, M, N, K, ALPHA, Api, LDA, Bpr, LDB, ZERO, Cpi, LDC);
//                 if( mxIsComplex(B) ) {                    // (complex matrix) * (complex matrix)
//                     alpha = -ai * bi;
//                     xGEMM(PTRANSA, PTRANSB, M, N, K, ALPHA, Api, LDA, Bpi, LDB,  ONE, Cpr, LDC);
//                     alpha = bi;
//                     xGEMM(PTRANSA, PTRANSB, M, N, K, ALPHA, Apr, LDA, Bpi, LDB,  ONE, Cpi, LDC);
//                 } else {                                  // (complex matrix) * (real matrix)
//                     // already done
//                 }
//             } else {
//                 if( mxIsComplex(B) ) {                    // (real matrix) * (complex matrix)
//                     alpha = bi;
//                     xGEMM(PTRANSA, PTRANSB, M, N, K, ALPHA, Apr, LDA, Bpi, LDB, ZERO, Cpi, LDC);
//                 } else {                                  // (real matrix) * (real matrix)
//                     // already done
//                 }
//             }
            
        }
    }
    
//----------------------------------------------------------------------------
// End Outer Loop to process all of the individual matrix multiplies. Increment
// the matrix pointers to point to the next pair of matrices to be multiplied.
//----------------------------------------------------------------------------

    if( ip < p-1 ) {
        j = 2;
        while( ++Cindx[j] == Cdims[j] ) Cindx[j++] = 0;
        Ap = Bp = 0;
        Ablock = Asize;
        Bblock = Bsize;
        for( j=2; j<Cndim; j++ ) {
            if( Cindx[j] < Adimz[j] ) Ap += Cindx[j] * Ablock;
            Ablock *= Adimz[j];
            if( Cindx[j] < Bdimz[j] ) Bp += Cindx[j] * Bblock;
            Bblock *= Bdimz[j];
        }
        Apr = Apr0 + Ap;
        if( Api ) Api = Api0 + Ap;
        Bpr = Bpr0 + Bp;
        if( Bpi ) Bpi = Bpi0 + Bp;
        Cpr += Csize;
        if( Cpi ) Cpi += Csize;
    }
        
    }

    mxFree(Cindx);
    mxFree(Cdims);
    mxFree(Adimz);
    mxFree(Bdimz);
    
//---------------------------------------------------------------------------------
// If the imaginary part is all zero, then free it and set the pointer to NULL.
//---------------------------------------------------------------------------------
    
    Cpi = mxGetPi(C);
    if( AllRealZero(Cpi, m*n*p) ) {
        mxFree(Cpi);
        mxSetImagData(C, NULL);
    }
    
//---------------------------------------------------------------------------------
// Done.
//---------------------------------------------------------------------------------
    
    return result;

}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------
// Returns 1 (true) if all of the elements are zero. Returns 0 (false) if at least one
// of the elements is non-zero, or if the pointer to the data is NULL.
//--------------------------------------------------------------------------------------

int AllRealZero(RealKind *x, mwSignedIndex n)
{
    register mwSignedIndex i;
    if( x == NULL ) return 0;
    for(i=0; i<n; i++) {
        if( x[i] != zero ) return 0;
    }
    return 1;
}

//--------------------------------------------------------------------------------------
// C = (1 + 0*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqP1P0TimesRealKindN(RealKind *Cpr, RealKind *Cpi, 
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i];
            Cpi[i] = Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 + 0*i) * (Bpr - Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqP1P0TimesRealKindG(RealKind *Cpr, RealKind *Cpi, 
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] =  Bpr[i];
            Cpi[i] = -Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 + 0*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqP1P0TimesRealKindT(RealKind *Cpr, RealKind *Cpi, 
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    *Cpi++ = *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 + 0*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------

void RealKindEqP1P0TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    *Cpi++ = -(*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 + 1*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqP1P1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i];
            Cpi[i] = Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i] - Bpi[i];
            Cpi[i] = Bpr[i] + Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 + 1*i) * (Bpr - Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqP1P1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i];
            Cpi[i] = Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i] + Bpi[i];
            Cpi[i] = Bpr[i] - Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 + 1*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqP1P1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    *Cpi++ = *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr - *bpi;
                    *Cpi++ = *bpr + *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 + 1*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------

void RealKindEqP1P1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    *Cpi++ = *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr + *bpi;
                    *Cpi++ = *bpr - *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 - 1*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqP1M1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =  Bpr[i];
            Cpi[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i] + Bpi[i];
            Cpi[i] = Bpi[i] - Bpr[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 - 1*i) * (Bpr - Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqP1M1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =  Bpr[i];
            Cpi[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] =  Bpr[i] - Bpi[i];
            Cpi[i] = -Bpi[i] - Bpr[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 - 1*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqP1M1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =   *bpr;
                    *Cpi++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr + *bpi;
                    *Cpi++ = *bpi - *bpr;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 - 1*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------

void RealKindEqP1M1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =   *bpr;
                    *Cpi++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =   *bpr - *bpi;
                    *Cpi++ = - *bpr - *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 + ai*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------
void RealKindEqP1PxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =      Bpr[i];
            Cpi[i] = ai * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i] - ai * Bpi[i];
            Cpi[i] = ai * Bpr[i] + Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 + ai*i) * (Bpr - Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqP1PxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =      Bpr[i];
            Cpi[i] = ai * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i] + ai * Bpi[i];
            Cpi[i] = ai * Bpr[i] - Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 + ai*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqP1PxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =       *bpr;
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr - ai * (*bpi);
                    *Cpi++ = ai * (*bpr) + *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (1 + ai*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------

void RealKindEqP1PxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =       *bpr;
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr + ai * (*bpi);
                    *Cpi++ = ai * (*bpr) - *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + 1*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqM1P1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] =  Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] - Bpi[i];
            Cpi[i] =  Bpr[i] - Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + 1*i) * (Bpr - Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqM1P1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] =  Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] + Bpi[i];
            Cpi[i] =  Bpr[i] + Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + 1*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqM1P1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ =   *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) - (*bpi);
                    *Cpi++ =   *bpr  - (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + 1*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------

void RealKindEqM1P1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ =   *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) + (*bpi);
                    *Cpi++ =   *bpr  + (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 - 1*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqM1M1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] + Bpi[i];
            Cpi[i] = -Bpi[i] - Bpr[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 - 1*i) * (Bpr - Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqM1M1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] - Bpi[i];
            Cpi[i] =  Bpi[i] - Bpr[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 - 1*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqM1M1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) + (*bpi);
                    *Cpi++ = -(*bpr) - (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 - 1*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------
void RealKindEqM1M1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) - (*bpi);
                    *Cpi++ = -(*bpr) + (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + 0*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqM1P0TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] = -Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + 0*i) * (Bpr - Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqM1P0TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] =  Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + 0*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqM1P0TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ = -(*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + 0*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------

void RealKindEqM1P0TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ =   *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + ai*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqM1PxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =    - Bpr[i];
            Cpi[i] = ai * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] - ai * Bpi[i];
            Cpi[i] =  ai * Bpr[i] - Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + ai*i) * (Bpr - Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqM1PxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =    - Bpr[i];
            Cpi[i] = ai * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] + ai * Bpi[i];
            Cpi[i] =  ai * Bpr[i] + Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + ai*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqM1PxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =    - (*bpr);
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) - ai * (*bpi);
                    *Cpi++ = ai * (*bpr) - *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (-1 + ai*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------

void RealKindEqM1PxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =    - (*bpr);
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) + ai * (*bpi);
                    *Cpi++ = ai * (*bpr) + *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + 1*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqPxP1TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
            Cpi[i] =      Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i] - Bpi[i];
            Cpi[i] = Bpr[i] + ar * Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + 1*i) * (Bpr - Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqPxP1TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
            Cpi[i] =      Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i] + Bpi[i];
            Cpi[i] = Bpr[i] - ar * Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + 1*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqPxP1TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ =       *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr) - (*bpi);
                    *Cpi++ = (*bpr) + ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + 1*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------

void RealKindEqPxP1TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ =       *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr) + (*bpi);
                    *Cpi++ = (*bpr) - ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar - 1*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqPxM1TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
            Cpi[i] =    - Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] =  ar * Bpr[i] + Bpi[i];
            Cpi[i] = -Bpr[i] + ar * Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar - 1*i) * (Bpr - Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqPxM1TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
            Cpi[i] =    - Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] =  ar * Bpr[i] - Bpi[i];
            Cpi[i] = -Bpr[i] - ar * Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar - 1*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqPxM1TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ =    - (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =  ar * (*bpr) + (*bpi);
                    *Cpi++ = -(*bpr) + ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar - 1*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------

void RealKindEqPxM1TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ =    - (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =  ar * (*bpr) - (*bpi);
                    *Cpi++ = -(*bpr) - ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + 0*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqPxP0TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
            Cpi[i] = ar * Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + 0*i) * (Bpr - Bpi * i)
void RealKindEqPxP0TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] =   ar * Bpr[i];
            Cpi[i] =  -ar * Bpi[i];
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + 0*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqPxP0TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =  ar * (*bpr);
                    *Cpi++ =  ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + 0*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------

void RealKindEqPxP0TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =   ar * (*bpr);
                    *Cpi++ =  -ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + ai*i) * (Bpr + Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqPxPxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        if( ai == zero ) {
            for(i=0; i<n; i++) {
                Cpr[i] = ar * Bpr[i];
            }
        } else {
            for(i=0; i<n; i++) {
                Cpr[i] = ar * Bpr[i];
                Cpi[i] = ai * Bpr[i];
            }
        }
    } else {
        if( ai == zero ) {
            for(i=0; i<n; i++) {
                Cpr[i] =  ar * Bpr[i];
                Cpi[i] =  ar * Bpi[i];
            }
        } else {
            for(i=0; i<n; i++) {
                Cpr[i] =  ar * Bpr[i] - ai * Bpi[i];
                Cpi[i] =  ai * Bpr[i] + ar * Bpi[i];
            }
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + ai*i) * (Bpr - Bpi * i)
//--------------------------------------------------------------------------------------

void RealKindEqPxPxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;
    
    if( Bpi == NULL ) {
        if( ai == zero ) {
            for(i=0; i<n; i++) {
                Cpr[i] = ar * Bpr[i];
            }
        } else {
            for(i=0; i<n; i++) {
                Cpr[i] = ar * Bpr[i];
                Cpi[i] = ai * Bpr[i];
            }
        }
    } else {
        if( ai == zero ) {
            for(i=0; i<n; i++) {
                Cpr[i] =   ar * Bpr[i];
                Cpi[i] =  -ar * Bpi[i];
            }
        } else {
            for(i=0; i<n; i++) {
                Cpr[i] =  ar * Bpr[i] + ai * Bpi[i];
                Cpi[i] =  ai * Bpr[i] - ar * Bpi[i];
            }
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + ai*i) * (Bpr + Bpi * i)T
//--------------------------------------------------------------------------------------

void RealKindEqPxPxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =  ar * (*bpr) - ai * (*bpi);
                    *Cpi++ =  ai * (*bpr) + ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// C = (ar + ai*i) * (Bpr + Bpi * i)C
//--------------------------------------------------------------------------------------

void RealKindEqPxPxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =  ar * (*bpr) + ai * (*bpi);
                    *Cpi++ =  ai * (*bpr) - ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

//--------------------------------------------------------------------------------------
// Fill the upper triangle with contents of the lower triangle
//--------------------------------------------------------------------------------------

void xFILLPOS(RealKind *Cpr, mwSignedIndex n)
{
    RealKind *source, *target;
    register mwSignedIndex i, j;
    
    source = Cpr + 1;
    target = Cpr + n;
    for( i=1; i<n; i++ ) {
        for( j=i; j<n; j++ ) {
            *target = *source;
            target += n;
            source++;
        }
        source += i + 1;
        target = source + n - 1;
    }
}

//--------------------------------------------------------------------------------------
// Add the -Cpr transpose to the current array
//--------------------------------------------------------------------------------------

void xFILLNEG(RealKind *Cpr, mwSignedIndex n)
{
    RealKind *source, *target;
    register mwSignedIndex i, j;
    
    source = Cpr;
    target = Cpr;
    for( i=0; i<n; i++ ) {
        for( j=i; j<n; j++ ) {
            *target -= *source;
            if( i != j ) {
                *source = -(*target);
            }
            target += n;
            source++;
        }
        source += i + 1;
        target = source;
    }
}

//------------------------------------------------------------------------------------------
// Dot Product type calculation. For conjugate cases, where ai == -1 or bi == -1, use custom
// loops instead of doing the actual multply by ai or bi. Also, use loop unrolling to speed
// up the calculations and improve the accuracy. For PC WinXP, the balance between speed and
// accuracy seemed to be optimal at a blocksize of 10. If the MATLAB mode is set, then just
// duplicate the BLAS calls that MATLAB uses (slower since it accesses the variables twice).
//------------------------------------------------------------------------------------------

struct RealKindComplex RealKindDotProduct(mwSignedIndex k,
                                          RealKind *Apr, RealKind *Api, RealKind ai, 
                                          RealKind *Bpr, RealKind *Bpi, RealKind bi)
{
    double sr = 0.0, si = 0.0;
    struct RealKindComplex z;
    mwSignedIndex i, k10, inc = 1;
    
    if( matlab ) {
        sr = xDOT( &k, Apr, &inc, Bpr, &inc );
        if( Api != NULL ) {
            si = xDOT( &k, Api, &inc, Bpr, &inc ) * ai;
            if( Bpi != NULL ) {
                sr -= xDOT( &k, Api, &inc, Bpi, &inc ) * ai * bi;
                si += xDOT( &k, Apr, &inc, Bpi, &inc ) * bi;
            }
        } else if( Bpi != NULL ) {
            si = xDOT( &k, Apr, &inc, Bpi, &inc ) * bi;
        }
        z.r = (RealKind) sr;
        z.i = (RealKind) si;
        return z;
    }
    
    k10 = k % 10;
    if( Api != NULL ) {
        if( Bpi != NULL ) {                    // (complex vector) dot (complex vector)
            if( ai == one ) {
                if( bi == one ) {
                    for( i=0; i<k10; i++ ) {
                        sr += Apr[i] * Bpr[i] - Api[i] * Bpi[i];
                        si += Api[i] * Bpr[i] + Apr[i] * Bpi[i];
                    }
                    for( i=k10; i<k; i+=10 ) {
                        sr += Apr[i  ] * Bpr[i  ] - Api[i  ] * Bpi[i  ]
                           +  Apr[i+1] * Bpr[i+1] - Api[i+1] * Bpi[i+1]
                           +  Apr[i+2] * Bpr[i+2] - Api[i+2] * Bpi[i+2]
                           +  Apr[i+3] * Bpr[i+3] - Api[i+3] * Bpi[i+3]
                           +  Apr[i+4] * Bpr[i+4] - Api[i+4] * Bpi[i+4]
                           +  Apr[i+5] * Bpr[i+5] - Api[i+5] * Bpi[i+5]
                           +  Apr[i+6] * Bpr[i+6] - Api[i+6] * Bpi[i+6]
                           +  Apr[i+7] * Bpr[i+7] - Api[i+7] * Bpi[i+7]
                           +  Apr[i+8] * Bpr[i+8] - Api[i+8] * Bpi[i+8]
                           +  Apr[i+9] * Bpr[i+9] - Api[i+9] * Bpi[i+9];
                        si += Api[i  ] * Bpr[i  ] + Apr[i  ] * Bpi[i  ]
                           +  Api[i+1] * Bpr[i+1] + Apr[i+1] * Bpi[i+1]
                           +  Api[i+2] * Bpr[i+2] + Apr[i+2] * Bpi[i+2]
                           +  Api[i+3] * Bpr[i+3] + Apr[i+3] * Bpi[i+3]
                           +  Api[i+4] * Bpr[i+4] + Apr[i+4] * Bpi[i+4]
                           +  Api[i+5] * Bpr[i+5] + Apr[i+5] * Bpi[i+5]
                           +  Api[i+6] * Bpr[i+6] + Apr[i+6] * Bpi[i+6]
                           +  Api[i+7] * Bpr[i+7] + Apr[i+7] * Bpi[i+7]
                           +  Api[i+8] * Bpr[i+8] + Apr[i+8] * Bpi[i+8]
                           +  Api[i+9] * Bpr[i+9] + Apr[i+9] * Bpi[i+9];
                    }
                } else {
                    for( i=0; i<k10; i++ ) {
                        sr += Apr[i] * Bpr[i] + Api[i] * Bpi[i];
                        si += Api[i] * Bpr[i] - Apr[i] * Bpi[i];
                    }
                    for( i=k10; i<k; i+=10 ) {
                        sr += Apr[i  ] * Bpr[i  ] + Api[i  ] * Bpi[i  ]
                           +  Apr[i+1] * Bpr[i+1] + Api[i+1] * Bpi[i+1]
                           +  Apr[i+2] * Bpr[i+2] + Api[i+2] * Bpi[i+2]
                           +  Apr[i+3] * Bpr[i+3] + Api[i+3] * Bpi[i+3]
                           +  Apr[i+4] * Bpr[i+4] + Api[i+4] * Bpi[i+4]
                           +  Apr[i+5] * Bpr[i+5] + Api[i+5] * Bpi[i+5]
                           +  Apr[i+6] * Bpr[i+6] + Api[i+6] * Bpi[i+6]
                           +  Apr[i+7] * Bpr[i+7] + Api[i+7] * Bpi[i+7]
                           +  Apr[i+8] * Bpr[i+8] + Api[i+8] * Bpi[i+8]
                           +  Apr[i+9] * Bpr[i+9] + Api[i+9] * Bpi[i+9];
                        si += Api[i  ] * Bpr[i  ] - Apr[i  ] * Bpi[i  ]
                           +  Api[i+1] * Bpr[i+1] - Apr[i+1] * Bpi[i+1]
                           +  Api[i+2] * Bpr[i+2] - Apr[i+2] * Bpi[i+2]
                           +  Api[i+3] * Bpr[i+3] - Apr[i+3] * Bpi[i+3]
                           +  Api[i+4] * Bpr[i+4] - Apr[i+4] * Bpi[i+4]
                           +  Api[i+5] * Bpr[i+5] - Apr[i+5] * Bpi[i+5]
                           +  Api[i+6] * Bpr[i+6] - Apr[i+6] * Bpi[i+6]
                           +  Api[i+7] * Bpr[i+7] - Apr[i+7] * Bpi[i+7]
                           +  Api[i+8] * Bpr[i+8] - Apr[i+8] * Bpi[i+8]
                           +  Api[i+9] * Bpr[i+9] - Apr[i+9] * Bpi[i+9];
                    }
                }
            } else {
                if( bi == one ) {
                    for( i=0; i<k10; i++ ) {
                        sr += Apr[i] * Bpr[i] + Api[i] * Bpi[i];
                        si += Apr[i] * Bpi[i] - Api[i] * Bpr[i];
                    }
                    for( i=k10; i<k; i+=10 ) {
                        sr += Apr[i  ] * Bpr[i  ] + Api[i  ] * Bpi[i  ]
                           +  Apr[i+1] * Bpr[i+1] + Api[i+1] * Bpi[i+1]
                           +  Apr[i+2] * Bpr[i+2] + Api[i+2] * Bpi[i+2]
                           +  Apr[i+3] * Bpr[i+3] + Api[i+3] * Bpi[i+3]
                           +  Apr[i+4] * Bpr[i+4] + Api[i+4] * Bpi[i+4]
                           +  Apr[i+5] * Bpr[i+5] + Api[i+5] * Bpi[i+5]
                           +  Apr[i+6] * Bpr[i+6] + Api[i+6] * Bpi[i+6]
                           +  Apr[i+7] * Bpr[i+7] + Api[i+7] * Bpi[i+7]
                           +  Apr[i+8] * Bpr[i+8] + Api[i+8] * Bpi[i+8]
                           +  Apr[i+9] * Bpr[i+9] + Api[i+9] * Bpi[i+9];
                        si += Apr[i  ] * Bpi[i  ] - Api[i  ] * Bpr[i  ]
                           +  Apr[i+1] * Bpi[i+1] - Api[i+1] * Bpr[i+1]
                           +  Apr[i+2] * Bpi[i+2] - Api[i+2] * Bpr[i+2]
                           +  Apr[i+3] * Bpi[i+3] - Api[i+3] * Bpr[i+3]
                           +  Apr[i+4] * Bpi[i+4] - Api[i+4] * Bpr[i+4]
                           +  Apr[i+5] * Bpi[i+5] - Api[i+5] * Bpr[i+5]
                           +  Apr[i+6] * Bpi[i+6] - Api[i+6] * Bpr[i+6]
                           +  Apr[i+7] * Bpi[i+7] - Api[i+7] * Bpr[i+7]
                           +  Apr[i+8] * Bpi[i+8] - Api[i+8] * Bpr[i+8]
                           +  Apr[i+9] * Bpi[i+9] - Api[i+9] * Bpr[i+9];
                    }
                } else {
                    for( i=0; i<k10; i++ ) {
                        sr += Apr[i] * Bpr[i] - Api[i] * Bpi[i];
                        si -= Api[i] * Bpr[i] + Apr[i] * Bpi[i];
                    }
                    for( i=k10; i<k; i+=10 ) {
                        sr += Apr[i  ] * Bpr[i  ] - Api[i  ] * Bpi[i  ]
                           +  Apr[i+1] * Bpr[i+1] - Api[i+1] * Bpi[i+1]
                           +  Apr[i+2] * Bpr[i+2] - Api[i+2] * Bpi[i+2]
                           +  Apr[i+3] * Bpr[i+3] - Api[i+3] * Bpi[i+3]
                           +  Apr[i+4] * Bpr[i+4] - Api[i+4] * Bpi[i+4]
                           +  Apr[i+5] * Bpr[i+5] - Api[i+5] * Bpi[i+5]
                           +  Apr[i+6] * Bpr[i+6] - Api[i+6] * Bpi[i+6]
                           +  Apr[i+7] * Bpr[i+7] - Api[i+7] * Bpi[i+7]
                           +  Apr[i+8] * Bpr[i+8] - Api[i+8] * Bpi[i+8]
                           +  Apr[i+9] * Bpr[i+9] - Api[i+9] * Bpi[i+9];
                        si -= Api[i  ] * Bpr[i  ] + Apr[i  ] * Bpi[i  ]
                           +  Api[i+1] * Bpr[i+1] + Apr[i+1] * Bpi[i+1]
                           +  Api[i+2] * Bpr[i+2] + Apr[i+2] * Bpi[i+2]
                           +  Api[i+3] * Bpr[i+3] + Apr[i+3] * Bpi[i+3]
                           +  Api[i+4] * Bpr[i+4] + Apr[i+4] * Bpi[i+4]
                           +  Api[i+5] * Bpr[i+5] + Apr[i+5] * Bpi[i+5]
                           +  Api[i+6] * Bpr[i+6] + Apr[i+6] * Bpi[i+6]
                           +  Api[i+7] * Bpr[i+7] + Apr[i+7] * Bpi[i+7]
                           +  Api[i+8] * Bpr[i+8] + Apr[i+8] * Bpi[i+8]
                           +  Api[i+9] * Bpr[i+9] + Apr[i+9] * Bpi[i+9];
                    }
                }
            }
            z.i = (RealKind) si;
        } else {                                  // (complex vector) dot (real vector)
            for( i=0; i<k10; i++ ) {
                sr += Apr[i] * Bpr[i];
                si += Api[i] * Bpr[i];
            }
            for( i=k10; i<k; i+=10 ) {
                sr += Apr[i  ] * Bpr[i  ]
                   +  Apr[i+1] * Bpr[i+1]
                   +  Apr[i+2] * Bpr[i+2]
                   +  Apr[i+3] * Bpr[i+3]
                   +  Apr[i+4] * Bpr[i+4]
                   +  Apr[i+5] * Bpr[i+5]
                   +  Apr[i+6] * Bpr[i+6]
                   +  Apr[i+7] * Bpr[i+7]
                   +  Apr[i+8] * Bpr[i+8]
                   +  Apr[i+9] * Bpr[i+9];
                si += Api[i  ] * Bpr[i  ]
                   +  Api[i+1] * Bpr[i+1]
                   +  Api[i+2] * Bpr[i+2]
                   +  Api[i+3] * Bpr[i+3]
                   +  Api[i+4] * Bpr[i+4]
                   +  Api[i+5] * Bpr[i+5]
                   +  Api[i+6] * Bpr[i+6]
                   +  Api[i+7] * Bpr[i+7]
                   +  Api[i+8] * Bpr[i+8]
                   +  Api[i+9] * Bpr[i+9];
            }
            z.i = (RealKind) (si * ai);
        }
    } else {
        if( Bpi != NULL ) {                    // (real vector) dot (complex vector)
            for( i=0; i<k10; i++ ) {
                sr += Apr[i] * Bpr[i];
                si += Apr[i] * Bpi[i];
            }
            for( i=k10; i<k; i+=10 ) {
                sr += Apr[i  ] * Bpr[i  ]
                   +  Apr[i+1] * Bpr[i+1]
                   +  Apr[i+2] * Bpr[i+2]
                   +  Apr[i+3] * Bpr[i+3]
                   +  Apr[i+4] * Bpr[i+4]
                   +  Apr[i+5] * Bpr[i+5]
                   +  Apr[i+6] * Bpr[i+6]
                   +  Apr[i+7] * Bpr[i+7]
                   +  Apr[i+8] * Bpr[i+8]
                   +  Apr[i+9] * Bpr[i+9];
                si += Apr[i  ] * Bpi[i  ]
                   +  Apr[i+1] * Bpi[i+1]
                   +  Apr[i+2] * Bpi[i+2]
                   +  Apr[i+3] * Bpi[i+3]
                   +  Apr[i+4] * Bpi[i+4]
                   +  Apr[i+5] * Bpi[i+5]
                   +  Apr[i+6] * Bpi[i+6]
                   +  Apr[i+7] * Bpi[i+7]
                   +  Apr[i+8] * Bpi[i+8]
                   +  Apr[i+9] * Bpi[i+9];
            }
            z.i = (RealKind) (si * bi);
        } else {                                  // (real vector) dot (real vector)
            for( i=0; i<k10; i++ ) {
                sr += Apr[i] * Bpr[i];
            }
            for( i=k10; i<k; i+=10 ) {
                sr += Apr[i  ] * Bpr[i  ]
                   +  Apr[i+1] * Bpr[i+1]
                   +  Apr[i+2] * Bpr[i+2]
                   +  Apr[i+3] * Bpr[i+3]
                   +  Apr[i+4] * Bpr[i+4]
                   +  Apr[i+5] * Bpr[i+5]
                   +  Apr[i+6] * Bpr[i+6]
                   +  Apr[i+7] * Bpr[i+7]
                   +  Apr[i+8] * Bpr[i+8]
                   +  Apr[i+9] * Bpr[i+9];
            }
        }
    }
    z.r = (RealKind) sr;
    return z;
}

//----------------------------------------------------------------------------------------
// Outer Product calculation. Use custom loops for all of the special cases involving
// conjugates and transposes to minimize the total number of operations involved.
//----------------------------------------------------------------------------------------

void RealKindOuterProduct(mwSignedIndex m, mwSignedIndex n,
                          RealKind *Apr, RealKind *Api, char transa, 
                          RealKind *Bpr, RealKind *Bpi, char transb,
                          RealKind *Cpr, RealKind *Cpi)
{
    register mwSignedIndex i, j; 
    mwSignedIndex kk;

    kk = 0;
    if( Api != NULL ) {
        if( Bpi != NULL ) {
            if( (transa == 'C' || transa == 'G') && (transb == 'C' || transb == 'G') ) {
                for( j=0; j<n; j++ ) { // (ar + bi*i)(C or G) * (br + bi*i)(C or G)
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] =   Apr[i] * Bpr[j] - Api[i] * Bpi[j];
                        Cpi[kk] = - Apr[i] * Bpi[j] - Api[i] * Bpr[j];
                        kk++;
                    }
                }
            } else if( (transa == 'C' || transa == 'G') && (transb == 'N' || transb == 'T') ) {
                for( j=0; j<n; j++ ) { // (ar + bi*i)(C or G) * (br + bi*i)(N or T)
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] = Apr[i] * Bpr[j] + Api[i] * Bpi[j];
                        Cpi[kk] = Apr[i] * Bpi[j] - Api[i] * Bpr[j];
                        kk++;
                    }
                }
            } else if( (transa == 'N' || transa == 'T') && (transb == 'C' || transb == 'G') ) {
                for( j=0; j<n; j++ ) { // (ar + bi*i)(N or T) * (br + bi*i)(C or G)
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] = Apr[i] * Bpr[j] + Api[i] * Bpi[j];
                        Cpi[kk] = Api[i] * Bpr[j] - Apr[i] * Bpi[j];
                        kk++;
                    }
                }
            } else {                   // (ar + bi*i)(N or T) * (br + bi*i)(N or T)
                for( j=0; j<n; j++ ) {
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] = Apr[i] * Bpr[j] - Api[i] * Bpi[j];
                        Cpi[kk] = Apr[i] * Bpi[j] + Api[i] * Bpr[j];
                        kk++;
                    }
                }
            }
        } else {
            if( transa == 'C' || transa == 'G' ) {
                for( j=0; j<n; j++ ) { // (ar + bi*i)(C or G) * (br)(N or T or C or G)
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] =   Apr[i] * Bpr[j];
                        Cpi[kk] = - Api[i] * Bpr[j];
                        kk++;
                    }
                }
            } else {
                for( j=0; j<n; j++ ) { // (ar + bi*i)(N or T) * (br)(N or T or C or G)
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] = Apr[i] * Bpr[j];
                        Cpi[kk] = Api[i] * Bpr[j];
                        kk++;
                    }
                }
            }
        }
    } else {
        if( Bpi != NULL ) {
            if( transb == 'C' || transb == 'G' ) {
                for( j=0; j<n; j++ ) { // (ar)(N or T or C or G) * (br + bi*i)(C or G)
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] =   Apr[i] * Bpr[j];
                        Cpi[kk] = - Apr[i] * Bpi[j];
                        kk++;
                    }
                }
            } else {
                for( j=0; j<n; j++ ) { // (ar)(N or T or C or G) * (br + bi*i)(N or T)
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] = Apr[i] * Bpr[j];
                        Cpi[kk] = Apr[i] * Bpi[j];
                        kk++;
                    }
                }
            }
        } else {
            for( j=0; j<n; j++ ) { // (ar)(N or T or C or G) * (br)(N or T or C or G)
                for( i=0; i<m; i++ ) {
                    Cpr[kk++] = Apr[i] * Bpr[j];
                }
            }
        }
    }
}

//--------------------------------------------------------------------------------------

mxArray *RealScalarTimesReal(mxArray *A, char transa, mwSize m1, mwSize n1,
                             mxArray *B, char transb, mwSize m2, mwSize n2)
{
    mwSize nzmax, Bndim, Cndim, Cp, p, Bsize, Csize;
    mwSize *Bdims, *Cdims;
    mwSignedIndex mn, n, k;
    mxArray *C, *result, *Bt = NULL;
    RealKind *Apr, *Api, *Bpr, *Bpi, *Cpr, *Cpi;
    RealKind ar, ai, br, bi, sr, si;
    char trans;
    mxComplexity complexflag;
    mwIndex *jc;
//-----
    if( m2 == 1 && n2 == 1 ) {  // Make sure scalar is A and array is B
        if( m1 == 1 && n1 == 1 ) {  // Check for scalar * scalar
            Apr = mxGetData(A);
            ar = *Apr;
            Api = mxGetImagData(A);
            ai = Api ? *Api : zero;
            if( transa == 'C' || transa == 'G' ) {
                ai = -ai;
            }
            Bpr = mxGetData(B);
            br = *Bpr;
            Bpi = mxGetImagData(B);
            bi = Bpi ? *Bpi : zero;
            if( transb == 'C' || transb == 'G' ) {
                bi = -bi;
            }
            sr = (RealKind) (((double)ar) * ((double)br) - ((double)ai) * ((double)bi));
            si = (RealKind) (((double)ar) * ((double)bi) + ((double)ai) * ((double)br));
            complexflag = (si == zero) ? mxREAL : mxCOMPLEX;
            if( mxIsSparse(A) || mxIsSparse(B) ) {
                result = mxCreateSparse(1, 1, 1, complexflag);
                *mxGetIr(result) = 0;
                jc = mxGetJc(result);
                jc[0] = 0;
                jc[1] = 1;
            } else {
                result = mxCreateNumericMatrix(1, 1, MatlabReturnType, complexflag);
            }
            Cpr = mxGetData(result);
            *Cpr = (RealKind) sr;
            if( complexflag == mxCOMPLEX ) {
                Cpi = mxGetImagData(result);
                *Cpi = (RealKind) si;
            }
            return result;
        } else {
            C = A;
            A = B;
            B = C;
            trans  = transa;
            transa = transb;
            transb = trans;
        }
    }
    
//--------------------------------------------------------------------------------------
// Check for multiplying by 1 and no actual transpose or conjugate is involved. In this
// case just return a shared data copy, which is very fast since no data copying is done.
// Also, if there *is* a transpose but the first two dimensions are a vector and there
// is no actual conjugate involved, we can return a shared data copy with a slight
// modification of the dimensions ... just switch the first two.
//--------------------------------------------------------------------------------------
    
    Apr = mxGetData(A);
    Api = mxGetImagData(A);
    ar = *Apr;
    ai = Api ? ((transa == 'C' || transa == 'G') ? -(*Api) :*Api) : zero;
    
    Bpr = mxGetData(B);
    Bpi = mxGetImagData(B);    
    Bndim = mxGetNumberOfDimensions(B);
    Bdims = mxGetDimensions(B);
    
    if( ar == one && ai == zero ) {
        if( transb == 'N' || (transb == 'G' && Bpi == NULL) ) {
            result = mxCreateSharedDataCopy(B);
            return result;
        } else if( (Bdims[0] == 1 || Bdims[1] == 1) && (transb == 'T' || (transb == 'C' && Bpi == NULL)) ) {
            result = mxCreateSharedDataCopy(B);
            Cdims = mxMalloc( Bndim * sizeof(*Cdims) );
            for( Cp=2; Cp<Bndim; Cp++) {
                Cdims[Cp] = Bdims[Cp];
            }
            Cdims[0] = Bdims[1];
            Cdims[1] = Bdims[0];
            mxSetDimensions(result,Cdims,Bndim);
            mxFree(Cdims);
            return result;
        }
    }
    
//--------------------------------------------------------------------------------------
// For sparse matrix, do the transpose now and then do the scalar multiply. That way if
// B is real and A is complex we are only doing the transpose work on one array.
//--------------------------------------------------------------------------------------
    
    if( (transb == 'T' || transb == 'C') && mxIsSparse(B) ) {
        mexCallMATLAB(1, &Bt, 1, &B, "transpose");
        B = Bt;
        transb = (transb == 'T') ? 'N' : 'G';
    }
        
    if( mxIsSparse(B) ) {
        jc = mxGetJc(B);
        n2 = mxGetN(B);
        n = jc[n2];
    } else {
        n = mxGetNumberOfElements(B);
    }
    if( n < 0 ) {
        mexErrMsgTxt("Number of elements too large ... overflows a signed integer");
    }
        
    complexflag = (n > 0 && (mxIsComplex(A) || mxIsComplex(B))) ? mxCOMPLEX : mxREAL;
    
//-------------------------------------------------------------------------------
// Construct the dimensions of the result. Also use the p variable to keep track
// of the total number of individual matrix multiples that are involved. The
// first two dimensions are simply the result of a single matrix multiply, with
// accouting for the transb pre-operation. The remaining dimensions are simply
// copied from B.
//-------------------------------------------------------------------------------
    
    Cndim = Bndim;
    Cdims = mxMalloc( Cndim * sizeof(*Cdims) );
    if( transb == 'N' || transb == 'G' ) {
        Cdims[0] = Bdims[0];
        Cdims[1] = Bdims[1];
    } else {
        Cdims[0] = Bdims[1];
        Cdims[1] = Bdims[0];
    }
    p = 1;
    for( Cp=2; Cp<Cndim; Cp++) {
        p *= (Cdims[Cp] = Bdims[Cp]);
    }
    
//------------------------------------------------------------------------------
// Create output array
//------------------------------------------------------------------------------
    
    if( mxIsSparse(B) ) {
        result = mxCreateSparse(Cdims[0], Cdims[1], n, complexflag);
        memcpy(mxGetIr(result), mxGetIr(B), n * sizeof(mwIndex));
        memcpy(mxGetJc(result), mxGetJc(B), (n2+1) * sizeof(mwIndex));
    } else if( mxGetNumberOfElements(B) == 0 ) {
        result = mxCreateNumericArray(Cndim, Cdims, MatlabReturnType, mxREAL);
        mxFree(Cdims);
        return result;
    } else {
        result = mxCreateNumericArray(Cndim, Cdims, MatlabReturnType, complexflag);
    }
    mxFree(Cdims);
    C = result;
    Cpr = mxGetData(C);
    Cpi = mxGetImagData(C);
    
    m2 = Bdims[0];
    n2 = Bdims[1];
    mn = m2 * n2;

    if( n == 0 ) {  // If result is empty, just return right now
        return result;
    }
    C = result;
    Cpr = mxGetData(C);
    Cpi = mxGetImagData(C);
    Bpr = mxGetData(B);
    Bpi = mxGetImagData(B);

//------------------------------------------------------------------------------------------------
// Check for matlab flag. If set, then use a BLAS call for this function.
//------------------------------------------------------------------------------------------------
    
//    if( matlab ) {
        // Future upgrade ... insert BLAS calls here and return?
//    }

//------------------------------------------------------------------------------------------------
// If the matrix is really a vector, then there is no need to do an actual transpose. We can
// simply strip the transpose operation away and rely on the dimensions to do the transpose.
//------------------------------------------------------------------------------------------------
    
    if( m2 ==1 || n2 == 1 ) {
        if( transb == 'T' ) transb = 'N';
        if( transb == 'C' ) transb = 'G';
    }
    
//------------------------------------------------------------------------------------------------
// Some specialized cases, no need to multiply by +1 or -1, we can program that directly into the
// calculations without a multiply. We do need to multiply by zero, however, so that the sign of
// any -0 that might be present gets carried over into the result, and also any inf or NaN that is
// present gets a proper result. So no special code for multiplying by zero.
//------------------------------------------------------------------------------------------------
    
    if( ar == one ) {
        if( ai == one ) {
            if( transb == 'N' ) {
                RealKindEqP1P1TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); // C = (1 + 1*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqP1P1TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); // C = (1 + 1*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqP1P1TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (1 + 1*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqP1P1TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (1 + 1*i) * (Bpr + Bpi * i)C
            }
        } else if( ai == -one ) {
            if( transb == 'N' ) {
                RealKindEqP1M1TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); // C = (1 - 1*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqP1M1TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); // C = (1 - 1*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqP1M1TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (1 - 1*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqP1M1TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (1 - 1*i) * (Bpr + Bpi * i)C
            }
        } else if( ai == zero ) {
            if( transb == 'N' ) {  // this case never reached ... it is the shared data copy above
                RealKindEqP1P0TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); // C = (1 + 0*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqP1P0TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); // C = (1 + 0*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqP1P0TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (1 + 0*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqP1P0TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (1 + 0*i) * (Bpr + Bpi * i)C
            }
        } else {
            if( transb == 'N' ) {
                RealKindEqP1PxTimesRealKindN(Cpr, Cpi, ai, Bpr, Bpi, n); // C = (1 + ai*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqP1PxTimesRealKindG(Cpr, Cpi, ai, Bpr, Bpi, n); // C = (1 + ai*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqP1PxTimesRealKindT(Cpr, Cpi, ai, Bpr, Bpi, m2, n2, p); // C = (1 + ai*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqP1PxTimesRealKindC(Cpr, Cpi, ai, Bpr, Bpi, m2, n2, p); // C = (1 + ai*i) * (Bpr + Bpi * i)C
            }
        }
    } else if( ar == -one ) {
        if( ai == one ) {
            if( transb == 'N' ) {
                RealKindEqM1P1TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); // C = (-1 + 1*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqM1P1TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); // C = (-1 + 1*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqM1P1TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (-1 + 1*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqM1P1TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (-1 + 1*i) * (Bpr + Bpi * i)C
            }
        } else if( ai == -one ) {
            if( transb == 'N' ) {
                RealKindEqM1M1TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); // C = (-1 - 1*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqM1M1TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); // C = (-1 - 1*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqM1M1TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (-1 - 1*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqM1M1TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (-1 - 1*i) * (Bpr + Bpi * i)C
            }
        } else if( ai == zero ) {
            if( transb == 'N' ) {
                RealKindEqM1P0TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); // C = (-1 + 0*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqM1P0TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); // C = (-1 + 0*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqM1P0TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (-1 + 0*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqM1P0TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); // C = (-1 + 0*i) * (Bpr + Bpi * i)C
            }
        } else {
            if( transb == 'N' ) {
                RealKindEqM1PxTimesRealKindN(Cpr, Cpi, ai, Bpr, Bpi, n); // C = (-1 + ai*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqM1PxTimesRealKindG(Cpr, Cpi, ai, Bpr, Bpi, n); // C = (-1 + ai*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqM1PxTimesRealKindT(Cpr, Cpi, ai, Bpr, Bpi, m2, n2, p); // C = (-1 + ai*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqM1PxTimesRealKindC(Cpr, Cpi, ai, Bpr, Bpi, m2, n2, p); // C = (-1 + ai*i) * (Bpr + Bpi * i)C
            }
        }
    } else {  // ar != one && ar != -one
        if( ai == one ) {
            if( transb == 'N' ) {
                RealKindEqPxP1TimesRealKindN(Cpr, Cpi, ar, Bpr, Bpi, n); // C = (ar + 1*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqPxP1TimesRealKindG(Cpr, Cpi, ar, Bpr, Bpi, n); // C = (ar + 1*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqPxP1TimesRealKindT(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); // C = (ar + 1*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqPxP1TimesRealKindC(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); // C = (ar + 1*i) * (Bpr + Bpi * i)C
            }
        } else if( ai == -one ) {
            if( transb == 'N' ) {
                RealKindEqPxM1TimesRealKindN(Cpr, Cpi, ar, Bpr, Bpi, n); // C = (ar - 1*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqPxM1TimesRealKindG(Cpr, Cpi, ar, Bpr, Bpi, n); // C = (ar - 1*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqPxM1TimesRealKindT(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); // C = (ar - 1*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqPxM1TimesRealKindC(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); // C = (ar - 1*i) * (Bpr + Bpi * i)C
            }
        } else if( ai == zero ) {
            if( transb == 'N' ) {
                RealKindEqPxP0TimesRealKindN(Cpr, Cpi, ar, Bpr, Bpi, n); // C = (ar + 0*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqPxP0TimesRealKindG(Cpr, Cpi, ar, Bpr, Bpi, n); // C = (ar + 0*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqPxP0TimesRealKindT(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); // C = (ar + 0*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqPxP0TimesRealKindC(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); // C = (ar + 0*i) * (Bpr + Bpi * i)C
            }
        } else {
            if( transb == 'N' ) {
                RealKindEqPxPxTimesRealKindN(Cpr, Cpi, ar, ai, Bpr, Bpi, n); // C = (ar + ai*i) * (Bpr + Bpi * i)
            } else if( transb == 'G' ) {
                RealKindEqPxPxTimesRealKindG(Cpr, Cpi, ar, ai, Bpr, Bpi, n); // C = (ar + ai*i) * (Bpr - Bpi * i)
            } else if( transb == 'T' ) {
                RealKindEqPxPxTimesRealKindT(Cpr, Cpi, ar, ai, Bpr, Bpi, m2, n2, p); // C = (ar + ai*i) * (Bpr + Bpi * i)T
            } else { // if( transb == 'C' ) {
                RealKindEqPxPxTimesRealKindC(Cpr, Cpi, ar, ai, Bpr, Bpi, m2, n2, p); // C = (ar + ai*i) * (Bpr + Bpi * i)C
            }
        }
    }
    
//-------------------------------------------------------------------------------------------
// If the imaginary part is all zero, then free it and set the pointer to NULL.
//-------------------------------------------------------------------------------------------
    
    if( AllRealZero(Cpi, n) ) {
        mxFree(Cpi);
        Cpi = NULL;
        mxSetImagData(C, NULL);
    }
    
//-------------------------------------------------------------------------------------------
// Clean up sparse matrix and realloc if appropriate. Also do the transpose now if necessary.
//-------------------------------------------------------------------------------------------
    
    if( mxIsSparse(C) ) {
        nzmax = mxGetNzmax(C);
        k = spclean(C);
        if( nzmax - k > REALLOCTOL ) {
            mxSetPr(C, mxRealloc(Cpr, k * sizeof(RealKind)));
            mxSetIr(C, mxRealloc(mxGetIr(C), k * sizeof(mwIndex)));
            if( Cpi != NULL ) {
                mxSetPi(C, mxRealloc(Cpi, k * sizeof(RealKind)));
            }
            mxSetNzmax(C, k);
        }
    }
    if( Bt != NULL ) {
        mxDestroyArray(Bt);
    }
    return result;
}
