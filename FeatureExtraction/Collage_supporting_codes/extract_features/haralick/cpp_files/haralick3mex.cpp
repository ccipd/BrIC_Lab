#include "mex.h"
#include <math.h>
#include <cmath>
#include <cstdlib>

#ifdef _WIN32
	#include <windows.h>
#else
	#include <unistd.h>
#endif

#ifndef min
#define min(X, Y)  (X < Y ? X : Y)
#endif
#ifndef max
#define max(X, Y)  (X < Y ? Y : X)
#endif

/* haralick3mex -- Haralick for 3D volumes. Syntax:
 * haralickims = haralick3mex(double volume, double graylevels, double window_size, double dist, double background [optional]) 
 *
 *To compile, it is necessary to link against libut:
 *-WINDOWS (64-bit, visual studio):
 *	mex('-v','-largeArrayDims','haralick3mex.cpp',[matlabroot '\extern\lib\win64\microsoft\libut.lib'])
 *-WINDOWS (32-bit, visual studio):
 *  mex('-v','-largeArrayDims','haralick3mex.cpp',[matlabroot '\extern\lib\win32\microsoft\libut.lib'])
 *-WINDOWS (32-bit, lcc):
 *  mex('-v','-largeArrayDims','haralick3mex.cpp',[matlabroot '\extern\lib\win32\lcc\libut.lib'])
 *-LINUX: mex -v -largeArrayDims haralick3mex.cpp -lut
 */

inline bool greater (int i,int j) { return (j<i); }

inline double logb(double x, double b) { return log(x)/log(b); }

// prototype the break handling functions in libut (C library)
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
    extern "C" void utSetInterruptPending(bool);
#else
    extern bool utIsInterruptPending();
    extern void utSetInterruptPending(bool);
#endif

void graycomtx(double ***imagemat, double *comtx, int ws, int dist, 
        int graylevels, const int background, int rows, int cols, int slices, 
        int row, int col, int slice) {
    
    int h, i, j, k, l, m, center_value, hws, slicepixels;
    //int centerind, pixind;
    int d_start_row, d_start_col, d_start_slice, d_end_row, d_end_col, d_end_slice;
    int block_start_row, block_start_col, block_start_slice;
    int block_end_row, block_end_col, block_end_slice;
    
    for (i = 0; i < graylevels*graylevels; i++)
        comtx[i] = 0.0;

    slicepixels=rows*cols;
    
    hws=(int) floor((float) ws/2);
    block_start_row = max(0, row-hws);
    block_start_col = max(0, col-hws);
    block_start_slice = max(0, slice-hws);
    block_end_row = min(rows-1, row+hws);
    block_end_col = min(cols-1, col+hws);
    block_end_slice = min(slices-1, slice+hws);
    
    for (h = block_start_slice; h <= block_end_slice; h++)
    for (j = block_start_col; j <= block_end_col; j++)
        for (i = block_start_row; i <= block_end_row; i++) {
            //centerind=i+j*rows+h*slicepixels;
            center_value = (int) imagemat[i][j][h];
            if (center_value == background)
                continue;

            d_start_row = max((int) 0, i-dist);
            d_start_col = max((int) 0, j-dist);
            d_start_slice = max((int) 0, h-dist);
            d_end_row = min((int) rows-1, i+dist);
            d_end_col = min((int) cols-1, j+dist);
            d_end_slice = min((int) slices-1, h+dist);
            for (m = d_start_slice; m <= d_end_slice; m++)
            for (l = d_start_col; l <= d_end_col; l++)
                for (k = d_start_row; k <= d_end_row; k++) {
                    //pixind=k+l*rows+m*slicepixels;
                    //if (imagemat[k][l][m]!=background && pixind!=centerind)
                    if (imagemat[k][l][m]!=background)
                        comtx[center_value + (int) (imagemat[k][l][m]+0.5)*graylevels] += 1;
                }
            
            comtx[center_value + center_value*graylevels] -= 1;
            //if (comtx[center_value + center_value*graylevels]<0) mexErrMsgTxt("Crap.");
        }
    /*
    for(k = block_start_slice; k < block_end_slice; k++)
        for(j = block_start_col; j < block_end_col; j++)
            for(i = block_start_row; i < block_end_row; i++) {
                center_value = (int) imagemat[i][j][k];
                if (center_value!=background)
                    comtx[center_value + center_value*graylevels] -= 1;
            }
     */
}

void haralick3(double ***imagemat, double *haralicks, int ws, int dist, int graylevels, int background, int rows, int cols, int slices, int nharalicks) {
    
    int h, i, j, k, ii, jj, nbins, nzeros, nnonzeros, somepct, tenpct; //, pynzs, pxnzs;
    int volumepixelind, volumepixels;
    int *hi, *hj, *himhj, *hiphj;
    double *comtx, *p, *pnz, *nzcomtx, *px, *py, *pxplusy, *pxminusy;
    double entropyval, energyval, inertiaval, idmval, 
            correlationval, info1val, info2val, H1, H2,
            sigma_x, sigma_y, mu_x, mu_y, h_x, h_y, h_max,
            saval, svval, seval, daval, dvval, deval, cosum;
    
    volumepixels=rows*cols*slices;
    
    nbins=graylevels*graylevels;
    somepct = (int) floor(.025*rows*cols*slices-1);
    tenpct = (int) floor(.1*rows*cols*slices-1);
    
    comtx = (double *) mxMalloc(nbins*sizeof(double));
    nzcomtx = (double *) mxMalloc(nbins*sizeof(double));
    
    p = (double *) mxMalloc(nbins*sizeof(double));
    pnz = (double *) mxMalloc(nbins*sizeof(double));
    px = (double *) mxMalloc(graylevels*sizeof(double));
    py = (double *) mxMalloc(graylevels*sizeof(double));
    pxplusy = (double *) mxMalloc(2*graylevels*sizeof(double));
    pxminusy = (double *) mxMalloc(graylevels*sizeof(double));
    
    hi = (int *) mxMalloc(nbins*sizeof(int));
    hj = (int *) mxMalloc(nbins*sizeof(int));
    himhj = (int *) mxMalloc(nbins*sizeof(int));
    hiphj = (int *) mxMalloc(nbins*sizeof(int));
    
    for(k = 0; k < slices; k++)
        for(j = 0; j < cols; j++)
            for(i = 0; i < rows; i++)
                if(imagemat[i][j][k] >= graylevels && imagemat[i][j][k]!=background)
                    mexErrMsgTxt("Graylevels of image fall outside acceptable range.");
    
    for (h=0; h<slices; ++h) {
    for (j=0; j<cols; ++j) {
        for (i=0; i<rows; ++i) {
            if (imagemat[i][j][h]!=background) {
            /* Get co-occurrence matrix */
            graycomtx(imagemat, comtx, ws, dist, graylevels, background, rows, cols, slices, i, j, h);
            
            /* Initialize feature values */
            entropyval=0; energyval=0; inertiaval=0; idmval=0;
            correlationval=0; info1val=0; info2val=0;
            saval=0; svval=0; seval=0; daval=0; dvval=0; deval=0;
            H1=0; H2=0; h_x=0; h_y=0; h_max=0; mu_x=0; mu_y=0; sigma_x=0; sigma_y=0;
            cosum=0;
            
            /* Non-zero elements & locations in comtx and distribution */
            //nzeros=std::count(comtx,comtx+nbins,0);
            //nnonzeros=nbins-nzeros;
            for (k=0; k<nbins; k++) cosum+=comtx[k];
            if (cosum<2) continue;
            for (k=0, ii=0; k<nbins; k++) {
                if (comtx[k]>0) {
                    p[k]=comtx[k]/cosum;
                    pnz[ii]=p[k];
                    nzcomtx[ii]=comtx[k];
                    hi[ii]=k % graylevels;
                    hj[ii]=(int) floor((float) k/(float) graylevels);
                    himhj[ii]=hi[ii]-hj[ii];
                    hiphj[ii]=hi[ii]+hj[ii];
                    ii++;
                } else {
                    p[k]=0;
                }
            }
            nnonzeros=ii; nzeros=nbins-nnonzeros;
            
            /* Entropy, Energy, Inertial, Inv. Diff. Moment */
            for (k=0; k<nnonzeros; k++) {
                //pnz[k]=nzcomtx[k]/nbins;
                entropyval-=pnz[k]*logb(pnz[k],2.0);
                energyval+=pnz[k]*pnz[k];
                inertiaval+=himhj[k]*himhj[k]*pnz[k];
                idmval+=pnz[k]/(1.0+himhj[k]*himhj[k]);
            }
            
            /* Marginal distributions */
            for (ii=0; ii<graylevels; ii++) { px[ii]=0; py[ii]=0; }
            for (k=0, ii=0; ii<graylevels; ii++)
                for (jj=0; jj<graylevels; jj++, k++) {
                    py[ii]+=p[k];
                    px[jj]+=p[k];
                }
            /*
            for (ii=0, pynzs=0, pxnzs=0; ii<graylevels; ii++) {
                pynzs+=py[ii]>0;
                pxnzs+=px[ii]>0;
            }
            if (pynzs<2 || pxnzs<2) continue;
             */
            
            /* Correlation */
            for (ii=0; ii<graylevels; ii++) {
                h_x-=(px[ii]>0 ? px[ii]*logb(px[ii],2.0) : 0);
                h_y-=(py[ii]>0 ? py[ii]*logb(py[ii],2.0) : 0);
                mu_x+=ii*px[ii];
                mu_y+=ii*py[ii];
            }
            
            for (ii=0; ii<graylevels; ii++) {
                sigma_x+=pow(ii-mu_x,2) * px[ii];
                sigma_y+=pow(ii-mu_y,2) * py[ii];
            }
            
            if (sigma_x>(1e-4) && sigma_y>(1e-4)) {
                for (k=0; k<nnonzeros; k++)
                    correlationval+=(hi[k]-mu_x)*(hj[k]-mu_y)*pnz[k];
                correlationval/=sqrt(sigma_x*sigma_y);
            } else correlationval=0;
            
            /* Information measures of correlation */
             for (k=0, ii=0; ii<graylevels; ii++)
                for (jj=0; jj<graylevels; jj++, k++) {
                    H1-=(p[k]>0 && px[jj]>0 && py[ii]>0 ? p[k]*logb(px[jj]*py[ii],2.0) : 0);
                    H2-=(px[jj]>0 && py[ii]>0 ? px[jj]*py[ii]*logb(px[jj]*py[ii],2.0) : 0);
                }
            h_max=max(h_x,h_y);
            info1val=(h_max!=0 ? (entropyval-H1)/h_max : 0);
            info2val=sqrt(abs(1-exp(-2*(H2-entropyval))));
            
            /* Sum average, variance and entropy */
            for (k=0; k<(2*graylevels); k++)
                pxplusy[k]=0;
            for (k=0; k<nnonzeros; k++)
                pxplusy[hiphj[k]]+=pnz[k];
            
            for (k=0; k<(2*graylevels); k++) {
                saval+=k*pxplusy[k];
                seval-=(pxplusy[k]>0 ? pxplusy[k]*logb(pxplusy[k],2.0) : 0);
            }
            for (k=0; k<(2*graylevels); k++)
                svval+=pow(k-saval,2) * pxplusy[k];
                
            /* Difference average, variance and entropy */
            for (k=0; k<graylevels; k++)
                pxminusy[k]=0;
            for (k=0; k<nnonzeros; k++)
                pxminusy[abs(himhj[k])]+=pnz[k];
            
            for (k=0; k<graylevels; k++) {
                daval+=k*pxminusy[k];
                deval-=(pxminusy[k]>0 ? pxminusy[k]*logb(pxminusy[k],2.0) : 0);
            }
            for (k=0; k<graylevels; k++)
                dvval+=pow(k-daval,2) * pxminusy[k];
            
            /* Work on unsorted comtx */
            /*
            for (k=0; k<nbins; k++) {
                p[k]=comtx[k]/nbins;
                entropyval-=(p[k]>0 ? p[k]*log(p[k]) : 0);
                energyval+=p[k]*p[k];
            }
             */
            
            /* Sorted comtx */
            /*
            std::sort(comtx,comtx+nbins,greater);
            zeroloc=std::find(comtx,comtx+nbins,0);
            
            for (k=0; k<zeroloc; k++) {
                p[k]=comtx[k]/nbins;
                entropyval-=p[k]*log(p[k]);
                energyval+=p[k]*p[k];
            }
             */
            
            /* Put feature values in output volume */
            volumepixelind=i + j*rows + h*rows*cols;
            haralicks[volumepixelind+0*volumepixels]=entropyval;
            haralicks[volumepixelind+1*volumepixels]=energyval;
            haralicks[volumepixelind+2*volumepixels]=inertiaval;
            haralicks[volumepixelind+3*volumepixels]=idmval;
            haralicks[volumepixelind+4*volumepixels]=correlationval;
            haralicks[volumepixelind+5*volumepixels]=info1val;
            haralicks[volumepixelind+6*volumepixels]=info2val;
            haralicks[volumepixelind+7*volumepixels]=saval;
            haralicks[volumepixelind+8*volumepixels]=svval;
            haralicks[volumepixelind+9*volumepixels]=seval;
            haralicks[volumepixelind+10*volumepixels]=daval;
            haralicks[volumepixelind+11*volumepixels]=dvval;
            haralicks[volumepixelind+12*volumepixels]=deval;
            
            } else {
                volumepixelind=i + j*rows + h*rows*cols;
                for (k=0; k<nharalicks; k++) haralicks[volumepixelind+k*volumepixels]=0;
            }
            
            if ((volumepixelind % 2000) == 0 && utIsInterruptPending()) {
                utSetInterruptPending(false);
                mexErrMsgTxt("<user interrupted>\n");
                return;
            }
            if (((volumepixelind+1) % somepct)==0) {
                mexPrintf("."); mexEvalString("drawnow");
            } else if (((volumepixelind) % tenpct)==0) {
                mexPrintf("%d%%",(int) ceil((float) 100*(volumepixelind)/(volumepixels-1))); mexEvalString("drawnow");
            }
            
        }
    }
    }
    mexPrintf("\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *haralicks, *image, ***imagemat;
    int dist, rows, cols, slices, i, j, k, ndims, slicepixels;
    const mwSize *imdims;
    int graylevels, background, ws;
    mwSize hardims[4];
    int nharalicks=13;
    //unsigned short int *X;
    
    if(nrhs > 5 || nrhs < 4)
        mexErrMsgTxt("haralick3mex(volume,graylevels,ws,dist,[background])");
    
    if(!mxIsDouble(prhs[0]))
        mexErrMsgTxt("Input volume must be DOUBLE.");
    
    image = mxGetPr(prhs[0]);
    imdims=mxGetDimensions(prhs[0]);
    ndims=(int) mxGetNumberOfDimensions(prhs[0]);
    if (ndims!=3) mexErrMsgTxt("Input image must have 3 dimensions.");
    rows=(int) imdims[0]; cols=(int) imdims[1]; slices=(int) imdims[2];
    slicepixels=rows*cols;
    //rows = (int) mxGetM(prhs[0]);
    //cols = (int) mxGetN(prhs[0]);
    graylevels=(int) mxGetScalar(prhs[1]);
    ws = (int) mxGetScalar(prhs[2]);
    dist = (int) mxGetScalar(prhs[3]);
    if (nrhs==4)
        background = -1;
    else
        background = (int) mxGetScalar(prhs[4]);
    
    if(graylevels < 0 || graylevels > 65535)
        mexErrMsgTxt("GRAYLEVELS must be between 0 and 2^16-1.");
    
    imagemat = (double ***) mxMalloc(rows*sizeof(double **));
    for (i=0; i<rows; i++) {
        imagemat[i] = (double **) mxMalloc(cols*sizeof(double *));
        for (j=0; j<cols; j++)
            imagemat[i][j] = (double *) mxMalloc(slices*sizeof(double));
    }
    
    for (k=0; k<slices; k++)
        for (j=0; j<cols; j++)
            for (i=0; i<rows; i++)
                imagemat[i][j][k]=image[i+j*rows+k*slicepixels];
    
    hardims[0] = rows; hardims[1] = cols; hardims[2] = slices;
    hardims[3] = nharalicks;
    plhs[0] = mxCreateNumericArray(4, hardims, mxDOUBLE_CLASS, mxREAL);
    haralicks = mxGetPr(plhs[0]);
    
    haralick3(imagemat,haralicks,ws,dist,graylevels,background,rows,cols,slices,nharalicks);
}
