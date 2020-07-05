#include "mex.h"
#include <math.h>
#include <cmath>
#include "maxNumCompThreads.h"
#include <cstdlib>

/*   undef needed for LCC compiler  */
#undef EXTERN_C
#ifdef _WIN32
	#include <windows.h>
	#include <process.h>
#else
	#include <pthread.h>
#endif

#ifndef min
#define min(X, Y)  (X < Y ? X : Y)
#endif
#ifndef max
#define max(X, Y)  (X < Y ? Y : X)
#endif
#ifndef isnan 
#define isnan(x) ((x)!=(x)) 
#endif
/* haralick2mex -- Haralick for 2D images. Syntax:
 * haralickims = haralick2mex(double image, double graylevels, double window_size, double dist, double background [optional])
  *
 *To compile, it is necessary to link against libut:
 *-WINDOWS (64-bit, visual studio):
 *	mex('-v','-largeArrayDims','haralick2mexmt.cpp',[matlabroot '\extern\lib\win64\microsoft\libut.lib'])
 *-WINDOWS (32-bit, visual studio):
 *  mex('-v','-largeArrayDims','haralick2mexmt.cpp',[matlabroot '\extern\lib\win32\microsoft\libut.lib'])
 *-WINDOWS (32-bit, lcc):
 *  mex('-v','-largeArrayDims','haralick2mexmt.cpp',[matlabroot '\extern\lib\win32\lcc\libut.lib'])
 *-LINUX: mex -v -largeArrayDims haralick2mexmt.cpp -lut
 */

struct haralick2_threadargs{
    double *image;
    double *haralicks;
    int ws;
    int dist;
    int graylevels;
    int background;
    int rows;
    int cols;
    int nharalicks;
    size_t ThreadID;
    size_t Nthreads;
};

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

/* GLOBALS */
size_t THREADLIMIT=4;

// Variables used to detect if the threads are finished
// (volatile: reload the variable instead of using the value available in a register)
static volatile int WaitForThreads[256];

// Finished percentage of the total input numbers
static volatile double percentage;

// Mutex used to lock percentage variable, to allow only one thread to 
// read or write it at the same time.
#ifdef _WIN32
static HANDLE percentageMutex;
#else
static pthread_mutex_t percentageMutex;
#endif

void graycomtx(const double *image, double *comtx, int ws, int dist, 
        int graylevels, const int background, int rows, int cols, int row, int col) {
    
    int i, j, k, l, centerind, pixind, center_value, hws;
    int d_start_row, d_start_col, d_end_row, d_end_col;
    int block_start_row, block_start_col, block_end_row, block_end_col;
    
    for (i = 0; i < graylevels*graylevels; i++)
        comtx[i] = 0.0;
    
    hws=(int) floor((float) ws/2);
    block_start_row = max(0, row-hws);
    block_start_col = max(0, col-hws);
    block_end_row = min(rows-1, row+hws);
    block_end_col = min(cols-1, col+hws);
    
    for (j = block_start_col; j <= block_end_col; j++) {
        for (i = block_start_row; i <= block_end_row; i++) {
            centerind=i+j*rows;
			if (isnan(image[centerind]) || (image[centerind] == background))
                continue;

            center_value = (int) image[centerind];

            d_start_row = max((int) 0, i-dist);
            d_start_col = max((int) 0, j-dist);
            d_end_row = min((int) rows-1, i+dist);
            d_end_col = min((int) cols-1, j+dist);
            for (l = d_start_col; l <= d_end_col; l++) {
                for (k = d_start_row; k <= d_end_row; k++) {
                    pixind=k+l*rows;
					if (!isnan(image[pixind]) && (image[pixind]!=background)) {
						//if ((dist==0) || (pixind!=centerind)) //either dist=0 or exclude dist=0
                        comtx[center_value + (int) (image[pixind]+0.5)*graylevels] += 1;
                    }
                }
			}
		    
        }
	}
	//testing purposes
    //for (i = 0; i < graylevels*graylevels; i++) {
      //  if ((int) comtx[i] != 0)
      //      mexPrintf("\nWith window centered at: [%i][%i]: comtx[%i]=%i",row,col,i,(int) comtx[i]);
    //}
   // mexPrintf(" (all else zeros)\n");
    //
}

#ifdef _WIN32
unsigned __stdcall haralick2(void *ThreadArgsV) {
#else
void *haralick2(void *ThreadArgsV) {
#endif
    int i, j, k, ii, jj, nbins, nzeros, nnonzeros, onepct, onequarterpct, somepct, tenpct, ip;
    int imagepixelind, imagepixels;
    int *hi, *hj, *himhj, *hiphj;
    double *comtx, *p, *pnz, *nzcomtx, *px, *py, *pxplusy, *pxminusy;
    double entropyval, energyval, inertiaval, idmval, 
            correlationval, info1val, info2val, H1, H2,
            sigma_x, sigma_y, mu_x, mu_y, h_x, h_y, h_max,
            saval, svval, seval, daval, dvval, deval, cosum;
    
    haralick2_threadargs *ThreadArgs=(haralick2_threadargs *) ThreadArgsV;
    
    double *image, *haralicks;
    int ws, dist, graylevels, background, rows, cols, nharalicks, ThreadOffset, Nthreads;
    
    image=ThreadArgs->image;
    haralicks=ThreadArgs->haralicks;
    ws=ThreadArgs->ws;
    dist=ThreadArgs->dist;
    graylevels=ThreadArgs->graylevels;
    background=ThreadArgs->background;
    rows=ThreadArgs->rows;
    cols=ThreadArgs->cols;
    nharalicks=ThreadArgs->nharalicks;
    ThreadOffset=ThreadArgs->ThreadID;
    Nthreads=ThreadArgs->Nthreads;
    
    imagepixels=rows*cols;
    
    nbins=graylevels*graylevels;
    onepct = (int) floor(.01*imagepixels-1);
    onequarterpct = (int) floor(.0025*imagepixels-1);
    somepct = (int) floor(.025*imagepixels-1);
    tenpct = (int) floor(.1*imagepixels-1);
    
    for(j = 0; j < cols; j++)
        for(i = 0; i < rows; i++)
            if(image[i+j*rows] >= graylevels && image[i+j*rows]!=background) {
                    WaitForThreads[ThreadOffset]=0;
                    #ifdef _WIN32
                    _endthreadex( 0 );
                    return 0;
                    #else
                    pthread_exit(NULL);
                    #endif
                    //mexErrMsgTxt("Graylevels of image fall outside acceptable range.");
                }
    
    comtx = (double *) malloc(nbins*sizeof(double));
    nzcomtx = (double *) malloc(nbins*sizeof(double));
    
    p = (double *) malloc(nbins*sizeof(double));
    pnz = (double *) malloc(nbins*sizeof(double));
    px = (double *) malloc(graylevels*sizeof(double));
    py = (double *) malloc(graylevels*sizeof(double));
    pxplusy = (double *) malloc(2*graylevels*sizeof(double));
    pxminusy = (double *) malloc(graylevels*sizeof(double));
    
    hi = (int *) malloc(nbins*sizeof(int));
    hj = (int *) malloc(nbins*sizeof(int));
    himhj = (int *) malloc(nbins*sizeof(int));
    hiphj = (int *) malloc(nbins*sizeof(int));
    
    for (j=ThreadOffset, ip=1; j<cols; j+=Nthreads) {
        for (i=0; i<rows; i++, ip++) {
            if (image[i+j*rows]!=background) {
            /* Get co-occurrence matrix */
            graycomtx(image, comtx, ws, dist, graylevels, background, rows, cols, i, j);
            
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
            imagepixelind=i + j*rows;
            haralicks[imagepixelind+0*imagepixels]=entropyval;
            haralicks[imagepixelind+1*imagepixels]=energyval;
            haralicks[imagepixelind+2*imagepixels]=inertiaval;
            haralicks[imagepixelind+3*imagepixels]=idmval;
            haralicks[imagepixelind+4*imagepixels]=correlationval;
            haralicks[imagepixelind+5*imagepixels]=info1val;
            haralicks[imagepixelind+6*imagepixels]=info2val;
            haralicks[imagepixelind+7*imagepixels]=saval;
            haralicks[imagepixelind+8*imagepixels]=svval;
            haralicks[imagepixelind+9*imagepixels]=seval;
            haralicks[imagepixelind+10*imagepixels]=daval;
            haralicks[imagepixelind+11*imagepixels]=dvval;
            haralicks[imagepixelind+12*imagepixels]=deval;
            
            if ((imagepixelind % 2000) == 0 && utIsInterruptPending()) {
                //utSetInterruptPending(false);
                free(comtx); free(nzcomtx);
                free(p); free(pnz); free(px); free(py); free(pxplusy); free(pxminusy);
                free(hi); free(hj); free(himhj); free(hiphj);
                WaitForThreads[ThreadOffset]=0;
                #ifdef _WIN32
                _endthreadex( 0 );
                return 0;
                #else
                pthread_exit(NULL);
                #endif
            }
            
            if ((ip % onequarterpct)==0) {
                #ifdef _WIN32
                WaitForSingleObject(percentageMutex, INFINITE);
                #else
                pthread_mutex_lock(&percentageMutex);
                #endif
                percentage+=onequarterpct;//100/((double)DataSize[0]);
                #ifdef _WIN32
                ReleaseMutex(percentageMutex);
                #else
                pthread_mutex_unlock(&percentageMutex);
                #endif
            }
            
            } else {
                imagepixelind=i + j*rows;
                for (k=0; k<nharalicks; k++) haralicks[imagepixelind+k*imagepixels]=0;
            }
        } /* rows */
    } /* cols */
    
    free(comtx);
    free(nzcomtx);
    
    free(p);
    free(pnz);
    free(px);
    free(py);
    free(pxplusy);
    free(pxminusy);
    
    free(hi);
    free(hj);
    free(himhj);
    free(hiphj);
    
    WaitForThreads[ThreadOffset]=0;
    
    /* end thread */
    #ifdef _WIN32
	_endthreadex( 0 );
    return 0;
	#else
	pthread_exit(NULL);
	#endif
    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *haralicks, *image;
    int dist, rows, cols, i, lastpct=0, pct;
    int graylevels, background, ws;
    mwSize dims[3];
    int nharalicks=13;
    size_t nmaxthreads;
    bool waitingforthreads;
    double impixels;
    //unsigned short int *X;
    mexPrintf("\haralick2mexmt called.\n");

    if(nrhs > 6 || nrhs < 4)
        mexErrMsgTxt("haralick2mexmt(image,graylevels,ws,dist,[background],[nmaxthreads])");
    
    if(!mxIsDouble(prhs[0]))
        mexErrMsgTxt("Input image must be DOUBLE.");
    
    image = mxGetPr(prhs[0]);
    rows = (int) mxGetM(prhs[0]);
    cols = (int) mxGetN(prhs[0]);
    impixels=rows*cols;
    graylevels=(int) mxGetScalar(prhs[1]);
    ws = (int) mxGetScalar(prhs[2]);
    dist = (int) mxGetScalar(prhs[3]);
    if (nrhs>=5 && !mxIsEmpty(prhs[4]))
        background = (int) mxGetScalar(prhs[4]);
    else
        background = -1;
    
    if (nrhs>=6)
        nmaxthreads=(size_t) mxGetScalar(prhs[5]);
	else
        nmaxthreads=1024;
    
    if(graylevels < 0 || graylevels > 65535)
        mexErrMsgTxt("GRAYLEVELS must be between 0 and 2^16-1.");
    
    dims[0] = rows; dims[1] = cols; dims[2] = nharalicks;
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    haralicks = mxGetPr(plhs[0]);
    
    // Threading
    #ifdef _WIN32
    percentageMutex = CreateMutex(NULL, FALSE, NULL);
    #else
    pthread_mutex_init(&percentageMutex, NULL);
    #endif
    percentage=0;
    
    size_t *ThreadID;
    haralick2_threadargs **ThreadArgs;
    haralick2_threadargs *ThreadArgsi;
    
    /* Handles to the worker threads */
	#ifdef _WIN32
		HANDLE *ThreadList; 
    #else
		pthread_t *ThreadList;
	#endif
    
    size_t Nthreads;
    size_t maxthreads=(size_t) getNumCores();
    Nthreads=(maxthreads>THREADLIMIT ? THREADLIMIT : maxthreads); // hard thread limit
    if (Nthreads>nmaxthreads) Nthreads=nmaxthreads; // input thread limit
    //mexPrintf("%d\n",Nthreads); mexEvalString("drawnow");
    
    //WaitForThreads=(volatile int *) calloc(Nthreads,sizeof(int));
    
    /* Reserve room for handles of threads in ThreadList */
    #ifdef _WIN32
		ThreadList = (HANDLE *) malloc(Nthreads*sizeof(HANDLE));
    #else
		ThreadList = (pthread_t *) malloc(Nthreads*sizeof(pthread_t));
	#endif
	
	ThreadID = (size_t *) malloc(Nthreads*sizeof(size_t));
	ThreadArgs = (haralick2_threadargs **) malloc(Nthreads*sizeof(haralick2_threadargs *));
    
    // Find kNNs for each point by traversing kd tree
    for (i=0; i<Nthreads; i++) {
        
        ThreadID[i]=i;
        
        ThreadArgsi = (haralick2_threadargs *) malloc(sizeof(haralick2_threadargs));
        ThreadArgsi->image=image;
        ThreadArgsi->haralicks=haralicks;
        ThreadArgsi->ws=ws;
        ThreadArgsi->dist=dist;
        ThreadArgsi->graylevels=graylevels;
        ThreadArgsi->background=background;
        ThreadArgsi->rows=rows;
        ThreadArgsi->cols=cols;
        ThreadArgsi->nharalicks=nharalicks;
        ThreadArgsi->ThreadID=ThreadID[i];
        ThreadArgsi->Nthreads=Nthreads;
        
        /* Start thread  */
        ThreadArgs[i]=ThreadArgsi; // now we can overwrite ThreadArgsi for the next thread
        WaitForThreads[i]=1;
        
        #ifdef _WIN32
            ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &haralick2, ThreadArgs[i] , 0, NULL );
        #else
            pthread_create((pthread_t*)&ThreadList[i], NULL, &haralick2, ThreadArgs[i]);
        #endif
    }
    
    mexPrintf("0%%");
    mexEvalString("drawnow");
    waitingforthreads=true;
    while (waitingforthreads) {
        for (i=0, waitingforthreads=false; i<Nthreads; i++) waitingforthreads=waitingforthreads || WaitForThreads[i]==1;
        #ifdef _WIN32
        WaitForSingleObject(percentageMutex, INFINITE);
        #else
        pthread_mutex_lock(&percentageMutex);
        #endif
        pct=ceil(100. * percentage/impixels);
        if (pct>lastpct) {
            lastpct=pct;
            if (((int) pct % 10) == 0)
                mexPrintf("%d%%",(int) pct);
            else if (((int) pct % 2) == 0)
                mexPrintf(".");
            mexEvalString("drawnow");
        }
        #ifdef _WIN32
        ReleaseMutex(percentageMutex);
        Sleep( 1 );
        #else
        pthread_mutex_unlock(&percentageMutex);
        usleep(1e3);
        #endif
    }
    
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* threads don't clear it */
        mexPrintf("<user interrupted>");
    }
    
    mexPrintf("\n");
    
    #ifdef _WIN32
            for (i=0; i<Nthreads; i++) WaitForSingleObject(ThreadList[i], INFINITE);
            for (i=0; i<Nthreads; i++) CloseHandle( ThreadList[i] );
            CloseHandle(percentageMutex);
    #else
            for (i=0; i<Nthreads; i++) pthread_join(ThreadList[i], NULL);
            pthread_mutex_destroy(&percentageMutex);
    #endif
    
	for (i=0; i<Nthreads; i++) free(ThreadArgs[i]);
    free(ThreadArgs);
    free(ThreadID);
    free(ThreadList);
    
    //haralick2(image,haralicks,ws,dist,graylevels,background,rows,cols,nharalicks);
}
