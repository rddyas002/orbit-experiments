/*
 * File : aSGP4.c
 * Abstract:
 *       Implementation of NRossouw MSc Thesis
 *		 Using Sumbandila's orbit from 01/02/2010-07/02/2010
 *  
 * Reference for SGP4 source/headers https://celestrak.com/software/vallado-sw.asp
 */


#define S_FUNCTION_NAME  aSGP4
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>

#include "asgp4_util.h"
#include "mex.h"

#define BLOCK_SAMPLE    1

static double jdm_sim_start = 0;

/************************ Standard S-function interface ************************/
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 1);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    if (!ssSetNumInputPorts(S, 3)) return;
    ssSetInputPortWidth(S, 0, 1);				
    ssSetInputPortDirectFeedThrough(S, 0, 1);	 // time in seconds
    ssSetInputPortWidth(S, 1, 6);
    ssSetInputPortDirectFeedThrough(S, 1, 1);    // GPS: X,Y,Z,Vx,Vy,Vz
    ssSetInputPortWidth(S, 2, 1);
    ssSetInputPortDirectFeedThrough(S, 2, 1);    // enable aSGP4

    if (!ssSetNumOutputPorts(S,3)) return;
    ssSetOutputPortWidth(S, 0, 3);
    ssSetOutputPortWidth(S, 1, 3);
    ssSetOutputPortWidth(S, 2, 4);

    ssSetNumSampleTimes(S, 1);

    /* specify the sim state compliance to be same as a built-in block */
    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);

    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S,
                 SS_OPTION_WORKS_WITH_CODE_REUSE |
                 SS_OPTION_EXCEPTION_FREE_CODE |
                 SS_OPTION_USE_TLC_WITH_ACCELERATOR);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{   
    ssSetSampleTime(S, 0, BLOCK_SAMPLE);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S); 
    
    real_T year = ( mxGetPr(ssGetSFcnParam(S,0)) )[0];
    real_T month = ( mxGetPr(ssGetSFcnParam(S,0)) )[1];
    real_T day = ( mxGetPr(ssGetSFcnParam(S,0)) )[2];
    real_T hour = ( mxGetPr(ssGetSFcnParam(S,0)) )[3];
    real_T min = ( mxGetPr(ssGetSFcnParam(S,0)) )[4];
    real_T sec = ( mxGetPr(ssGetSFcnParam(S,0)) )[5];    
      
	jday((int) year, (int) month, (int) day, (int) hour, (int) min, (double) sec, jdm_sim_start);
    // convert jd to jdm
    jdm_sim_start *= 24*60;
    
    /* Initialise SGP4 with TLE */
    load_TLEs();
	init_filters();
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    double ro[3];
    double vo[3];
    
    int_T             i;
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType pntGPS = ssGetInputPortRealSignalPtrs(S,1);
    InputRealPtrsType enable_aSGP4 = ssGetInputPortRealSignalPtrs(S,2);
    real_T            *y0    = ssGetOutputPortRealSignal(S,0);
    real_T            *y1    = ssGetOutputPortRealSignal(S,1);
    real_T            *y2    = ssGetOutputPortRealSignal(S,2);
    int_T             width = ssGetOutputPortWidth(S,0);

    double minutes_from_start_relative = *uPtrs[0]/60;      // convert seconds to jdm
    double minutes_from_start_absolute = minutes_from_start_relative + jdm_sim_start;
    
    double GPS_r[3], GPS_v[3];
    GPS_r[0] = *pntGPS[0];
    GPS_r[1] = *pntGPS[1];
    GPS_r[2] = *pntGPS[2];
    GPS_v[0] = *pntGPS[3];
    GPS_v[1] = *pntGPS[4];
    GPS_v[2] = *pntGPS[5];
    double q_gps = *enable_aSGP4[0];
    double delta_kep[4] = {0};
    
	asgp4(minutes_from_start_absolute, GPS_r, GPS_v, q_gps, ro, vo, delta_kep);

    *y0++ = ro[0];
    *y0++ = ro[1];
    *y0 = ro[2];
    *y1++ = vo[0];
    *y1++ = vo[1];
    *y1 = vo[2];
        
    *y2++ = delta_kep[0];
    *y2++ = delta_kep[1];
    *y2++ = delta_kep[2];
    *y2++ = delta_kep[3];
}

static void mdlTerminate(SimStruct *S)
{
}
/************************ Standard S-function interface ************************/


#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
