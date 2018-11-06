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

#include "sgp4ext.h"
#include "sgp4unit.h"
#include "sgp4io.h"
#include "mex.h"

#define TLE_FILE    "REFERENCE.TLE"
#define TLE_DATA	"SUMBANDILA.TLE"
#define TLE_LOAD_MAX	25
#define BLOCK_SAMPLE    1

/* Global variables */
elsetrec satrec;
gravconsttype  whichconst;
static double jdm_sim_start = 0;
static char longstr1_tle[TLE_LOAD_MAX][130];
static char longstr2_tle[TLE_LOAD_MAX][130];
static int tles_loaded = 0;
static double filter_state[4] = {0};  

// main functions
static void modifyTleParameters(double delta_param[4]);
static void eci2kep(const double r[3], const double v[3], double *incl, double *raan, double *e, double *w);
static void load_TLEs();

// helper functions
static real_T jday(real_T yr, real_T mon, real_T day, real_T hr, real_T b_min, real_T sec);
static void calc_eci2ntw_dcm(const double r[3], const double v[3], double ECI2NTW_DCM[3][3]);
static void map_eci2ntw(const double ECI2NTW_DCM[3][3], const double u[3], double y[3]);
static void crossProduct3(const double a[3], const double b[3], double c[3]);
static void vectorMult3d(const double a[3], const double b, double c[3]);
static void vectorDiv3d(const double a[3], const double b, double c[3]);
static void vectorAdd3d(const double a[3], const double b[3], double c[3]);
static void vectorSubtract3d(const double a[3], const double b[3], double c[3]);
static double norm3(const double a[3]);
static double dotProduct3(const double a[3], const double b[3]);
static void normaliseVector3(double v[3]);
static void copyVec3(const double a[3], double b[3]);

/************************ Standard S-function interface ************************/
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 1);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    if (!ssSetNumInputPorts(S, 2)) return;
    ssSetInputPortWidth(S, 0, 1);				
    ssSetInputPortDirectFeedThrough(S, 0, 1);	 // time in seconds
    ssSetInputPortWidth(S, 1, 6);
    ssSetInputPortDirectFeedThrough(S, 1, 1);    // GPS: X,Y,Z,Vx,Vy,Vz

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
    
    const real_T year = ( mxGetPr(ssGetSFcnParam(S,0)) )[0];
    const real_T month = ( mxGetPr(ssGetSFcnParam(S,0)) )[1];
    const real_T day = ( mxGetPr(ssGetSFcnParam(S,0)) )[2];
    const real_T hour = ( mxGetPr(ssGetSFcnParam(S,0)) )[3];
    const real_T min = ( mxGetPr(ssGetSFcnParam(S,0)) )[4];
    const real_T sec = ( mxGetPr(ssGetSFcnParam(S,0)) )[5];    
      
    jdm_sim_start = jday(year, month, day, hour, min, sec)*24*60;
    
    /* Initialise SGP4 with TLE */
    load_TLEs();
	
	filter_state[0] = filter_state[1] = filter_state[2] = filter_state[3] = 0;
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
	const double alpha[4] = {0.0001,0.0001,0.0001,0.0001};
    double ro[3];
    double vo[3];
    
    int_T             i;
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType pntGPS = ssGetInputPortRealSignalPtrs(S,1);
    real_T            *y0    = ssGetOutputPortRealSignal(S,0);
    real_T            *y1    = ssGetOutputPortRealSignal(S,1);
    real_T            *y2    = ssGetOutputPortRealSignal(S,2);
    int_T             width = ssGetOutputPortWidth(S,0);

    double minutes_from_start_relative = *uPtrs[0]/60;
    double minutes_from_start_absolute = minutes_from_start_relative + jdm_sim_start;
    // compute julian day minutes from epoch to request time
    double delta_t_pred_mjd = minutes_from_start_absolute - satrec.jdsatepoch*24*60;
    sgp4(whichconst, satrec,  delta_t_pred_mjd, ro,  vo);
    
	// compute position error in ECI
	double error_r_eci[3], error_r_ntw[3];
	double GPS_r[3], GPS_v[3];
	GPS_r[0] = *pntGPS[0];
	GPS_r[1] = *pntGPS[1];
	GPS_r[2] = *pntGPS[2];
	GPS_v[0] = *pntGPS[3];
	GPS_v[1] = *pntGPS[4];
	GPS_v[2] = *pntGPS[5];
	// compute instaneous kepler elements
	double GPS_incl=0, GPS_RAAN=0, GPS_e=0, GPS_w=0;
	double SGP4_incl=0, SGP4_RAAN=0, SGP4_e=0, SGP4_w=0;
	double delta_incl, delta_RAAN, delta_e, delta_w;
	eci2kep(GPS_r, GPS_v, &GPS_incl, &GPS_RAAN, &GPS_e, &GPS_w);
	eci2kep(ro, vo, &SGP4_incl, &SGP4_RAAN, &SGP4_e, &SGP4_w);
	
	delta_incl = GPS_incl - SGP4_incl;
	delta_RAAN = GPS_RAAN - SGP4_RAAN;
	delta_e = GPS_e - SGP4_e;
	delta_w = GPS_w - SGP4_w;
	
	filter_state[0] = delta_incl*alpha[0] + (1 - alpha[0])*filter_state[0];
	filter_state[1] = delta_RAAN*alpha[1] + (1 - alpha[1])*filter_state[1];
	filter_state[2] = delta_e*alpha[2] + (1 - alpha[2])*filter_state[2];
	filter_state[3] = delta_w*alpha[3] + (1 - alpha[3])*filter_state[3];
	
	modifyTleParameters(filter_state);	
	
	// error mapping
	error_r_eci[0] = GPS_r[0] - ro[0];
	error_r_eci[1] = GPS_r[1] - ro[1];
	error_r_eci[2] = GPS_r[2] - ro[2];
	// compute ECI2NTW DCM
	double ECI2NTW_DCM[3][3];
	calc_eci2ntw_dcm(ro, vo, ECI2NTW_DCM);
	// map error into NTW
	map_eci2ntw(ECI2NTW_DCM, error_r_eci, error_r_ntw);
	
	*y2++ = filter_state[0];
	*y2++ = filter_state[1];
	*y2++ = filter_state[2];
	*y2 = filter_state[3];
	
    *y0++ = ro[0];
    *y0++ = ro[1];
    *y0 = ro[2];
    *y1++ = vo[0];
    *y1++ = vo[1];
    *y1 = vo[2];
}
/************************ Standard S-function interface ************************/

static void load_TLEs(){
    FILE *infile;
    double sec,  jd, rad, tsince, startmfe, stopmfe, deltamin;
    char typerun, typeinput, opsmode;
    double ro[3];
    double vo[3];
    
    infile = fopen(TLE_FILE, "r");
    
    // Set mode variables
    opsmode = 'i';
    typerun = 'm';
    typeinput = 'v';
    whichconst = wgs84;
    
	tles_loaded = 0;
    // load lines into buffer
	while(infile != NULL){
		fgets(&longstr1_tle[tles_loaded][0],130,infile);
		fgets(&longstr2_tle[tles_loaded][0],130,infile);	
		tles_loaded++;
	}

    fclose(infile);
    twoline2rv(&longstr1_tle[0][0], &longstr2_tle[0][0], typerun, typeinput, opsmode, whichconst, 
        startmfe, stopmfe, deltamin, satrec);
    
    // call the propagator to get the initial state vector value
    sgp4 (whichconst, satrec,  0.0, ro,  vo);    
}

static void modifyTleParameters(double delta_param[4]){
    // need to reinitialise sgp4 parameters
	elsetrec& satrec_tmp;
	int cardnumb_tmp;
    long revnum_tmp = 0;
	// edit longstr2_tle
	if (longstr2_tle[52] == ' ')
		sscanf(longstr2_tle,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %10lf %6ld \n",
			&cardnumb_tmp,&satrec_tmp.satnum, &satrec_tmp.inclo,
            &satrec_tmp.nodeo,&satrec_tmp.ecco, &satrec_tmp.argpo, &satrec_tmp.mo, &satrec_tmp.no,
            &revnum_tmp );
    else
		sscanf(longstr2_tle,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld \n",
			&cardnumb_tmp,&satrec_tmp.satnum, &satrec_tmp.inclo,
            &satrec_tmp.nodeo,&satrec_tmp.ecco, &satrec_tmp.argpo, &satrec_tmp.mo, &satrec_tmp.no,
            &revnum_tmp );
	
    //satrec.inclo += delta_param[0];
    //satrec.nodeo += delta_param[1];  
	twoline2rv(longstr1_tle, longstr2_tle, 'm', 'v', 'i', wgs84, 
        startmfe, stopmfe, deltamin, satrec);
}

static void eci2kep(const double r[3], const double v[3], double *incl, double *raan, double *e, double *w){
	double h[3];
	double mu = 3.986004418e5; //km^3/s^2
	crossProduct3(r, v, h);
	double h_norm = norm3(h);
	double v_norm = norm3(v);
	double khat[3] = {0,0,1};
	double n[3];
	crossProduct3(khat, h, n);
	double n_norm = norm3(n);
	double r_norm = norm3(r);
	double vh[3];
	crossProduct3(v, h, vh);
	vectorDiv3d(vh, mu, vh);
	double r_normalised[3];
	copyVec3(r,r_normalised);
	normaliseVector3(r_normalised);
	double e_vec[3];
	vectorSubtract3d(vh, r_normalised, e_vec);
	double e_norm = norm3(e_vec);
	*e = e_norm;
	
	*incl = acos(h[2]/h_norm);
	// compute specific energy
	double E = pow(v_norm,2)/2 - mu/r_norm;
	// semi-major axis
	double a = -mu/(2*E); 		// alternative	a = dot(h,h)/(mu*(1-e^2));
	*raan = acos(n[0]/n_norm);	//	khat = [0;0;1];	n = cross(khat,h);	W = acos(n(1)/norm(n));
	if (n[1] < 0)
		*raan = 2*M_PI - *raan;
	*w = acos(dotProduct3(n,e_vec)/(n_norm*e_norm));
	if (e_vec[2] < 0)
		*w = 2*M_PI - *w;
}

static void map_eci2ntw(const double ECI2NTW_DCM[3][3], const double u[3], double y[3]){
    y[0] = ECI2NTW_DCM[0][0]*u[0] + ECI2NTW_DCM[0][1]*u[1] + ECI2NTW_DCM[0][2]*u[2];
    y[1] = ECI2NTW_DCM[1][0]*u[0] + ECI2NTW_DCM[1][1]*u[1] + ECI2NTW_DCM[1][2]*u[2];
    y[2] = ECI2NTW_DCM[2][0]*u[0] + ECI2NTW_DCM[2][1]*u[1] + ECI2NTW_DCM[2][2]*u[2];
}

static void calc_eci2ntw_dcm(const double r[3], const double v[3], double ECI2NTW_DCM[3][3]){
    double T[3];	// in-line with velocity vector
    double W[3];	//
    double N[3];	// opp nadir
      
    T[0] = v[0];
    T[1] = v[1];
    T[2] = v[2];
    normaliseVector3(T);
    
    crossProduct3(r, v, W);
    normaliseVector3(W);
    
    crossProduct3(T, W, N);
    
    ECI2NTW_DCM[0][0] = N[0]; ECI2NTW_DCM[0][1] = N[1]; ECI2NTW_DCM[0][2] = N[2];
    ECI2NTW_DCM[1][0] = T[0]; ECI2NTW_DCM[1][1] = T[1]; ECI2NTW_DCM[1][2] = T[2];
    ECI2NTW_DCM[2][0] = W[0]; ECI2NTW_DCM[2][1] = W[1]; ECI2NTW_DCM[2][2] = W[2];
}

static void crossProduct3(const double a[3], const double b[3], double c[3]) {
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

static double dotProduct3(const double a[3], const double b[3]) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void vectorMult3d(const double a[3], const double b, double c[3]) {
	c[0] = a[0]*b;
	c[1] = a[1]*b;
	c[2] = a[2]*b;
}

static void vectorDiv3d(const double a[3], const double b, double c[3]) {
	c[0] = a[0]/b;
	c[1] = a[1]/b;
	c[2] = a[2]/b;
}

static void vectorAdd3d(const double a[3], const double b[3], double c[3]){
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
}

static void vectorSubtract3d(const double a[3], const double b[3], double c[3]){
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

static double norm3(const double a[3]) {
	return sqrt(dotProduct3(a, a));
}

static void copyVec3(const double a[3], double b[3]){
	b[0] = a[0];
	b[1] = a[1];
	b[2] = a[2];
}

static void normaliseVector3(double v[3]) {
    double norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v[0] /= norm;
    v[1] /= norm;
    v[2] /= norm;
}

static real_T jday(real_T yr, real_T mon, real_T day, real_T hr, real_T b_min, real_T sec){
  real_T jd;

  /*  Computes the Julian day given YY:MM:DD-hh:mm:ss format */
  /*   references    : */
  /*     vallado       2007, 189, alg 14, ex 3-14 */
  /*  ----------------------------------------------------------------------------- */
  jd = ((((367.0 * yr - floor(7.0 * (yr + floor((mon + 9.0) / 12.0)) * 0.25)) +
          floor(275.0 * mon / 9.0)) + day) + 1.7210135E+6) + ((sec / 60.0 +
    b_min) / 60.0 + hr) / 24.0;

  /*   - 0.5 * sign(100.0 * yr + mon - 190002.5) + 0.5; */
  return jd;
}

/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S)
{
}


#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
