#include "asgp4_util.h"

/* Global variables */
elsetrec satrec;
gravconsttype  whichconst;
static char longstr1_tle[130];
static char longstr2_tle[130];
static double filter_state[5] = {0}; 
static double prev_q_gps = 0;
static double time_correction = 0;
static double time_error = 0;

/*
 *  Takes in gps position and velocity and a gps quality factor.
 *  q_gps = -1 implies no gps data
 */
void asgp4(const double jdm_abs, const double GPS_r[3], const double GPS_v[3], const double q_gps, double ro[3], double vo[3], double delta_kep[4]){
    const double alpha[5] = {0.01,0.01,0.01,0.01,0.1};
    static double last_GPS_r[3] = {0}, last_GPS_v[3] = {0};
    static double last_SGP4_r[3] = {0}, last_SGP4_v[3] = {0};
    
    if (prev_q_gps != q_gps){
        if (q_gps != 1.0){
            // gps switched off
            // propagate last SGP4 sample
            //double delta_t_pred_mjd = jdm_abs - satrec.jdsatepoch*24*60 + time_correction/60 - 1/60;
            //sgp4(whichconst, satrec,  delta_t_pred_mjd, ro,  vo);            
            //time_error += -TIME_CORRECTION_GAIN*time_correction_calculation(last_GPS_r, last_GPS_v, ro, vo);
            //time_correction = time_error;
            time_correction += time_error;//filter_state[4];
        }
    }    
    
    // compute julian day minutes from epoch to request time
    double delta_t_pred_mjd = jdm_abs - satrec.jdsatepoch*24*60 + time_correction/60;
    sgp4(whichconst, satrec,  delta_t_pred_mjd, ro,  vo);
       
    if (q_gps != 0){
        // compute position error in ECI
        double error_r_eci[3], error_r_ntw[3];
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
        
        // handle RAAN wrapping
        if (delta_RAAN > pi)
            delta_RAAN -= 2*pi;
	
        filter_state[0] += delta_incl*alpha[0];// + (1 - alpha[0])*filter_state[0];
        filter_state[1] += delta_RAAN*alpha[1];// + (1 - alpha[1])*filter_state[1];
        filter_state[2] += delta_e*alpha[2];// + (1 - alpha[2])*filter_state[2];
        filter_state[3] += delta_w*alpha[3];// + (1 - alpha[3])*filter_state[3];
                
        modifyTleParameters(filter_state);
        time_error = -time_correction_calculation(GPS_r, GPS_v, ro, vo);
        filter_state[4] = time_error*alpha[4] + (1 - alpha[4])*filter_state[4];         
        copyVec3(GPS_r, last_GPS_r);
        copyVec3(GPS_v, last_GPS_v);
    }
    prev_q_gps = q_gps;
    
    // send to calling function
    delta_kep[0] = filter_state[0];
    delta_kep[1] = filter_state[1];
    delta_kep[2] = filter_state[2];
    delta_kep[3] = filter_state[3];
}

double time_correction_calculation(const double GPS_r[3], const double GPS_v[3], const double SGP4_r[3], const double SGP4_v[3]){
    // compute position error in ECI
    double error_r_eci[3], error_r_ntw[3];
    // SGP4 error w.r.t. GPS
    vectorSubtract3d(SGP4_r, GPS_r, error_r_eci);
    // compute ECI2NTW DCM
    double ECI2NTW_DCM[3][3];
    // NTW frame constructed w.r.t. GPS coordinates
    calc_eci2ntw_dcm(GPS_r, GPS_v, ECI2NTW_DCM); 
    // map error into NTW
    map_eci2ntw(ECI2NTW_DCM, error_r_eci, error_r_ntw);
    
    return error_r_ntw[1]/norm3(GPS_v);
}

void init_filters(void){
    prev_q_gps = 0;
    time_correction = 0;
    time_error = 0;
	filter_state[0] = filter_state[1] = filter_state[2] = filter_state[3] = filter_state[4] = 0;
}

void load_TLEs(){
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
    
	fgets(&longstr1_tle[0],130,infile);
	fgets(&longstr2_tle[0],130,infile);	
    fclose(infile);
    
    char longstr1_tle_sgp4[130], longstr2_tle_sgp4[130];
    memcpy(&longstr1_tle_sgp4[0], &longstr1_tle[0], 130);
    memcpy(&longstr2_tle_sgp4[0], &longstr2_tle[0], 130);

    twoline2rv(&longstr1_tle_sgp4[0], &longstr2_tle_sgp4[0], typerun, typeinput, opsmode, whichconst, 
        startmfe, stopmfe, deltamin, satrec); 
}

void modifyTleParameters(double delta_param[4]){
    double startmfe, stopmfe, deltamin;
    char typerun, typeinput, opsmode;    
    opsmode = 'i';
    typerun = 'm';
    typeinput = 'v';
    whichconst = wgs84;
    
    char longstr1_tle_sgp4[130], longstr2_tle_sgp4[130];
    memcpy(&longstr1_tle_sgp4[0], &longstr1_tle[0], 130);
    memcpy(&longstr2_tle_sgp4[0], &longstr2_tle[0], 130);    
    
    // need to read inclination, RAAN, eccentricity and AOP from TLE
    longstr2_tle_sgp4[25] = '.';            // add decimal point for e
    for (int j = 26; j <= 32; j++)          // ensure spaces are 0 for e
        if (longstr2_tle_sgp4[j] == ' ')
               longstr2_tle_sgp4[j] = '0';
    
    int cardnumb; long revnum = 0;
    long int  satnum;
    double inclo,nodeo,ecco,argpo,mo,no;    
    sscanf(longstr2_tle_sgp4,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld \n",
        &cardnumb,&satnum,&inclo,&nodeo,&ecco,&argpo,&mo,&no,&revnum);
    // edit inclination
    inclo += delta_param[0]*180/M_PI;
    sprintf(&longstr2_tle_sgp4[8], "%8.4f", inclo);
    longstr2_tle_sgp4[16] = ' ';
    // edit RAAN
    nodeo += delta_param[1]*180/M_PI;
    sprintf(&longstr2_tle_sgp4[17], "%8.4f", nodeo);
    longstr2_tle_sgp4[25] = ' ';
    // edit eccentricity
    ecco += delta_param[2];
    char temp_ecc[12];
    sprintf(&temp_ecc[0], "%9.7f", ecco);
    memcpy(&longstr2_tle_sgp4[26], &temp_ecc[2], 7);
    // edit AOP
    argpo += delta_param[3]*180/M_PI;
    sprintf(&longstr2_tle_sgp4[34], "%8.4f", argpo);
    longstr2_tle_sgp4[42] = ' ';
    //edit MA
    mo -= delta_param[3]*180/M_PI;
    sprintf(&longstr2_tle_sgp4[43], "%8.4f", mo);
    longstr2_tle_sgp4[51] = ' ';
    // let's try to edit the drag term too
    char temp_B[11];
    memcpy(temp_B, &longstr1_tle_sgp4[52], 10);
    if (temp_B[1] != ' ')
        temp_B[0] = temp_B[1];
    temp_B[1] = '.';
    temp_B[9] = '\n';
    int ibexp = 0;
    double bstar = 0;
    sscanf(temp_B,"%7lf%2d\n", &bstar, &ibexp);
    bstar = bstar * pow(10.0, ibexp);
    
    twoline2rv(&longstr1_tle_sgp4[0], &longstr2_tle_sgp4[0], typerun, typeinput, opsmode, whichconst, 
        startmfe, stopmfe, deltamin, satrec);  
}

void eci2kep(const double r[3], const double v[3], double *incl, double *raan, double *e, double *w){
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
