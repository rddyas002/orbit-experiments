#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>

#include "sgp4ext.h"
#include "sgp4unit.h"
#include "sgp4io.h"

#define TLE_FILE    "REFERENCE.TLE"
#define TLE_DATA	"SUMBANDILA.TLE"

// main functions
void asgp4(const double jdm_abs, const double GPS_r[3], const double GPS_v[3], const double q_gps, double ro[3], double vo[3], double delta_kep[4]);
void modifyTleParameters(double delta_param[4]);
double time_correction_calculation(const double GPS_r[3], const double GPS_v[3], const double SGP4_r[3], const double SGP4_v[3]);
void eci2kep(const double r[3], const double v[3], double *incl, double *raan, double *e, double *w);
void load_TLEs();
void init_filters();

// helper functions
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