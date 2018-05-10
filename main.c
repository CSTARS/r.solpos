
/****************************************************************************
 *
 * MODULE:       r.solpos
 * AUTHOR(S):    Quinn Hart
 * PURPOSE:      Borrowed heavily from r.sunhours, this application computes, NRELs solpos on 
 *               a raster, and can provide solar azimuth, angle, sunshine hours, sunrise/sunset time,
 *               and an instantaneous sun angle.
 *               Uses NREL SOLPOS
 * COPYRIGHT:    (C) 2018 by the GRASS Development Team
 *
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with GRASS
 *               for details.
 *
 *****************************************************************************/

#define MAIN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/gprojects.h>
#include <grass/glocale.h>
#include "solpos00.h"

void set_solpos_time(struct posdata *, int, int, int, int, int, int);
void set_solpos_longitude(struct posdata *, double );
int roundoff(double *);

int main(int argc, char *argv[])
{
    struct GModule *module;
    struct
    {
      struct Option *elev, *azimuth, *sretr, *ssetr, *hrang, *sunhours, *year,
	    *month, *day, *hour, *minutes, *seconds;
    } parm;
    struct Cell_head window;
    FCELL *elevbuf, *azimuthbuf, *sretrbuf, *ssetrbuf, *hrangbuf,  *sunhourbuf;
    struct History hist;

    /* projection information of input map */
    struct Key_Value *in_proj_info, *in_unit_info;
    struct pj_info iproj;	/* input map proj parameters  */
    struct pj_info oproj;	/* output map proj parameters  */

    char *elev_name, *azimuth_name, *sretr_name, *ssetr_name, *hrang_name, *sunhour_name;
    int elev_fd, azimuth_fd, sretr_fd, ssetr_fd, hrang_fd, sunhour_fd;
    double east, east_ll, north, north_ll;
    double north_gc, north_gc_sin, north_gc_cos;  /* geocentric latitude */
    double ba2;
    int year, month, day, hour, minutes, seconds;
    int doy;   					/* day of year */
    int row, col, nrows, ncols;
    int do_reproj = 0;
    struct posdata pd;

    G_gisinit(argv[0]);

    module = G_define_module();
    G_add_keyword(_("raster"));
    G_add_keyword(_("solar"));
    G_add_keyword(_("sun energy"));
    G_add_keyword(_("sun position"));
    module->label = _("Calculates NREL's SOLPOS parameters; solar elevation, solar azimuth, sun hours sunrise/sunset, hour angle.");
    module->description = _("Solar elevation: the angle between the direction of the geometric center "
			    "of the sun's apparent disk and the (idealized) horizon. "
			    "Solar azimuth: the angle from due north in clockwise direction.");
    
    parm.elev = G_define_standard_option(G_OPT_R_OUTPUT);
    parm.elev->key = "elevation";
    parm.elev->label = _("Output raster map with solar elevation angle");
    parm.elev->required = NO;

    parm.azimuth = G_define_standard_option(G_OPT_R_OUTPUT);
    parm.azimuth->key = "azimuth";
    parm.azimuth->label = _("Output raster map with solar azimuth angle");
    parm.azimuth->required = NO;
    
    parm.sretr = G_define_standard_option(G_OPT_R_OUTPUT);
    parm.sretr->key = "sretr";
    parm.sretr->label = _("Output raster map with sunrise time, minutes from midnight");
    parm.sretr->description = _("Sunrise hour requires SOLPOS use and Greenwich standard time");
    parm.sretr->required = NO;

    parm.ssetr = G_define_standard_option(G_OPT_R_OUTPUT);
    parm.ssetr->key = "ssetr";
    parm.ssetr->label = _("Output raster map with sunset time, minutes from midnight");
    parm.ssetr->description = _("Sunset hour requires SOLPOS use and Greenwich standard time");
    parm.ssetr->required = NO;

    parm.hrang = G_define_standard_option(G_OPT_R_OUTPUT);
    parm.hrang->key = "hrang";
    parm.hrang->label = _("Output raster map with hour angle--hour of sun from solar noon, degrees WEST");
    parm.hrang->description = _("Hour angle requires SOLPOS use and Greenwich standard time");
    parm.hrang->required = NO;

    parm.sunhours = G_define_standard_option(G_OPT_R_OUTPUT);
    parm.sunhours->key = "sunhour";
    parm.sunhours->label = _("Output raster map with sunshine hours");
    parm.sunhours->description = _("Sunshine hours require SOLPOS use and Greenwich standard time");
    parm.sunhours->required = NO;

    parm.year = G_define_option();
    parm.year->key = "year";
    parm.year->type = TYPE_INTEGER;
    parm.year->required = YES;
    parm.year->description = _("Year");
    parm.year->options = "1950-2050";
    parm.year->guisection = _("Time");

    parm.month = G_define_option();
    parm.month->key = "month";
    parm.month->type = TYPE_INTEGER;
    parm.month->required = NO;
    parm.month->label = _("Month");
    parm.month->description = _("If not given, day is interpreted as day of the year");
    parm.month->options = "1-12";
    parm.month->guisection = _("Time");

    parm.day = G_define_option();
    parm.day->key = "day";
    parm.day->type = TYPE_INTEGER;
    parm.day->required = YES;
    parm.day->description = _("Day");
    parm.day->options = "1-366";
    parm.day->guisection = _("Time");

    parm.hour = G_define_option();
    parm.hour->key = "hour";
    parm.hour->type = TYPE_INTEGER;
    parm.hour->required = NO;
    parm.hour->description = _("Hour");
    parm.hour->options = "0-24";
    parm.hour->answer = "12";
    parm.hour->guisection = _("Time");

    parm.minutes = G_define_option();
    parm.minutes->key = "minute";
    parm.minutes->type = TYPE_INTEGER;
    parm.minutes->required = NO;
    parm.minutes->description = _("Minutes");
    parm.minutes->options = "0-60";
    parm.minutes->answer = "0";
    parm.minutes->guisection = _("Time");

    parm.seconds = G_define_option();
    parm.seconds->key = "second";
    parm.seconds->type = TYPE_INTEGER;
    parm.seconds->required = NO;
    parm.seconds->description = _("Seconds");
    parm.seconds->options = "0-60";
    parm.seconds->answer = "0";
    parm.seconds->guisection = _("Time");

    if (G_parser(argc, argv))
	exit(EXIT_FAILURE);

    G_get_window(&window);

    /* require at least one output */
    elev_name = parm.elev->answer;
    azimuth_name = parm.azimuth->answer;
    sretr_name = parm.sretr->answer;
    ssetr_name = parm.ssetr->answer;
    hrang_name = parm.hrang->answer;
    sunhour_name = parm.sunhours->answer;
    if (!elev_name && !azimuth_name && !sretr_name && !ssetr_name && !sunhour_name)
	G_fatal_error(_("No output requested, exiting."));

    year = atoi(parm.year->answer);
    if (parm.month->answer)
	month = atoi(parm.month->answer);
    else
	month = -1;

    day = atoi(parm.day->answer);
    hour = atoi(parm.hour->answer);
    minutes = atoi(parm.minutes->answer);
    seconds = atoi(parm.seconds->answer);

    /* init variables */
    north_gc_cos = 0;
    north_gc_sin = 1;

    if ((G_projection() != PROJECTION_LL)) {
	if (window.proj == 0)
	    G_fatal_error(_("Current projection is x,y (undefined)."));

	do_reproj = 1;

	/* read current projection info */
	if ((in_proj_info = G_get_projinfo()) == NULL)
	    G_fatal_error(_("Cannot get projection info of current location"));

	if ((in_unit_info = G_get_projunits()) == NULL)
	    G_fatal_error(_("Cannot get projection units of current location"));

	if (pj_get_kv(&iproj, in_proj_info, in_unit_info) < 0)
	    G_fatal_error(_("Cannot get projection key values of current location"));

	G_free_key_value(in_proj_info);
	G_free_key_value(in_unit_info);

	/*  output projection to lat/long w/ same ellipsoid as input */
	oproj.zone = 0;
	oproj.meters = 1.;
	sprintf(oproj.proj, "ll");
	if ((oproj.pj = pj_latlong_from_proj(iproj.pj)) == NULL)
	    G_fatal_error(_("Unable to update lat/long projection parameters"));
    }

    /* always init pd */
    S_init(&pd);
    pd.function = S_GEOM;
    pd.function = S_ZENETR;
    if (azimuth_name)
      pd.function = S_SOLAZM;
    if (sunhour_name || sretr_name || ssetr_name || hrang_name )
      pd.function |= S_SRSS;

    if (month == -1)
	doy = day;
    else
	doy = dom2doy2(year, month, day);
    
    set_solpos_time(&pd, year, 1, doy, hour, minutes, seconds);
    set_solpos_longitude(&pd, 0);
    pd.latitude = 0;
    S_solpos(&pd);

    if (elev_name) {
	if ((elev_fd = Rast_open_new(elev_name, FCELL_TYPE)) < 0)
	    G_fatal_error(_("Unable to create raster map <%s>"), elev_name);

	elevbuf = Rast_allocate_f_buf();
    }
    else {
	elevbuf = NULL;
	elev_fd = -1;
    }

    if (azimuth_name) {
	if ((azimuth_fd = Rast_open_new(azimuth_name, FCELL_TYPE)) < 0)
	    G_fatal_error(_("Unable to create raster map <%s>"), azimuth_name);

	azimuthbuf = Rast_allocate_f_buf();
    }
    else {
	azimuthbuf = NULL;
	azimuth_fd = -1;
    }

    if (sretr_name) {
	if ((sretr_fd = Rast_open_new(sretr_name, FCELL_TYPE)) < 0)
	    G_fatal_error(_("Unable to create raster map <%s>"), sretr_name);

	sretrbuf = Rast_allocate_f_buf();
    }
    else {
	sretrbuf = NULL;
	sretr_fd = -1;
    }

    if (ssetr_name) {
	if ((ssetr_fd = Rast_open_new(ssetr_name, FCELL_TYPE)) < 0)
	    G_fatal_error(_("Unable to create raster map <%s>"), ssetr_name);

	ssetrbuf = Rast_allocate_f_buf();
    }
    else {
	ssetrbuf = NULL;
	ssetr_fd = -1;
    }

    if (hrang_name) {
	if ((hrang_fd = Rast_open_new(hrang_name, FCELL_TYPE)) < 0)
	    G_fatal_error(_("Unable to create raster map <%s>"), hrang_name);

	hrangbuf = Rast_allocate_f_buf();
    }
    else {
	hrangbuf = NULL;
	hrang_fd = -1;
    }

    if (sunhour_name) {
      if ((sunhour_fd = Rast_open_new(sunhour_name, FCELL_TYPE)) < 0)
	    G_fatal_error(_("Unable to create raster map <%s>"), sunhour_name);

	sunhourbuf = Rast_allocate_f_buf();
    }
    else {
	sunhourbuf = NULL;
	sunhour_fd = -1;
    }

    if (elev_name && azimuth_name) {
	G_message(_("Calculating solar elevation and azimuth..."));
    }
    else if (elev_name) {
	G_message(_("Calculating solar elevation..."));
    }
    else if (azimuth_name) {
	G_message(_("Calculating solar azimuth..."));
    }

    nrows = Rast_window_rows();
    ncols = Rast_window_cols();

    ba2 = 6356752.3142 / 6378137.0;
    ba2 = ba2 * ba2;

    for (row = 0; row < nrows; row++) {

	G_percent(row, nrows, 2);
	
	/* get cell center northing */
	north = window.north - (row + 0.5) * window.ns_res;
	north_ll = north;

	for (col = 0; col < ncols; col++) {
	    long int retval;

	    /* get cell center easting */
	    east = window.west + (col + 0.5) * window.ew_res;
	    east_ll = east;

	    if (do_reproj) {
		north_ll = north;

		if (pj_do_proj(&east_ll, &north_ll, &iproj, &oproj) < 0)
		    G_fatal_error(_("Error in pj_do_proj (projection of input coordinate pair)"));
	    }

	    /* geocentric latitude */
	    north_gc = atan(ba2 * tan(DEG2RAD * north_ll));
	    north_gc_sin = sin(north_gc);
	    roundoff(&north_gc_sin);
	    north_gc_cos = cos(north_gc);
	    roundoff(&north_gc_cos);

	    set_solpos_longitude(&pd, east_ll);
	    pd.latitude = north_gc * RAD2DEG;
	    retval = S_solpos(&pd);
	    S_decode(retval, &pd);
	    G_debug(3, "solpos hour angle: %.5f", pd.hrang);

	    if (elev_name)
		elevbuf[col] = pd.elevetr;

	    if (azimuth_name) {
	      azimuthbuf[col] = pd.azim; 
	    }

	    if (sretr_name) {
	      sretrbuf[col] = pd.sretr;
	    }

	    if (ssetr_name) {
	      ssetrbuf[col] = pd.ssetr;
	    }

	    if (hrang_name) {
	      ssetrbuf[col] = pd.hrang;
	    }

	    if (sunhour_name) {
		sunhourbuf[col] = (pd.ssetr - pd.sretr) / 60.;
		if (sunhourbuf[col] > 24.)
		    sunhourbuf[col] = 24.;
		if (sunhourbuf[col] < 0.)
		    sunhourbuf[col] = 0.;
	    }
	}
	if (elev_name)
	    Rast_put_f_row(elev_fd, elevbuf);
	if (azimuth_name)
	    Rast_put_f_row(azimuth_fd, azimuthbuf);
	if (sretr_name)
	    Rast_put_f_row(sretr_fd, sretrbuf);
	if (ssetr_name)
	    Rast_put_f_row(ssetr_fd, ssetrbuf);
	if (hrang_name)
	    Rast_put_f_row(hrang_fd, hrangbuf);
	if (sunhour_name)
	    Rast_put_f_row(sunhour_fd, sunhourbuf);
    }
    G_percent(1, 1, 2);

    if (elev_name) {
	Rast_close(elev_fd);
	/* writing history file */
	Rast_short_history(elev_name, "raster", &hist);
	Rast_command_history(&hist);
	Rast_write_history(elev_name, &hist);
    }
    if (azimuth_name) {
	Rast_close(azimuth_fd);
	/* writing history file */
	Rast_short_history(azimuth_name, "raster", &hist);
	Rast_command_history(&hist);
	Rast_write_history(azimuth_name, &hist);
    }
    if (sretr_name) {
	Rast_close(sretr_fd);
	/* writing history file */
	Rast_short_history(sretr_name, "raster", &hist);
	Rast_command_history(&hist);
	Rast_write_history(sretr_name, &hist);
    }
    if (ssetr_name) {
	Rast_close(ssetr_fd);
	/* writing history file */
	Rast_short_history(ssetr_name, "raster", &hist);
	Rast_command_history(&hist);
	Rast_write_history(ssetr_name, &hist);
    }
    if (hrang_name) {
	Rast_close(hrang_fd);
	/* writing history file */
	Rast_short_history(hrang_name, "raster", &hist);
	Rast_command_history(&hist);
	Rast_write_history(hrang_name, &hist);
    }
    if (sunhour_name) {
	Rast_close(sunhour_fd);
	/* writing history file */
	Rast_short_history(sunhour_name, "raster", &hist);
	Rast_command_history(&hist);
	Rast_write_history(sunhour_name, &hist);
    }

    G_done_msg(" ");

    exit(EXIT_SUCCESS);
}

void set_solpos_time(struct posdata *pdat, int year, int month, int day,
                    int hour, int minute, int second)
{
    pdat->year = year; 
    pdat->month = month; 
    pdat->day = day; 
    pdat->daynum = day; 
    pdat->hour = hour; 
    pdat->minute = minute; 
    pdat->second = second;
    pdat->timezone = 0;

    pdat->time_updated = 1; 
    pdat->longitude_updated = 1;
}

void set_solpos_longitude(struct posdata *pdat, double longitude)
{
    pdat->longitude = longitude;

    pdat->longitude_updated = 1;
}

int roundoff(double *x)
{
    /* watch out for the roundoff errors */
    if (fabs(*x) > 1.0) {
	if (*x > 0.0)
	    *x = 1.0;
	else
	    *x = -1.0;

	return 1;
    }

    return 0;
}

