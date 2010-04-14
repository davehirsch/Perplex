/*
	graphic_output.c

	A library of subroutines to generate PDF/PS/SVG files and draw graphics primitives.
	Designed to be a drop-in substitute for Perplex's pslib.f postscript generation routines.
	(Hence all the psxxxx names instead of more reasonable names).  Uses Cairo to do the
	heavy lifting.

	Created by David Hirsch on 1/30/10.
	Copyright 2010 Western Washington University. All rights reserved.

	Implementation notes:
	Coordinates: Cairo includes its own transformation matrix.  I use this instead of
	the pslib transformation and scaling features.  The user space coordinate system
	is a 100.0 x 100.0 square.  This gets translated to be centered within the page
	and scaled to allow a 1-inch margin on all sides.
	
	Version History:
	
	0.8 - Beta testing release
	0.9 - 26 March 2010
		- Added clearing of rotation before line drawing

	0.91 - 28 March 2010
		- Fixed problem with spline function.  Wasn't calling getControlPoints correctly.
		- Removed code from psublk_. Trimming now gets done in pstext (psublk wasn't actually doing anything anyway)
		- Added code to adjust rotation of text based upon plot aspect ratio
	
	0.92 - 28 March 2010
		- Added function to get single plot option value that may have spaces in it, in order
			to be able to specify font names that contain spaces
 		
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "cairo.h"
#include "cairo-pdf.h"
#include "cairo-svg.h"
#include "cairo-ps.h"
#include "graphic_output.h"

/*==================================================================================
	GLOBALS - required for working correctly with calling Perplex routines.  (i.e., 
	we cannot change the setup to return and keep track of any variables).  Prefixed
	with 'dmh_' to avoid name conflicts (hopefully).
 ==================================================================================*/

cairo_surface_t *dmh_surf;
cairo_t *dmh_cr;
cairo_matrix_t dmh_unrotatedMatrix;
char dmh_inRotatedTransform = 0;
char dmh_debug = 0;
char dmh_track_min_x = 0;
double dmh_min_tracked_x = DBL_MAX;
double dmh_xscale = 1;	/* these are set in psssc2 or psssc1 */
double dmh_yscale = 1;
double dmh_xoffset = 1;
double dmh_yoffset = 1;
double dmh_aspectRatio = 1;
double dmh_inRelProcess = 0;	/* boolean describing whether we have begun a set of relative operations and not yet ended them.  This is required, because we cannot stroke each relative-move line separately if they are to be followed by another relative operation.  We only know when we are done with a set of relative operations by beginning an absolute operation. */
int dmh_pageWidth;
int dmh_pageHeight;
char *dmh_fontFace;
char *dmh_fileNameRoot;
int dmh_outputFileType;

/*==================================================================================
	MACROS
==================================================================================*/
#define DEBUGPRINT(x) if (dmh_debug) {printf x;};
#define MYMAX(x,y) ((x) > (y)) ? (x) : (y)

/*==================================================================================
	DEFINED PARAMETERS 
==================================================================================*/
#define TOPMARGIN 72			/* 1 inch * 72 pts/inch */
#define BOTTOMMARGIN 72			/* 1 inch * 72 pts/inch */
#define LEFTMARGIN 72			/* 1 inch * 72 pts/inch */
#define LEFTAXISMARGIN 72		/* Space between margin and vertical axis: 1 inch * 72 pts/inch */
#define RIGHTMARGIN 72			/* 1 inch * 72 pts/inch */
#define EXTRASPACEPTS 4			/* Extra space between left/right/top/bottom-aligned text and reference point */
#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif
#define MINCIRCLESIZE 2
#define REFERENCE_TO_TRACKED_MIN_X 32700
#define MAXPLOTSPERPAGE 6
#define MULTIPLOTGUTTER 24

/*==================================================================================
	psopen - opens a document for writing
==================================================================================*/
void psopen_(char *fname, int fnamelen) {
	
	char *outFileName;
	char outputType[255];
	char pageWidthString[255];
	char pageHeightString[255];
	
	/* Set debug status based upon presence of file named 'debug_yes' in directory */
	FILE *debugFile = fopen("debug_yes", "r");
	if (debugFile == NULL) {
		dmh_debug = 0;
	} else {
		dmh_debug = 1;
		fclose(debugFile);
	}

	
	fname[fnamelen]='\0';
	fname = trim(fname);
	outFileName = malloc((strlen(fname) + 50) * sizeof(char));
	strcpy(outFileName, fname);
	DEBUGPRINT(("Found file name:%s of length: %lu\n", fname, strlen(fname)));
	dmh_fileNameRoot = malloc(255*sizeof(char));
	strcpy(dmh_fileNameRoot, fname);
	
	/* Look for plot options file to specify what sort of output we want: SVG/PS/PDF
	 Could use the Perplex routine to do this, but calling into C from
	 fortran is bad enough. */
	
	/* Default Values - If these get changed, then need to change the display in
		pscom_new.f as well */
	if (!getOptionForKey("plot_output_type", outputType, 255)) {
		strcpy(outputType, "PDF");
	}
	
	if (get2OptionsForKey("page_size", pageWidthString, pageHeightString, 255)) {
		DEBUGPRINT(("Found page size strings: (%s, %s)\n", pageWidthString, pageHeightString));
		sscanf(pageWidthString, "%i", &dmh_pageWidth);
		sscanf(pageHeightString, "%i", &dmh_pageHeight);
		DEBUGPRINT(("Found page size values: (%i, %i)\n", dmh_pageWidth, dmh_pageHeight));
	} else {
		dmh_pageWidth = 612;
		dmh_pageHeight = 792;
	}
	
	dmh_fontFace = malloc(255*sizeof(char));
	if (!getCompleteOptionForKey("new_font", dmh_fontFace, 255)) {
		strcpy(dmh_fontFace, "Arial");
	}
	
	
	if (!strcmp(outputType, "SVG") || !strcmp(outputType, "svg")) {
		dmh_outputFileType = SVGTYPE;
		strcat(outFileName, ".svg");
		dmh_surf = cairo_svg_surface_create (outFileName, dmh_pageWidth, dmh_pageHeight);
	} else if (!strcmp(outputType, "PS") || !strcmp(outputType, "ps")) {
		dmh_outputFileType = PSTYPE;
		strcat(outFileName, ".ps");
		dmh_surf = cairo_ps_surface_create (outFileName, dmh_pageWidth, dmh_pageHeight);
		cairo_ps_surface_set_eps (dmh_surf, 1);
	} else {
		/* PDF is the default if nothing is specified in the plot options file */
		dmh_outputFileType = PDFTYPE;
		strcat(outFileName, ".pdf");
		dmh_surf = cairo_pdf_surface_create (outFileName, dmh_pageWidth, dmh_pageHeight);
	}
	
	dmh_cr = cairo_create (dmh_surf);
	cairo_identity_matrix(dmh_cr);
	cairo_set_line_join(dmh_cr, CAIRO_LINE_JOIN_ROUND);
	
	dmh_min_tracked_x = DBL_MAX;
	
	free(outFileName);
}

/*==================================================================================
 psclos - closes a pdf document after writing
 ==================================================================================*/
void psclos_() {
	
	DEBUGPRINT(("In psclos.\n"));
/* DEBUGGING TESTS:
	
	cairo_identity_matrix(dmh_cr);	// remove transformations
	cairo_move_to(dmh_cr, 100, 200);
	cairo_rectangle(dmh_cr, 100, 200, 10, 10);
	cairo_set_source_rgb (dmh_cr, 0, 1, 1);
	cairo_fill(dmh_cr);

	cairo_rotate(dmh_cr, -90.0*M_PI/180);
	double x=100, y=200;
	cairo_device_to_user(dmh_cr, &x, &y);
	cairo_move_to(dmh_cr, x, y);
	
	cairo_text_extents_t te;
	cairo_set_font_size(dmh_cr, 44);
	cairo_text_extents(dmh_cr,"test",&te);
	cairo_rel_move_to(dmh_cr, -te.width * 0.5, 0);
	cairo_set_source_rgb (dmh_cr, 0, 1, 0);
	cairo_show_text (dmh_cr, "test");
	DEBUGPRINT(("Status:%s\n",cairo_status_to_string(cairo_status(dmh_cr))));
*/
	
	closeSurface();
	
	free(dmh_fontFace);
	free(dmh_fileNameRoot);
}


/*==================================================================================
 closeSurface - cleans up and finalizes the Cairo surface.  Called from psclos or
 pspltrgn.
 ==================================================================================*/
void closeSurface() {
	cairo_surface_flush(dmh_surf);
	cairo_surface_finish(dmh_surf);
	cairo_destroy (dmh_cr);
	cairo_surface_destroy (dmh_surf);
}

/*==================================================================================
 completeRelativeOperation - checks to see if we have begun a relative operation that
 therefore needs to be completed and stroked.  Handling things this way means that the
 line and fill settings for the last relative operation call will take precedence over
 earlier ones.  It has the advantage of keeping all the elements (e.g., tick marks) as
 a single group.
 ==================================================================================*/
void completeRelativeOperation () {
	if (dmh_inRelProcess) {
		DEBUGPRINT(("Completing relative operations.\n"));
		cairo_stroke_preserve(dmh_cr);
		cairo_fill(dmh_cr);
		dmh_inRelProcess = 0;
	}
}

/*==================================================================================
 getOptionForKey - examines the plot options file for the value of a keyword. Returns
 the value in the supplied string variable, which should be pre-allocated.
 The value may have spaces in it.
 ==================================================================================*/
char getCompleteOptionForKey(char *keyword, char *value, int valueSize) {
	char optionFileName[655] = "perplex_plot_option.dat";
	char *oneLine = malloc(655 * sizeof(char));
	char currentKeyword[655];
	char currentValue[655];
	FILE *optionFile;
	char foundValue = 0;
	int argsAssigned;
	
	optionFile = fopen(optionFileName, "r");
	
	if (optionFile) {
		while (!feof(optionFile) && !foundValue) {
			oneLine = fgets(oneLine, 655, optionFile);
			
			argsAssigned = sscanf(oneLine, "%654[^ 	|] %654[^|]", currentKeyword, currentValue);
			
			if (!strcmp(currentKeyword, keyword)) {
				foundValue = 1;
				trim(currentValue);
				DEBUGPRINT(("Found %i values: '%s' for keyword '%s'.\n", argsAssigned, currentValue, currentKeyword));
				trim(currentValue);
				strncpy(value, currentValue, valueSize-1);	/* copy only as many characters as will fit in value */
			}
		}
		fclose(optionFile);
	}
	free(oneLine);
	return foundValue;
}

/*==================================================================================
 getOptionForKey - examines the plot options file for the value of a keyword. Returns
 the value in the supplied string variable, which should be pre-allocated.
 ==================================================================================*/
char getOptionForKey(char *keyword, char *value, int valueSize) {
	char optionFileName[655] = "perplex_plot_option.dat";
	char *oneLine = malloc(655 * sizeof(char));
	char currentKeyword[655];
	char currentValue[655];
	FILE *optionFile;
	char foundValue = 0;
	int argsAssigned;
	
	optionFile = fopen(optionFileName, "r");
	
	if (optionFile) {
		
		while (!feof(optionFile) && !foundValue) {
			oneLine = fgets(oneLine, 655, optionFile);
			
			argsAssigned = sscanf(oneLine, "%654[^ 	|] %654[^ 	|]", currentKeyword, currentValue);
			
			if (!strcmp(currentKeyword, keyword)) {
				foundValue = 1;
				DEBUGPRINT(("Found %i values: '%s' for keyword '%s'.\n", argsAssigned, currentValue, currentKeyword));
				strncpy(value, currentValue, valueSize-1);	/* copy only as many characters as will fit in value */
			}
		}
		fclose(optionFile);
	}
	free(oneLine);
	return foundValue;
}

/*==================================================================================
 get2OptionsForKey - examines the plot options file for the value of a keyword. Returns
 the value in the supplied string variable, which should be pre-allocated.
 (2 parameters)
 ==================================================================================*/
char get2OptionsForKey(char *keyword, char *value1, char *value2, int valueSize) {
	char optionFileName[655] = "perplex_plot_option.dat";
	char *oneLine = malloc(655 * sizeof(char));
	char currentKeyword[655];
	char currentValue1[655];
	char currentValue2[655];
	FILE *optionFile;
	char foundValue = 0;
	int argsAssigned;
	
	optionFile = fopen(optionFileName, "r");
	
	if (optionFile) {
		while (!feof(optionFile) && !foundValue) {
			oneLine = fgets(oneLine, 655, optionFile);
			
			argsAssigned = sscanf(oneLine, "%654[^ 	|] %654[^ 	|] %654[^ 	|]", currentKeyword, currentValue1, currentValue2);
			
			if (!strcmp(currentKeyword, keyword)) {
				foundValue = 1;
				DEBUGPRINT(("Found %i values: '%s', '%s' for keyword '%s'.\n", argsAssigned, currentValue1, currentValue2, currentKeyword));
				strncpy(value1, currentValue1, valueSize-1);	/* copy only as many characters as will fit in value */
				strncpy(value2, currentValue2, valueSize-1);	/* copy only as many characters as will fit in value */
			}
		}
		fclose(optionFile);
	}
	free(oneLine);
	return foundValue;
}


/*==================================================================================
 pselip - generates an ellipse
	if dx or dy = 0, then they are set to minimum values in device space
 ==================================================================================*/
void pselip_ (double *xor, double *yor, double *dx, double *dy, double *rline, double *width, int *ifill) {
	double devX = deviceX(*xor);
	double devY = deviceY(*yor);
	double devW = deviceW(*dx);
	double devH = deviceH(*dy);
	
	DEBUGPRINT(("In pselip.  Origin: (%f, %f), Size: (%f, %f), DevOrigin: (%f, %f), DevSize: (%f, %f).\n", 
				*xor, *yor, *dx, *dy, devX, devY, devW, devH));
	completeRelativeOperation();
	cairo_new_path(dmh_cr);
	cairo_save (dmh_cr);
	
	if (devW == 0 || devH == 0) {
		devW = MINCIRCLESIZE;
		devH = MINCIRCLESIZE;
	}

	cairo_translate (dmh_cr, devX + devW/ 2., devY + devH / 2.);
	cairo_scale (dmh_cr, devW / 2., devH / 2.);
	cairo_arc (dmh_cr, 0., 0., 1., 0., 2 * M_PI);
	cairo_restore (dmh_cr);
	
	setLineProperties((int) *rline, *width);
	cairo_set_source_rgb (dmh_cr, 0, 0, 0);	/* black */
	cairo_stroke_preserve(dmh_cr);
	setFillType (*ifill);
	cairo_fill(dmh_cr);
}

/*==================================================================================
 getControlPoints - a routine to calculate control points based on the vertices
	in a Bezier spline.
		From http://www.antigrain.com/research/bezier_interpolation/
		Assume we need to calculate the control
		points between (x1,y1) and (x2,y2).
		Then x0,y0 - the previous vertex,
		x3,y3 - the next one. 
 ==================================================================================*/
void getControlPoints(double x0, double y0,
					  double x1, double y1,
					  double x2, double y2,
					  double x3, double y3,
					  double *ctrl1_x, double *ctrl1_y,
					  double *ctrl2_x, double *ctrl2_y) {
    double xc1 = (x0 + x1) / 2.0;
    double yc1 = (y0 + y1) / 2.0;
    double xc2 = (x1 + x2) / 2.0;
    double yc2 = (y1 + y2) / 2.0;
    double xc3 = (x2 + x3) / 2.0;
    double yc3 = (y2 + y3) / 2.0;

    double len1 = sqrt((x1-x0) * (x1-x0) + (y1-y0) * (y1-y0));
    double len2 = sqrt((x2-x1) * (x2-x1) + (y2-y1) * (y2-y1));
    double len3 = sqrt((x3-x2) * (x3-x2) + (y3-y2) * (y3-y2));

    double k1 = len1 / (len1 + len2);
    double k2 = len2 / (len2 + len3);

    double xm1 = xc1 + (xc2 - xc1) * k1;
    double ym1 = yc1 + (yc2 - yc1) * k1;

    double xm2 = xc2 + (xc3 - xc2) * k2;
    double ym2 = yc2 + (yc3 - yc2) * k2;

	/* Resulting control points. Here smooth_value is mentioned
	 above coefficient K whose value should be in range [0...1]. */
	double smooth_value = 0.8;
	*ctrl1_x = xm1 + (xc2 - xm1) * smooth_value + x1 - xm1;
	*ctrl1_y = ym1 + (yc2 - ym1) * smooth_value + y1 - ym1;
	
	*ctrl2_x = xm2 + (xc2 - xm2) * smooth_value + x2 - xm2;
	*ctrl2_y = ym2 + (yc2 - ym2) * smooth_value + y2 - ym2;
}


/*==================================================================================
 psbspl - subroutine to generate open bsplines
	We receive only the vertices from perplex, so we must figure out good control
	points to use for the spline.
 ==================================================================================*/
void psbspl_ (double *x, double *y, int *npts, double *rline, double *width, int *ifill) {
	int i, N;
	double ctrl1_x, ctrl1_y, ctrl2_x, ctrl2_y;
	
	N = *npts;
	
	DEBUGPRINT(("In psbspl. With %i points, rline=%f, wdith=%f, fill=%i.\n", N, *rline, *width, *ifill));
	completeRelativeOperation();

	if (dmh_inRotatedTransform) {
		cairo_set_matrix(dmh_cr, &dmh_unrotatedMatrix);	/* remove rotation */
		dmh_inRotatedTransform = 0;
		DEBUGPRINT(("Removing rotation.\n"));
	}

	if (N >= 4) {
		cairo_move_to(dmh_cr, deviceX(x[0]), deviceY(y[0]));
		getControlPoints(deviceX(x[0]), deviceY(y[0]), deviceX(x[0]), deviceY(y[0]), 
						 deviceX(x[1]), deviceY(y[1]), deviceX(x[2]), deviceY(y[2]), 
						 &ctrl1_x, &ctrl1_y, &ctrl2_x, &ctrl2_y);
		cairo_curve_to(dmh_cr, ctrl1_x, ctrl1_y, ctrl2_x, ctrl2_y, deviceX(x[1]), deviceY(y[1]));
//		cairo_line_to(dmh_cr, deviceX(x[1]), deviceY(y[1]));
		for (i=1; i <= N-3; i++) {
			getControlPoints(deviceX(x[i-1]), deviceY(y[i-1]), deviceX(x[i]), deviceY(y[i]), 
							 deviceX(x[i+1]), deviceY(y[i+1]), deviceX(x[i+2]), deviceY(y[i+2]), 
							 &ctrl1_x, &ctrl1_y, &ctrl2_x, &ctrl2_y);
			cairo_curve_to(dmh_cr, ctrl1_x, ctrl1_y, ctrl2_x, ctrl2_y, deviceX(x[i+1]), deviceY(y[i+1]));
		}
		
		/* We've now drawn all the way to N-2.  We need to draw through N-1.
			We can get one more control point by working backwards from the end: */
		getControlPoints(deviceX(x[N-3]), deviceY(y[N-3]), deviceX(x[N-2]), deviceY(y[N-2]), 
						 deviceX(x[N-1]), deviceY(y[N-1]), deviceX(x[N-1]), deviceY(y[N-1]), 
						 &ctrl1_x, &ctrl1_y, &ctrl2_x, &ctrl2_y);
		cairo_curve_to(dmh_cr, ctrl1_x, ctrl1_y, ctrl2_x, ctrl2_y, deviceX(x[N-1]), deviceY(y[N-1]));
//		cairo_line_to(dmh_cr, deviceX(x[N-1]), deviceY(y[N-1]));
	} else if (N == 3) {
		/* for three points, we'll pretend for the purposes of control point calculation that
		 it's a closed loop */
		cairo_move_to(dmh_cr, deviceX(x[0]), deviceY(y[0]));
		getControlPoints(deviceX(x[2]), deviceY(y[2]), deviceX(x[0]), deviceY(y[0]), 
						 deviceX(x[1]), deviceY(y[1]), deviceX(x[2]), deviceY(y[2]), 
						 &ctrl1_x, &ctrl1_y, &ctrl2_x, &ctrl2_y);
		cairo_curve_to(dmh_cr, ctrl1_x, ctrl1_y, ctrl2_x, ctrl2_y, deviceX(x[1]), deviceY(y[1]));
		getControlPoints(deviceX(x[0]), deviceY(y[0]), deviceX(x[1]), deviceY(y[1]), 
						 deviceX(x[2]), deviceY(y[2]), deviceX(x[0]), deviceY(y[0]), 
						 &ctrl1_x, &ctrl1_y, &ctrl2_x, &ctrl2_y);
		cairo_curve_to(dmh_cr, ctrl1_x, ctrl1_y, ctrl2_x, ctrl2_y, deviceX(x[2]), deviceY(y[2]));
	} else if (N == 2) {
		/* for two points, it's a line segment.  Duh. */
		cairo_move_to(dmh_cr, deviceX(x[0]), deviceY(y[0]));
		cairo_line_to(dmh_cr, deviceX(x[1]), deviceY(y[1]));
	}
	setLineProperties((int) *rline, *width);
	cairo_set_source_rgb (dmh_cr, 0, 0, 0);	/* black */
	cairo_stroke_preserve(dmh_cr);
	setFillType (*ifill);
	cairo_fill(dmh_cr);
}
/*==================================================================================
 psrpgn - subroutine to generate closed polygons using rel. coordinates for all after
	the first
==================================================================================*/
void psrpgn_ (double *x1, double *y1, double *rx, double *ry, int *npts, double *rline, double *width,int *ifill) {

	int i;

	DEBUGPRINT(("In psrpgn\n"));
	completeRelativeOperation();
	cairo_move_to(dmh_cr, *x1, *y1);
	for (i = 0; i <= *npts-1; i++) {
		cairo_rel_line_to(dmh_cr, deviceW(rx[i]), deviceH(ry[i]));
	}
	cairo_close_path(dmh_cr);

	setLineProperties((int) *rline, *width);
	cairo_set_source_rgb (dmh_cr, 0, 0, 0);	/* black */
	cairo_stroke_preserve(dmh_cr);
	setFillType (*ifill);
	cairo_fill(dmh_cr);
}	

/*==================================================================================
 setLineProperties - subroutine to set line qualities used by cairo_stroke
 width - width of line in points (I think).
 dashtype  - line style indicator:
	0 - no line
	1 - solid line
	2 - short dashed line
	3 - uneven dashed line
	4 - very sparse dotted line
	5 - long even dashed line
	6 - dashed line
	7 - very short dashed line
	8 - sparse dotted
	9 - heavily dotted
	10 - ultra sparse
 ==================================================================================*/
void setLineProperties (int dashtype, double width) {
	double *dashes = malloc(4 * sizeof(double));
	
	if (width > 0) {
		cairo_set_line_width(dmh_cr, width);
	}
	
	switch (dashtype) {
		case 0:
			cairo_set_line_width(dmh_cr, 0);
			cairo_set_dash (dmh_cr, dashes, 0, 0);
			break;
		case 1:
			cairo_set_dash (dmh_cr, dashes, 0, 0);
			break;
		case 2:
			dashes[0] = 4.0;
			cairo_set_dash (dmh_cr, dashes, 1, 0);
			cairo_set_line_cap (dmh_cr, CAIRO_LINE_CAP_BUTT);
			break;
		case 3:
			dashes[0] = 7.0;	
			dashes[1] = 3.0;	
			dashes[2] = 3.0;	
			dashes[3] = 3.0;	
			cairo_set_dash (dmh_cr, dashes, 4, 0);
			cairo_set_line_cap (dmh_cr, CAIRO_LINE_CAP_ROUND);
			break;
		case 4:
			dashes[0] = 0.0;	/* Zero dash length = dot (if line cap type = round) */
			dashes[1] = 10.0;	
			cairo_set_dash (dmh_cr, dashes, 2, 0);
			cairo_set_line_cap (dmh_cr, CAIRO_LINE_CAP_ROUND);
			break;
		case 5:
			dashes[0] = 12.0;
			dashes[1] = 4.0;
			cairo_set_dash (dmh_cr, dashes, 2, 0);
			cairo_set_line_cap (dmh_cr, CAIRO_LINE_CAP_BUTT);
			break;
		case 6:
			dashes[0] = 8.0;
			cairo_set_dash (dmh_cr, dashes, 1, 0);
			cairo_set_line_cap (dmh_cr, CAIRO_LINE_CAP_BUTT);
			break;
		case 7:
			dashes[0] = 2.0;
			cairo_set_dash (dmh_cr, dashes, 1, 0);
			cairo_set_line_cap (dmh_cr, CAIRO_LINE_CAP_BUTT);
			break;
		case 8:
			dashes[0] = 0.0;	
			dashes[1] = 5.0;
			cairo_set_dash (dmh_cr, dashes, 2, 0);
			cairo_set_line_cap (dmh_cr, CAIRO_LINE_CAP_ROUND);
			break;
		case 9:
			dashes[0] = 0.0;
			dashes[1] = 2.0;
			cairo_set_dash (dmh_cr, dashes, 2, 0);
			cairo_set_line_cap (dmh_cr, CAIRO_LINE_CAP_ROUND);
			break;
		case 10:
			dashes[0] = 0.0;
			dashes[1] = 15.0;
			cairo_set_dash (dmh_cr, dashes, 2, 0);
			cairo_set_line_cap (dmh_cr, CAIRO_LINE_CAP_ROUND);
			break;
	}
}


/*==================================================================================
 setFillType - subroutine to set fill type, depending on the integer ifill:
	0 - transparent
	1 - white
	2 - v. light grey
	3 - light grey
	4 - med-light grey
	5 - med grey
	6 - dark grey
	7 - black
	(8-15 should be patterns, but for now, they will be colors.  Making a vector pattern
		is a bit challenging in Cairo)
	8 - light red
	9 - light blue
	10 - light green
	11 - light orange
	12 - yellow
	13 - light purple
	14 - light cyan
	15 - light brown
 ==================================================================================*/
void setFillType (int ifill) {
	
	switch (ifill) {
		case 0:
			cairo_set_source_rgba (dmh_cr, 0, 0, 0, 0);	/* transparent black */
			break;
		case 1:
			cairo_set_source_rgb (dmh_cr, 1, 1, 1);	/* white */
			break;
		case 2:
			cairo_set_source_rgb (dmh_cr, 0.95, 0.95, 0.95);	/* v. lt grey */
			break;
		case 3:
			cairo_set_source_rgb (dmh_cr, 0.85, 0.85, 0.85);	/* lt grey */
			break;
		case 4:
			cairo_set_source_rgb (dmh_cr, 0.75, 0.75, 0.75);	/* med-lt grey */
			break;
		case 5:
			cairo_set_source_rgb (dmh_cr, 0.5, 0.5, 0.5);	/* med grey */
			break;
		case 6:
			cairo_set_source_rgb (dmh_cr, 0.25, 0.25, 0.25);	/* dk grey */
			break;
		case 7:
			cairo_set_source_rgb (dmh_cr, 0, 0, 0);		/*  black */
			break;
		case 8:
			cairo_set_source_rgb (dmh_cr, 1, 0.6, 0.6);	/*  light red */
			break;
		case 9:
			cairo_set_source_rgb (dmh_cr, 0.6, 0.6, 1);	/*  light blue */
			break;
		case 10:
			cairo_set_source_rgb (dmh_cr, 0.6, 1, 0.6);	/*  light green */
			break;
		case 11:
			cairo_set_source_rgb (dmh_cr, 1, 0.5, 0.2);	/*  light orange */
			break;
		case 12:
			cairo_set_source_rgb (dmh_cr, 1, 1, 0.4);	/*  yellow */
			break;
		case 13:
			cairo_set_source_rgb (dmh_cr, 1, 0.5, 1);	/*  light purple */
			break;
		case 14:
			cairo_set_source_rgb (dmh_cr, 0.5, 1, 1);	/*  light cyan */
			break;
		case 15:
			cairo_set_source_rgb (dmh_cr, 0.7, 0.7, 0.4);	/*  light brown */
			break;
		default:
			break;
	}
}
/*      implicit none
 
 integer ifill 
 
 integer nps
 double precision xscale,yscale,xmn,ymn
 common/ scales /xscale,yscale,xmn,ymn,nps
 
 character*30 fill(15)
 
 save fill
 //
 //                               first seven fills vary in intensity
 //                               from white to black
 data fill/'1','0.95','0.85','0.75','0.5','0.25','0',
 //                               next 8 fills are patterned
 *          '< 11 22 44 88 11 22 44 88 > -1',
 *          '< 88 44 22 11 88 44 22 11 > -1',
 *          '< ff 00 00 00 ff 00 00 00 > -1',
 *          '< 88 88 88 88 88 88 88 88 > -1',
 *          '< ff 88 88 88 ff 88 88 88 > -1',
 *          '< 88 55 22 55 88 55 22 55 > -1',
 *          '< cc cc 33 33 cc cc 33 33 > -1',
 *          '< 77 bb ee dd 77 bb ee dd > -1'/
 
 if (ifill.eq.0) then
 write (nps,1000) 
 else if (ifill.le.15) then
 write (nps,1010) fill(ifill)
 else
 write (*,*) 'invalid fill choice'
 stop
 end if
 
 1000  format ('none SetP %I p n')
 1010  format ('%I p',/,a30,' SetP')
 
 end
 */


/*==================================================================================
 pspygn - subroutine to generate closed polygons, abs. coordinates
 ==================================================================================*/
void pspygn_ (double *x, double *y, int *npts, double *rline, double *width, int *ifill) {

	int i;
	
	DEBUGPRINT(("In pspygn.  Received %i points.  First=(%f, %f); Dev=(%f,%f). rline=%f, width=%f, fill=%i\n", *npts, x[0], y[0], deviceX(x[0]), deviceY(y[0]), *rline, *width, *ifill));
	completeRelativeOperation();
	cairo_move_to(dmh_cr, deviceX(x[0]), deviceY(y[0]));
	for (i = 1; i <= *npts-1; i++) {
		cairo_line_to(dmh_cr, deviceX(x[i]), deviceY(y[i]));
	}
	cairo_close_path(dmh_cr);
	
	setLineProperties((int) *rline, *width);
	cairo_set_source_rgb (dmh_cr, 0, 0, 0);	/* black */
	cairo_stroke_preserve(dmh_cr);
	setFillType (*ifill);
	cairo_fill(dmh_cr);
}

/*==================================================================================
 pspyln - subroutine to generate open polylines
 ==================================================================================*/
void pspyln_ (double *x, double *y, int *npts, double *rline, double *width, int *ifill) {
	int i;
	
	DEBUGPRINT(("In pspyln. Received %i points, with rline=%f, width=%f, fill=%i\n", *npts, *rline, *width, *ifill));
	completeRelativeOperation();

	if (dmh_inRotatedTransform) {
		cairo_set_matrix(dmh_cr, &dmh_unrotatedMatrix);	/* remove rotation */
		dmh_inRotatedTransform = 0;
		DEBUGPRINT(("Removing rotation.\n"));
	}

	cairo_move_to(dmh_cr, deviceX(x[0]), deviceY(y[0]));
	for (i = 1; i <= *npts-1; i++) {
		DEBUGPRINT(("In pspyln, lineto: (%f, %f).\n", deviceX(x[i]), deviceY(y[i])));
		cairo_line_to(dmh_cr, deviceX(x[i]), deviceY(y[i]));
	}
	
	setLineProperties((int) *rline, *width);
	cairo_set_source_rgb (dmh_cr, 0, 0, 0);	/* black */
	cairo_stroke_preserve(dmh_cr);
	setFillType (*ifill);
	cairo_fill(dmh_cr);
}

/*==================================================================================
 psstrn - subroutine to set transformation matrix
 	Obsolete, but easier to make it a stub than to remove it altogether.
 ==================================================================================*/
void psstrn_ (double *xs, double *ys, double *xt, double *yt, double *theta) {
	DEBUGPRINT(("In psstrn\n"));
}

/*==================================================================================
 pspltrgn - subroutine to set which small plot region to use
 	This routine is called by pschem in psvdraw_new in order to plot multiple
 	ternary chemographies on a single page.  Each page can hold MAXPLOTSPERPAGE.  If
 	another is requested, then a new document will be created with a name similar to 
 	that of the first, and the process will begin again.
 	
 	This overrides any settings made in psssc2, and assumes that x and y (real-unit) 
 	bounds are 0.0-1.0.  Note that if this routine is called, then plot_aspect_ratio will be ignored.  
 	
 	Note that plotnum is a zero-based index
 ==================================================================================*/
void pspltrgn_ (int *plotnum) {
	char *outFileName = malloc((strlen(dmh_fileNameRoot)+50) * sizeof(char));
			
	int plotPosition = *plotnum % MAXPLOTSPERPAGE;
	int pageNum = *plotnum / MAXPLOTSPERPAGE;
	int boxEdge, rowNum, colNum;
	
	DEBUGPRINT(("In pspltrgn. Plotnum = %i, plotPosition = %i, pageNum = %i\n", *plotnum, plotPosition, pageNum));
	
	if (plotPosition == 0 && pageNum > 0) {
		/* we need to start a new page, so let's close the current one, and start a new one
			with a related name */
		DEBUGPRINT(("In pspltrgn. Closing current page and starting a new one with\n page number=%i and file type=%i\n", pageNum, dmh_outputFileType));
		closeSurface();	// close existing surface
		switch(dmh_outputFileType) {
			case PDFTYPE:
				sprintf(outFileName, "%s_%i.%s", dmh_fileNameRoot, pageNum, "pdf");
				dmh_surf = cairo_pdf_surface_create (outFileName, dmh_pageWidth, dmh_pageHeight);
				break;
			case PSTYPE:
				sprintf(outFileName, "%s_%i.%s", dmh_fileNameRoot, pageNum, "ps");
				dmh_surf = cairo_ps_surface_create (outFileName, dmh_pageWidth, dmh_pageHeight);
				cairo_ps_surface_set_eps (dmh_surf, 1);
				break;
			case SVGTYPE:
				sprintf(outFileName, "%s_%i.%s", dmh_fileNameRoot, pageNum, "svg");
				dmh_surf = cairo_svg_surface_create (outFileName, dmh_pageWidth, dmh_pageHeight);
				break;
		}
		
		dmh_cr = cairo_create (dmh_surf);
		cairo_identity_matrix(dmh_cr);
		cairo_set_line_join(dmh_cr, CAIRO_LINE_JOIN_ROUND);
		
		dmh_min_tracked_x = DBL_MAX;
	}
	
	/* Set the location on the page for the small plot */
	dmh_aspectRatio = 1.0;	/* Ignores plot_aspect_ratio.  */
	boxEdge = (dmh_pageWidth - LEFTMARGIN - RIGHTMARGIN - MULTIPLOTGUTTER) / 2;
	rowNum = plotPosition / 2;
	colNum = plotPosition % 2;
	dmh_xoffset = LEFTMARGIN + ((boxEdge + MULTIPLOTGUTTER) * colNum);
	dmh_yoffset = (dmh_pageHeight * 0.5) + (boxEdge * 1.5) + MULTIPLOTGUTTER - (boxEdge * (rowNum+1)) - (MULTIPLOTGUTTER * rowNum); 
	dmh_xscale = boxEdge;
	dmh_yscale = boxEdge;
	DEBUGPRINT(("End pspltrgn.  boxEdge=%i; r,c=(%i,%i); scale=(%f, %f); offset=(%f, %f); DevPtRange=(%f,%f)-(%f,%f).\n", boxEdge, rowNum, colNum, dmh_xscale, dmh_yscale, dmh_xoffset, dmh_yoffset,deviceX(0), deviceY(0), deviceX(1), deviceY(1)));
}


/*==================================================================================
 psssc2 - subroutine to set scaling of axes (independently)
	This routine is called once per drawing session to tell the graphic system
	what the bounds in real units are for x and y.  All other calls provide x and y
	in real units (e.g., bars, Kelvin, mole fraction, etc.).
 
	Here we set the scaling and offset factors to translate these units into page
	coordinates.  We cannot use the cairo transformation calls, because the line
	widths get screwed up royally if we do that. We dictate that the graph occupies
	a square region on the page, centered vertically and right-justified, leaving 1 inch
	between the graph edge and the right page edge, and 2 inches between the graph edge (of
	the box) and the left page edge.  This should allow enough space for the final x tick label
	on the right, and the vertical axis labels on the left.
 ==================================================================================*/
void psssc2_ (double *xmin, double *xmax, double *ymin, double *ymax) {
	double deviceAxisX, deviceAxisY;
	double scaledXmin, scaledYmin;
	double deviceOriginX, deviceOriginY;
	char plotAspectRatioString[255];
	float plotAspectRatio;	/* x axis/y axis */
	
	DEBUGPRINT(("In psssc2.  min=(%f, %f); max=(%f, %f).\n", *xmin, *ymin, *xmax, *ymax));
	if (!getOptionForKey("plot_aspect_ratio", plotAspectRatioString, 255)) {
		plotAspectRatio = 1.0;
	}
	sscanf(plotAspectRatioString, "%f", &plotAspectRatio);
	DEBUGPRINT(("In psssc2.  plotAspectRatioString=%s, and plotAspectRatio=%f.\n", plotAspectRatioString, plotAspectRatio));
	
	dmh_aspectRatio = plotAspectRatio;	/* save this value for use later */
	deviceAxisX = dmh_pageWidth - LEFTMARGIN - LEFTAXISMARGIN - RIGHTMARGIN;
	deviceAxisY = deviceAxisX / plotAspectRatio;
	deviceOriginX = LEFTMARGIN + LEFTAXISMARGIN;
	deviceOriginY = (dmh_pageHeight * 0.5) - (deviceAxisY * 0.5); 
	dmh_xscale = deviceAxisX / (*xmax - *xmin);
	dmh_yscale = deviceAxisY / (*ymax - *ymin);
	scaledXmin = dmh_xscale * *xmin;
	scaledYmin = dmh_yscale * *ymin;
	dmh_xoffset = deviceOriginX - scaledXmin;
	dmh_yoffset = deviceOriginY - scaledYmin;
	DEBUGPRINT(("End psssc2.  scale=(%f, %f); offset=(%f, %f).\n", dmh_xscale, dmh_yscale, dmh_xoffset, dmh_yoffset));
}
/*==================================================================================
 psrect - subroutine to output a rectangle, with integer fill code
 ==================================================================================*/
void psrect_ (double *x1, double *x2, double *y1, double *y2, double *rline, double *width, int *ifill) {
	DEBUGPRINT(("In psrect. minPt(1)=(%f, %f); maxPt(2)=(%f, %f), rline=%f, width=%f, fill=%i\n", *x1, *y1, *x2, *y2, *rline, *width, *ifill));
	DEBUGPRINT(("In psrect. Device Coords: minPt=(%f, %f); Size=(%f, %f)\n", deviceX(*x1), 
				deviceY(*y2),
				deviceW(*x2 - *x1), 
				deviceH(fabs(*y2 - *y1))));
	completeRelativeOperation();
	cairo_rectangle(dmh_cr, deviceX(*x1), 
					deviceY(*y2),	/* this is y2, not y1, to deal with Cairo's flipped coordinate system */
					deviceW(*x2 - *x1),
					deviceH(*y2 - *y1));
	setLineProperties((int) *rline, *width);
	cairo_set_source_rgb (dmh_cr, 0, 0, 0);	/* black */
	cairo_stroke_preserve(dmh_cr);
	setFillType (*ifill);
	cairo_fill(dmh_cr);
}

/*==================================================================================
 device_ - these four routines translate X, Y, Width, and Height values from 
	real (P, T, X) coordinates into device coordinates (points).  Note that because
	cairo's coordinate system is flipped (origin at upper left, not lower left), 
	deviceY is a bit more complex.
 ==================================================================================*/

double deviceX (double inX) {
	return (inX * dmh_xscale) + dmh_xoffset;
}

double deviceY (double inY) {
	return dmh_pageHeight - ((inY * dmh_yscale) + dmh_yoffset);
}

double deviceW (double inWidth) {
	return inWidth * dmh_xscale;
}

double deviceH (double inHeight) {
	return inHeight * dmh_yscale;
}

/*==================================================================================
 psrecr - subroutine to output a rectangle, with real fill
 ==================================================================================*/
void psrecr_ (double *x1, double *x2, double *y1, double *y2, double *rline, double *width, double *rfill) {
	DEBUGPRINT(("In psrecr. minPt(1)=(%f, %f); maxPt(2)=(%f, %f), rline=%f, width=%f, fill=%f\n", *x1, *y1, *x2, *y2, *rline, *width, *rfill));
	DEBUGPRINT(("In psrecr. Device Coords: minPt=(%f, %f); Size=(%f, %f)\n", deviceX(*x1), 
				deviceY(*y2),
				deviceW(*x2 - *x1), 
				deviceH(*y2 - *y1)));
	completeRelativeOperation();
	cairo_rectangle(dmh_cr, deviceX(*x1), 
					deviceY(*y2),	/* this is y2, not y1, to deal with Cairo's flipped coordinate system */
					deviceW(*x2 - *x1),
					deviceH(*y2 - *y1));
	setLineProperties((int) *rline, *width);
	cairo_set_source_rgb (dmh_cr, 0, 0, 0);
	cairo_stroke_preserve(dmh_cr);
	cairo_set_source_rgb (dmh_cr, *rfill, *rfill, *rfill);
	cairo_fill(dmh_cr);
}

/*==================================================================================
 psline - subroutine to output a line (absolute).
 ==================================================================================*/
void psline_ (double *x1, double *y1, double *x2, double *y2, double *rline, double *width) {
	DEBUGPRINT(("In psline from (%f, %f) to (%f, %f), with rline=%f and width=%f.\n", *x1, *y1, *x2, *y2, *rline, *width));
	completeRelativeOperation();
	cairo_move_to(dmh_cr, deviceX(*x1), deviceY(*y1));
	cairo_line_to(dmh_cr, deviceX(*x2), deviceY(*y2));
	setLineProperties(*rline, *width);
	cairo_set_source_rgb (dmh_cr, 0, 0, 0);	/* black */
	cairo_stroke(dmh_cr);
}
/*==================================================================================
 psmove - subroutine to move the current point
 ==================================================================================*/
void psmove_ (double *x1, double *y1) {
	DEBUGPRINT(("In psmove to (%f, %f)\n", deviceX(*x1), deviceY(*y1)));
	completeRelativeOperation();
	cairo_move_to(dmh_cr, deviceX(*x1), deviceY(*y1));
}

/*==================================================================================
 psrmov - subroutine to make a relative move
 ==================================================================================*/
void psrmov_ (double *dx, double *dy) {
	DEBUGPRINT(("In psrmov. delta= (%f, %f), device=(%f, %f)\n", *dx, *dy,deviceW(*dx), deviceH(*dy)));
	dmh_inRelProcess = 1;
	cairo_rel_move_to(dmh_cr, deviceW(*dx), -deviceH(*dy));	/* negative is due to flipped Cairo coordinates */
}

/*==================================================================================
 psrlin - subroutine to draw a relative line
 ==================================================================================*/
void psrlin_ (double *dx, double *dy, double *rline, double *width) {
	DEBUGPRINT(("In psrlin. delta= (%f, %f), device=(%f, %f) rline=%f, width=%f\n", *dx, *dy, deviceW(*dx), deviceH(*dy), *rline, *width));
	dmh_inRelProcess = 1;
	cairo_rel_line_to(dmh_cr, deviceW(*dx), -deviceH(*dy));	/* negative is due to flipped Cairo coordinates */
	setLineProperties((int) *rline, *width);
	cairo_set_source_rgb (dmh_cr, 0, 0, 0);	/* black */
}

/*==================================================================================
 pswtod - 
 ==================================================================================
void pswtod (x1,y1,x2,y2) {
 
// pswtod - subroutine to do something i can't remember

      implicit none

      double precision x1,y1,x2,y2,x0,y0

      double precision a,b,c,d,xt,yt
      common/ trans /a,b,c,d,xt,yt

      integer nps
      double precision xscale, yscale, xmn, ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      x0 = ((x1-xmn) * xscale)
      y0 = ((y1-ymn) * yscale)
 
      x2 = x0 * a + y0 * c + xt
      y2 = x0 * b + y0 * d + yt
 
      end

 /*==================================================================================
 pssctr - subroutine to set text properties.
 kfont: code for font size, weight, and face to use (face is currently overridden by plot option)
 scale: scaling for font (max of xscale, yscale)
 theta: angle of rotation, given as positive value for ccw rotation.
 ==================================================================================*/
void pssctr_ (int *kfont, double *xscale, double *yscale, double *theta) {
	cairo_font_slant_t slant;
	cairo_font_weight_t weight;

	char *fontFace;
	int	fontSize;
	double scaling = 1;
	double inTheta = *theta;	/* make copy for internal use */
	int myKFont = *kfont;

	DEBUGPRINT(("In pssctr. kfont=%i, Scaling=(%f, %f), Theta=%f\n", *kfont, *xscale, *yscale, *theta));
	switch (myKFont) {
		case 1:
			fontFace = "serif";
			weight = CAIRO_FONT_WEIGHT_NORMAL;
			slant = CAIRO_FONT_SLANT_ITALIC;
			fontSize = 14;
			break;
		case 2:
			fontFace = "serif";
			weight = CAIRO_FONT_WEIGHT_BOLD;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 14;
			break;
		case 3:
			fontFace = "serif";
			weight = CAIRO_FONT_WEIGHT_NORMAL;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 14;
			break;
		case 4:
			fontFace = "serif";
			weight = CAIRO_FONT_WEIGHT_NORMAL;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 12;
			break;
		case 5:
			fontFace = "sans";
			weight = CAIRO_FONT_WEIGHT_NORMAL;
			slant = CAIRO_FONT_SLANT_OBLIQUE;
			fontSize = 14;
			break;
		case 6:
			fontFace = "sans";
			weight = CAIRO_FONT_WEIGHT_BOLD;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 14;
			break;
		case 7:
			fontFace = "sans";
			weight = CAIRO_FONT_WEIGHT_NORMAL;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 14;
			break;
		case 8:
			fontFace = "sans";
			weight = CAIRO_FONT_WEIGHT_NORMAL;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 12;
			break;
		case 9:
			fontFace = "monospace";
			weight = CAIRO_FONT_WEIGHT_BOLD;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 12;
			break;
		case 10:
			fontFace = "monospace";
			weight = CAIRO_FONT_WEIGHT_NORMAL;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 10;
			break;
		case 11:
			fontFace = "monospace";
			weight = CAIRO_FONT_WEIGHT_NORMAL;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 8;
			break;
		case 12:
			fontFace = "monospace";
			weight = CAIRO_FONT_WEIGHT_NORMAL;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 24;
			break;
		case 13:
			fontFace = "monospace";
			weight = CAIRO_FONT_WEIGHT_NORMAL;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 30;
			break;
		default:
			fontFace = "sans";
			weight = CAIRO_FONT_WEIGHT_NORMAL;
			slant = CAIRO_FONT_SLANT_NORMAL;
			fontSize = 14;
			break;
	}
	
	scaling = MYMAX(*xscale, *yscale);
	fontSize *= scaling;
	
	if (dmh_inRotatedTransform) {
		cairo_set_matrix(dmh_cr, &dmh_unrotatedMatrix);	/* remove rotation */
		dmh_inRotatedTransform = 0;
		DEBUGPRINT(("Removing rotation.\n"));
	}
	
	if (inTheta != 0.0) {
		cairo_get_matrix(dmh_cr, &dmh_unrotatedMatrix);
		if (dmh_aspectRatio != 1.0) {
  		    /* the supplied Theta value does not take the plot aspect ratio into account. Need to compensate. */
			DEBUGPRINT(("Adjusting theta for aspect ratio.  Initial Theta=%f, AR=%f\n", inTheta, dmh_aspectRatio));

			inTheta = inTheta * M_PI / 180.0;
			inTheta = atan(tan(inTheta) / dmh_aspectRatio);
			inTheta = inTheta * 180.0 / M_PI;
			DEBUGPRINT(("Adjusted theta for aspect ratio.  Fixed Theta=%f\n", inTheta));
		}
		cairo_rotate(dmh_cr, -(inTheta)*M_PI/180); /* cairo's rotation is positive for ccw; perplex's is opposite */
		dmh_inRotatedTransform = 1;
	}
	
	DEBUGPRINT(("Setting toy font.\n"));
	cairo_select_font_face (dmh_cr, dmh_fontFace, slant, weight);

	DEBUGPRINT(("FontFace Status:%s\n",cairo_status_to_string(cairo_status(dmh_cr))));
	cairo_set_font_size(dmh_cr, (double) fontSize);
}

/*==================================================================================
	pstext - subroutine to output text strings.

		x, y - coordinates of the label's reference point, in real coordinates (but see valign and halign, below)
		text - character string to be output
		jchar - length of meaningful character string, 0 if unknown.
		textlen - automatic argument supplied by fortran, gives allocated length of string
		valign: integer alignment code for vertical (from possibly rotated character standpoint) alignment.
			0 = top-aligned
			1 = middle-aligned
			2 = bottom-aligned
			10-12 = y should be considered a number of lines from top margin of page, then alignment=valign-10
			100-902 = valign is given by the final digit.  The y-coordinate should be calculated
				based upon the given y, but then add a number of lines to move down the page from there
				given by this code/100 (e.g., 1-9 lines)
			1000-9002 = valign is given by the final digit.  The y-coordinate should be calculated
			based upon the given y, but then add a number of lines to move up the page from there
			given by this code/1000 (e.g., 1-9 lines)
			(e.g., 202 would mean that the text should be bottom-aligned to the reference point, and 
				the y coordinate is two lines below the supplied y coordinate
		halign: integer alignment code for horizontal (from possibly rotated character standpoint) alignment.
			0 = left-aligned
			1 = center-aligned
			2 = right-aligned
			10-12 = x should be considered a number of lines from left margin of page, then alignment=halign-10
			100-9002 = halign is given by the final digit.  The x coordinate should be calculated based upon
			    the given x, but add a number of points given by this code / 100 (e.g., 1-90).
			REFERENCE_TO_TRACKED_MIN_X - REFERENCE_TO_TRACKED_MIN_X+2 = halign is given by the final digit. The
				x coordinate to be used is the tracked minimum x.
 
==================================================================================*/
void pstext_ (double *x, double *y, char *text, int *jchar, int *valign, int *halign, int textlen) {
	cairo_text_extents_t te;
	cairo_font_extents_t fe;
	cairo_matrix_t curMatrix;
	
	double userX, userY;
	double tempX, tempY;
	int valignCode, halignCode;
	int linesToMove = 0;
	double pointsToMove = 0;
	double xOffset = 0;
	double yOffset = 0;
	char *temptext=malloc(255 * sizeof(char));
	char *mytext;
	
	strncpy(temptext, text, 254);	/* copy only as many characters as will fit in value */

	if (*jchar == 0) {
		temptext[textlen]='\0';
	} else {
		temptext[*jchar]='\0';
	}
	
	mytext = trim(temptext);
	compressSpaces(mytext);
	completeRelativeOperation();

	DEBUGPRINT(("In pstext.  Pt=(%f, %f), DevPoint=(%f, %f), Text=%s, jchar=%i, textlen=%i, valign=%i, halign=%i.\n", *x, *y, deviceX(*x), deviceY(*y), mytext, *jchar, textlen, *valign, *halign));
	cairo_text_extents(dmh_cr,mytext,&te);
	cairo_font_extents (dmh_cr, &fe);
	cairo_set_source_rgb (dmh_cr, 0, 0, 0);	/* black */

	userX = deviceX(*x);	/* change coordinates from real to device */
	userY = deviceY(*y);

	if (dmh_debug) {
		tempX = userX;
		tempY = userY;
		cairo_device_to_user(dmh_cr, &tempX, &tempY);	/* need this to deal with possibly rotated transform */
		cairo_rectangle(dmh_cr, tempX-4, tempY-4, 8, 8);
		cairo_set_source_rgb (dmh_cr, 1, 0.5, 0);
		cairo_fill(dmh_cr);
		cairo_set_source_rgb (dmh_cr, 0, 0, 0);
		cairo_move_to(dmh_cr, tempX, tempY);
	}
	
	DEBUGPRINT(("In pstext.  UserPt=(%f, %f) tempPt=(%f, %f).\n", userX, userY, tempX, tempY));
	if (*valign > 9 && *valign < 100) {
		/* y should be considered a number of lines from top of page.  Note that cairo thinks the origin is at top-left */
		userY = TOPMARGIN + ((int)*y * fe.height);
	} else if (*valign > 99 && *valign < 1000) {
		linesToMove = *valign / 100;
	} else if (*valign > 990 && *valign < 10000) {
		linesToMove = -(*valign / 1000);
	}
	DEBUGPRINT(("In pstext.  LinesToMove=%i.\n", linesToMove));
	if (*halign > 9 && *halign < 100) {
		/* x should be considered a number of lines from top of page.  Note that cairo thinks the origin is at top-left */
		userX = LEFTMARGIN + ((int)*x * fe.height);
	} else if (*halign > 99 && *halign < 10000) {
		pointsToMove = (*halign)/100;
	} else if (*halign >= REFERENCE_TO_TRACKED_MIN_X && *halign <= REFERENCE_TO_TRACKED_MIN_X+2) {
		userX = dmh_min_tracked_x;
		dmh_min_tracked_x = DBL_MAX;	/* Now turn off tracking */
		dmh_track_min_x = 0;
	}

	DEBUGPRINT(("In pstext.  UserPt=(%f, %f).\n", userX, userY));
	cairo_device_to_user(dmh_cr, &userX, &userY);	/* need this to deal with possibly rotated transform */
	cairo_move_to(dmh_cr, userX, userY);
	if (linesToMove != 0) {
		cairo_rel_move_to(dmh_cr, 0, linesToMove * fe.height);
	}
	
	if (pointsToMove != 0) {
		cairo_rel_move_to(dmh_cr, pointsToMove, 0);
	}
	
	if (dmh_debug) {
		cairo_get_current_point(dmh_cr, &tempX, &tempY);
		cairo_rectangle(dmh_cr, tempX-2, tempY-2, 4, 4);
		cairo_set_source_rgb (dmh_cr, 0, 0.5, 1);
		cairo_fill(dmh_cr);
		cairo_set_source_rgb (dmh_cr, 0, 0, 0);
		cairo_move_to(dmh_cr, tempX, tempY);
	}

	valignCode = *valign % 10;
	halignCode = *halign % 10;
	switch (valignCode) {
		case 0:	/* Top-aligned */
			yOffset = te.height + EXTRASPACEPTS;
			break;
		case 1:	/* Middle-aligned */
			yOffset = te.height * 0.5;
			break;
		case 2:	/* Bottom-aligned */
			yOffset = -EXTRASPACEPTS;
			break;
	}

	switch (halignCode) {
		case 0:	/* Left-aligned */
			xOffset = EXTRASPACEPTS;
			break;
		case 1:	/* Center-aligned */
			xOffset = -te.width * 0.5;
			break;
		case 2:	/* Right-aligned */
			xOffset = -te.width - EXTRASPACEPTS;
			break;
	}
	DEBUGPRINT(("In pstext.  Offset=(%f, %f).\n", xOffset, yOffset));
	cairo_rel_move_to(dmh_cr, xOffset, yOffset);

	if (dmh_debug) {
		cairo_get_current_point(dmh_cr, &tempX, &tempY);
		cairo_rectangle(dmh_cr, tempX-2, tempY-2, 4, 4);
		cairo_set_source_rgb (dmh_cr, 1, 0, 0.5);
		cairo_fill(dmh_cr);
		cairo_set_source_rgb (dmh_cr, 0, 0, 0);
		cairo_move_to(dmh_cr, tempX, tempY);
	}

	cairo_get_matrix(dmh_cr, &curMatrix);
	DEBUGPRINT(("Current matrix: %f, %f, %f, %f, %f, %f.\n", curMatrix.xx, curMatrix.yx, curMatrix.xy, curMatrix.yy, curMatrix.x0, curMatrix.y0));

	if (dmh_track_min_x) {
		cairo_get_current_point(dmh_cr, &tempX, &tempY);
		if (tempX < dmh_min_tracked_x) {
			dmh_min_tracked_x = tempX;
		}
	}
	cairo_show_text (dmh_cr, mytext);

	DEBUGPRINT(("Status:%s\n",cairo_status_to_string(cairo_status(dmh_cr))));

	free(temptext);
}

/*==================================================================================
 trkmnx - subroutine to turn on the tracking of the minimum x-coordinate of all text 
 	drawn until the correct code is given in pstext
 ==================================================================================*/
void trkmnx_ () {
	dmh_track_min_x = 1;
}

/*==================================================================================
 psublk - subroutine to remove double and leading blanks from text
 
 //     text - character string 
 //     jchar - length of unblanked character string, 0 if unknown.
 
 ==================================================================================*/
void psublk_ (char *text, int *jchar, int textlen) {
	// This function iis just a stub.  All the trimming gets done now from pstext
}

/*==================================================================================
 trim - trims spaces from ends of string p
 ==================================================================================*/
char *trim( register char *str)
{
  char *end;

  // Trim leading space
  while(isspace(*str)) str++;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  // Write new null terminator
  *(end+1) = 0;

  return str;
}

/*==================================================================================
 compressSpaces - replaces multiple spaces within string p with a single space
 ==================================================================================*/
void compressSpaces(char *str)
{
	
	/* Internet version: */
	char *dst = str;
	
	for (; *str; ++str) {
		*dst++ = *str;
		if (*str == ' ') {
			do ++str; while (*str == ' ');
			--str;
		}
	}
	*dst = 0;
}


