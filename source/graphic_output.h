/*
 *  graphic_output.h
 *  
 *  A library of subroutines to generate PDF/PS/SVG files and draw graphics primitives.
 *	Designed to be a drop-in substitute for Perplex's pslib.f postscript generation routines.
 *	(Hence all the psxxxx names instead of more reasonable names).  Uses Cairo to do the
 *	heavy lifting.
 *
 *  Created by David Hirsch on 1/30/10.
 *  Copyright 2010 Western Washington University. All rights reserved.
 *
 */

char *trim( register char *p);
void psopen_(char *fname, int fnamelen);
void psclos_();
void closeSurface();
void pselip_ (double *xor, double *yor, double *dx, double *dy, double *rline, double *width, int *ifill);
void pssctr_ (int *kfont, double *xs, double *ys, double *theta);
void psbspl_ (double *x, double *y, int *npts, double *rline, double *width, int *ifill);
void psrpgn_ (double *x1, double *y1, double *rx, double *ry, int *npts, double *rline, double *width,int *ifill);
void setLineProperties (int dashtype, double width);
void setFillType (int ifill);
void compressSpaces(char *str);
char getOptionForKey(char *keyword, char *value, int valueSize);
char getCompleteOptionForKey(char *keyword, char *value, int valueSize);
char get2OptionsForKey(char *keyword, char *value1, char *value2, int valueSize);
void pspygn_ (double *x, double *y, int *npts, double *rline, double *width, int *ifill);
void pspyln_ (double *x, double *y, int *npts, double *rline, double *width, int *ifill);
void psstrn_ (double *xs, double *ys, double *xt, double *yt, double *theta);
void pspltrgn_ (int *plotnum);
void psssc2_ (double *xmin, double *xmax, double *ymin, double *ymax);
void psrect_ (double *x1, double *x2, double *y1, double *y2, double *rline, double *width, int *ifill);
void psrecr_ (double *x1, double *x2, double *y1, double *y2, double *rline, double *width, double *rfill);
void psline_ (double *x1, double *y1, double *x2, double *y2, double *rline, double *width);
void psmove_ (double *x1, double *y1);
void psrmov_ (double *dx, double *dy);
void psrlin_ (double *dx, double *dy, double *rline, double *width);
void pstext_ (double *x, double *y, char *text, int *jchar, int *valign, int *halign, int textlen);
void psublk_ (char *text, int *jchar, int textlen);
void getControlPoints(double x0, double y0,
					  double x1, double y1,
					  double x2, double y2,
					  double x3, double y3,
					  double *ctrl1_x, double *ctrl1_y,
					  double *ctrl2_x, double *ctrl2_y);
double deviceX (double inX);
double deviceY (double inY);
double deviceW (double inWidth);
double deviceH (double inHeight);
void completeRelativeOperation ();

enum OUTPUTTYPE {PDFTYPE, PSTYPE, SVGTYPE };