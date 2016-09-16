/*************************************************************************\
* Copyright (c) 2016 UChicago Argonne, LLC,
*               as Operator of Argonne National Laboratory.
* This file is distributed subject to a Software License Agreement
* found in file LICENSE that is included with this distribution. 
\*************************************************************************/


/*
  Written by Dohn A. Arms, Argonne National Laboratory
  Send comments to dohnarms@anl.gov
  
  0.1   -- July 2005
  0.1.1 -- December 2006
           Added support for files that have more than 32k points
  1.0.0 -- November 2009
           Basically gutted the old code, getting rid of mda-load library
           in order to access data directly
  1.0.1 -- August 2010
           Show actual offset of the scan and PV Extras as encountered
  1.1   -- November 2010
  1.1.1 -- March 2011
           Have counted strings immediately freed after printing them.
  1.2   -- March 2011
           Fixed integer issues by tying short to int16_t, long to int32_t,
           and char to int8_t.  Changed %li to %i in printf's.  For MacOS
           Darwin, add fix to use xdr_char instead of xdr_int8_t.
  1.2.1 -- January 2012
           Minor build tweak
  1.2.2 -- June 2012
  1.3.0 -- February 2013
           Use printf correctly
  1.3.1 -- February 2014
           Apply XDR hack to file
  1.4.0 -- July 2016

 */


/*****************  mda_dump.c  *****************/

#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>

#include <unistd.h>

#define VERSION       "1.4.1 (August 2016)"
#define YEAR          "2016"
#define VERSIONNUMBER "1.4.1"


#ifdef XDR_HACK
  #include "xdr_hack.h"
#else
  #include <rpc/types.h>
  #include <rpc/xdr.h>

  static bool_t xdr_counted_string( XDR *xdrs, char **p)
  {
    int mode;
    int32_t length;

    mode = (xdrs->x_op == XDR_DECODE);

    /* If writing, obtain the length */
    if( !mode)
      length = strlen(*p);

    /* Transfer the string length */
    if( !xdr_int32_t(xdrs, &length))
      return 0;

    /* If reading, obtain room for the string */
    if (mode)
      {
        *p = (char *) malloc( (length + 1) * sizeof(char) );
        (*p)[length] = '\0'; /* Null termination */
      }

    /* If the string length is nonzero, transfer it */
    return(length ? xdr_string(xdrs, p, length) : 1);
  }
#endif


/*
   This function was set up to fflush after every printf in case of a crash
   but that was too slow, but I left it in case someone wants to reenable .
*/
void print( char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vprintf( fmt, args);
  va_end(args);

  /*  fflush( stdout);  */
}


int32_t li_print( XDR *xdrs, char *fmt )
{
  int32_t ll;

  if( !xdr_int32_t(xdrs, &ll) )
    {
      fflush( stdout);
      exit(1);
    }

  printf( fmt, ll);

  return ll;
}


int16_t si_print( XDR *xdrs, char *fmt )
{
  int16_t ss;

  if( !xdr_int16_t(xdrs, &ss) )
    {
      fflush( stdout);
      exit(1);
    }

  printf( fmt, ss);

  return ss;
}


float f_print( XDR *xdrs, char *fmt )
{
  float ff;

  if( !xdr_float(xdrs, &ff) )
    {
      fflush( stdout);
      exit(1);
    }

  printf( fmt, ff);

  return ff;
}

double d_print( XDR *xdrs, char *fmt )
{
  double dd;

  if( !xdr_double(xdrs, &dd) )
    {
      fflush( stdout);
      exit(1);
    }

  printf( fmt, dd);

  return dd;
}

// does NOT return string, but simply frees it as nothing cares
// about the actual value, and this eats memory stupidly
void cs_print( XDR *xdrs, char *fmt )
{
  char *cs = NULL;

  if( !xdr_counted_string( xdrs, &cs ) )
    {
      fflush( stdout);
      exit(1);
    }

  printf( fmt, cs);
  free( cs);
}


int mda_dump_header( XDR *xdrs)
{
  int16_t rank;
  
  int32_t extrapvs;

  int i;

  f_print( xdrs, "Version = %g\n" );
  li_print( xdrs, "Scan number = %i\n");
  rank = si_print( xdrs, "Data rank = %i\n");

  if( rank <= 0)
    {
      fflush(stdout);
      exit(1);
    }
  
  print("Largest scan bounds = ");
  for( i = 0; i < rank; i++)
    {
      if( i)
	print(" x ");
      li_print( xdrs, "%i");
    }
  print("\n");

  li_print( xdrs, "Regularity = %i\n");

  extrapvs = li_print( xdrs, "File offset to extra PV's = %i bytes\n");

  if( extrapvs < 0)
    {
      fflush(stdout);
      exit(1);
    }
  
  print("\n\n");
  
  if( extrapvs)
    return 1;
  else
    return 0;
}


void mda_dump_scan( XDR *xdrs)
{
  int i, j;
  
  int16_t rank, num_pos, num_det, num_trg;
  
  int32_t req_pts, cmp_pts;
  int32_t *offsets;
  
  offsets = NULL;
  
  
  rank = si_print( xdrs, "This scan's rank = %i\n");
  req_pts = li_print( xdrs, "Number of requested points = %i\n");
  cmp_pts = li_print( xdrs, "Last completed point = %i\n");
  
  if( (rank <= 0) || (req_pts <= 0) || (cmp_pts < 0) )
    {
      fflush(stdout);
      exit(1);
    }
  
  if( rank > 1)
    {
      offsets = (int32_t *) malloc( sizeof( int32_t) * req_pts);
      
      print("Offsets to lower scans = ");
      for( i = 0; i < req_pts; i++)
	{
	  if( i)
	    print(", ");
	  offsets[i] = li_print( xdrs, "%i");
	}
      print("\n");
    }
  
  cs_print( xdrs, "Scan name = \"%s\"\n");
  cs_print( xdrs, "Scan time = \"%s\"\n");
  
  num_pos = si_print( xdrs, "Number of positioners = %i\n");
  num_det = si_print( xdrs, "Number of detectors = %i\n");
  num_trg = si_print( xdrs, "Number of triggers = %i\n");

  if( num_pos < 0)
    {
      fflush(stdout);
      exit(1);
    }

  for( i = 0; i < num_pos; i++)
    {
      print( "\nPositioner Data Set #%i\n", i+1);

      si_print( xdrs, "              Positioner: %i\n");
      cs_print( xdrs, "                    Name: %s\n");
      cs_print( xdrs, "             Description: %s\n");
      cs_print( xdrs, "               Step Mode: %s\n");
      cs_print( xdrs, "                    Unit: %s\n");
      cs_print( xdrs, "           Readback Name: %s\n");
      cs_print( xdrs, "    Readback Description: %s\n");
      cs_print( xdrs, "           Readback Unit: %s\n");
    }

  if( num_det < 0)
    {
      fflush(stdout);
      exit(1);
    }

  for( i = 0; i < num_det; i++)
    {
      print( "\nDetector Data Set #%i\n", i+1);
      
      si_print( xdrs, "        Detector: %2i\n");
      cs_print( xdrs, "            Name: %s\n");
      cs_print( xdrs, "     Description: %s\n");
      cs_print( xdrs, "            Unit: %s\n");
    }

  if( num_trg < 0)
    {
      fflush(stdout);
      exit(1);
    }

  for( i = 0; i < num_trg; i++)
    {
      print( "\nTrigger #%i\n", i+1);

      si_print( xdrs, "    Trigger: %i\n");
      cs_print( xdrs, "       Name: %s\n");
      f_print(  xdrs, "    Command: %g\n");
    }

  for( i = 0; i < num_pos; i++)
    {
      print( "\nPositioner #%i data:\n", i+1);
      for( j = 0; j < req_pts; j++)
	d_print( xdrs, "%.9g ");
      print( "\n");
    }

  for( i = 0; i < num_det; i++)
    {
      print( "\nDetector #%i data:\n", i+1);
      for( j = 0; j < req_pts; j++)
	f_print( xdrs, "%.9g ");
      print( "\n");
    }

  if( rank > 1)
    for( i = 0; i < req_pts; i++)
      {
      if( offsets[i] == 0)
          break;
        print( "\n\n%i-D Subscan #%i\n", rank - 1, i+1);
        print( "Offset = %i\n\n", xdr_getpos( xdrs) );
        mda_dump_scan( xdrs);
      }

}


void mda_dump_extra( XDR *xdrs)
{
  enum PV_TYPES { EXTRA_PV_STRING=0, EXTRA_PV_INT8=32,  EXTRA_PV_INT16=29,
		  EXTRA_PV_INT32=33, EXTRA_PV_FLOAT=30, EXTRA_PV_DOUBLE=34 };

  int i, j;
  
  int16_t pvs, type, count;

  count = 0;

  print( "\n\nExtra PV Offset = %i", xdr_getpos( xdrs) );

  pvs = si_print( xdrs, "\n\nNumber of Extra PV's = %i.\n");
  if( pvs < 0)
    {
      fflush(stdout);
      exit(1);
    }


  for( i = 0 ; i < pvs; i++)
    {
      print( "\nExtra PV #%i:\n", i+1);
      
      cs_print( xdrs, "    Name = %s\n");
      cs_print( xdrs, "    Description = %s\n");

      if( !xdr_int16_t(xdrs, &type) )
	return;

      if( (type != EXTRA_PV_STRING) && (type != EXTRA_PV_INT8) &&
	  (type != EXTRA_PV_INT16) &&  (type != EXTRA_PV_INT32) &&
	  (type != EXTRA_PV_FLOAT) && (type != EXTRA_PV_DOUBLE))
	{
	  print( "    Type = %i (UNKNOWN)\n", type);
	  print( "\nExiting......\n");
	  exit(2);
	}

      if( type != EXTRA_PV_STRING)
	{
	  count = si_print( xdrs, "    Count = %i\n");
	  cs_print( xdrs, "    Unit = %s\n");
	}

      switch(type)
	{
	case EXTRA_PV_STRING:
	  print( "    Type = %i (STRING)\n", type);
	  cs_print( xdrs, "    Value = \"%s\"\n");
	  break;
	case EXTRA_PV_INT8:
	  {
	    int8_t *bytes;

	    print( "    Type = %i (INT8)\n", type);
	    print( "    Value%s = ", (count == 1) ? "" : "s");

	    bytes = (int8_t *) malloc( count * sizeof(int8_t));
#ifndef DARWIN
	    if( !xdr_vector( xdrs, (char *) bytes, count,
                             sizeof( int8_t), (xdrproc_t) xdr_int8_t))
#else
              // MacOS Darwin is missing xdr_int8_t,
              //have to fake it with xdr_char
	    if( !xdr_vector( xdrs, (char *) bytes, count,
			     sizeof( int8_t), (xdrproc_t) xdr_char))
#endif
	      return;

	    for( j = 0; j < count; j++)
	      {
		if( j)
		  print( ", ");
		print( "%i", bytes[j]);
	      }
	    print( "\n");
	  }
	  break;
	case EXTRA_PV_INT16:
	  {
	    int16_t *shorts;

	    print( "    Type = %i (INT16)\n", type);
	    print( "    Value%s = ", (count == 1) ? "" : "s");

	    shorts = (int16_t *) malloc( count * sizeof(int16_t));
	    if( !xdr_vector( xdrs, (char *) shorts, count,
			     sizeof( int16_t), (xdrproc_t) xdr_int16_t))
	      return;

	    for( j = 0; j < count; j++)
	      {
		if( j)
		  print( ", ");
		print( "%i", shorts[j]);
	      }
	    print( "\n");

	    free (shorts);
	  }
	  break;
	case EXTRA_PV_INT32:
	  {
	    int32_t *longs;
	    
	    print( "    Type = %i (INT32)\n", type);
	    print( "    Value%s = ", (count == 1) ? "" : "s");

	    longs = (int32_t *) malloc( count * sizeof(int32_t));
	    if( !xdr_vector( xdrs, (char *) longs, count,
			     sizeof( int32_t), (xdrproc_t) xdr_int32_t))
	      return ;

	    for( j = 0; j < count; j++)
	      {
		if( j)
		  print( ", ");
		print( "%i", longs[j]);
	      }
	    print( "\n");

	    free( longs);
	  }
	  break;
	case EXTRA_PV_FLOAT:
	  {
	    float *floats;

	    print( "    Type = %i (FLOAT)\n", type);
	    print( "    Value%s = ", (count == 1) ? "" : "s");

	    floats = (float *) malloc( count * sizeof(float));
	    if( !xdr_vector( xdrs, (char *) floats, count,
			     sizeof( float), (xdrproc_t) xdr_float))
	      return ;

	    for( j = 0; j < count; j++)
	      {
		if( j)
		  print( ", ");
		print( "%.9g", floats[j]);
	      }
	    print( "\n");

	    free( floats);
	  }
	  break;
	case EXTRA_PV_DOUBLE:
	  {
	    double *doubles;

	    print( "    Type = %i (DOUBLE)\n", type);
	    print( "    Value%s = ", (count == 1) ? "" : "s");

	    doubles = (double *) malloc( count * sizeof(double));
	    if( !xdr_vector( xdrs, (char *) doubles, count,
			     sizeof( double), (xdrproc_t) xdr_double))
	      return;

	    for( j = 0; j < count; j++)
	      {
		if( j)
		  print( ", ");
		print( "%.9g", doubles[j]);
	      }
	    print( "\n");

	    free( doubles);
	  }
	  break;
	}
    }
}



int mda_dump( char *file)
{
#ifndef XDR_HACK
  XDR xdrs;
#endif
  XDR *xdrstream;

  FILE *input;

  int extraflag;


  if( (input = fopen( file, "rb")) == NULL)
    {
      fprintf(stderr, "Can't open file!\n");
      return 1;
    }


#ifdef XDR_HACK
  xdrstream = input;
#else
  xdrstream = &xdrs;
  xdrstdio_create(xdrstream, input, XDR_DECODE);
#endif

  extraflag = mda_dump_header( xdrstream);
  mda_dump_scan( xdrstream);
  if( extraflag)
    mda_dump_extra( xdrstream);


#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif

  fclose(input);

  return 0;
}


void help(void)
{
  printf("Usage:  mda-dump [-hv] FILE\n"
         "Does a dump of everything in the EPICS MDA file, FILE.\n"
         "\n"
         "-h  This help text.\n"
         "-v  Show version information.\n"
         "\n"
         "The output format of the information exactly follows the binary "
         "format\n"
         "in the MDA file.\n"
);

}

void version(void)
{
  printf("mda-dump %s\n"
         "\n"
         "Copyright (c) %s UChicago Argonne, LLC,\n"
         "as Operator of Argonne National Laboratory.\n"
         "\n"
         "Written by Dohn Arms, dohnarms@anl.gov.\n",
         VERSION, YEAR);
}


int main( int argc, char *argv[])
{
  int opt;

  while((opt = getopt( argc, argv, "hv")) != -1)
    {
      switch(opt)
        {
        case 'h':
          help();
          return 0;
          break;
        case 'v':
          version();
          return 0;
          break;
        case ':':
          // option normally resides in 'optarg'
          printf("Error: option missing its value!\n");
          return -1;
          break;
        }
    }

  if( ((argc - optind) == 0) || ((argc - optind) > 1) )
    {
      printf("For help, type: mda-dump -h\n");
      return 0;
    }

  printf("********** mda-dump %s generated output **********\n\n\n", 
         VERSIONNUMBER);

  if( mda_dump( argv[optind]) )
    {
      fprintf(stderr, "Error!\n");
      return 1;
    }
  
  return 0;
}

