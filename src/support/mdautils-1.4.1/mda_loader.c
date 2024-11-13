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
           Added support for files that have more than 32k points.
  0.2.0 -- October 2007
           Added several new functions for accessing scans or the extra PV's
           without loading the entire file (never publically released).
  0.2.1 -- March 2009
           Removed several memory leaks. Made trace structure consistent.
  1.0.0 -- October 2009
           Release with vastly changed mda-utils.
           Renamed structures to give them "mda_" prefix.
  1.0.1 -- June 2010
           Added checking in header code to make sure there were no screwy
           dimensions, by looking for 0xFFFFFFFF in one of them.
  1.0.2 -- August 2010
           Use offsets to load scans, as bad files sometimes have them out 
           of order.
  1.1   -- November 2010
  1.1.1 -- February 2011
           Realloc memory in screwy xdr counted strings.
           Added <stdio.h> to remove MacOS warning. (From J. Lewis Muir)
  1.2   -- March 2011
           Fixed integer issues by tying short to int16_t, long to int32_t,
           and char to int8_t.
           Renamed DBR_* type variables to EXTRA_PV_* type variables, to
           make explicit the type of integer used (which broke compatibility
           with the EPICS DBR type names).
  1.2.1 -- January 2012
  1.2.2 -- June 2012
           Fixed major bug with INT8 Extra PVs.
  1.3.0 -- February 2013
           Don't load files that have nonsensical scan dimensions.
  1.3.1 -- February 2014
           Added a check for bad zero values in the scan offsets.
           Added support for XDR hack code.
  1.4.0 -- July 2016
           Added better error handling: checking for malloc() returning zero,
           cleaning up allocated memory when there is an error in loading.
           It checks MDA file version in case of future format change.
           Added mda_update(), for loading an MDA file, but using an earlier
           loaded version's data to reduce loading.
           Added mda_test() function that does a version of mda_load to use
           for all the error checks, but doesn't allocate memory for detectors
           or positioners to keep the memory usage from being crazy.
 */


/* 
   comments around else blocks that assign NULL values to pointers
   are there to show that it is acknowledged that this needs to happen,
   but the calls to calloc() have taken care of NULL assignment
*/


/****************  mda_loader.c  **********************/


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "mda-load.h"

enum test_option {NOTEST, TEST};
enum recurse_option {NORECURSE, RECURSE};

#ifdef XDR_HACK
  #include "xdr_hack.h"
#else
  #include <rpc/types.h>
  #include <rpc/xdr.h>

  static bool_t xdr_counted_string( XDR *xdrs, char **p)
  {
    int mode;
    int32_t length;

    bool_t ret_bool;
 
    mode = (xdrs->x_op == XDR_DECODE);

    /* If writing, obtain the length */
    if( !mode) 
      length = strlen(*p);
 
    /* Transfer the string length */
    /* had to change this to int32_t, as values would sometimes turn negative
       resizing is done below */
    if( !xdr_int32_t(xdrs, &length)) 
      return 0;

    /* If reading, obtain room for the string */
    if (mode)
      {
        *p = malloc( (length + 1) * sizeof(char) );
        if( *p == NULL)
          return 0;
        (*p)[length] = '\0'; /* Nul termination */
      }

    /* If the string length is nonzero, transfer it */
    if( length)
      {
        ret_bool = xdr_string(xdrs, p, length);

        // this fix is for when the size is not set correctly,
        // reducing the memory allocated; it only checks if size is large;
        // sometimes 4MB is allocated for an 8 character string....
        if( (length > 4095) && ret_bool && mode )
          { 
            int32_t l;

            l = strlen(*p);
            if( l < length)
              {
                char *q;

                q = realloc( *p, l+1);
                if( q != NULL)  // if NULL, just use the original
                  *p = q;
              }
          }
      }
    else
      ret_bool = 1;


    return ret_bool;
  }
#endif


static struct mda_header *header_read( XDR *xdrs)
{
  struct mda_header *header;

  int ver;
  int i;

  header = calloc( 1, sizeof(struct mda_header));
  if( header == NULL)
    return NULL;

  if( !xdr_float(xdrs, &(header->version)) )
    return mda_header_unload( header), NULL;
  // because version floats weren't precisely the correct values
  ver = (int) (100.0*header->version + 0.5);
  if( ((ver != 120) && (ver != 130) && (ver != 140)) ||
      !xdr_int32_t(xdrs, &(header->scan_number)) ||
      !xdr_int16_t(xdrs, &(header->data_rank)) ||
      (header->data_rank < 1) )
    return mda_header_unload( header), NULL;

  header->dimensions = malloc( header->data_rank * sizeof(int32_t));
  if( (header->dimensions == NULL) ||
      !xdr_vector( xdrs, (char *) header->dimensions, header->data_rank, 
		   sizeof( int32_t), (xdrproc_t) xdr_int32_t) ||
      !xdr_int16_t(xdrs, &(header->regular)) ||
      !xdr_int32_t(xdrs, &(header->extra_pvs_offset)) ) 
    return mda_header_unload( header), NULL;

  for( i = 0; i < header->data_rank; i++)
    if(header->dimensions[i] < 1)
      return mda_header_unload( header), NULL;
      
  return header;
}


static struct mda_positioner *positioner_read(XDR *xdrs)
{
  struct mda_positioner *positioner;


  positioner = calloc( 1, sizeof(struct mda_positioner));
  if( positioner == NULL)
    return NULL;

  if( !xdr_int16_t(xdrs, &(positioner->number) ) ||
      !xdr_counted_string( xdrs, &(positioner->name) ) ||
      !xdr_counted_string( xdrs, &(positioner->description) ) ||
      !xdr_counted_string( xdrs, &(positioner->step_mode) ) ||
      !xdr_counted_string( xdrs, &(positioner->unit) ) ||
      !xdr_counted_string( xdrs, &(positioner->readback_name) ) ||
      !xdr_counted_string( xdrs, &(positioner->readback_description) ) ||
      !xdr_counted_string( xdrs, &(positioner->readback_unit) ) )
    {
      free( positioner->readback_description);
      free( positioner->readback_name);
      free( positioner->unit);
      free( positioner->step_mode);
      free( positioner->description);
      free( positioner->name);
      free( positioner);
      return NULL;
    }

  return positioner;
}



static struct mda_detector *detector_read(XDR *xdrs)
{
  struct mda_detector *detector;


  detector = calloc( 1, sizeof(struct mda_detector));
  if( detector == NULL)
    return NULL;

  if( !xdr_int16_t(xdrs, &(detector->number) ) ||
      !xdr_counted_string( xdrs, &(detector->name) ) ||
      !xdr_counted_string( xdrs, &(detector->description) ) ||
      !xdr_counted_string( xdrs, &(detector->unit) ) )
    {
      free( detector->description);
      free( detector->name);
      free( detector);
      return NULL;
    }

  return detector;
}




static struct mda_trigger *trigger_read(XDR *xdrs)
{
  struct mda_trigger *trigger;


  trigger = calloc( 1, sizeof(struct mda_trigger));
  if( trigger == NULL)
    return NULL;

  if( !xdr_int16_t(xdrs, &(trigger->number) ) ||
      !xdr_counted_string( xdrs, &(trigger->name ) )||
      !xdr_float(xdrs, &(trigger->command) ) )
    {
      free( trigger->name);
      free( trigger);
      return NULL;
    }

  return trigger;
}



/* this function is recursive, due to the nature of the file format */
/* it can be turned off by making recursive 0, the rank says what */
/* the file structure thinks this scan rank is to do error checking, */
/* test_flag switches to not allocating memory for data */
static struct mda_scan *scan_read(XDR *xdrs, enum recurse_option recurse_flag, 
                                  int rank, enum test_option test_flag)
{
  struct mda_scan *scan;
 
  int i;

  scan = calloc( 1, sizeof(struct mda_scan));
  if( scan == NULL)
    return NULL;


  if( !xdr_int16_t(xdrs, &(scan->scan_rank)) ||
      ( scan->scan_rank != rank) ||
      !xdr_int32_t(xdrs, &(scan->requested_points)) ||
      ( scan->requested_points < 1) ||
      !xdr_int32_t(xdrs, &(scan->last_point)) ||
      ( scan->last_point < 0) ||
      ( scan->last_point > scan->requested_points) )
    return mda_scan_unload( scan), NULL;

  if( scan->scan_rank > 1)
    {
      scan->offsets = malloc( scan->requested_points * sizeof(int32_t));
      if( (scan->offsets == NULL) ||
          !xdr_vector( xdrs, (char *) scan->offsets, scan->requested_points, 
		       sizeof( int32_t), (xdrproc_t) xdr_int32_t) )
        return mda_scan_unload( scan), NULL;
      
      for( i = 0; i < scan->requested_points; i++)
        if( scan->offsets[i] < 0)
          return mda_scan_unload( scan), NULL;
      // there can be no zero offsets for the first "last_point" values
      for (i = 0; i < scan->last_point; i++)
      {
          if (scan->offsets[i] == 0)
          {
              if (i > 0)
              {
                  scan->last_point = i - 1;
                  break;
              }
              return mda_scan_unload(scan), NULL;
          }
      }
    }
  /* else */
  /*   scan->scan_rank = NULL; */

  if( !xdr_counted_string( xdrs, &(scan->name)) ||
      !xdr_counted_string( xdrs, &(scan->time)) ||
      !xdr_int16_t(xdrs, &(scan->number_positioners)) ||
      (scan->number_positioners < 0) ||
      !xdr_int16_t(xdrs, &(scan->number_detectors)) ||
      (scan->number_detectors < 0) ||
      !xdr_int16_t(xdrs, &(scan->number_triggers)) ||
      (scan->number_triggers < 0) )
    return mda_scan_unload( scan), NULL;

  scan->positioners =
    calloc( scan->number_positioners, sizeof(struct mda_positioner *));
  if( scan->positioners == NULL)
    return mda_scan_unload( scan), NULL;

  for( i = 0; i < scan->number_positioners; i++)
    {
      if( (scan->positioners[i] = positioner_read( xdrs)) == NULL )
        return mda_scan_unload( scan), NULL;
    }

  scan->detectors =  
    calloc( scan->number_detectors, sizeof(struct mda_detector *));
  if( scan->detectors == NULL)
    return mda_scan_unload( scan), NULL;
  for( i = 0; i < scan->number_detectors; i++)
    {
      if( (scan->detectors[i] = detector_read( xdrs)) == NULL )
        return mda_scan_unload( scan), NULL;
    }

  scan->triggers =  
    calloc( scan->number_triggers, sizeof(struct mda_trigger *));
  if( scan->triggers == NULL)
    return mda_scan_unload( scan), NULL;
  for( i = 0; i < scan->number_triggers; i++)
    {
      if( (scan->triggers[i] = trigger_read( xdrs)) == NULL )
        return mda_scan_unload( scan), NULL;
    }

  if( test_flag == NOTEST)
    {
      scan->positioners_data = 
        calloc( scan->number_positioners, sizeof(double *));
      if( scan->positioners_data == NULL)
        return mda_scan_unload( scan), NULL;
      for( i = 0 ; i < scan->number_positioners; i++)
        {
          scan->positioners_data[i] =
            malloc( scan->requested_points * sizeof(double));
          if( (scan->positioners_data[i] == NULL) || 
              !xdr_vector( xdrs, (char *) scan->positioners_data[i], 
                           scan->requested_points, 
                           sizeof( double), (xdrproc_t) xdr_double))
            return mda_scan_unload( scan), NULL;
        }
    }
  else
    {
      // I could change the file position to skip reading these values,
      // but I'm not actually supposed to know how XDR is formatted 
      double *pos;

      /* scan->positioners_data = NULL; */
      if( (pos = malloc( scan->requested_points * sizeof(double))) == NULL)
        return mda_scan_unload( scan), NULL;
      for( i = 0 ; i < scan->number_positioners; i++)
        {
          if( !xdr_vector( xdrs, (char *) pos, scan->requested_points, 
                           sizeof( double), (xdrproc_t) xdr_double))
            return mda_scan_unload( scan), NULL;
        }
      free(pos);
    }

  if( test_flag == NOTEST)
    {
      scan->detectors_data = 
        calloc( scan->number_detectors, sizeof(float *));
      if( scan->detectors_data == NULL)
        return mda_scan_unload( scan), NULL;
      for( i = 0 ; i < scan->number_detectors; i++)
        {
          scan->detectors_data[i] = 
            malloc( scan->requested_points * sizeof(float));
          if( (scan->detectors_data[i] == NULL) ||
              !xdr_vector( xdrs, (char *) scan->detectors_data[i], 
                           scan->requested_points, 
                           sizeof( float), (xdrproc_t) xdr_float))
            return mda_scan_unload( scan), NULL;
        }
    }
  else
    {
      float *det;

      /* scan->detectors_data = NULL; */
      if( (det = malloc( scan->requested_points * sizeof(float))) == NULL)
        return mda_scan_unload( scan), NULL;
      for( i = 0 ; i < scan->number_detectors; i++)
        {
          if( !xdr_vector( xdrs, (char *) det, scan->requested_points, 
                           sizeof( float), (xdrproc_t) xdr_float))
            return mda_scan_unload( scan), NULL;
        }
      free(det);
    }

  if( (scan->scan_rank > 1) && (recurse_flag == RECURSE))
    {
      scan->sub_scans = 
	calloc( scan->requested_points, sizeof( struct mda_scan *) );
      if( scan->sub_scans == NULL)
        return mda_scan_unload( scan), NULL;
      // scan->sub_scans[i] can remain NULL for i > scan->last_point
      for( i = 0 ; (i < scan->requested_points) && 
	     (scan->offsets[i] != 0); i++)
	/* it recurses here */
        {
          if( xdr_getpos( xdrs) != scan->offsets[i])
            {
              if( !xdr_setpos( xdrs,scan->offsets[i] ) )
                return mda_scan_unload( scan), NULL;
            }
          scan->sub_scans[i] = scan_read(xdrs, recurse_flag, 
                                         scan->scan_rank - 1, test_flag);
          // if a subscan beyond last point is read, allow to fail
          // as it is considered an unfinished scan
          if( (i < scan->last_point) && (scan->sub_scans[i] == NULL) )
            return mda_scan_unload( scan), NULL;
        }
    }
  /* else */
  /*   scan->sub_scans = NULL; */

  return scan;
}




static struct mda_pv *pv_read(XDR *xdrs)
{
  struct mda_pv *pv;

  pv = calloc( 1, sizeof(struct mda_pv));
  if( pv == NULL)
    return NULL;

  if( !xdr_counted_string( xdrs, &(pv->name)) ||
      !xdr_counted_string( xdrs, &(pv->description)) ||
      !xdr_int16_t(xdrs, &(pv->type)) )
    goto PV_Error;

  if( (pv->type != EXTRA_PV_STRING) && (pv->type != EXTRA_PV_INT8) && 
      (pv->type != EXTRA_PV_INT16) && (pv->type != EXTRA_PV_INT32) && 
      (pv->type != EXTRA_PV_FLOAT) && (pv->type != EXTRA_PV_DOUBLE) )
    goto PV_Error;
  
  if( pv->type != EXTRA_PV_STRING)
    {
      if( !xdr_int16_t(xdrs, &(pv->count)) ||
          (pv->count < 1) ||
          !xdr_counted_string( xdrs, &(pv->unit)) )
        goto PV_Error;
    }
  /* else */
  /*   { */
  /*     pv->count = 0; */
  /*     pv->unit = NULL; */
  /*   } */

  switch( pv->type)
    {
    case EXTRA_PV_STRING:
      if( !xdr_counted_string( xdrs, &(pv->values) ))
        goto PV_Error;
      break;
    case EXTRA_PV_INT8:
      pv->values = malloc( pv->count * sizeof(int8_t));
      if( (pv->values == NULL) ||
#ifndef DARWIN
          !xdr_vector( xdrs, pv->values, pv->count, 
		       sizeof( int8_t), (xdrproc_t) xdr_int8_t))
#else
        // MacOS Darwin is missing xdr_int8_t, have to fake it with xdr_char
          !xdr_vector( xdrs, pv->values, pv->count, 
		       sizeof( int8_t), (xdrproc_t) xdr_char))
#endif
        goto PV_Error;
      break;
    case EXTRA_PV_INT16:
      pv->values = malloc( pv->count * sizeof(int16_t));
      if( (pv->values == NULL) ||
          !xdr_vector( xdrs, pv->values, pv->count, sizeof( int16_t), 
                       (xdrproc_t) xdr_int16_t))
        goto PV_Error;
      break;
    case EXTRA_PV_INT32:
      pv->values = malloc( pv->count * sizeof(int32_t));
      if( (pv->values == NULL) || 
          !xdr_vector( xdrs, pv->values, pv->count, 
		       sizeof( int32_t), (xdrproc_t) xdr_int32_t))
        goto PV_Error;
      break;
    case EXTRA_PV_FLOAT:
      pv->values = malloc( pv->count * sizeof(float));
      if( (pv->values == NULL) || 
          !xdr_vector( xdrs, pv->values, pv->count, 
		       sizeof( float), (xdrproc_t) xdr_float))
        goto PV_Error;
      break;
    case EXTRA_PV_DOUBLE:
      pv->values = malloc( pv->count * sizeof(double));
      if( (pv->values == NULL) ||
          !xdr_vector( xdrs, pv->values, pv->count, 
		       sizeof( double), (xdrproc_t) xdr_double))
        goto PV_Error;
      break;
    }

  return pv;

 PV_Error:
  free( pv->name);
  free( pv->description);
  free( pv->unit);
  free( pv->values);
  free( pv);
  return NULL;
}




static struct mda_extra *extra_read(XDR *xdrs)
{
  struct mda_extra *extra;

  int i;

  extra = calloc( 1, sizeof(struct mda_extra));
  if( extra == NULL)
    return NULL;

  if( !xdr_int16_t(xdrs, &(extra->number_pvs) ) ||
      (extra->number_pvs < 0) )
    return mda_extra_unload( extra), NULL;

  if( extra->number_pvs)
    {
      extra->pvs = calloc( extra->number_pvs, sizeof( struct mda_pv *) );
      if( extra->pvs == NULL)
        return mda_extra_unload( extra), NULL;

      for( i = 0 ; i < extra->number_pvs; i++)
        { 
          if( (extra->pvs[i] = pv_read(xdrs)) == NULL )
            return mda_extra_unload( extra), NULL;
        } 
    }
  /* else */
  /*   { */
  /*     extra->pvs = NULL; */
  /*   } */

  return extra;
}


/////////////////////////////////////////////////////////////




static struct mda_file *mda_load_full( FILE *fptr, enum test_option test_flag)
{
#ifndef XDR_HACK
  XDR xdrs;
#endif
  XDR *xdrstream;

  struct mda_file *mda;
  struct mda_scan *scan;

  int i;

  rewind( fptr);

#ifdef XDR_HACK
  xdrstream = fptr;
#else
  xdrstream = &xdrs;
  xdrstdio_create(xdrstream, fptr, XDR_DECODE);
#endif

  if( (mda = calloc( 1, sizeof(struct mda_file))) == NULL)
    goto Load_Error;

  if( ((mda->header = header_read( xdrstream)) == NULL) ||
      ((mda->scan = scan_read( xdrstream, RECURSE, mda->header->data_rank, 
                               test_flag)) == NULL) )
        goto Load_Error;
  for( scan = mda->scan, i = 0; i < (mda->header->data_rank - 1); i++)
    {
      if( scan == NULL)
        goto Load_Error;
      scan = scan->sub_scans[0];
    }
  if( mda->header->extra_pvs_offset)
    {
      xdr_setpos(xdrstream, mda->header->extra_pvs_offset);
      mda->extra = extra_read(xdrstream);
    }
  /* else */
  /*   mda->extra = NULL; */
  
#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif

  return mda;

 Load_Error:
#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif
  mda_unload(mda);
  return NULL;
}

struct mda_file *mda_load( FILE *fptr)
{
  return mda_load_full( fptr, NOTEST);
}

int mda_test( FILE *fptr)
{
  struct mda_file *mda;

  mda = mda_load_full( fptr, TEST);
  if( mda == NULL)
    return 1;
  else
    mda_unload( mda);
  
  return 0;
}


struct mda_header *mda_header_load( FILE *fptr)
{
#ifndef XDR_HACK
  XDR xdrs;
#endif
  XDR *xdrstream;

  struct mda_header *header;

  rewind( fptr);

#ifdef XDR_HACK
  xdrstream = fptr;
#else
  xdrstream = &xdrs;
  xdrstdio_create(xdrstream, fptr, XDR_DECODE);
#endif

  header = header_read( xdrstream);
#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif
  if( header == NULL)
    return NULL;
  
  return header;
}



struct mda_scan *mda_scan_load( FILE *fptr)
{
  return mda_subscan_load( fptr, 0, NULL, 1);
}


struct mda_scan *mda_subscan_load( FILE *fptr, int depth, int *indices, 
                                   int recursive)
{
  struct mda_scan *scan;

#ifndef XDR_HACK
  XDR xdrs;
#endif
  XDR *xdrstream;

  struct mda_header *header;

  rewind( fptr);

#ifdef XDR_HACK
  xdrstream = fptr;
#else
  xdrstream = &xdrs;
  xdrstdio_create(xdrstream, fptr, XDR_DECODE);
#endif

  if( (header = header_read( xdrstream)) == NULL)
    goto Subscan_Error;

  if( (depth < 0) || (depth >= header->data_rank ) )
    goto Subscan_Error;

  if( depth)
    {
      int i;

      int16_t scan_rank;
      int32_t  requested_points;
      int32_t  last_point;
      int32_t *offsets;

      for( i = 0; i < depth; i++)
	{
	  if( indices[i] >= header->dimensions[i])
            goto Subscan_Error;
	}

      for( i = 0; i < depth; i++)
	{
	  if( !xdr_int16_t(xdrstream, &scan_rank) ||
              (scan_rank <= 1) ||
              !xdr_int32_t(xdrstream, &requested_points) ||
              (requested_points < 1) ||
              !xdr_int32_t(xdrstream, &last_point) ||
              (last_point < 0) ||
              (last_point > requested_points) ||
              (indices[i] >= last_point) )
            goto Subscan_Error;

	  offsets = malloc( requested_points * sizeof(int32_t));
          if( offsets == NULL)
            goto Subscan_Error;
          if( !xdr_vector( xdrstream, (char *) offsets, requested_points, 
                           sizeof( int32_t), (xdrproc_t) xdr_int32_t))
            {
              free(offsets);
              goto Subscan_Error;
            }

	  if( (offsets[indices[i]] == 0) || 
              !xdr_setpos( xdrstream, offsets[indices[i]]) )
            {
              free( offsets);
              goto Subscan_Error;
            }

	  free( offsets);
	}
    }

  if( (scan = scan_read( xdrstream, recursive ? RECURSE : NORECURSE, 
                         header->data_rank - depth, NOTEST)) == NULL)
    return NULL;

#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif

  mda_header_unload(header);

  return scan;

 Subscan_Error:
#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif
  mda_header_unload(header);
  return NULL;
}


struct mda_extra *mda_extra_load( FILE *fptr)
{
  struct mda_extra *extra;

#ifndef XDR_HACK
  XDR xdrs;
#endif
  XDR *xdrstream;

  struct mda_header *header;

  rewind( fptr);

#ifdef XDR_HACK
  xdrstream = fptr;
#else
  xdrstream = &xdrs;
  xdrstdio_create(xdrstream, fptr, XDR_DECODE);
#endif

  if( (header = header_read( xdrstream)) == NULL)
    goto Extra_Error;

  extra = NULL;
  if( header->extra_pvs_offset)
    {
      if( !xdr_setpos( xdrstream, header->extra_pvs_offset) ||
          ((extra = extra_read( xdrstream)) == NULL) )
        goto Extra_Error;
    }

#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif

  mda_header_unload(header);

  return extra;

 Extra_Error:
#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif
  mda_header_unload(header);
  return NULL;
}

//////////////////////////////////////////////////////////////////
// updater

// These are fairly complicated functions, to allow for reading only 
// the part of the MDA file that has changed since the last load or update.
// It will back up and utilize the standard loaders if situation call for it.
static struct mda_scan *scan_check( XDR *xdrs, int rank,
                                    struct mda_scan *orig_scan)
{
  struct mda_scan *scan, *temp_scan;
  int org_pos;
  
  int i;

  if( (orig_scan == NULL) || (rank == 1) )
    return scan_read( xdrs, RECURSE, rank, NOTEST);

  org_pos = xdr_getpos(xdrs);

  if( (scan = scan_read( xdrs, NORECURSE, rank, NOTEST)) == NULL)
    return NULL;
  if( orig_scan->sub_scans == NULL) // nonrecursive read of original scan
    goto Scan_And_Go;

  // Do the easy checks to see that the scans are the same
  // if not, then reread the whole thing
  if( (scan->scan_rank != orig_scan->scan_rank) ||
      (scan->requested_points != orig_scan->requested_points) ||
      (scan->scan_rank != orig_scan->scan_rank) ||
      strcmp( scan->name, orig_scan->name) ||
      strcmp( scan->time, orig_scan->time) ||
      (scan->number_positioners != orig_scan->number_positioners) ||
      (scan->number_detectors != orig_scan->number_detectors) ||
      (scan->number_triggers != orig_scan->number_triggers) )
    goto Scan_And_Go;

  // steal the scans from original scan
  scan->sub_scans = orig_scan->sub_scans;
  orig_scan->sub_scans = NULL;
  for( i = orig_scan->last_point; i < scan->last_point; i++)
    {
      if( !xdr_setpos( xdrs, scan->offsets[i]) )
        goto Scan_Update_Error;
      // if temp_scan is NULL, taken care of in scan_check
      temp_scan = scan->sub_scans[i];
      scan->sub_scans[i] = scan_check( xdrs, rank - 1, temp_scan);
      mda_scan_unload( temp_scan);
      if( scan->sub_scans[i] == NULL)
        goto Scan_Update_Error;
    }
  if( (scan->last_point < scan->requested_points) && (scan->offsets[i] != 0) )
    {
      if( !xdr_setpos( xdrs, scan->offsets[i])  )
        goto Scan_Update_Error;
      temp_scan = scan->sub_scans[i];
      // allow returned value to be NULL
      scan->sub_scans[i] = scan_check( xdrs, rank - 1, temp_scan);
      mda_scan_unload( temp_scan);
    }

  return scan;

 Scan_Update_Error:
  mda_scan_unload( scan);
  return NULL;

 Scan_And_Go:
  mda_scan_unload( scan);
  if( !xdr_setpos( xdrs, org_pos) )
    return NULL; 
  return scan_read( xdrs, RECURSE, rank, NOTEST);
}

struct mda_file *mda_update( FILE *fptr, struct mda_file *previous_mda)
{
#ifndef XDR_HACK
  XDR xdrs;
#endif
  XDR *xdrstream;

  struct mda_file *mda;

  int i;

  if( previous_mda == NULL)
    return mda_load( fptr);

  rewind( fptr);

#ifdef XDR_HACK
  xdrstream = fptr;
#else
  xdrstream = &xdrs;
  xdrstdio_create(xdrstream, fptr, XDR_DECODE);
#endif

  if( (mda = calloc( 1, sizeof(struct mda_file))) == NULL)
    goto Update_Error;

  if((mda->header = header_read( xdrstream)) == NULL)
    goto Update_Error;

  if( (mda->header->version != previous_mda->header->version) ||
      (mda->header->scan_number != previous_mda->header->scan_number) ||
      (mda->header->data_rank != previous_mda->header->data_rank) ||
      (mda->header->regular != previous_mda->header->regular) )
    goto Reload_And_Go;
  for( i = 0; i < mda->header->data_rank; i++)
    if( mda->header->dimensions[i] != previous_mda->header->dimensions[i])
      goto Reload_And_Go;
      
  if( (mda->scan = scan_check( xdrstream, mda->header->data_rank, 
                               previous_mda->scan)) == NULL) 
    goto Update_Error;

  // always reread extra PVs
  if( mda->header->extra_pvs_offset)
    {
      if( !xdr_setpos( xdrstream, mda->header->extra_pvs_offset ) ||
          ( (mda->extra = extra_read( xdrstream)) == NULL) )
        goto Update_Error;
    }
  /* else */
  /*   mda->extra = NULL; */
  
#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif

  mda_unload(previous_mda);
  return mda;

Reload_And_Go: // the file is different, just reload the thing from scratch
#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif
  mda_unload(previous_mda); // get rid of old data
  mda_unload(mda);  // get rid of what was started
return mda_load(fptr);

  Update_Error:
#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif
  mda_unload(previous_mda);
  mda_unload(mda);
  return NULL;
}



///////////////////////////////////////////////////////////////////
// unloaders


void mda_header_unload(struct mda_header *header)
{
  if( header == NULL)
    return;

  free( header->dimensions);
  free( header);
}


/* this function is recursive */
void mda_scan_unload( struct mda_scan *scan)
{
  int i;

  if( scan == NULL)
    return;

  if( (scan->scan_rank > 1) && (scan->sub_scans != NULL))
    {
      for( i = 0; (i < scan->requested_points) && (scan->sub_scans[i] != NULL);
           i++)
	mda_scan_unload( scan->sub_scans[i]);
    }
  free( scan->sub_scans);
  
  free( scan->offsets);
  free( scan->name);
  free( scan->time);

  if( scan->positioners != NULL)
    {
      for( i = 0; i < scan->number_positioners; i++)
        {
          if( scan->positioners[i] == NULL)
            continue;
          free(scan->positioners[i]->name);
          free(scan->positioners[i]->description);
          free(scan->positioners[i]->step_mode);
          free(scan->positioners[i]->unit);
          free(scan->positioners[i]->readback_name);
          free(scan->positioners[i]->readback_description);
          free(scan->positioners[i]->readback_unit);
          free(scan->positioners[i]);
        }
      free( scan->positioners);
    }

  if( scan->detectors != NULL)
    {
      for( i = 0; i < scan->number_detectors; i++)
        {
          if( scan->detectors[i] == NULL)
            continue;
          free(scan->detectors[i]->name);
          free(scan->detectors[i]->description);
          free(scan->detectors[i]->unit);
          free(scan->detectors[i]);
        }
      free( scan->detectors);
    }

  if( scan->triggers != NULL)
    {
      for( i = 0; i < scan->number_triggers; i++)
        {
          if( scan->triggers[i] == NULL)
            continue;
          free(scan->triggers[i]->name);
          free(scan->triggers[i]);
        }      
      free( scan->triggers);
    }

  if( scan->positioners_data != NULL)
    {
      for( i = 0 ; i < scan->number_positioners; i++)
        free( scan->positioners_data[i] );
      free( scan->positioners_data );
    }

  if( scan->detectors_data != NULL)
    {
      for( i = 0 ; i < scan->number_detectors; i++)
        free( scan->detectors_data[i] );
      free( scan->detectors_data );
    }

  free( scan);
}


void mda_extra_unload(struct mda_extra *extra)
{
  int i;

  if( extra == NULL)
    return;

  if( extra->pvs != NULL)
    {
      for( i = 0; i < extra->number_pvs; i++)
        {
          if( extra->pvs[i] == NULL)
            continue;
          free( extra->pvs[i]->name);
          free( extra->pvs[i]->description);
          free( extra->pvs[i]->unit);
          free( extra->pvs[i]->values);
          free( extra->pvs[i]);
        }
      free( extra->pvs);
    }
  free( extra);
}


/*  deallocates all the memory used for the mda file loading  */
void mda_unload( struct mda_file *mda)
{
  if( mda == NULL)
    return;

  // these functions check to make sure not NULL
  mda_header_unload(mda->header);
  mda_scan_unload(mda->scan);
  mda_extra_unload(mda->extra);
  
  free( mda);
}


//////////////////////////////////////////////////////////////////////////


struct mda_fileinfo *mda_info_load( FILE *fptr)
{
  struct mda_fileinfo *fileinfo;

#ifndef XDR_HACK
  XDR xdrs;
#endif
  XDR *xdrstream;

  int32_t  last_point;
  int32_t first_offset;
  char *time;

  int ver;
  int i, j;

  int32_t t;
  
 
  rewind( fptr);

#ifdef XDR_HACK
  xdrstream = fptr;
#else
  xdrstream = &xdrs;
  xdrstdio_create(xdrstream, fptr, XDR_DECODE);
#endif

  fileinfo = calloc( 1, sizeof(struct mda_fileinfo));
  if( fileinfo == NULL)
    goto Info_Error;

  if( !xdr_float(xdrstream, &(fileinfo->version)) )
    goto Info_Error;
  // because version floats weren't precisely the correct values
  ver = (int) (100.0*fileinfo->version + 0.5);
  if( ((ver != 120) && (ver != 130) && (ver != 140)) ||
      !xdr_int32_t(xdrstream, &(fileinfo->scan_number)) ||
      !xdr_int16_t(xdrstream, &(fileinfo->data_rank)) ||
      (fileinfo->data_rank <= 0) )
    goto Info_Error;

  fileinfo->dimensions = malloc( fileinfo->data_rank * sizeof(int32_t));
  if( (fileinfo->dimensions == NULL) ||
      !xdr_vector( xdrstream, (char *) fileinfo->dimensions, 
                   fileinfo->data_rank, sizeof( int32_t), 
                   (xdrproc_t) xdr_int32_t))
    goto Info_Error;

  for( i = 0; i < fileinfo->data_rank; i++)
    // -1 is it was int16_t, not int32_t
    if( fileinfo->dimensions[i] < 1)
      goto Info_Error;

  if( !xdr_int16_t(xdrstream, &(fileinfo->regular)) ||
      !xdr_int32_t(xdrstream, &t )) 
    goto Info_Error;

  // This double pointer business is a bit of overkill, but to be consistent, 
  // I'll do it here as well
  fileinfo->scaninfos = 
    calloc( fileinfo->data_rank, sizeof(struct mda_scaninfo *));
  if( fileinfo->scaninfos == NULL)
        goto Info_Error;


  for( i = 0; i < fileinfo->data_rank; i++)
    {
      fileinfo->scaninfos[i] = calloc( 1, sizeof(struct mda_scaninfo ));

      if( !xdr_int16_t(xdrstream, &(fileinfo->scaninfos[i]->scan_rank)) ||
          (fileinfo->scaninfos[i]->scan_rank != (fileinfo->data_rank - i)) ||
          !xdr_int32_t(xdrstream, 
                       &(fileinfo->scaninfos[i]->requested_points)) ||
          (fileinfo->scaninfos[i]->requested_points < 1) ||
          !xdr_int32_t(xdrstream, &last_point ) ||
          (last_point < 0) || 
          (last_point > fileinfo->scaninfos[i]->requested_points) )
        goto Info_Error;

      if( fileinfo->scaninfos[i]->scan_rank > 1)
	{
          int32_t *offsets;

	  offsets = malloc( fileinfo->scaninfos[i]->requested_points 
				     * sizeof(int32_t));
	  if( !xdr_vector( xdrstream, (char *) offsets, 
			   fileinfo->scaninfos[i]->requested_points, 
			   sizeof( int32_t), (xdrproc_t) xdr_int32_t) ||
              (offsets[0] == 0) )
            goto Info_Error;
          first_offset = offsets[0];
          free(offsets);
	}
      else
        first_offset = 0;

      if( !xdr_counted_string( xdrstream, &(fileinfo->scaninfos[i]->name)) ||
          !xdr_counted_string( xdrstream, &time) ||
          !xdr_int16_t(xdrstream, 
                       &(fileinfo->scaninfos[i]->number_positioners)) ||
          (fileinfo->scaninfos[i]->number_positioners < 0) ||
          !xdr_int16_t(xdrstream, 
                       &(fileinfo->scaninfos[i]->number_detectors)) ||
          (fileinfo->scaninfos[i]->number_detectors < 0) ||
          !xdr_int16_t(xdrstream, &(fileinfo->scaninfos[i]->number_triggers)) ||
          (fileinfo->scaninfos[i]->number_triggers < 0) )
        goto Info_Error;

      // only want this stuff for outer loop
      if( !i)
	{
	  fileinfo->last_topdim_point = last_point;
	  fileinfo->time = time;
	}
      else
	free(time);
     
      if( (fileinfo->scaninfos[i]->positioners =
           calloc( fileinfo->scaninfos[i]->number_positioners,
                   sizeof(struct mda_positioner *))) == NULL)
        goto Info_Error;
      for( j = 0; j < fileinfo->scaninfos[i]->number_positioners; j++)
	{
	  if( (fileinfo->scaninfos[i]->positioners[j] = 
	       positioner_read(xdrstream)) == NULL )
            goto Info_Error;
	}

      if( (fileinfo->scaninfos[i]->detectors =
           calloc( fileinfo->scaninfos[i]->number_detectors,
                   sizeof(struct mda_detector *))) == NULL)
        goto Info_Error;
      for( j = 0; j < fileinfo->scaninfos[i]->number_detectors; j++)
	{
	  if( (fileinfo->scaninfos[i]->detectors[j] = 
	       detector_read( xdrstream)) == NULL )
            goto Info_Error;
	}

      if( (fileinfo->scaninfos[i]->triggers =
           malloc( fileinfo->scaninfos[i]->number_triggers * 
                   sizeof(struct mda_trigger *))) == NULL)
        goto Info_Error;
      for( j = 0; j < fileinfo->scaninfos[i]->number_triggers; j++)
	{
	  if( (fileinfo->scaninfos[i]->triggers[j] = 
	       trigger_read( xdrstream)) == NULL )
            goto Info_Error;
	}

      if( first_offset != 0)
	{
          if( !xdr_setpos( xdrstream, first_offset) )
            goto Info_Error;
	}
    }

#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif

  return fileinfo;

 Info_Error:
#ifndef XDR_HACK
  xdr_destroy( xdrstream);
#endif
  mda_info_unload( fileinfo);
  return NULL;
}


void mda_info_unload( struct mda_fileinfo *fileinfo)
{
  int i, j;

  struct mda_scaninfo *scaninfo;

  if( fileinfo == NULL)
    return;

  if( fileinfo->scaninfos != NULL)
    for( j = 0; j < fileinfo->data_rank; j++)
      {
        scaninfo = fileinfo->scaninfos[j];
        if( scaninfo == NULL)
          continue;

        free( scaninfo->name);
    
        if( scaninfo->positioners != NULL)
          {
            for( i = 0; i < scaninfo->number_positioners; i++)
              {
                if( scaninfo->positioners[i] == NULL)
                  continue;
                free(scaninfo->positioners[i]->name);
                free(scaninfo->positioners[i]->description);
                free(scaninfo->positioners[i]->step_mode);
                free(scaninfo->positioners[i]->unit);
                free(scaninfo->positioners[i]->readback_name);
                free(scaninfo->positioners[i]->readback_description);
                free(scaninfo->positioners[i]->readback_unit);
                free(scaninfo->positioners[i]);
              }
            free( scaninfo->positioners);
          }

        if( scaninfo->detectors != NULL)
          {
            for( i = 0; i < scaninfo->number_detectors; i++)
              {
                if( scaninfo->detectors[i] == NULL)
                  continue;
                free(scaninfo->detectors[i]->name);
                free(scaninfo->detectors[i]->description);
                free(scaninfo->detectors[i]->unit);
                free(scaninfo->detectors[i]);
              }
            free( scaninfo->detectors);
          }

        if( scaninfo->triggers != NULL)
          {
            for( i = 0; i < scaninfo->number_triggers; i++)
              {
                if( scaninfo->triggers[i] == NULL)
                  continue;
                free(scaninfo->triggers[i]->name);
                free(scaninfo->triggers[i]);
              }      
            free( scaninfo->triggers);
          }

        free( scaninfo);
      }
  free( fileinfo->scaninfos);

  free( fileinfo->dimensions);
  free( fileinfo->time);
  
  free( fileinfo);
}
