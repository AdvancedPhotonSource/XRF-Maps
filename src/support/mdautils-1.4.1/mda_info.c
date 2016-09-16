/*************************************************************************\
* Copyright (c) 2016 UChicago Argonne, LLC,
*               as Operator of Argonne National Laboratory.
* This file is distributed subject to a Software License Agreement
* found in file LICENSE that is included with this distribution. 
\*************************************************************************/


/*
  Written by Dohn A. Arms, Argonne National Laboratory
  Send comments to dohnarms@anl.gov
  
  0.1   -- August 2005
  0.1.1 -- December 2006
           Added support for files that have more than 32k points.
  0.2.0 -- November 2007
           Use new fileinfo routines to get information.
  1.0.0 -- October 2009
           Overhauled the output format, adding information about 
           detectors and triggers.
  1.0.1 -- May 2010
           Added showing numbering of detectors, positioners, 
           and triggers as done by saveData.
  1.1   -- November 2010
  1.1.1 -- March 2011
  1.2   -- March 2011
           Fixed integer issues by tying short to int16_t, long to int32_t,
           and char to int8_t.  Changed %li to %i in printf's.
  1.2.1 -- January 2012
           Cleaned up the overuse of pointer dereferencing, hopefully
           making it faster as well as easier to understand
  1.2.2 -- June 2012
  1.3.0 -- February 2013
  1.3.1 -- February 2014
  1.4.0 -- July 2016
           New version of load library is used, with better error checking.
           By default, try to load the data file completely using mda_test(),
           in order to find any data error not found by info.
           This can be turned off with new -s switch to speed it up.
  1.4.1 -- August 2016
           Changed regularity reporting to simply state as true or false.
 */



/****************  mda_info.c  **********************/


#include <stdio.h>
#include <string.h>

#include <unistd.h>

#include "mda-load.h"


#define VERSION "1.4.1 (August 2016)"
#define YEAR "2016"



/////////////////////////////////////////////////////////////////////////////

int test_print(char *format, char *string)
{
  if( string[0] != '\0')
    return printf( format, string);

  return 0;
}


int information( struct mda_fileinfo *fileinfo)
{
  struct mda_scaninfo   *scinf;
  struct mda_positioner *pos;  
  struct mda_detector   *det;  
  struct mda_trigger    *trig;  

  int i, j;


  printf("MDA file version = %g\n", fileinfo->version);
  printf("     Scan number = %i\n", fileinfo->scan_number);
  printf("  Dimensionality = %i\n", fileinfo->data_rank);
  printf("       Scan size = ");
  if( fileinfo->last_topdim_point != fileinfo->dimensions[0] )
    printf("(%i)", fileinfo->last_topdim_point );
  for( i = 0; i < fileinfo->data_rank; i++)
    {
      if( i)
	printf("x");
      printf( "%i", fileinfo->dimensions[i]);
    }
  printf("\n");
  printf(" Scan start time = %s\n", fileinfo->time);
  if( fileinfo->regular)
    printf("      Regularity = TRUE\n"); 
  else
    printf("      Regularity = FALSE\n"); 

  for( i = 0; i < fileinfo->data_rank; i++)
    {
      scinf = fileinfo->scaninfos[i];

      printf("\n%i-D scanner name = %s\n", scinf->scan_rank, scinf->name);
      printf("%i-D triggers:%s\n", scinf->scan_rank,
             !scinf->number_triggers ? "  None" : "" );
      for( j = 0; j < scinf->number_triggers; j++)
	{
          trig = scinf->triggers[j];
	  printf( "    %i (T%d) = %s\n", j + 1, trig->number + 1, trig->name);
        }
      printf("%i-D positioners:%s\n", scinf->scan_rank,
             !scinf->number_positioners ? "  None" : "");
      for( j = 0; j < scinf->number_positioners; j++)
	{
          pos = scinf->positioners[j];

	  printf("    %i (P%i) = %s", j + 1, pos->number + 1, pos->name);
          test_print(" (%s)", pos->description);
          test_print(" [%s]", pos->unit);
          printf("\n");
          if( pos->readback_name[0] != '\0')
            {
              printf("             %s", pos->readback_name );
              test_print(" (%s)", pos->readback_description);
              test_print(" [%s]", pos->readback_unit);
              printf("\n");
            }
	}
      
      printf("%i-D detectors:%s\n", scinf->scan_rank,
             !scinf->number_detectors ? "  None" : "");
      for( j = 0; j < fileinfo->scaninfos[i]->number_detectors; j++)
	{
          det = scinf->detectors[j];

	  printf("  %02i (D%02i) = %s", j + 1, det->number + 1, det->name);
          test_print(" (%s)", det->description);
          test_print(" [%s]", det->unit);
          printf("\n");
	}
    }
  
  return 0;
}


void help(void)
{
  printf("Usage: mda-info [-hvs] FILE\n"
         "Prints the basic scan information about the EPICS MDA file, FILE.\n"
         "\n"
         "-h  This help text.\n"
         "-v  Show version information.\n"
         "-s  Only do simple data file checks, instead of a full load test.\n"
         "\n"
         "Information such as dimensionality and time of scan start are shown,\n"
         "as well as all the positioners, detectors, and triggers for each dimension.\n"
         "For each trigger, positioner, and detector there is both a number and\n"
         "(in parentheses) the PV field number used in saveData.\n"
         );
}

void version(void)
{
  printf("mda-info %s\n"
         "\n"
         "Copyright (c) %s UChicago Argonne, LLC,\n"
         "as Operator of Argonne National Laboratory.\n"
         "\n"
         "Written by Dohn Arms, dohnarms@anl.gov.\n",
         VERSION, YEAR);
}

int main( int argc, char *argv[])
{
  FILE *input;
  struct mda_fileinfo *fileinfo;

  int opt;
  int check_flag = 1;

  while((opt = getopt( argc, argv, "hvs")) != -1)
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
        case 's':
          check_flag = 0;
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
      printf("For help, type: mda-info -h\n");
      return 0;
    }

  if( (input = fopen( argv[optind], "rb")) == NULL)
    {
      fprintf(stderr, "Can't open file \"%s\"!\n", argv[optind]);
      return 1;
    }

  if( check_flag)
    {
      if( mda_test( input) )
        {
          fprintf(stderr, "Loading file \"%s\" failed!\n", argv[optind]);
          return 1;
        }
    }

  if( (fileinfo = mda_info_load( input)) == NULL )
    {
      fprintf(stderr, "Reading information from file \"%s\" failed!\n", 
              argv[optind]);
      return 1;
    }

  fclose(input);

  information(fileinfo);

  mda_info_unload(fileinfo);

  return 0;
}


