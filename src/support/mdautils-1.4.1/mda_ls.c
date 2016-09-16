/*************************************************************************\
* Copyright (c) 2016 UChicago Argonne, LLC,
*               as Operator of Argonne National Laboratory.
* This file is distributed subject to a Software License Agreement
* found in file LICENSE that is included with this distribution. 
\*************************************************************************/


/*
  Written by Dohn A. Arms, Argonne National Laboratory
  Send comments to dohnarms@anl.gov
  
  0.1.0 -- May 2009
  1.0.0 -- October 2009
           Added Search capabilities
  1.0.1 -- November 2009
           Redid directory scanning code to not use scandir()
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
           Used printf better, removed formatting strings
  1.3.1 -- February 2014
           If there is an unopenable file, it's ignored instead of halting.
  1.4.0 -- July 2016
           New version of load library is used, with better error checking.
           By default, try to load the data file completely using mda_test(),
           in order to find any data error not found by info.
           This can be turned off with new -s switch to speed it up.
 */



#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <unistd.h>
#include <dirent.h>


#include "mda-load.h"


#define VERSION "1.4.1 (August 2016)"
#define YEAR "2016"

// this function relies too much on the input format not changing
void time_reformat( char *original, char *new)
{
  int i;

  char *months[12] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };


  for( i = 0; i < 4; i++)
    new[i] = original[8+i];

  new[4] = '-';

  for( i = 0; i < 12; i++)
    {
      if( !strncasecmp( months[i], original, 3) )
        {
          new[5] = '0' + (i+1)/10;
          new[6] = '0' + (i+1)%10;
          break;
        }
    }
  if( i == 12)
    {
      new[0] = '\0';
      return;
    }

  new[7] = '-';

  for( i = 0; i < 2; i++)
    new[8+i] = original[4+i];

  new[10] = ' ';

  strncpy( &new[11], &original[13], 8);
  
  new[19] = '\0';
}




void helper(void)
{
  printf( "Usage: mda-ls [-hvfs] [-p POSITIONER] [-d DETECTOR] [-t TRIGGER] [DIRECTORY]\n"
          "Shows scan information for all MDA files in a directory.\n"
          "\n"
          "-h  This help text.\n"
          "-v  Show version information.\n"
	  "-f  Show a fuller listing, including the time of the scan start.\n"
          "-s  Only do simple data file checks, instead of a full load test.\n"
          "-p  Include only listings with POSITIONER as a scan positioner.\n"
          "-d  Include only listings with DETECTOR as a scan detector.\n"
          "-t  Include only listings with TRIGGER as a scan trigger.\n"
          "\n"
          "The directory searched is specified with DIRECTORY else the current directory\n"
          "is used. The format of the listing is: filename, [date and time,] scan size,\n"
          "and positioners. The scansize shows the number of points for each dimension;\n"
          "an incomplete highest level scan has the intended number shown in parentheses.\n"
          "The positioner entries show: scan level, positioner PV, and description.\n"
          );

}


void version(void)
{
  printf( "mda-ls %s\n"
          "\n"
          "Copyright (c) %s UChicago Argonne, LLC,\n"
          "as Operator of Argonne National Laboratory.\n"
          "\n"
          "Written by Dohn Arms, dohnarms@anl.gov.\n", VERSION, YEAR);
}


int filter(const struct dirent *a)
{
  int len;

  len = strlen( a->d_name);

  if( len < 5)
    return 0;
  if( strcmp( ".mda", &a->d_name[len-4]) )
    return 0;

  return 1;
}

int sort(const char **a, const char **b)
{
  return(strcmp(*a, *b));
}

int mda_dir(const char *dirname, char ***namelist )
{
  // shove entries on linked list stack, and pop them back onto array
  struct link
  {
    char *name;
    struct link *next;
  };

  struct link *head, *entry;
  int count;
  DIR *dir;
  struct dirent *direntry;
  int i;

  head = NULL;

  if ((dir = opendir(dirname)) == NULL)
     return -1;
  count = 0;
  while( (direntry = readdir(dir)) != NULL)
    {
      if( !filter(direntry)) 
        continue;
      entry = (struct link *) malloc( sizeof( struct link));
      if( entry == NULL)
        return -1;

      entry->name = strdup( direntry->d_name);
      if( entry->name == NULL)
        return -1;
      entry->next = head;
      head = entry;
      count++;
    }
  if( closedir(dir)) 
    return -1;
  if( !count) 
    return 0;

  *namelist = (char **) malloc( count * sizeof( char *) );
  if( *namelist == NULL)
    return -1;
  for( i = count - 1; i >= 0 ; i--)
    { 
      entry = head;
      (*namelist)[i] = entry->name;
      head = entry->next;
      free( entry);
    }

  qsort((void *) (*namelist), (size_t) count, sizeof( char *), 
        (const void *) &sort);
  
  return count;
}


int main( int argc, char *argv[])
{
  char **filelist;
  int dir_number;

  struct mda_fileinfo **fileinfos;

  // these are used to reduce pointer dereferencing
  struct mda_fileinfo   *finf;
  struct mda_scaninfo   *scinf;
  struct mda_positioner *pos;  

  FILE *fptr;

  char *dir;


#define STRING_SIZE (4096)

  char string[STRING_SIZE];

  int max_namelen, max_dimlen, max_timelen, offset;

  int full_flag = 0;
  int check_flag = 1;
  int search_flag = 0;
  int positioner_flag = 0;
  char *positioner_term = NULL;
  int detector_flag = 0;
  char *detector_term = NULL;
  int trigger_flag = 0;
  char *trigger_term = NULL;

  // values are to shut up gcc
  int   allow_count = 0;
  char *allow_list = NULL;

  int opt;

  int skip_line = 0;

  int i, j, k, m, n;


  while((opt = getopt( argc, argv, "hvfsp:d:t:")) != -1)
    {
      switch(opt)
        {
        case 'h':
          helper();
          return 0;
          break;
        case 'v':
          version();
          return 0;
          break;
        case 'f':
          full_flag = 1;
          break;
        case 's':
          check_flag = 0;
          break;
	case 'p':
	  positioner_flag = 1;
          search_flag = 1;
	  positioner_term = strdup( optarg);
	  break;
	case 'd':
	  detector_flag = 1;
          search_flag = 1;
	  detector_term = strdup( optarg);
	  break;
	case 't':
	  trigger_flag = 1;
          search_flag = 1;
	  trigger_term = strdup( optarg);
	  break;
        case ':':
          // option normally resides in 'optarg'
          printf("Error: option missing its value!\n");  
          return -1;
          break;
        }
    }

  if( (argc - optind) > 1)
    {
      helper();
      return 1;
    }

  if( (argc - optind) == 0)
    {
      dir = NULL;
      dir_number = mda_dir( ".", &filelist);
    }
  else
    {
      dir = strdup( argv[argc - 1]);
      dir_number = mda_dir( dir, &filelist);
    }

  if (dir_number <= 0)
    {
      printf("No MDA files in the directory.\n");
      return 1;
    }

  if( (dir != NULL) &&  chdir( dir) )
    {
      printf("Can't change to that directory.\n");
      return 1;
    }

  allow_count = 0;
  allow_list = (char *) malloc( dir_number * sizeof(char) );
  for( i = 0; i < dir_number; i++)
    allow_list[i] = 0;
  
  fileinfos = (struct mda_fileinfo **) 
    malloc( dir_number * sizeof(struct mda_fileinfo *) );
  for( i = 0; i < dir_number; i++)
    {
      if( (fptr = fopen( filelist[i], "rb")) == NULL)
        {
          printf("Can't open file \"%s\", skipping.\n", filelist[i] );
          fileinfos[i] = NULL;
          skip_line = 1;
        }
      else
        {
          if( check_flag && mda_test( fptr) )
            fileinfos[i] = NULL;
          else
            fileinfos[i] = mda_info_load( fptr);
          fclose(fptr);
          allow_list[i] = 1;
          allow_count++;
        }
    }

  if( skip_line)
    printf("\n");

  if( search_flag)
    {
      // resets the search, as the first one filtered out unreadable files.
      allow_count = 0;

      for( i = 0; i < dir_number; i++)
        {
          if( !allow_list[i])
            continue;
          else
            allow_list[i] = 0;

          for( k = 0; k < fileinfos[i]->data_rank; k++)
            {
              scinf = fileinfos[i]->scaninfos[k];
              
              if( positioner_flag && scinf->number_positioners)
                for( j = 0; j < scinf->number_positioners; j++)
                  if( !strcmp( positioner_term, (scinf->positioners[j])->name ))
                    {
                      allow_list[i] = 1;
                      allow_count++;
                      goto Shortcut;
                    }
              if( detector_flag && scinf->number_detectors)
                for( j = 0; j < scinf->number_detectors; j++)
                  if( !strcmp( detector_term, (scinf->detectors[j])->name ))
                    {
                      allow_list[i] = 1;
                      allow_count++;
                      goto Shortcut;
                    }
              if( trigger_flag && scinf->number_triggers)
                for( j = 0; j < scinf->number_triggers; j++)
                  if( !strcmp( trigger_term, (scinf->triggers[j])->name ))
                    {
                      allow_list[i] = 1;
                      allow_count++;
                      goto Shortcut;
                    }
            }
        Shortcut:
          ;
        }
    }
      
  if( !allow_count)
    {
      printf("No MDA files in the directory fit your search.\n");
      return 1;
    }

  max_namelen = 0;
  max_dimlen = 0;
  for( i = 0; i < dir_number; i++)
    {
      if( !allow_list[i] )
        continue;

      j = strlen( filelist[i]);
      if( j > max_namelen)
        max_namelen = j;

      finf = fileinfos[i];

      if( finf == NULL)
        {
          //          printf("Invalid\n");
          if( 7 > max_dimlen)
            max_dimlen = 7;
        }
      else
        {
          if( finf->dimensions[0] != finf->last_topdim_point)
            j = snprintf( string, STRING_SIZE, "%i(%i)", 
                          finf->last_topdim_point, finf->dimensions[0]);
          else
            j = snprintf( string, STRING_SIZE, "%i", finf->dimensions[0]);

          for( k = 1; k < finf->data_rank; k++)
            j += snprintf( &string[j], STRING_SIZE - j, "x%i", 
                           finf->dimensions[k]);

          if( j > max_dimlen)
            max_dimlen = j;
        }

    }

  max_timelen = 19;

  offset = max_namelen + max_dimlen + 4;
  if( full_flag)
    offset += (2 + max_timelen);

  for( i = 0; i < dir_number; i++)
    {
      if( !allow_list[i] )
        continue;

      printf("%*s  ", max_namelen, filelist[i]);

      finf = fileinfos[i];

      if( finf == NULL)
        {
          printf("Invalid\n");
          continue;
        }

      if( full_flag)
        {
          time_reformat( finf->time, string);
          printf( "%s  ", string );
        }

      if( finf->dimensions[0] != finf->last_topdim_point)
        j = snprintf( string, STRING_SIZE, "%i(%i)", 
                      finf->last_topdim_point, finf->dimensions[0]);
      else
        j = snprintf( string, STRING_SIZE, "%i", finf->dimensions[0]);
      for( k = 1; k < finf->data_rank; k++)
        j += snprintf( &string[j], STRING_SIZE - j, "x%i", 
                       finf->dimensions[k]);
      printf("%-*s -", max_dimlen, string);

      m = 0;
      for( k = 0; k < finf->data_rank; k++)
        {
          scinf = finf->scaninfos[k];
          if( scinf->number_positioners)
            for( j = 0; j < scinf->number_positioners; j++)
              {
                pos = scinf->positioners[j];
                if( m)
                  for( n = 0; n < offset; n++)
                    printf(" ");
                printf(" %dD %s", scinf->scan_rank, pos->name );
                if( pos->description[0] != '\0' )
                  printf(" (%s)", pos->description );
                printf("\n");
                m = 1;
              }
        }
      if( !m)
        printf("\n");
    }

  for( i = 0; i < dir_number; i++)
    {
      if( fileinfos[i] != NULL) 
        mda_info_unload(fileinfos[i]);
    }
  free( filelist);
  free( fileinfos);

  if( search_flag)
    free(allow_list );

  return 0;
}


