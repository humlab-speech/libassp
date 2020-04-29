/***********************************************************************
*                                                                      *
* This file is part of the Advanced Speech Signal Processor library.   *
*                                                                      *
* Copyright (C) 1989 - 2010  Michel Scheffers                          *
*                            IPdS, CAU Kiel                            *
*                            Leibnizstr. 10                            *
*                            24118 Kiel, Germany                       *
*                            ms@ipds.uni-kiel.de                       *
*                                                                      *
* This library is free software: you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* This library is distributed in the hope that it will be useful,      *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with this library. If not, see <http://www.gnu.org/licenses/>. *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* File:     acfana.c                                                   *
* Contents: Command line interface to autocorrelation analysis.        *
* Author:   Michel T.M. Scheffers                                      *
*                                                                      *
***********************************************************************/
/* $Id: acfana.c,v 1.19 2010/07/14 13:44:30 mtms Exp $ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>    /* [f/s]printf() remove() stderr FILE NULL */
#include <stddef.h>   /* size_t */
#include <stdlib.h>   /* strtod() strtol() */
#include <string.h>   /* str... */
#include <ctype.h>    /* isdigit() isgraph() */

#include <miscdefs.h> /* TRUE FALSE LOCAL */
#include <mylimits.h> /* PATH_MAX NAME_MAX SUFF_MAX */
#include <misc.h>     /* mybasename() parsepath() fgetl() */
#include <acf.h>      /* processing parameters & ACF functions */
#include <assp.h>     /* message & trace handler */
#include <asspfio.h>  /* file handler */
#include <asspdsp.h>  /* wfunc_e wfType() listWFs() wfShortList */
#include <asspana.h>  /* AOPTS */

/*
 * constants
 */
#define MAXFILES 2000 /* multiple input files */

/*
 * global arrays and variables
 */
char  progName[NAME_MAX+1];
char *argList[MAXFILES], *batList, *fixOutFile, *outDir, *acfTrack;
char  smpFile[PATH_MAX+1], acfFile[PATH_MAX+1], outExt[SUFF_MAX+1];
int   numFiles;
int   X_OPTS;
FILE *batFP;
AOPTS anaOpts;

/*
 * prototypes of local functions
 */
LOCAL void initGlobals(void);
LOCAL int  evalArgs(int argc, char *argv[]);
LOCAL int  optError(char *opt);
LOCAL void usage(void);
LOCAL int  nextFile(void);

/**********************************************************************/

int main(int argc, char *argv[])
{
  char  info[128];
  int   err;
  DOBJ *smpDOp, *acfDOp;

  strcpy(progName, mybasename(argv[0]));
  strcpy(traceFile, progName);
  strcat(traceFile, TRACE_SUFFIX);
  initGlobals();
  if((err=evalArgs(argc, argv)) < 0)
    exit(1);
  if(err)
    exit(0);
  if(TRACE['F'] || TRACE['f']) {
    sprintf(info, "%s version %d.%d", progName, ACF_MAJOR, ACF_MINOR);
    openTrace(info);
  }
  else openTrace(NULL);
/*
 * loop over input files
 */
  smpDOp = acfDOp = NULL;
  while((err=nextFile()) > 0) {
    smpDOp = asspFOpen(smpFile, AFO_READ, NULL);
    if(smpDOp == NULL) {
      err = prtAsspMsg(NULL);
      break;
    }
    acfDOp = createACF(smpDOp, &anaOpts);
    if(acfDOp == NULL) {
      err = prtAsspMsg(NULL);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(asspFOpen(acfFile, AFO_WRITE, acfDOp) == NULL) {
      err = prtAsspMsg(NULL);
      acfDOp = freeDObj(acfDOp);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(computeACF(smpDOp, &anaOpts, acfDOp) == NULL) {
      err = prtAsspMsg(NULL);
    }
    asspFClose(smpDOp, AFC_FREE);
    smpDOp = NULL;
    asspFClose(acfDOp, AFC_FREE);
    acfDOp = NULL;
    if(err < 0) {
      remove(acfFile); /* contents invalid */
      break;
    }
    if(asspWarning)
      prtAsspMsg(NULL);
  } /* END loop over input files */
  if(batFP != NULL) {
    fclose(batFP);
    batFP = NULL;
  }
  closeTrace();
  if(err < 0)
    exit(1);
  exit(0);
}
/***********************************************************************
* initialize global parameters, pointers and analysis options          *
***********************************************************************/
LOCAL void initGlobals(void)
{
  KDTAB *entry;

  clrAsspMsg();
  batList = fixOutFile = outDir = acfTrack = NULL;
  strcpy(outExt, ACF_DEF_SUFFIX);
  X_OPTS = FALSE;
  batFP = NULL;
  setACFdefaults(&anaOpts);
  entry = dtype2entry(DT_ACF, KDT_SSFF);
  if(entry != NULL && entry->keyword != NULL)
    acfTrack = entry->keyword;
  else {
    setAsspMsg(AEB_ERR_TRACK, "for data type ACF");
    exit(prtAsspMsg(NULL));
  }
  return;
}
/***********************************************************************
* check program options and arguments                                  *
***********************************************************************/
LOCAL int evalArgs(int argc, char *argv[])
{
  char   *cPtr;
  int     i;
  wfunc_e winFunc;

  numFiles = 0;
  for(i = 1; i < argc; i++) {
    cPtr = argv[i];
    if(*cPtr == '-') {      /* OPTION */
      cPtr++;
      do {
	switch(*cPtr) {
	case '-':           /* long option */
	  cPtr++;
	  switch(*cPtr) {
	  case 'H':                                           /* HELP */
	  case 'h':
	    return(optError(NULL));
	    break;
	  default:
	    return(optError(cPtr-2));
	  }
	  break;            /* end long option */
	case 'h':                                             /* help */
	  return(optError(NULL));
	  break;
	case 'b':                                       /* begin time */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	    cPtr++;
	    anaOpts.beginTime = strtod(cPtr, &cPtr);
	    anaOpts.options &= ~AOPT_USE_CTIME;
	  }
	  else return(optError("-b"));
	  break;
	case 'c':                                      /* centre time */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	    cPtr++;
	    anaOpts.centreTime = strtod(cPtr, &cPtr);
	    anaOpts.options |= AOPT_USE_CTIME;
	  }
	  else return(optError("-c"));
	  break;
	case 'e':                                         /* end time */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	    cPtr++;
	    anaOpts.endTime = strtod(cPtr, &cPtr);
	    anaOpts.options &= ~AOPT_USE_CTIME;
	  }
	  else return(optError("-e"));
	  break;
	case 'm':                                   /* analysis order */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.order = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.order < 1)
	      return(optError("-m"));
	  }
	  else return(optError("-m"));
	  break;
	case 'n':                /* energy-normalized autocorrelation */
	  anaOpts.options |= ACF_OPT_NORM;
	  cPtr++;
	  break;
	case 'N':                /* length-normalized autocorrelation */
	  anaOpts.options |= ACF_OPT_MEAN;
	  cPtr++;
	  break;
	case 's':                                      /* frame shift */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	    cPtr++;
	    anaOpts.msShift = strtod(cPtr, &cPtr);
	  }
	  else return(optError("-s"));
	  break;
	case 'L':                                 /* effective length */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	    cPtr++;
	    anaOpts.msSize = strtod(cPtr, &cPtr);
	    anaOpts.options |= AOPT_EFFECTIVE;
	  }
	  else return(optError("-L"));
	  break;
	case 'S':                                 /* exact frame size */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	    cPtr++;
	    anaOpts.msSize = strtod(cPtr, &cPtr);
	    anaOpts.options &= ~AOPT_EFFECTIVE;
	  }
	  else return(optError("-S"));
	  break;
	case 'C':                                   /* channel number */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.channel = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.channel < 1 || anaOpts.channel > ACF_I_CHANS)
	      return(optError("-C"));
	  }
	  else return(optError("-C"));
	  break;
	case 'o':           /* output modifiers */
	  cPtr++;
	  switch(*cPtr) {
	  case 'A':                                          /* ASCII */
	    X_OPTS = TRUE;                         /* extended option */
	    strcpy(anaOpts.format, ACF_ASC_FORMAT);
	    cPtr++;
	    break;
	  case 'd':                                      /* directory */
	    cPtr++;
	    if(*cPtr == '=') {
	      cPtr++;
	      if(*cPtr != EOS) {
		outDir = cPtr;
		cPtr += strlen(cPtr);
		anaOpts.options &= ~AOPT_IN_DIR;         /* overruled */
	      }
	      else return(optError("-od"));
	    }
	    else
	      anaOpts.options |= AOPT_IN_DIR;  /* use input directory */
	    break;
	  case 'f':                                      /* file name */
	    X_OPTS = TRUE;                         /* extended option */
	    cPtr++;
	    if(*cPtr == '=' && cPtr[1] != EOS) {
	      cPtr++;
	      fixOutFile = cPtr;
	      cPtr += strlen(cPtr);
	    }
	    else return(optError("-of"));
	    break;
	  case 'a':                                       /* accuracy */
	    X_OPTS = TRUE;                         /* extended option */
	    cPtr++;
	    if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	      cPtr++;
	      anaOpts.accuracy = (int)strtol(cPtr, &cPtr, 10);
	    }
	    else return(optError("-oa"));
	    break;
	  case 'x':                                      /* extension */
	    cPtr++;
	    if(*cPtr == '=') {
	      cPtr++;
	      if(*cPtr != EOS) {
		if(*cPtr != '.') {
		  strcpy(outExt, ".");
		  strcat(outExt, cPtr);
		}
		else
		  strcpy(outExt, cPtr);
		cPtr += strlen(cPtr);
	      }
	      else return(optError("-ox"));
	    }
	    else return(optError("-ox"));
	    break;
	  default:
	    return(optError(--cPtr));
	  }
	  break;            /* end output modifiers */
	case 'w':                                  /* window function */
	  cPtr++;
	  if(*cPtr == '=') {
	    cPtr++;
	    winFunc = wfType(cPtr);
	    if(winFunc <= WF_NONE || winFunc >= WF_NUM_FIX)
	      return(optError("-w"));
	    strcpy(anaOpts.winFunc, cPtr);
	    cPtr += strlen(cPtr);
	  }
	  else return(optError("-w"));
	  break;
	case 'W':                               /* window information */
	  listWFs(wfShortList, NULL);
	  return(1);
	  break;
	case 'z':                                       /* batch mode */
	  cPtr++;
	  if(*cPtr == '=' && cPtr[1] != EOS) {
	    anaOpts.options |= AOPT_BATCH;
	    cPtr++;
	    batList = cPtr;
	    cPtr += strlen(cPtr);
	  }
	  else return(optError("-z"));
	  break;
	case 'X':                            /* show extended options */
	  X_OPTS = TRUE;
	  return(optError(NULL));
	  break;
	case 'Y':                                    /* trace options */
	  TRACE[0] = TRUE;
	  cPtr++;
	  while(*cPtr != EOS && isgraph((int)*cPtr)) {
	    TRACE[(int)*cPtr] = TRUE;
	    cPtr++;
	  }
	  break;
	default:                                         /* not found */
	  return(optError(cPtr));
        }
      } while(*cPtr != EOS);
    }                       /* END OPTION */
    else {                  /* ARGUMENT */
      if(numFiles < MAXFILES) {
	argList[numFiles] = cPtr;
	numFiles++;
      }
      else {
	fprintf(stderr, "\nERROR: more than %d input files\n", MAXFILES);
	return(-1);
      }
    }                       /* END ARGUMENT */
  }

  if(numFiles == 0 && !(anaOpts.options & AOPT_BATCH))  /* naked call */
    return(optError(NULL));
  if(anaOpts.options & AOPT_BATCH) {
    numFiles = 0;                                        /* overruled */
    fixOutFile = NULL;                                   /* overruled */
  }
  else if(numFiles > 1)
    fixOutFile = NULL;                                   /* overruled */
  if(fixOutFile != NULL) {
    outDir = NULL;                                       /* overruled */
    anaOpts.options &= ~AOPT_IN_DIR;                     /* overruled */
  }
  return(0);
}
/***********************************************************************
* print option error message                                           *
***********************************************************************/
LOCAL int optError(char *opt)
{ 
  usage();
  if(opt == NULL || *opt == EOS) return(1);
  if(*opt == '-')
    fprintf(stderr, "ERROR: %s option: specification missing or unacceptable\n", opt);
  else
    fprintf(stderr, "ERROR: unknown option '-%s'\n", opt);
  return(-1);
}
/***********************************************************************
* print program syntax                                                 *
***********************************************************************/
LOCAL void usage(void)
{
  printf("\n");
  printf("Syntax  : %s [<opts>] <file> [<opts>] {<file> [<opts>]}\n",\
	 progName);
  printf("Release : %d.%d (%s)\n", ACF_MAJOR, ACF_MINOR, __DATE__);
  printf("Function: Analysis of short-term autocorrelation function of\n");
  printf("          the signal in <file>.\n");
  printf("          Analysis results will be written to a file with the\n");
  printf("          base name of the input file and extension '%s'.\n",\
	 ACF_DEF_SUFFIX);
  printf("          Default output is in SSFF binary format (track '%s').\n",\
	 acfTrack);
  printf("Options:\n");
  printf(" -h/--help  print this text\n");
  printf(" -b=<time>  set begin of analysis interval to <time> seconds\n");
  printf("            (default: begin of file)\n");
  printf(" -e=<time>  set end of analysis interval to <time> seconds\n");
  printf("            (default: end of file)\n");
  printf(" -s=<dur>   set analysis window shift to <dur> ms\n");
  printf("            (default: %.1f)\n", ACF_DEF_SHIFT);
  printf(" -L=<dur>   set effective length of analysis window to <dur> ms\n");
  printf("            (default: %.1f)\n", ACF_DEF_SIZE);
  printf(" -m=<num>   set analysis order to <num>\n");
  printf("            (default: %s)\n", DFLT_ORDER_STR);
  printf(" -n         calculate energy-normalized autocorrelation\n");
  printf(" -N         calculate length-normalized autocorrelation\n");
  printf(" -w=<type>  set analysis window function to <type>\n");
  printf("            (default: %s)\n", ACF_DEF_WINDOW);
  printf(" -W         print summary of window functions\n");
  printf(" -od        store output files in directory of input files\n");
  printf("            (default: current working directory)\n");
  printf(" -od=<dir>  store output files in directory <dir>\n");
  printf(" -ox=<ext>  set extension of output file names to <ext>\n");
  printf(" -z=<list>  batch mode: take input file names from file <list>;\n");
  printf("            overrules arguments <file>\n");
  printf(" -X         show extended options\n");
if(TRACE[0] || X_OPTS) {
  printf("         extended options\n");
  printf(" -C=<num>   for multi-channel input files: analyse channel <num>\n");
  printf("            with 1 <= <num> <= %d  (default: channel %d)\n",\
	 ACF_I_CHANS, ACF_DEF_CHANNEL);
  printf(" -c=<time>  set single-frame analysis with the analysis window\n");
  printf("            centred at <time> seconds; overrules -b, -e and -s\n");
  printf("            options\n");
  printf(" -S=<dur>   set exact analysis window size to <dur> ms;\n");
  printf("            overrules -L option\n");
  printf(" -oA        output in plain ASCII format\n");
  printf(" -oa=<num>  set accuracy of ASCII output to <num> digits\n");
  printf("            (default: %d)\n", ACF_DEF_DIGITS);
  printf(" -of=<file> set name of output file to <file>\n");
  printf("            (only in single-file mode; overrules -od option)\n");
}
if(TRACE[0]) {
  printf("         debugging options\n");
  printf(" -Y<abc>    output trace info to 'stderr' for options <abc>\n");
  printf("            the following options are available:\n");
  printf("   'F'   append trace output to file '%s'\n", traceFile);
  printf("   'f'   write trace output to file '%s'\n", traceFile);
  printf("   'A'   analysis parameters\n");
}
  printf("NOTE:\n");
  printf(" o   Existing output files will be overwritten without notice.\n");
  printf(" o   The analysis interval will be rounded to the nearest frame\n");
  printf("     boundary (integral multiple of the window shift).\n");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf("This program is free software under the GNU General Public License.\n");
  printf("It comes with ABSOLUTELY NO WARRANTY.\n");
  printf("For details see the file COPYING provided with this program or visit\n");
  printf("http://www.gnu.org/licenses/\n");
  printf("\n");
  return;
}
/***********************************************************************
* get next input file name; if there is one: construct path to output  *
* file and return 1, otherwise return 0; return -1 upon error          *
***********************************************************************/
LOCAL int nextFile(void)
{
  LOCAL int fileNum=0;
  char *dPtr, *bPtr;
  int   n;
  
  clrAsspMsg();
  if(anaOpts.options & AOPT_BATCH) {
    if(batFP == NULL) {
      if((batFP=fopen(batList, "r")) == NULL) {
	setAsspMsg(AEF_ERR_OPEN, batList);
	return(prtAsspMsg(NULL));                        /* break off */
      }
    }
    while((n=fgetl(smpFile, PATH_MAX + 1, batFP, NULL)) == 0) NIX;
    if(n == EOF) {                                   /* we're through */
      fclose(batFP);
      batFP = NULL;
      return(0);
    }
  }
  else {
    if(fileNum >= numFiles)
      return(0);                                     /* we're through */
    strcpy(smpFile, argList[fileNum]);
    fileNum++;
  }
  /*
   * construct name of output file
   */
  if(fixOutFile != NULL) {
    strcpy(acfFile, fixOutFile);
    fixOutFile = NULL;
  }
  else {
    parsepath(smpFile, &dPtr, &bPtr, NULL);
    if(anaOpts.options & AOPT_IN_DIR)
      strcpy(acfFile, dPtr);       /* store output in input directory */
    else if(outDir != NULL) {
      strcpy(acfFile, outDir);
      n = strlen(outDir);
      if(n > 0 && outDir[n-1] != DIR_SEP_CHR)
	strcat(acfFile, DIR_SEP_STR);
    }
    else strcpy(acfFile, "");           /* store in current directory */
    strcat(acfFile, bPtr);   /* append base name */
    strcat(acfFile, outExt); /* append extension */
  }
  if(strcmp(smpFile, acfFile) == 0) {
    setAsspMsg(AEC_IO_CLASH, acfFile);
    return(prtAsspMsg(NULL));
  }
  if(TRACE['F'] || TRACE['f'])
    fprintf(traceFP, "Input file: %s\nOutput file: %s\n",\
	    smpFile, acfFile);
  return(1); /* we got one */
}
