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
* File:     rfcana.c                                                   *
* Contents: Command line interface to Linear Prediction analysis.      *
* Author:   Michel T.M. Scheffers                                      *
*                                                                      *
* ToDo: add calculation of residual signal                             *
*                                                                      *
***********************************************************************/
/* $Id: rfcana.c,v 1.14 2010/07/14 13:44:30 mtms Exp $ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>    /* [f/s]printf() remove() stderr FILE NULL */
#include <stdlib.h>   /* strtod() strtol() */
#include <stddef.h>   /* size_t */
#include <string.h>   /* str..() */
#include <ctype.h>    /* isspace() isdigit() isgraph() */

#include <miscdefs.h> /* TRUE FALSE LOCAL */
#include <mylimits.h> /* PATH_MAX NAME_MAX SUFF_MAX */
#include <misc.h>     /* mybasename() parsepath() fgetl() */
#include <rfc.h>      /* processing parameters & LP functions */
#include <assp.h>     /* message & trace handler */
#include <asspfio.h>  /* file handler */
#include <asspdsp.h>  /* wfunc_e wfType() wfShortList listWFs() */
#include <asspana.h>  /* AOPTS */

/*
 * constants
 */
#define MAXFILES 2000 /* multiple input files */

/*
 * global arrays and variables
 */
char  progName[NAME_MAX+1];
char *argList[MAXFILES], *batList, *fixOutFile, *outDir;
char *rmsTrack, *gainTrack, *lpTrack;
char  smpFile[PATH_MAX+1], lpFile[PATH_MAX+1], outExt[SUFF_MAX+1];
int   numFiles, X_OPTS;
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
  DOBJ *smpDOp, *lpDOp;

  strcpy(progName, mybasename(argv[0]));
  strcpy(traceFile, progName);
  strcat(traceFile, TRACE_SUFFIX);
  initGlobals();
  if((err=evalArgs(argc, argv)) < 0)
    exit(1);
  if(err != 0)
    exit(0);
  if(TRACE['F'] || TRACE['f']) {
    sprintf(info, "%s version %d.%d", progName, RFC_MAJOR, RFC_MINOR);
    openTrace(info);
  }
  else openTrace(NULL);
/*
 * loop over input files
 */
  smpDOp = lpDOp = NULL;
  while((err=nextFile()) > 0) {
    smpDOp = asspFOpen(smpFile, AFO_READ, NULL);
    if(smpDOp == NULL) {
      err = prtAsspMsg(NULL);
      break;
    }
    lpDOp = createLP(smpDOp, &anaOpts);
    if(lpDOp == NULL) {
      err = prtAsspMsg(NULL);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(asspFOpen(lpFile, AFO_WRITE, lpDOp) == NULL) {
      err = prtAsspMsg(NULL);
      lpDOp = freeDObj(lpDOp);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(computeLP(smpDOp, &anaOpts, lpDOp) == NULL) {
      err = prtAsspMsg(NULL);
    }
    asspFClose(smpDOp, AFC_FREE);
    smpDOp = NULL;
    asspFClose(lpDOp, AFC_FREE);
    lpDOp = NULL;
    if(err < 0) {
      remove(lpFile); /* contents invalid */
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
  batList = fixOutFile = outDir = NULL;
  rmsTrack = gainTrack = lpTrack = NULL;
  strcpy(outExt, "");
  X_OPTS = FALSE;
  batFP = NULL;
  setLPdefaults(&anaOpts);
  entry = dtype2entry(DT_RMS, KDT_SSFF);
  if(entry != NULL && entry->keyword != NULL)
    rmsTrack = entry->keyword;
  else {
    setAsspMsg(AEB_ERR_TRACK, "for data type RMS");
    exit(prtAsspMsg(NULL));
  }
  entry = dtype2entry(DT_GAIN, KDT_SSFF);
  if(entry != NULL && entry->keyword != NULL)
    gainTrack = entry->keyword;
  else {
    setAsspMsg(AEB_ERR_TRACK, "for data type GAIN");
    exit(prtAsspMsg(NULL));
  }
  entry = keyword2entry(anaOpts.type, KDT_SSFF);
  if(entry != NULL && entry->keyword != NULL)
    lpTrack = entry->keyword;
  else {
    setAsspMsg(AEB_ERR_TRACK, NULL);
    sprintf(applMessage, "for data type %s", anaOpts.type);
    exit(prtAsspMsg(NULL));
  }
  return;
}
/***********************************************************************
* check program options and arguments                                  *
***********************************************************************/
LOCAL int evalArgs(int argc, char *argv[])
{
  char    *cPtr;
  int      i;
  wfunc_e  winFunc;
  LP_TYPE *lPtr;
  KDTAB   *entry;
  
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
	  break;           /* end long option */
	case 'h':                                             /* help */
	  return(optError(NULL));
	  break;
	case 'b':                                       /* begin time */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.' )) {
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
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.' )) {
	    cPtr++;
	    anaOpts.endTime = strtod(cPtr, &cPtr);
	    anaOpts.options &= ~AOPT_USE_CTIME;
	  }
	  else return(optError("-e"));
	  break;
	case 'C':                                   /* channel number */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.channel = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.channel < 1 || anaOpts.channel > LP_I_CHANS)
	      return(optError("-C"));
	  }
	  else return(optError("-C"));
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
	case 'm':                                         /* LP order */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.order = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.order < 1)
	      return(optError("-m"));
	  }
	  else return(optError("-m"));
	  break;
	case 'o':           /* output modifiers */
	  cPtr++;
	  switch(*cPtr) {
	  case 'A':                                          /* ASCII */
	    X_OPTS = TRUE;                         /* extended option */
	    strcpy(anaOpts.format, LP_ASC_FORMAT);
	    cPtr++;
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
	  case 'p':                                      /* precision */
	    X_OPTS = TRUE;                         /* extended option */
	    cPtr++;
	    if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	      cPtr++;
	      anaOpts.precision = (int)strtol(cPtr, &cPtr, 10);
	    }
	    else return(optError("-op"));
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
	case 'p':                                     /* pre-emphasis */
	  cPtr++;
	  if(*cPtr == '=') {
	    cPtr++;
	    if(isdigit((int)*cPtr) || *cPtr == '.' || *cPtr == '-') {
	      anaOpts.preEmph = strtod(cPtr, &cPtr);
	      if(anaOpts.preEmph < -1.0 || anaOpts.preEmph > 1.0)
		return(optError("-p"));
	    }
	    else return(optError("-p"));
	  }
	  else return(optError("-p"));
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
	case 't':                                   /* parameter type */
	  cPtr++;
	  if(*cPtr == '=') {
	    cPtr++;
	    for(lPtr = lpType; lPtr->ident != NULL; lPtr++) {
	      if(strnxcmp(cPtr, lPtr->ident, 2) == 0)
		break;
	    }
	    if(lPtr->ident == NULL)
	      return(optError("-t"));
	    strcpy(anaOpts.type, lPtr->ident);
	    cPtr += strlen(cPtr);
	  }
	  else return(optError("-t"));
	  break;
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
	case 'X':                            /* show extended options */
	  X_OPTS = TRUE;
	  return(optError(NULL));
	  break;
	case 'Y':                                    /* trace options */
	  TRACE[0] = TRUE;
	  cPtr++;
	  while(*cPtr && isgraph((int)*cPtr)) {
	    TRACE[(int)*cPtr] = TRUE;
	    cPtr++;
	  }
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
	default:                                         /* not found */
	  return(optError(cPtr));
        }
      } while(*cPtr);
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
  for(lPtr = lpType; lPtr->ident != NULL; lPtr++) {
    if(strnxcmp(anaOpts.type, lPtr->ident, 2) == 0)
      break;
  }
  if(lPtr->ident == NULL) {
    setAsspMsg(AEG_ERR_BUG, "main: unknown LP parameter type");
    return(prtAsspMsg(NULL));
  }
  entry = dtype2entry(lPtr->type, KDT_SSFF);
  if(entry != NULL && entry->keyword != NULL)
    lpTrack = entry->keyword;
  else {
    setAsspMsg(AEB_ERR_TRACK, "for LP coefficients");
    return(prtAsspMsg(NULL));
  }
  if(strlen(outExt) == 0)
    strcpy(outExt, lPtr->ext);

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
  printf("Release : %d.%d (%s)\n", RFC_MAJOR, RFC_MINOR, __DATE__);
  printf("Function: Linear Prediction analysis of <file>(s) using the\n");
  printf("          autocorrelation method and the Durbin recursion.\n");
  printf("          This program calculates the RMS amplitudes of the input\n");
  printf("          and residual signal in dB and, per default, reflection\n");
  printf("          coefficients (see '-t' option).\n");
  printf("          Analysis results will be written to a file with the\n");
  printf("          base name of the input file and the parameter type in\n");
  printf("          lower case as extension (e.g. '.rfc').\n");
  printf("          Default output is in SSFF binary format (tracks '%s',\n"\
	 "          '%s' and the LP type in lower case).\n",\
	 rmsTrack, gainTrack);
  printf("Options :\n");
  printf(" -h/--help  print this text\n");
  printf(" -b=<time>  set begin of analysis interval to <time> seconds\n");
  printf("            (default: begin of file)\n");
  printf(" -e=<time>  set end of analysis interval to <time> seconds\n");
  printf("            (default: end of file)\n");
  printf(" -s=<dur>   set analysis window shift to <dur> ms\n");
  printf("            (default: %.1f)\n", LP_DEF_SHIFT);
  printf(" -L=<dur>   set effective length of analysis window to <dur> ms\n");
  printf("            (default: %.1f)\n", LP_DEF_SIZE);
  printf(" -m=<num>   set prediction order to <num>\n");
  printf("            (default: %s)\n", DFLT_ORDER_STR);
  printf(" -p=<val>   set pre-emphasis factor to <val>\n");
  printf("            (default: %.2f)\n", LP_DEF_PREEMPH);
  printf(" -t=<type>  calculate <type> LP parameters; <type> may be:\n");
  printf("            \"ARF\": area function\n");
  printf("            \"LAR\": log area ratios\n");
  printf("            \"LPC\": linear prediction filter coefficients\n");
  printf("            \"RFC\": reflection coefficients (default)\n");
  printf(" -w=<type>  set analysis window function to <type>\n");
  printf("            (default: %s)\n", LP_DEF_WINDOW);
  printf(" -W         print summary of window functions\n");
  printf(" -od        store output file(s) in directory of input file(s)\n");
  printf("            (default: current working directory)\n");
  printf(" -od=<dir>  store output file(s) in directory <dir>\n");
  printf(" -ox=<ext>  set extension of output file names to <ext>\n");
  printf(" -z=<list>  batch mode: take input file names from file <list>;\n");
  printf("            overrules arguments <file>\n");
  printf(" -X         show extended options\n");
if(TRACE[0] || X_OPTS) {
  printf("         extended options\n");
  printf(" -C=<num>   for multi-channel input files: analyse channel <num>\n");
  printf("            with 1 <= <num> <= %d  (default: channel %d)\n",\
	 LP_I_CHANS, LP_DEF_CHANNEL);
  printf(" -c=<time>  set single-frame analysis with the analysis window\n");
  printf("            centred at <time> seconds; overrules -b, -e and -s\n");
  printf("            options\n");
  printf(" -S=<dur>   set exact analysis window size to <dur> ms;\n");
  printf("            overrules -L option\n");
  printf(" -oA        output in plain ASCII\n");
  printf(" -oa=<num>  set accuracy of coefficients in ASCII to <num> digits\n");
  printf("            (default: %d)\n", LP_DEF_DIGITSA);
  printf(" -op=<num>  set precision of RMS and gain in ASCII to <num> digits\n");
  printf("            (default: %d)\n", LP_DEF_DIGITSP);
  printf(" -of=<file> set name of output file to <file>\n");
  printf("            (only in single file mode; overrules -od option)\n");
}
if(TRACE[0]) {
  printf("         debugging options\n");
  printf(" -Y<abc>    output trace info to 'stderr' for options <abc>\n");
  printf("            the following options are available:\n");
  printf("   'F'   append trace output to file '%s'\n", traceFile);
  printf("   'f'   write trace output to file '%s'\n", traceFile);
  printf("   'A'   analysis parameters\n");
  printf("   'N'   use length-normalized ACF (instable)\n");
}
  printf("NOTE:\n");
  printf(" o   An existing output file will be overwritten without notice.\n");
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
    strcpy(lpFile, fixOutFile);
    fixOutFile = NULL;
  }
  else {
    parsepath(smpFile, &dPtr, &bPtr, NULL);
    if(anaOpts.options & AOPT_IN_DIR)
      strcpy(lpFile, dPtr);        /* store output in input directory */
    else if(outDir != NULL) {
      strcpy(lpFile, outDir);
      n = strlen(outDir);
      if(n > 0 && outDir[n-1] != DIR_SEP_CHR)
	strcat(lpFile, DIR_SEP_STR);
    }
    else strcpy(lpFile, "");            /* store in current directory */
    strcat(lpFile, bPtr);   /* append base name */
    strcat(lpFile, outExt); /* append extension */
  }
  if(strcmp(smpFile, lpFile) == 0) {
    setAsspMsg(AEC_IO_CLASH, lpFile);
    return(prtAsspMsg(NULL));
  }
  if(TRACE['F'] || TRACE['f'])
    fprintf(traceFP, "Input file: %s\nOutput file: %s\n",\
	    smpFile, lpFile);
  return(1); /* we got one */
}
