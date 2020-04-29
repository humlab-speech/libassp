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
* File:     forest.c                                                   *
* Contents: Command line interface to the formant estimator 'forest'.  *
* Author:   Michel T.M. Scheffers                                      *
*                                                                      *
***********************************************************************/
/* $Id: forest.c,v 1.10 2010/07/14 13:44:30 mtms Exp $ */

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
#include <fmt.h>      /* processing parameters & public functions */
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
char *lp1Track, *ffrTrack, *fbwTrack;
char  smpFile[PATH_MAX+1], fmtFile[PATH_MAX+1], outExt[SUFF_MAX+1];
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
  int   err, n;
  DOBJ *smpDOp, *fmtDOp;

  strcpy(progName, mybasename(argv[0]));
  strcpy(traceFile, progName);
  strcat(traceFile, TRACE_SUFFIX);
  initGlobals();
  if((err=evalArgs(argc, argv)) < 0)
    exit(1);
  if(err)
    exit(0);
  if(TRACE['F'] || TRACE['f']) {
    sprintf(info, "%s version %d.%d", progName, FMT_MAJOR, FMT_MINOR);
    openTrace(info);
  }
  else openTrace(NULL);
  if(TRACE['s']) {
    initFMTstats();
    if(TRACE['H'] || TRACE['h'])
      statInclHist(&statPQ, 0.0, 5.0, 120);
  }
/*
 * loop over input files
 */
  smpDOp = fmtDOp = NULL;
  while((err=nextFile()) > 0) {
    if(TRACE['s'])
      totFMTfiles++;
    smpDOp = asspFOpen(smpFile, AFO_READ, NULL);
    if(smpDOp == NULL) {
      err = prtAsspMsg(NULL);
      break;
    }
    fmtDOp = createFMT(smpDOp, &anaOpts);
    if(fmtDOp == NULL) {
      err = prtAsspMsg(NULL);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(asspFOpen(fmtFile, AFO_WRITE, fmtDOp) == NULL) {
      err = prtAsspMsg(NULL);
      fmtDOp = freeDObj(fmtDOp);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(computeFMT(smpDOp, &anaOpts, fmtDOp) == NULL) {
      err = prtAsspMsg(NULL);
    }
    asspFClose(smpDOp, AFC_FREE);
    smpDOp = NULL;
    asspFClose(fmtDOp, AFC_FREE);
    fmtDOp = NULL;
    if(err < 0) {
      remove(fmtFile); /* contents invalid */
      break;
    }
    if(asspWarning)
      prtAsspMsg(NULL);
  } /* END loop over input files */
  if(batFP != NULL) {
    fclose(batFP);
    batFP = NULL;
  }
  if(TRACE['s']) {
    fprintf(traceFP, "%lu files, %lu frames processed (%lu frames non-silence)\n",\
	    totFMTfiles, totFMTframes, totFMTframes-totFMTsilent);
    fprintf(traceFP, "Iteration statistics [ ");
    if(anaOpts.options & FMT_NOT_TRACK_PQ) {
      if(anaOpts.options & FMT_NOT_USE_PF)
	fprintf(traceFP, "NOM_PQ ");
      else
	fprintf(traceFP, "USE_PF ");
    }
    else {
      fprintf(traceFP, "TRACK_PQ ");
      if(!(anaOpts.options & FMT_NOT_SORT_PQ))
	fprintf(traceFP, "SORT_PQ ");
      if(!(anaOpts.options & FMT_NOT_RETRY_PQ))
	fprintf(traceFP, "RETRY_PQ ");
    }
    if(anaOpts.options & FMT_OPT_PE_ADAPT)
      fprintf(traceFP, "PE_ADAPT ");
    fprintf(traceFP, "]\n");
    fprintf(traceFP, "  root solving: n = %u (%lu convergence failures)\n",\
	    statPQ.numX, totFMTfail);
    fprintf(traceFP, "%15s min = %.0f  mean = %.1f  max = %.0f\n", " ",\
	    statPQ.minX, statGetMean(&statPQ), statPQ.maxX);
    if(TRACE['H'] || TRACE['h'])
      fprintf(traceFP, "%15s Q1 = %.1f  median = %.1f  Q3 = %.1f\n", " ",\
	      statEstQuantile(&statPQ, 25), statEstQuantile(&statPQ, 50),\
	      statEstQuantile(&statPQ, 75) );
    if(statPF.numX > 0)
      fprintf(traceFP, "  P-frequencies: n = %u  min = %.0f  mean = %.1f"\
	      "  max = %.0f\n", statPF.numX, statPF.minX,\
	      statGetMean(&statPF), statPF.maxX);
    fprintf(traceFP, "Frequency statistics\n");
    for(n = 0; n < FMT_MAX_BUF; n++) {
      if(n < FMT_NUM_PSTATS && statP[n].numX > 0)
	fprintf(traceFP, "  p%i: n = %u  min = %4.0f  mean = %6.1f  "\
		"max = %4.0f  sd = %.2f\n", n+1, statP[n].numX,\
		statP[n].minX, statGetMean(&statP[n]), statP[n].maxX,\
		statGetSD(&statP[n]));
      if(n < FMT_NUM_FSTATS && statF[n].numX > 0)
	fprintf(traceFP, "  F%i: n = %u  min = %4.0f  mean = %6.1f  "\
		"max = %4.0f  sd = %.2f\n", n+1, statF[n].numX,\
		statF[n].minX, statGetMean(&statF[n]), statF[n].maxX,\
		statGetSD(&statF[n]));
    }
    if(TRACE['H']) {
      fprintf(traceFP, "Iteration histogram\n");
      statPrintHist(&statPQ, traceFP);
    }
    freeFMTstats();
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
  lp1Track = ffrTrack = fbwTrack = NULL;
  strcpy(outExt, FMT_DEF_SUFFIX);
  X_OPTS = FALSE;
  batFP = NULL;
  setFMTdefaults(&anaOpts);
  entry = dtype2entry(DT_LP1, KDT_SSFF);
  if(entry != NULL && entry->keyword != NULL)
    lp1Track = entry->keyword;
  else {
    setAsspMsg(AEB_ERR_TRACK, "for data type LP1");
    exit(prtAsspMsg(NULL));
  }
  entry = dtype2entry(DT_FFR, KDT_SSFF);
  if(entry != NULL && entry->keyword != NULL)
    ffrTrack = entry->keyword;
  else {
    setAsspMsg(AEB_ERR_TRACK, "for data type FFR");
    exit(prtAsspMsg(NULL));
  }
  entry = dtype2entry(DT_FBW, KDT_SSFF);
  if(entry != NULL && entry->keyword != NULL)
    fbwTrack = entry->keyword;
  else {
    setAsspMsg(AEB_ERR_TRACK, "for data type FBW");
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
	  case 'R':                                     /* references */
	  case 'r':
	    printFMTrefs();
	    return(1);
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
	case 'C':                                   /* channel number */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.channel = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.channel < 1 || anaOpts.channel > FMT_I_CHANS)
	      return(optError("-C"));
	  }
	  else return(optError("-C"));
	  break;
	case 'F':                                       /* nominal F1 */
	  cPtr++;
	  if(*cPtr == '1') {
	    cPtr++;
	    if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	      cPtr++;
	      anaOpts.nomF1 = strtod(cPtr, &cPtr);
	    }
	    else return(optError("-F1"));
	  }
	  else return(optError(--cPtr));
	  break;
	case 'f':                        /* defaults for female voice */
	  setFMTgenderDefaults(&anaOpts, 'f');
	  cPtr++;
	  break;
	case 'i':                 /* insert rough frequency estimates */
	  anaOpts.options |= FMT_OPT_INS_ESTS;
	  cPtr++;
	  break;
	case 'm':                                         /* LP order */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.order = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.order < 1)
	      return(optError("-m"));
	    if(anaOpts.order > MAXLPORDER) {
	      fprintf(stderr, "\nERROR: -m option: prediction order "\
		      "too high (maximally %d)\n", MAXLPORDER);
	      return(-1);
	    }
	    if(anaOpts.order < FMT_MIN_ORDER) {
	      fprintf(stderr, "\nERROR: -m option: prediction order "\
		      "too low (minimally %d)\n", FMT_MIN_ORDER);
	      return(-1);
	    }
	    if(ODD(anaOpts.order)) {
	      fprintf(stderr, "\nERROR: -m option: prediction order "\
	     	      "must be even\n");
	      /* otherwise zillions of convergence failures */
	      return(-1);
	    }
	    anaOpts.options |= FMT_OPT_LPO_FIXED;
	  }
	  else if(*cPtr == '+' || *cPtr == '-') {     /* in-/decrease */
	    while(*cPtr == '+' || *cPtr == '-') {       /* order by 2 */
	      if(*cPtr == '+')
		(anaOpts.increment)++;
	      else
		(anaOpts.increment)--;
	      cPtr++;
	    }
	    anaOpts.options &= ~FMT_OPT_LPO_FIXED;     /* re. default */
	  }
	  else return(optError("-m"));
	  break;
	case 'n':                        /* number of output formants */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.numFormants = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.numFormants > FMT_MAX_OUT) {
	      fprintf(stderr, "\nERROR: -n option: number of formants "\
		      "too large (maximally %d)\n", FMT_MAX_OUT);
	      return(-1);
	    }
	  }
	  else return(optError("-n"));
	  break;
	case 'p':                                     /* pre-emphasis */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.' || cPtr[1] == '-')) {
	    cPtr++;
	    anaOpts.preEmph = strtod(cPtr, & cPtr);
	    if(anaOpts.preEmph < -1.0 || anaOpts.preEmph > 0.0)
	      return(optError("-p"));
	    anaOpts.options |= FMT_OPT_PE_FIXED;
	    anaOpts.options &= ~FMT_OPT_PE_ADAPT;      /* clear flag */
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
	case 't':                                /* silence threshold */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.' || cPtr[1] == '-')) {
	    cPtr++;
	    anaOpts.threshold = strtod(cPtr, &cPtr);
	  }
	  else return(optError("-S"));
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
	case 'o':           /* output modifiers */
	  cPtr++;
	  switch(*cPtr) {
	  case 'A':                                          /* ASCII */
	    strcpy(anaOpts.format, FMT_ASC_FORMAT);
	    cPtr++;
	    break;
	  case 'a':                                       /* accuracy */
	    X_OPTS = TRUE;                         /* extended option */
	    cPtr++;
	    if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	      cPtr++;
	      anaOpts.accuracy = (int)strtol(cPtr, &cPtr, 10);
	    }
	    else return(optError("-op"));
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
  if(TRACE['a'])
    anaOpts.options |= FMT_OPT_PE_ADAPT;  /* signal-adaptive emphasis */
  if(TRACE['P'])
    anaOpts.options |= FMT_NOT_USE_PF;                  /* switch off */
  if(TRACE['R'])
    anaOpts.options |= FMT_NOT_RETRY_PQ;
  if(TRACE['S'])
    anaOpts.options |= FMT_NOT_SORT_PQ;
  if(TRACE['T']) {
    anaOpts.options |= FMT_NOT_TRACK_PQ;
    anaOpts.options |= FMT_NOT_RETRY_PQ;            /* waste of time */
  }
  if(TRACE['H'] || TRACE['h'])   /* histogram of bairstow iterations */
    TRACE['s'] = TRUE;                         /* implies statistics */

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
    fprintf(stderr, "ERROR: %s option: specification missing "\
	    "or unacceptable\n", opt);
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
  printf("Release : %d.%d (%s)\n", FMT_MAJOR, FMT_MINOR, __DATE__);
  printf("Function: Formant estimation of the speech signal in <file>.\n");
  printf("          Raw resonance frequency and bandwidth values are\n");
  printf("          obtained by root-solving of the Linear Prediction\n");
  printf("          polynomial from the autocorrelation method with the\n");
  printf("          Split-Levinson-Algorithm (SLA). Resonances are then\n");
  printf("          classified as formants using the so-called Pisarenko\n");
  printf("          frequencies (by-product of the SLA) and a formant\n");
  printf("          frequeny range table derived from the nominal F1\n");
  printf("          frequency. The latter should be increased by about\n");
  printf("          12%% for female voices (see -f and -F1 options).\n");
  printf("          Formant estimates will be written to a file with the\n");
  printf("          base name of the input file and extension '%s'.\n",\
	 FMT_DEF_SUFFIX);
  printf("          Default output is in SSFF binary format (tracks '%s'\n"\
	 "          and '%s')\n", ffrTrack, fbwTrack);
  printf("Options:\n");
  printf(" -h/--help  print this text\n");
  printf(" --refs     print references for the algorithm\n");
  printf(" -b=<time>  set begin of analysis interval to <time> seconds\n");
  printf("            (default: begin of file)\n");
  printf(" -e=<time>  set end of analysis interval to <time> seconds\n");
  printf("            (default: end of file)\n");
  printf(" -s=<dur>   set analysis window shift to <dur> ms\n");
  printf("            (default: %.1f)\n", FMT_DEF_SHIFT);
  printf(" -L=<dur>   set effective length of analysis window to <dur> ms\n");
  printf("            (default: %.1f)\n", FMT_DEF_EFFLENm);
  printf(" -F1=<freq> set nominal F1 frequency to <freq> Hz\n");
  printf("            (default: %.1f Hz)\n", FMT_DEF_NOMF1m);
  printf(" -f         use default settings for female voice\n");
  printf("            (eff. window length = %.1f ms  nominal F1 = %.1f Hz)\n",\
	 FMT_DEF_EFFLENf, FMT_DEF_NOMF1f);
  printf(" -i         insert rough frequency estimates of missing formants\n");
  printf("            (default: frequency and bandwidth set to zero)\n");
  printf(" -m=<num>   set Linear Prediction order to <num>\n");
  printf("            (default: smallest even value larger than the sample\n");
  printf("             rate divided by twice the nominal F1 frequency)\n");
  printf(" -m-        decrease default order by 2 (one resonance less)\n");
  printf(" -m+        increase default order by 2 (one resonance more)\n");
  printf(" -n=<num>   set number of output formants to <num>\n");
  printf("            (default: %d;  maximum: %d or half the LP order)\n",\
	  FMT_DEF_OUT, FMT_MAX_OUT);
  printf(" -w=<type>  set analysis window function to <type>\n");
  printf("            (default: %s)\n", FMT_DEF_WINDOW);
  printf(" -W         print summary of available window functions\n");
  printf(" -oA        output in plain ASCII format\n");
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
	 FMT_I_CHANS, FMT_DEF_CHANNEL);
  printf(" -c=<time>  set single-frame analysis with the analysis window\n");
  printf("            centred at <time> seconds; overrules -b, -e and -s\n");
  printf("            options\n");
  printf(" -p=<val>   set pre-emphasis factor to <val> (-1 <= val <= 0)\n");
  printf("            (default: dependent on sample rate and nominal F1)\n");
  printf(" -S=<dur>   set exact analysis window size to <dur> ms;\n");
  printf("            overrules -L option\n");
  printf(" -t=<num>   set silence threshold (no analysis) to <num> dB\n");
  printf("            (default: %.1f dB)\n", FMT_DEF_SILENCE);
  printf(" -of=<file> set name of output file to <file>\n");
  printf("            (only in single-file mode; overrules -od option)\n");
  printf(" -oa=<num>  set accuracy of LP1 in ASCII output to <num> digits\n");
  printf("            (default: %d)\n", FMT_DEF_DIGITSA);
}
if(TRACE[0]) {
  printf("         debugging options\n");
  printf(" -Y<abc>    output trace info to 'stderr' for options <abc>\n");
  printf("            the following options are available:\n");
  printf("   'F'   append trace output to file '%s'\n", traceFile);
  printf("   'f'   write trace output to file '%s'\n", traceFile);
  printf("   'l'   formant limits table\n");
  printf("   'A'   analysis parameters\n");
  printf("   'c'   convergence failures of bairstow()\n");
  printf("   '0'   raw analysis data\n");
  printf("   '1'   data after removal of non-resonances\n");
  printf("   '2'   data after resonance (re-)classification\n");
  printf("   '3'   data after enforced removal of double assignments\n");
  printf("   '4'   data after final classification\n");
  printf("   '5'   data with estimates for missing formants\n");
  printf("   'p'   probabilities of formant mapping (with 2 and 3)\n");
  printf("   's'   statistics of iterations and lower frequencies\n");
  printf("   'h'   quartiles of bairstow iterations from histogram\n");
  printf("   'H'   histogram of bairstow iterations\n");
  printf("   'a'   use signal-adaptive pre-emphasis\n");
  printf("   'P'   don't use Pisarenko frequencies\n");
  printf("   'R'   don't reset and try again on convergence failure\n");
  printf("   'S'   don't sort PQ parameters\n");
  printf("   'T'   don't track PQ parameters\n");
}
  printf("NOTE:\n");
  printf(" o   Existing output files will be overwritten without notice.\n");
  printf(" o   The analysis interval will be rounded to the nearest frame\n");
  printf("     boundary (integral multiple of the window shift).\n");
  printf(" o   The -f option overrules earlier specifications using the\n");
  printf("     -F1, -L and/or -S options.\n");
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
  static int fileNum=0;
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
    strcpy(fmtFile, fixOutFile);
    fixOutFile = NULL;
  }
  else {
    parsepath(smpFile, &dPtr, &bPtr, NULL);
    if(anaOpts.options & AOPT_IN_DIR)
      strcpy(fmtFile, dPtr);       /* store output in input directory */
    else if(outDir != NULL) {
      strcpy(fmtFile, outDir);
      n = strlen(outDir);
      if(n > 0 && outDir[n-1] != DIR_SEP_CHR)
	strcat(fmtFile, DIR_SEP_STR);
    }
    else strcpy(fmtFile, "");           /* store in current directory */
    strcat(fmtFile, bPtr);   /* append base name */
    strcat(fmtFile, outExt); /* append extension */
  }
  if(strcmp(smpFile, fmtFile) == 0) {
    setAsspMsg(AEC_IO_CLASH, fmtFile);
    return(prtAsspMsg(NULL));
  }
  if(TRACE['F'] || TRACE['f'])
    fprintf(traceFP, "Input file: %s\nOutput file: %s\n",\
	    smpFile, fmtFile);
  return(1); /* we got one */
}
