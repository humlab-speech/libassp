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
* File:     spectrum.c                                                 *
* Contents: Command line interface to spectral/cepstral analyses.      *
* Author:   Michel T.M. Scheffers                                      *
*                                                                      *
***********************************************************************/
/* $Id: spectrum.c,v 1.2 2010/07/22 12:12:36 mtms Exp $ */

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
#include <spectra.h>  /* processing parameters & SPECT functions */
#include <assp.h>     /* message & trace handler */
#include <asspfio.h>  /* file & data handler */
#include <asspdsp.h>  /* wfunc_e wfType() listWFs() wfShortList */
#include <asspana.h>  /* AOPTS */
#include <dataobj.h>  /* DOBJ DT_... */

/*
 * constants
 */
#define MAXFILES 2000 /* multiple input files */

/*
 * global arrays and variables
 */
char  progName[NAME_MAX+1];
char *argList[MAXFILES], *batList, *fixOutFile, *outDir;
char  smpFile[PATH_MAX+1], spectFile[PATH_MAX+1], outExt[SUFF_MAX+1];
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
  DOBJ *smpDOp, *spectDOp;
  
  strcpy(progName, mybasename(argv[0]));
  strcpy(traceFile, progName);
  strcat(traceFile, TRACE_SUFFIX);
  initGlobals();
  if((err=evalArgs(argc, argv)) < 0)
    exit(1);
  if(err != 0)
    exit(0);
  if(TRACE['F'] || TRACE['f']) {
    sprintf(info, "%s version %d.%d", progName, SPECT_MAJOR, SPECT_MINOR);
    openTrace(info);
  }
  else openTrace(NULL);
  /* loop over input files */
  smpDOp = spectDOp = NULL;
  while((err=nextFile()) > 0) {
    smpDOp = asspFOpen(smpFile, AFO_READ, NULL);
    if(smpDOp == NULL) {
      err = prtAsspMsg(NULL);
      break;
    }
    spectDOp = createSPECT(smpDOp, &anaOpts);
    if(spectDOp == NULL) {
      err = prtAsspMsg(NULL);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(asspFOpen(spectFile, AFO_WRITE, spectDOp) == NULL) {
      err = prtAsspMsg(NULL);
      spectDOp = freeDObj(spectDOp);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(computeSPECT(smpDOp, &anaOpts, spectDOp) == NULL) {
      err = prtAsspMsg(NULL);
    }
    asspFClose(smpDOp, AFC_FREE);
    smpDOp = NULL;
    asspFClose(spectDOp, AFC_FREE);
    spectDOp = NULL;
    if(err < 0) {
      remove(spectFile); /* contents invalid */
      break;
    }
    if(asspWarning)
      prtAsspMsg(NULL);
  } /* END loop over input files */
  closeTrace();
  if(batFP != NULL) {
    fclose(batFP);
    batFP = NULL;
  }
  if(err < 0)
    exit(1);
  exit(0);
}
/***********************************************************************
* initialize global parameters, pointers and file descriptors          *
***********************************************************************/
LOCAL void initGlobals(void)
{
  clrAsspMsg();
  batList = fixOutFile = outDir = NULL;
  batFP = NULL;
  strcpy(outExt, "");
  X_OPTS = FALSE;
  setSPECTdefaults(&anaOpts);
  return;
}
/***********************************************************************
* check program options and arguments                                  *
***********************************************************************/
LOCAL int evalArgs(int argc, char *argv[])
{
  char   *cPtr;
  int     i;
  dtype_e spectType;
  wfunc_e winFunc;

  numFiles = 0;
  for(i = 1; i < argc; i++) {
    cPtr = argv[i];
    if(*cPtr == '-') {          /* option */
      cPtr++;
      do {
	switch(*cPtr) {
	case '-':               /* long option */
	  cPtr++;
	  switch(*cPtr) {
	  case 'H':                                           /* HELP */
	  case 'h':
	    return(optError(NULL));
	    break;
	  default:
	    return(optError(cPtr-2));
	  }
	  break;
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
	case 'B':                              /* effective bandwidth */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.' )) {
	    cPtr++;
	    anaOpts.bandwidth = strtod(cPtr, &cPtr);
	    anaOpts.options |= AOPT_USE_ENBW;
	  }
	  else return(optError("-B"));
	  break;
	case 'C':                                   /* channel number */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.channel = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.channel < 1 || anaOpts.channel > SPECT_I_CHANS)
	      return(optError("-C"));
	  }
	  else return(optError("-C"));
	  break;
	case 'm':                                /* LP order/ CS lags */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.order = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.order < 1)
	      return(optError("-m"));
	  }
	  else return(optError("-m"));
	  break;
	case 'n':                             /* number of FFT points */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.FFTLen = strtol(cPtr, &cPtr, 10);
	    if(anaOpts.FFTLen < MIN_NFFT)
	      return(optError("-n"));
	    anaOpts.resolution = 0.0;
	  }
	  else return(optError("-n"));
	  break;
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
	case 'D':                                 /* omit de-emphasis */
	  anaOpts.options &= ~LPS_OPT_DEEMPH;
	  cPtr++;
	  break;
	case 'r':                              /* spectral resolution */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.' )) {
	    cPtr++;
	    anaOpts.resolution = strtod(cPtr, &cPtr);
	    if(anaOpts.resolution < SPECT_MIN_RES)
	      return(optError("-r"));
	    anaOpts.FFTLen = 0;                          /* overruled */
	  }
	  else return(optError("-r"));
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
	    anaOpts.options &= ~AOPT_USE_ENBW;
	    anaOpts.options |= AOPT_EFFECTIVE;
	  }
	  else return(optError("-L"));
	  break;
	case 'S':                                /* exact window size */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	    cPtr++;
	    anaOpts.msSize = strtod(cPtr, &cPtr);
	    anaOpts.options &= ~(AOPT_EFFECTIVE | AOPT_USE_ENBW);
	  }
	  else return(optError("-S"));
	  break;
	case 't':                                    /* spectrum type */
	  cPtr++;
	  if(*cPtr == '=') {
	    cPtr++;
	    spectType = getSPECTtype(cPtr, NULL);
	    if(spectType != DT_ERROR) {
	      strcpy(anaOpts.type, cPtr);
	      cPtr += strlen(cPtr);
	      switch(spectType) {            /* adjust default values */
	      case DT_FTPOW:
		setDFTdefaults(&anaOpts);
		break;
	      case DT_FTLPS:
		setLPSdefaults(&anaOpts);
		break;
	      case DT_FTCSS:
		setCSSdefaults(&anaOpts);
		break;
	      case DT_FTCEP:
		setCEPdefaults(&anaOpts);
		break;
	      default:
		return(optError("-t"));
	      }
	    }
	    else return(optError("-t"));
	  }
	  else return(optError("-t"));
	  break;
	case 'o':               /* output modifiers */
	  cPtr++;
	  switch(*cPtr) {
	  case 'A':                                          /* ASCII */
	    X_OPTS = TRUE;                         /* extended option */
	    strcpy(anaOpts.format, SPECT_ASC_FORMAT);
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
	    if(*cPtr == '=' && cPtr[1]) {
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
	  while(*cPtr && isgraph((int)*cPtr)) {
	    TRACE[(int)*cPtr] = TRUE;
	    cPtr++;
	  }
	  break;
	default:                                         /* not found */
	  return(optError(cPtr));
        }
      } while(*cPtr);
    }
    else {                      /* ARGUMENT */
      if(numFiles < MAXFILES) {
	argList[numFiles] = cPtr;
	numFiles++;
      }
      else {
	fprintf(stderr, "\nERROR: more than %d input files\n", MAXFILES);
	return(-1);
      }
    }
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
  if(strlen(outExt) == 0)
    spectType = getSPECTtype(anaOpts.type, outExt);
  else
    spectType = getSPECTtype(anaOpts.type, NULL);
  if(spectType == DT_ERROR)
    return(prtAsspMsg(NULL));

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
  printf("Release : %d.%d (%s)\n", SPECT_MAJOR, SPECT_MINOR, __DATE__);
  printf("Function: Short-term spectral analysis of the signal in <file>\n");
  printf("          using the Fast Fourier Transform. The default is to\n");
  printf("          calculate unsmoothed narrow-band spectra with the size\n");
  printf("          of the analysis window equal to the length of the FFT.\n");
  printf("          The output from the FFT will be converted to a power\n");
  printf("          spectrum in dB from 0 Hz up to and including the\n");
  printf("          Nyquist rate.\n");
  printf("          Alternatively, the program can calculate smoothed\n");
  printf("          spectra or cepstra. In the latter case the number of\n"); 
  printf("          coefficients per output record will also equal the\n");
  printf("          FFT length / 2 + 1 (i.e. be non-mirrored).\n"); 
  printf("          Analysis results will be written to a file with the\n");
  printf("          base name of the input file and the spectrum type in\n");
  printf("          lower case as extension (e.g. '.dft').\n");
  printf("          Default output is in SSFF binary format with the\n");
  printf("          spectrum type in lower case as track name.\n");
  printf("Options :\n");
  printf(" -h/---help print this text\n");
  printf(" -t=<type>  set spectrum type; <type> may be:\n");
  printf("            \"DFT\": unsmoothed spectrum (default)\n");
  printf("            \"LPS\": Linear Prediction smoothed spectrum\n");
  printf("            \"CSS\": cepstral smoothed spectrum\n");
  printf("            \"CEP\": cepstrum\n");
  printf(" -b=<time>  set begin of analysis interval to <time> seconds\n");
  printf("            (default: begin of file)\n");
  printf(" -e=<time>  set end of analysis interval to <time> seconds\n");
  printf("            (default: end of file)\n");
  printf(" -r=<freq>  set FFT length to the smallest value which results\n");
  printf("            in a frequency resolution of <freq> Hz or better\n");
  printf("            (default: %.1f)\n", SPECT_DEF_RES);
  printf(" -n=<num>   set FFT length to <num> points\n");
  printf("            (overrules default settings and '-r' option)\n");
  printf(" -s=<dur>   set analysis window shift to <dur> ms\n");
  printf("            (default: %.1f)\n", SPECT_DEF_SHIFT);
  printf(" -w=<type>  set analysis window function to <type>\n");
  printf("            (default: %s)\n", SPECT_DEF_WINDOW);
  printf(" -W         print summary of window functions\n");
  printf(" -od        store output file(s) in directory of input file(s)\n");
  printf("            (default: current working directory)\n");
  printf(" -od=<dir>  store output file(s) in directory <dir>\n");
  printf(" -ox=<ext>  set extension of output file name(s) to <ext>\n");
  printf(" -X         show extended options\n");
  printf("         options for DFT spectrum:\n");
  printf(" -B=<freq>  set the effective analysis bandwidth to <freq> Hz\n");
  printf("            (default: 0, yielding the smallest possible value\n");
  printf("             given the length of the FFT)\n");
  printf("         options for LP smoothed spectrum:\n");
  printf(" -L=<dur>   set effective length of analysis window to <dur> ms\n");
  printf("            (default: %.1f)\n", LPS_DEF_SIZE);
  printf(" -m=<num>   set prediction order to <num>\n");
  printf("            (default: %s)\n", DFLT_ORDER_STR);
  printf(" -p=<val>   set pre-emphasis factor to <val>\n");
  printf("            (default: %.2f)\n", LPS_DEF_PREEMPH);
  printf(" -D         omit de-emphasis (default: undo spectral tilt due\n");
  printf("            to pre-emphasis used in LP analysis)\n");
  printf("         options for cepstral smoothed spectrum:\n");
  printf(" -m=<num>   set number of cepstral coefficients used to <num>\n");
  printf("            (default: %s)\n", CSS_DFLT_LAGS_STR);
if(TRACE[0] || X_OPTS) {
  printf("         extended options\n");
  printf(" -C=<num>   for multi-channel input files: analyse channel <num>\n");
  printf("            with 1 <= <num> <= %d  (default: channel %d)\n",\
	 SPECT_I_CHANS, SPECT_DEF_CHANNEL);
  printf(" -c=<time>  set single-frame analysis with the analysis window\n");
  printf("            centred at <time> seconds; overrules -b, -e and -s\n");
  printf("            options\n");
  printf(" -S=<dur>   set exact analysis window size to <dur> ms;\n");
  printf("            overrules -B and -L options, respectively\n");
  printf(" -oA        output in plain ASCII\n");
  printf(" -of=<file> set name of output file to <file>\n");
  printf("            (only in single file mode; overrules -od option)\n");
  printf(" -op=<num>  set precision of ASCII output in dB to <num> digits\n");
  printf("            (default: %d)\n", SPECT_DEF_DIGITSP);
}
if(TRACE[0]) {
  printf(" -Y<abc>    output trace info to 'stderr' for options <abc>\n");
  printf("            the following options are available:\n");
  printf("   'F'   append trace output to file '%s'\n", traceFile);
  printf("   'f'   write trace output to file '%s'\n", traceFile);
  printf("   'A'   analysis parameters\n");
  printf("   'c'   cepstral coefficients\n");
  printf("   'C'   cepstral coefficients after zeroing\n");
  printf("   'd'   rounding errors in Durbin recursion\n");
}
  printf("NOTE:\n");
  printf(" o   The '-t' option should be the first option on the command\n");
  printf("     line because it adjusts some default values.\n");
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
    strcpy(spectFile, fixOutFile);
    fixOutFile = NULL;
  }
  else {
    parsepath(smpFile, &dPtr, &bPtr, NULL);
    if(anaOpts.options & AOPT_IN_DIR)
      strcpy(spectFile, dPtr);       /* store output in input directory */
    else if(outDir != NULL) {
      strcpy(spectFile, outDir);
      n = strlen(outDir);
      if(n > 0 && outDir[n-1] != DIR_SEP_CHR)
	strcat(spectFile, DIR_SEP_STR);
    }
    else strcpy(spectFile, "");           /* store in current directory */
    strcat(spectFile, bPtr);   /* append base name */
    strcat(spectFile, outExt); /* append extension */
  }
  if(strcmp(smpFile, spectFile) == 0) {
    setAsspMsg(AEC_IO_CLASH, spectFile);
    return(prtAsspMsg(NULL));
  }
  if(TRACE['F'] || TRACE['f'])
    fprintf(traceFP, "Input file: %s\nOutput file: %s\n",\
	    smpFile, spectFile);
  return(1); /* we got one */
}
