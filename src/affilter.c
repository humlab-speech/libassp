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
* File:     affilter.c                                                 *
* Contents: Command line interface to audio filtering functions.       *
* Author:   Michel T.M. Scheffers                                      *
*                                                                      *
***********************************************************************/
/* $Id: affilter.c,v 1.4 2010/03/30 14:12:31 mtms Exp $ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>    /* [f/s]printf() remove() stderr FILE NULL */
#include <stdlib.h>   /* strtod() strtol() */
#include <stddef.h>   /* size_t */
#include <string.h>   /* str..() */
#include <ctype.h>    /* isdigit() isgraph() */
#include <math.h>     /* fabs() */

#include <miscdefs.h> /* TRUE FALSE LOCAL */
#include <misc.h>     /* mybasename() parsepath() fgetl() */
#include <mylimits.h> /* PATH_MAX NAME_MAX SUFF_MAX */
#include <filter.h>   /* processing parameters & FILT functions */
#include <assp.h>     /* message & trace handler */
#include <asspfio.h>  /* file handler */
#include <asspdsp.h>  /* AF...GAIN FILTER FIR IIR2 */
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
char  inpFile[PATH_MAX+1], outFile[PATH_MAX+1], outExt[SUFF_MAX+1];
int   numFiles;
int   X_OPTS, filterType;
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
  DOBJ *inpDOp, *filtDOp, outDObj, *outDOp;

  strcpy(progName, mybasename(argv[0]));
  strcpy(traceFile, progName);
  strcat(traceFile, TRACE_SUFFIX);
  initGlobals();
  if((err=evalArgs(argc, argv)) < 0)
    exit(1);
  if(err)
    exit(0);
  if(TRACE['F'] || TRACE['f']) {
    sprintf(info, "%s version %d.%d", progName, FILT_MAJOR, FILT_MINOR);
    openTrace(info);
  }
  else openTrace(NULL);
  if(TRACE['p']) {
    fprintf(traceFP, "Filter parameters\n");
    switch(filterType) {
    case FILTER_LP:
      fprintf(traceFP, "  low-pass to %.1f Hz\n", anaOpts.lpCutOff);
      break;
    case FILTER_HP:
      fprintf(traceFP, "  high-pass from %.1f Hz\n", anaOpts.hpCutOff);
      break;
    case FILTER_BP:
      fprintf(traceFP, "  band-pass from %.1f to %.1f Hz\n",	\
	      anaOpts.hpCutOff, anaOpts.lpCutOff);
      break;
    case FILTER_BS:
      fprintf(traceFP, "  band-stop between %.1f and %.1f Hz\n",	\
	      anaOpts.lpCutOff, anaOpts.hpCutOff);
      break;
    }
    if(anaOpts.options & FILT_OPT_USE_IIR) {
      fprintf(traceFP, "  number of IIR sections: %i\n",	\
	      anaOpts.order);
    }
    else {
      fprintf(traceFP, "  width of transition band: %.1f Hz\n",	\
	      anaOpts.tbWidth);
      fprintf(traceFP, "  stop-band attenuation: %.1f dB\n",	\
	      anaOpts.stopDB);
    }
  }
/*
 * loop over input files
 */
  inpDOp = filtDOp = NULL;
  outDOp = &outDObj;
  while((err=nextFile()) > 0) {
    inpDOp = asspFOpen(inpFile, AFO_READ, NULL);
    if(inpDOp == NULL) {
      err = prtAsspMsg(NULL);
      break;
    }
    filtDOp = createFilter(inpDOp, &anaOpts);
    if(filtDOp == NULL) {
      err = prtAsspMsg(NULL);
      asspFClose(inpDOp, AFC_FREE);
      inpDOp = NULL;
      break;
    }
    err = copyDObj(outDOp, inpDOp);
    if(err != 0) {
      prtAsspMsg(NULL);
      if(err < 0) {
	asspFClose(inpDOp, AFC_FREE);
	inpDOp = NULL;
	filtDOp = destroyFilter(filtDOp);
	break;
      }
      clrAsspMsg();
      err = 0;
    }
    if(outDOp->ddl.numFields > FILT_O_CHANS) {
      outDOp->ddl.numFields = FILT_O_CHANS;
      setRecordSize(outDOp);              /* needs to be recalculated */
    }
    if(asspFOpen(outFile, AFO_WRITE, outDOp) == NULL) {
      err = prtAsspMsg(NULL);
      asspFClose(inpDOp, AFC_FREE);
      inpDOp = NULL;
      filtDOp = destroyFilter(filtDOp);
      clearDObj(outDOp);
      break;
    }
    if(filterSignal(inpDOp, filtDOp, outDOp) == NULL) {
      err = prtAsspMsg(NULL);
    }
    asspFClose(inpDOp, AFC_FREE);
    inpDOp = NULL;
    filtDOp = destroyFilter(filtDOp);
    asspFClose(outDOp, AFC_CLEAR); /* BEWARE: outDObj fixed! */
    if(err < 0) {
      remove(outFile); /* contents invalid */
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
* initialize global parameters, pointers and file descriptors          *
***********************************************************************/
LOCAL void initGlobals(void)
{
  clrAsspMsg();
  batList = fixOutFile = outDir = NULL;
  batFP = NULL;
  strcpy(outExt, "");
  X_OPTS = FALSE;
  setFILTdefaults(&anaOpts);
  filterType = FILTER_NONE;
  return;
}
/***********************************************************************
* check program options and arguments                                  *
***********************************************************************/
LOCAL int evalArgs(int argc, char *argv[])
{
  char  *cPtr;
  int    i;
  double temp;

  numFiles = 0;
  for(i = 1; i < argc; i++) {
    cPtr = argv[i];
    if(*cPtr == '-') {          /* OPTION */
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
	  break;               /* end long option */
	  /* the range options are omitted because we don't know */
	  /* what exactly the user wants with them: filter only the */
	  /* range and leave the rest unfiltered or filter and */
	  /* excise (both will produce clicks!) */
/* 	case 'b':                                        begin time */
/* 	  cPtr++; */
/* 	  if(*cPtr == '=' && (isdigit((int)cPtr[1]) || cPtr[1] == '.' )) { */
/* 	    cPtr++; */
/* 	    anaOpts.beginTime = strtod(cPtr, &cPtr); */
/* 	  } */
/* 	  else return(optError("-b")); */
/* 	  break; */
/* 	case 'e':                                          end time */
/* 	  cPtr++; */
/* 	  if(*cPtr == '=' && (isdigit((int)cPtr[1]) || cPtr[1] == '.' )) { */
/* 	    cPtr++; */
/* 	    anaOpts.endTime = strtod(cPtr, &cPtr); */
/* 	  } */
/* 	  else return(optError("-e")); */
/* 	  break; */
	case 'h':                                  /* help/hp cut-off */
	  cPtr++;
	  if(*cPtr == EOS || *cPtr != 'p')                    /* help */
	    return(optError(NULL));
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	    cPtr++;
	    anaOpts.hpCutOff = strtod(cPtr, &cPtr);
	  }
	  else return(optError("-hp"));
	  break;
	case 'l':                                       /* lp cut-off */
	  cPtr++;
	  if(*cPtr == 'p') {
	    cPtr++;
	    if(*cPtr == '=' &&
	       (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	      cPtr++;
	      anaOpts.lpCutOff = strtod(cPtr, &cPtr);
	    }
	    else return(optError("-lp"));
	  }
	  else return(optError(--cPtr));
	  break;
	case 'a':                            /* stop-band attenuation */
	  anaOpts.options &= ~FILT_OPT_USE_IIR;
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	    cPtr++;
	    anaOpts.stopDB = strtod(cPtr, &cPtr);
	    if(anaOpts.stopDB <= 0.0)
	      return(optError("-a"));
	  }
	  else return(optError("-a"));
	  break;
	case 't':                         /* width of transition band */
	  anaOpts.options &= ~FILT_OPT_USE_IIR;
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	    cPtr++;
	    anaOpts.tbWidth = strtod(cPtr, &cPtr);
	    if(anaOpts.tbWidth <= 2.0)
	      return(optError("-t"));
	  }
	  else return(optError("-t"));
	  break;
	case 'I':                                   /* use IIR filter */
	  anaOpts.options |= FILT_OPT_USE_IIR;
	  cPtr++;
	  if(*cPtr == '=') {          /* number of 2nd-order sections */
	    cPtr++;
	    if(isdigit((int)(*cPtr))) {
	      anaOpts.order = (int)strtol(cPtr, &cPtr, 10);
	      if(anaOpts.order < 1)
		return(optError("-I"));
	    }
	    else return(optError("-I"));
	  }
	  else {
	    anaOpts.order = FILT_DEF_SECTS;
	  }
	  break;
	case 'C':                                   /* channel number */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.channel = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.channel < 1 || anaOpts.channel > FILT_I_CHANS)
	      return(optError("-C"));
	  }
	  else return(optError("-C"));
	  break;
	case 'd':                             /* switch dithering off */
	  X_OPTS = TRUE;
	  anaOpts.options |= FILT_NOPT_DITHER;
	  cPtr++;
	  break;
	case 'g':                                        /* auto-gain */
	  X_OPTS = TRUE;
	  anaOpts.options |= FILT_OPT_AUTOGAIN;
	  cPtr++;
	  if(*cPtr == '=') {
	    cPtr++;
	    if(*cPtr != EOS) {
	      temp = strtod(cPtr, &cPtr);
	      if(temp < AF_MIN_GAIN || temp > AF_MAX_GAIN)
		return(optError("-g"));
	      anaOpts.gain = temp;
	    }
	    else return(optError("-g"));
	  }
	  else {
	    anaOpts.gain = AF_DEF_GAIN;
	  }
	  break;
	case 'o':               /* output modifiers */
	  cPtr++;
	  switch(*cPtr) {
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
	    if(*cPtr == '=' && (cPtr[1] != EOS)) {
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
	  break;                /* END output modifiers */
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
  if(numFiles > 1)
    fixOutFile = NULL;
  if(fixOutFile != NULL) {
    outDir = NULL;                                       /* overruled */
    anaOpts.options &= ~AOPT_IN_DIR;                     /* overruled */
  }
  if(strlen(outExt) == 0)
    filterType = getFILTtype(&anaOpts, outExt);
  else
    filterType = getFILTtype(&anaOpts, NULL);
  if(filterType <= 0)
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
  printf("Syntax  : %s [<opts>] <file> [<opts>] {<file> [<opts>]}\n",	\
	 progName);
  printf("Release : %d.%d (%s)\n", FILT_MAJOR, FILT_MINOR, __DATE__);
  printf("Function: Filters the audio signal in <file>.\n");
  printf("          By specifying the high-pass and/or low-pass cut-off [1]\n");
  printf("          frequency one of four filter characteristics may be\n");
  printf("          selected as shown in the table below.\n\n");
  printf("            hp     lp     filter characteristic         extension\n");
  printf("          -------------------------------------------------------\n");
  printf("           > 0     [0]    high-pass from hp             '.hpf'\n");
  printf("            [0]   > 0     low-pass up to lp             '.lpf'\n");
  printf("           > 0    >= hp   band-pass from hp to lp       '.bpf'\n");
  printf("           > lp   > 0     band-stop between lp and hp   '.bsf'\n\n");
  printf("          Per default the window design method (Kaiser window) will\n");
  printf("          be used to compute the coefficients of a linear-phase FIR\n");
  printf("          filter with unity gain pass-band and adjustable stop-band\n");
  printf("          attenuation and width of the transition band.\n");
  printf("          Since steep FIR filters are computationally expensive [2]\n");
  printf("          you may optionally use a recursive (IIR) filter with a\n");
  printf("          Butterworth characteristic. Such a filter works very much\n");
  printf("          faster but its non-linear phase may distort the waveform\n");
  printf("          of e.g. an EGG signal beyond recognition.\n");
  printf("          The filtered signal will be written to a file with the\n");
  printf("          base name of the input file and an extension corresponding\n");
  printf("          to the filter characteristic (see table). The format of\n");
  printf("          the output file will be the same as that of the input file.\n");
  printf("Options :\n");
  printf(" -h/--help  print this text\n");
/*   printf(" -b=<time>  set begin of processing interval to <time> seconds\n"); */
/*   printf("            (default: begin of file)\n"); */
/*   printf(" -e=<time>  set end of processing interval to <time> seconds\n"); */
/*   printf("            (default: end of file)\n"); */
  printf(" -hp=<num>  set the high-pass cut-off frequency to <num> Hz\n");
  printf("            (default: 0, no high-pass filtering)\n");
  printf(" -lp=<num>  set the low-pass cut-off frequency to <num> Hz\n");
  printf("            (default: 0, no low-pass filtering)\n");
  printf(" -a=<num>   FIR: set the stop-band attenuation to <num> dB\n");
  printf("            (default: %.1f dB, minimum: %.1f dB)\n",	\
	 FILT_DEF_ATTEN, FILT_MIN_ATTEN);
  printf(" -t=<num>   FIR: set the width of the transition band to <num> Hz\n");
  printf("            (default: %.1f Hz)\n", FILT_DEF_WIDTH);
  printf(" -I[=<num>] use an IIR filter with <num> 2nd order sections [3];\n");
  printf("            if <num> is not specified %d sections will be used\n", \
	 FILT_DEF_SECTS);
  printf(" -od        store output file(s) in directory of input file(s)\n");
  printf("            (default: current working directory)\n");
  printf(" -od=<dir>  store output file(s) in directory <dir>\n");
  printf(" -ox=<ext>  set extension of output file names to <ext>\n");
  printf(" -z=<list>  batch mode: take input file names from file <list>;\n");
  printf("            overrules arguments <file>\n");
  printf(" -X         show extended options\n");
if(TRACE[0] || X_OPTS) {
  printf("         extended options\n");
  printf(" -C=<num>   for multi-channel input files: extract and filter\n");
  printf("            channel <num> (1 <= <num> <= %d  default: channel %d)\n",	\
	 FILT_I_CHANS, FILT_DEF_CHANNEL);
  printf(" -d         do not add dithering noise\n");
  printf("            (default: dither signals in 16-bit integer)\n");
  printf(" -g[=<num>] gain output signal to <num>%% of the full scale with\n");
  printf("            %.1f <= num <= %.1f; if <num> is not specified %.0f%%\n",	\
	 AF_MIN_GAIN, AF_MAX_GAIN, AF_DEF_GAIN);
  printf("            will be used\n");
  printf(" -of=<file> set name of output file to <file>\n");
  printf("            (only in single file mode; overrules -od option)\n");
}
if(TRACE[0]) {
  printf(" -Y<abc>    output trace info to 'stderr' for options <abc>\n");
  printf("            the following options are available:\n");
  printf("   'F'   append trace output to file '%s'\n", traceFile);
  printf("   'f'   write trace output to file '%s'\n", traceFile);
  printf("   'p'   filter parameters\n");
  printf("   'c'   filter coefficients\n");
  printf("   'm'   in- and output magnitudes etc.\n");
}
  printf("NOTE:\n");
  printf(" o   An existing output file will be overwritten without notice.\n");
  printf("[1]  For FIR filters the cut-off frequency specifies the endpoint\n");
  printf("     of the pass-band (in 0.X versions of this program it was the\n");
  printf("     -6 dB point in the centre of the transition band).\n");
  printf("     For IIR filters it specifies the -3 dB point.\n");
  printf("[2]  Processing time linearly increases with decreasing width of\n");
  printf("     the transition band and increasing stop-band attenuation.\n");
  printf("[3]  Each section adds 12 dB/oct to the slope of the filter.\n");
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
    while((n=fgetl(inpFile, PATH_MAX + 1, batFP, NULL)) == 0) NIX;
    if(n == EOF) {                                   /* we're through */
      fclose(batFP);
      batFP = NULL;
      return(0);
    }
  }
  else {
    if(fileNum >= numFiles)
      return(0);                                     /* we're through */
    strcpy(inpFile, argList[fileNum]);
    fileNum++;
  }
  /*
   * construct name of output file
   */
  if(fixOutFile != NULL) {
    strcpy(outFile, fixOutFile);
    fixOutFile = NULL;
  }
  else {
    parsepath(inpFile, &dPtr, &bPtr, NULL);
    if(anaOpts.options & AOPT_IN_DIR)
      strcpy(outFile, dPtr);       /* store output in input directory */
    else if(outDir != NULL) {
      strcpy(outFile, outDir);
      n = strlen(outDir);
      if(n > 0 && outDir[n-1] != DIR_SEP_CHR)
	strcat(outFile, DIR_SEP_STR);
    }
    else strcpy(outFile, "");           /* store in current directory */
    strcat(outFile, bPtr);   /* append base name */
    strcat(outFile, outExt); /* append extension */
  }
  if(strcmp(inpFile, outFile) == 0) {
    setAsspMsg(AEC_IO_CLASH, outFile);
    return(prtAsspMsg(NULL));
  }
  if(TRACE['F'] || TRACE['f'])
    fprintf(traceFP, "Input file: %s\nOutput file: %s\n",	\
	    inpFile, outFile);
  return(1); /* we got one */
}

