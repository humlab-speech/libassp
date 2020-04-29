/***********************************************************************
*                                                                      *
* This file is part of the Advanced Speech Signal Processor library.   *
*                                                                      *
* Copyright (C) 2003 - 2010  Michel Scheffers                          *
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
* File:     afdiff.c                                                   *
* Contents: Command line interface to audio differentiation functions. *
* Author:   Michel T.M. Scheffers                                      *
*                                                                      *
***********************************************************************/
/* $Id: afdiff.c,v 1.3 2010/03/30 14:38:05 mtms Exp $ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>    /* [f/s]printf() remove() stderr FILE NULL */
#include <stdlib.h>   /* strtod() strtol() */
#include <string.h>   /* str..() */
#include <ctype.h>    /* isdigit() isgraph() */
#include <math.h>     /* fabs() */

#include <miscdefs.h> /* TRUE FALSE LOCAL */
#include <misc.h>     /* mybasename() parsepath() fgetl() */
#include <mylimits.h> /* PATH_MAX NAME_MAX SUFF_MAX */
#include <diff.h>     /* processing parameters & DIFF functions */
#include <assp.h>     /* message & trace handler */
#include <asspfio.h>  /* file handler */
#include <asspdsp.h>  /* AF...GAIN */
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
  DOBJ *inpDOp, *outDOp;

  strcpy(progName, mybasename(argv[0]));
  strcpy(traceFile, progName);
  strcat(traceFile, TRACE_SUFFIX);
  initGlobals();
  if((err=evalArgs(argc, argv)) < 0)
    exit(1);
  if(err)
    exit(0);
  if(TRACE['F'] || TRACE['f']) {
    sprintf(info, "%s version %d.%d", progName, DIFF_MAJOR, DIFF_MINOR);
    openTrace(info);
  }
  else openTrace(NULL);
/*
 * loop over input files
 */
  inpDOp = outDOp = NULL;
  while((err=nextFile()) > 0) {
    inpDOp = asspFOpen(inpFile, AFO_READ, NULL);
    if(inpDOp == NULL) {
      err = prtAsspMsg(NULL);
      break;
    }
    outDOp = createDiff(inpDOp, &anaOpts);
    if(outDOp == NULL) {
      err = prtAsspMsg(NULL);
      asspFClose(inpDOp, AFC_FREE);
      inpDOp = NULL;
      break;
    }
    if(asspFOpen(outFile, AFO_WRITE, outDOp) == NULL) {
      err = prtAsspMsg(NULL);
      asspFClose(inpDOp, AFC_FREE);
      inpDOp = NULL;
      outDOp = freeDObj(outDOp);
      break;
    }
    if(diffSignal(inpDOp, &anaOpts, outDOp) < 0) {
      err = prtAsspMsg(NULL);
    }
    asspFClose(inpDOp, AFC_FREE);
    inpDOp = NULL;
    asspFClose(outDOp, AFC_FREE);
    outDOp = NULL;
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
  setDiffDefaults(&anaOpts);
  return;
}
/***********************************************************************
* check program options and arguments                                  *
***********************************************************************/
LOCAL int evalArgs(int argc, char *argv[])
{
  char *cPtr;
  int   i;

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
	case 'h':                                             /* help */
	  return(optError(NULL));
	  break;
	case 'b':                              /* backward difference */
	  anaOpts.options &= ~DIFF_OPT_CENTRAL;
	  anaOpts.options |= DIFF_OPT_BACKWARD;
	  cPtr++;
	  break;
	case 'c':                               /* central difference */
	  anaOpts.options &= ~DIFF_OPT_BACKWARD;
	  anaOpts.options |= DIFF_OPT_CENTRAL;
	  cPtr++;
	  break;
	case 'C':                                   /* channel number */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.channel = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.channel < 1 || anaOpts.channel > DIFF_I_CHANS)
	      return(optError("-C"));
	  }
	  else return(optError("-C"));
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
  printf("Release : %d.%d (%s)\n", DIFF_MAJOR, DIFF_MINOR, __DATE__);
  printf("Function: Computes the first difference of the signal in the audio-\n");
  printf("          formatted file(s) <file>. The differentiated signal will\n");
  printf("          be written to a file with the base name of the input file\n");
  printf("          and an extension consisting of '.d', followed by the\n");
  printf("          extension of the input file. The format of the output file\n");
  printf("          will be the same as that of the input file.\n");
  printf("          Differentiation can improve results an F0 analysis of e.g.\n");
  printf("          EGG signals because it removes a DC offset, attenuates\n");
  printf("          very low frequency components - and in the case of central\n");
  printf("          differentiation also very high ones - and enhances the\n");
  printf("          moment of glottal closure.\n");
  printf("Options :\n");
  printf(" -h/--help  print this text\n");
  printf(" -b         compute backward difference (s'[n] = s[n] - s[n-1])\n");
  printf("            (default: forward difference s'[n] = s[n+1] - s[n])\n");
  printf(" -c         compute central/interpolated/3-point difference\n");
  printf(" -C=<num>   for multi-channel input files: extract and differentiate\n");
  printf("            channel <num> (1 <= <num> <= %d  default: channel %d)\n",\
	 DIFF_I_CHANS, DIFF_DEF_CHANNEL);
  printf(" -od        store output file(s) in directory of input file(s)\n");
  printf("            (default: current working directory)\n");
  printf(" -od=<dir>  store output file(s) in directory <dir>\n");
  printf(" -ox=<ext>  set extension of output file names to <ext>\n");
  printf(" -z=<list>  batch mode: take input file names from file <list>;\n");
  printf("            overrules arguments <file>\n");
  printf(" -X         show extended options\n");
if(TRACE[0] || X_OPTS) {
  printf("         extended options\n");
  printf(" -of=<file> set name of output file to <file>\n");
  printf("            (only in single file mode; overrules -od option)\n");
}
if(TRACE[0]) {
  printf("         debugging options\n");
  printf(" -Y<abc>    output trace info to 'stderr' for options <abc>\n");
  printf("            the following options are available:\n");
  printf("   'F'   append trace output to file '%s'\n", traceFile);
  printf("   'f'   write trace output to file '%s'\n", traceFile);
  printf("   'p'   processing parameters\n");
}
  printf("NOTE:\n");
  printf(" o   An existing output file will be overwritten without notice.\n");
  printf(" o   Differentiation may result in numerical overflow. If that\n");
  printf("     would occur the output signal will automatically be compressed\n");
  printf("     to %.0f%% of the full scale.\n", AF_MAX_GAIN);
  printf(" o   Unless central differentiation is used the output will be\n");
  printf("     shifted by half a sample with respect to the input. For higher-\n");
  printf("     order differentiation one should therefore alternate forward\n");
  printf("     and backward differentiation.\n");
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
  char *dPtr, *bPtr, *xPtr;
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
    parsepath(inpFile, &dPtr, &bPtr, &xPtr);
    if(anaOpts.options & AOPT_IN_DIR)
      strcpy(outFile, dPtr);       /* store output in input directory */
    else if(outDir != NULL) {
      strcpy(outFile, outDir);
      n = strlen(outDir);
      if(n > 0 && outDir[n-1] != DIR_SEP_CHR)
	strcat(outFile, DIR_SEP_STR);
    }
    else strcpy(outFile, "");           /* store in current directory */
    strcat(outFile, bPtr);     /* append base name */
    if(strlen(outExt) > 0)
      strcat(outFile, outExt); /* append extension */
    else {
      strcat(outFile, ".d");
      xPtr++;                  /* skip period */
      strcat(outFile, xPtr);   /* append original extension */
    }
  }
  if(strcmp(inpFile, outFile) == 0) {
    setAsspMsg(AEC_IO_CLASH, outFile);
    return(prtAsspMsg(NULL));
  }
  if(TRACE['F'] || TRACE['f'])
    fprintf(traceFP, "Input file: %s\nOutput file: %s\n",\
	    inpFile, outFile);
  return(1); /* we got one */
}
