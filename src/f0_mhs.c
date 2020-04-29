/***********************************************************************
*                                                                      *
* This file is part of the Advanced Speech Signal Processor library.   *
*                                                                      *
* Copyright (C) 1999 - 2010  Michel Scheffers                          *
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
* File:     f0_mhs.c                                                   *
* Contents: Command line interface to Michel's/Modified Harmonic Sieve *
*           pitch analysis.                                            *
* Author:   Michel T.M. Scheffers                                      *
* Note:     See printMHSrefs() in 'mhs.c' for references.              *
*                                                                      *
***********************************************************************/
/* $Id: f0_mhs.c,v 1.13 2010/06/04 09:08:55 mtms Exp $ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>    /* [f/s]printf() remove() stderr FILE NULL */
#include <stddef.h>   /* size_t */
#include <stdlib.h>   /* strtod() strtol() */
#include <string.h>   /* str... */
#include <ctype.h>    /* isdigit() isgraph() */
#include <math.h>     /* log() pow() */

#include <miscdefs.h> /* TRUE FALSE LOCAL */
#include <mylimits.h> /* PATH_MAX NAME_MAX SUFF_MAX */
#include <misc.h>     /* mybasename() parsepath() fgetl() */
#include <mhs.h>      /* processing parameters & MHS functions */
#include <assp.h>     /* message & trace handler */
#include <asspfio.h>  /* file handler */
#include <asspdsp.h>  /* wfunc_e wfType() rfft() */
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
char  smpFile[PATH_MAX+1], pitFile[PATH_MAX+1], outExt[SUFF_MAX+1];
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
  DOBJ *smpDOp, *pitDOp;

  strcpy(progName, mybasename(argv[0]));
  strcpy(traceFile, progName);
  strcat(traceFile, TRACE_SUFFIX);
  initGlobals();
  if((err=evalArgs(argc, argv)) < 0)
    exit(1);
  if(err)
    exit(0);
  if(TRACE['F'] || TRACE['f']) {
    sprintf(info, "%s version %d.%d", progName, MHS_MAJOR, MHS_MINOR);
    openTrace(info);
  }
  else openTrace(NULL);
/*
 * loop over input files
 */
  smpDOp = pitDOp = NULL;
  while((err=nextFile()) > 0) {
    smpDOp = asspFOpen(smpFile, AFO_READ, NULL);
    if(smpDOp == NULL) {
      err = prtAsspMsg(NULL);
      break;
    }
    pitDOp = createMHS(smpDOp, &anaOpts);
    if(pitDOp == NULL) {
      err = prtAsspMsg(NULL);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(asspFOpen(pitFile, AFO_WRITE, pitDOp) == NULL) {
      err = prtAsspMsg(NULL);
      pitDOp = freeDObj(pitDOp);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(computeMHS(smpDOp, &anaOpts, pitDOp) == NULL) {
      err = prtAsspMsg(NULL);
    }
    asspFClose(smpDOp, AFC_FREE);
    smpDOp = NULL;
    asspFClose(pitDOp, AFC_FREE);
    pitDOp = NULL;
    if(err < 0) {
      remove(pitFile); /* contents invalid */
      break;
    }
    if(asspWarning)
      prtAsspMsg(NULL);
  }  /* END loop over input files */
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
  clrAsspMsg();
  batList = fixOutFile = outDir = NULL;
  strcpy(outExt, MHS_DEF_SUFFIX);
  X_OPTS = FALSE;
  batFP = NULL;
  setMHSdefaults(&anaOpts);
  return;
}
/***********************************************************************
* check program options and arguments                                  *
***********************************************************************/
int evalArgs(int argc, char *argv[])
{
  char *cPtr;
  int   i;
  
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
	    printMHSrefs();
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
	  }
	  else return(optError("-b"));
	  break;
	case 'e':                                         /* end time */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit((int)cPtr[1]) || cPtr[1] == '.')) {
	    cPtr++;
	    anaOpts.endTime = strtod(cPtr, &cPtr);
	  }
	  else return(optError("-e"));
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
	case 'C':                                   /* channel number */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	    cPtr++;
	    anaOpts.channel = (int)strtol(cPtr, &cPtr, 10);
	    if(anaOpts.channel < 1 || anaOpts.channel > MHS_I_CHANS)
	      return(optError("-C"));
	  }
	  else return(optError("-C"));
	  break;
	case 'g':                                           /* gender */
	  cPtr++;
	  if(*cPtr == '=' && cPtr[1] != EOS) {
	    cPtr++;
	    if(setMHSgenderDefaults(&anaOpts, *cPtr) < 0)
	      return(optError("-g"));
	    cPtr += strlen(cPtr);
	  }
	  else return(optError("-g"));
	  break;
	case 'M':                                       /* maximum F0 */
	  cPtr++;
	  if(*cPtr == '=' &&\
	     (isdigit((int)cPtr[1]) ||\
	      (cPtr[1] == '.' && isdigit((int)cPtr[2]))) ) {
	    cPtr++;
	    anaOpts.maxF = strtod(cPtr, &cPtr);
	  }
	  else return(optError("-M"));
	  break;
	case 'm':                                       /* minimum F0 */
	  cPtr++;
	  if(*cPtr == '=' &&\
	     (isdigit((int)cPtr[1]) ||\
	      (cPtr[1] == '.' && isdigit((int)cPtr[2]))) ) {
	    cPtr++;
	    anaOpts.minF = strtod(cPtr, &cPtr);
	    if(anaOpts.minF < MHS_ABSMIN_F0)
	      return(optError("-m"));
	  }
	  else return(optError("-m"));
	  break;
	case 'p':                             /* plain power spectrum */
	  cPtr++;
	  anaOpts.options |= MHS_OPT_POWER;
	  break;
	case 'o':           /* output modifiers */
	  cPtr++;
	  switch(*cPtr) {
	  case 'A':                                          /* ASCII */
	    strcpy(anaOpts.format, MHS_ASC_FORMAT);
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
	case 'v':           /* voicing parameters */
	  X_OPTS = TRUE;                          /* extended options */
	  cPtr++;
	  switch(*cPtr) {
	  case 'a':                              /* minimum magnitude */
	    cPtr++;
	    if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	      cPtr++;
	      anaOpts.voiMag = strtod(cPtr, &cPtr);
	    }
	    else return(optError("-va"));
	    break;
	  case 'c':                    /* 1st correlation coefficient */
	    cPtr++;
	    if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	      cPtr++;
	      anaOpts.voiAC1 = strtod(cPtr, &cPtr);
	    }
	    else return(optError("-vc"));
	    break;
	  case 'q':                                /* minimum Q value */
	    cPtr++;
	    if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	      cPtr++;
	      anaOpts.voiProb = strtod(cPtr, &cPtr);
	    }
	    else return(optError("-vq"));
	    break;
	  case 'r':                          /* minimum RMS amplitude */
	    cPtr++;
	    if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	      cPtr++;
	      anaOpts.voiRMS = strtod(cPtr, &cPtr);
	    }
	    else return(optError("-vr"));
	    break;
	  case 'z':                     /* maximum zero crossing rate */
	    cPtr++;
	    if(*cPtr == '=' && isdigit((int)cPtr[1])) {
	      cPtr++;
	      anaOpts.voiZCR = strtod(cPtr, &cPtr);
	    }
	    else return(optError("-vz"));
	    break;
	  default:
	    return(optError(--cPtr));
	  }
	  break;            /* end voicing parameters */
	case 'X':                            /* show extended options */
	  X_OPTS = TRUE;
	  return(optError(NULL));
	  break;
	case 'Y':                                    /* trace options */
	  TRACE[0] = TRUE;
	  cPtr++;
	  if(*cPtr == '=')
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
  if(anaOpts.maxF <= anaOpts.minF) {
    fprintf(stderr, "\nERROR: maximum pitch (%.1f Hz) too low\n",\
	    anaOpts.maxF);
    return(-1);
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
    fprintf(stderr, "ERROR: %s option: specification missing or invalid\n", opt);
  else
    fprintf(stderr, "ERROR: unknown option '-%s'\n", opt);
  return(-1);
}
/***********************************************************************
* print program syntax                                                 *
***********************************************************************/
void usage(void)
{
  printf("\n");
  printf("Syntax  : %s [<opts>] <file> [<opts>] {<file> [<opts>]}\n",\
	 progName);
  printf("Release : %d.%d (%s)\n", MHS_MAJOR, MHS_MINOR, __DATE__);
  printf("Function: Pitch analysis of the speech signal in <file> using\n");
  printf("          Michel's/Modified Harmonic Sieve algorithm.\n");
  printf("          Analysis results will be written to a file with the\n");
  printf("          base name of the input file and extension '%s'.\n",\
	 MHS_DEF_SUFFIX);
  printf("          Default output is in SSFF binary format (track '%s').\n",\
         MHS_SSFF_ID);
  printf("Options :\n");
  printf(" -h/--help  print this text\n");
  printf(" --refs     print references for the algorithm\n");
  printf(" -b=<time>  set begin of analysis interval to <time> seconds\n");
  printf("            (default: begin of file)\n");
  printf(" -e=<time>  set end of analysis interval to <time> seconds\n");
  printf("            (default: end of file)\n");
  printf(" -s=<dur>   set analysis window shift to <dur> ms\n");
  printf("            (default: %.1f)\n", MHS_DEF_SHIFT);
  printf(" -g=<code>  set gender-specific pitch ranges; <code> may be:\n");
  printf("            \"f[emale]\" (%.1f - %.1f Hz)\n",\
	 MHS_MINF0_f, MHS_MAXF0_f);
  printf("            \"m[ale]\" (%.1f - %.1f Hz)\n",\
	 MHS_MINF0_m, MHS_MAXF0_m);
  printf("            \"u[nknown]\" (default; %.1f - %.1f Hz)\n",\
	 MHS_MINF0_u, MHS_MAXF0_u);
  printf(" -M=<freq>  set maximum pitch value to <freq> Hz\n");
  printf(" -m=<freq>  set minimum pitch value to <freq> Hz\n");
  printf("            (absolute minimum: %.1f)\n", MHS_ABSMIN_F0);
  printf(" -p         use plain rather than masked power spectrum\n");
  printf(" -oA        output in XASSP ASCII format\n");
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
	 MHS_I_CHANS, MHS_DEF_CHANNEL);
  printf(" -of=<file> set name of output file to <file>\n");
  printf("            (only in single-file mode; overrules -od option)\n");
  printf(" -op=<num>  set precision of ASCII output to <num> digits\n");
  printf("            (default: %d)\n", MHS_DEF_DIGITS);
  printf("         thresholds for voicing decision\n");
  printf(" -va=<amp>  minimum signal amplitude (default: %.0f)\n",\
	 MHS_DEF_VOIMAG);
  printf(" -vc=<num>  minimum 1st correlation coefficient (default: %.3f)\n",\
	 MHS_DEF_VOIAC1);
  printf(" -vr=<num>  minimum RMS amplitude in dB (default: %.1f)\n",\
	 MHS_DEF_VOIRMS);
  printf(" -vz=<freq> maximum zero crossing rate in Hz (default: %.0f)\n",\
	 MHS_DEF_VOIZCR);
  printf(" -vq=<num>  minimum quality value of F0 fit (default: %.3f)\n",\
	 MHS_DEF_MINQVAL);
}
if(TRACE[0]) {
  printf("         debugging options\n");
  printf(" -Y<abc>    output trace info to 'stderr' for options <abc>\n");
  printf("            the following options are available:\n");
  printf("   'F'   append trace output to file '%s'\n", traceFile);
  printf("   'f'   write trace output to file '%s'\n", traceFile);
  printf("   'A'   analysis parameters\n");
  printf("   'v'   values used in voicing decision\n");
  printf("   'P'   values of detected peaks\n");
  printf("   'c'   F0 candidates and their quality values\n");
  printf("   't'   End parameters of tracks\n");
}
  printf("NOTE:\n");
  printf(" o   Existing output files will be overwritten without notice.\n");
  printf(" o   The analysis interval will be rounded to the nearest frame\n");
  printf("     boundary (integral multiple of the window shift).\n");
  printf(" o   The -g option overrules earlier specification of the pitch\n");
  printf("     range using the -M and/or -m option.\n");
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
    strcpy(pitFile, fixOutFile);
    fixOutFile = NULL;
  }
  else {
    parsepath(smpFile, &dPtr, &bPtr, NULL);
    if(anaOpts.options & AOPT_IN_DIR)
      strcpy(pitFile, dPtr);       /* store output in input directory */
    else if(outDir != NULL) {
      strcpy(pitFile, outDir);
      n = strlen(outDir);
      if(n > 0 && outDir[n-1] != DIR_SEP_CHR)
	strcat(pitFile, DIR_SEP_STR);
    }
    else strcpy(pitFile, "");           /* store in current directory */
    strcat(pitFile, bPtr);   /* append base name */
    strcat(pitFile, outExt); /* append extension */
  }
  if(strcmp(smpFile, pitFile) == 0) {
    setAsspMsg(AEC_IO_CLASH, pitFile);
    return(prtAsspMsg(NULL));
  }
  if(TRACE['F'] || TRACE['f'])
    fprintf(traceFP, "Input file: %s\nOutput file: %s\n",\
	    smpFile, pitFile);
  return(1); /* we got one */
}
