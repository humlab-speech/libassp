/***********************************************************************
*                                                                      *
* This file is part of the Advanced Speech Signal Processor library.   *
*                                                                      *
* Copyright (C) 1989 - 2011  Michel Scheffers                          *
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
* File:     f0_ksv.c                                                   *
* Contents: Command line interface to Kurt Schäfer-Vincent's           *
*           Periodicity Detection Algorithm.                           *
* Author:   Michel T.M. Scheffers                                      *
* Note:     See printKSVrefs() in 'ksv.c' for references.              *
*                                                                      *
***********************************************************************/
/* $Id: f0_ksv.c,v 1.4 2011/01/06 10:47:59 mtms Exp $ */

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
#include <ksv.h>      /* processing parameters & KSV functions */
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
char  smpFile[PATH_MAX+1], f0File[PATH_MAX+1], outExt[SUFF_MAX+1];
char  prdFile[PATH_MAX+1], prdExt[SUFF_MAX+1];
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
  DOBJ *smpDOp, *f0DOp, *prdDOp;

  strcpy(progName, mybasename(argv[0]));
  strcpy(traceFile, progName);
  strcat(traceFile, TRACE_SUFFIX);
  initGlobals();
  if((err=evalArgs(argc, argv)) < 0)
    exit(1);
  if(err)
    exit(0);
  if(TRACE['F'] || TRACE['f']) {
    sprintf(info, "%s version %d.%d", progName, KSV_MAJOR, KSV_MINOR);
    openTrace(info);
  }
  else openTrace(NULL);
/*
 * loop over input files
 */
  smpDOp = f0DOp = prdDOp = NULL;
  while((err=nextFile()) > 0) {
    smpDOp = asspFOpen(smpFile, AFO_READ, NULL);
    if(smpDOp == NULL) {
      err = prtAsspMsg(NULL);
      break;
    }
    f0DOp = createKSV(smpDOp, &anaOpts);
    if(f0DOp == NULL) {
      err = prtAsspMsg(NULL);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(asspFOpen(f0File, AFO_WRITE, f0DOp) == NULL) {
      err = prtAsspMsg(NULL);
      f0DOp = freeDObj(f0DOp);
      asspFClose(smpDOp, AFC_FREE);
      smpDOp = NULL;
      break;
    }
    if(anaOpts.options & KSV_OPT_PRD_OUT) {
      prdDOp = createPRD(smpDOp, &anaOpts);
      if(prdDOp == NULL) {
	err = prtAsspMsg(NULL);
	asspFClose(smpDOp, AFC_FREE);
	smpDOp = NULL;
	asspFClose(f0DOp, AFC_FREE);
	f0DOp = NULL;
	remove(f0File); /* contents invalid */
	break;
      }
      if(asspFOpen(prdFile, AFO_WRITE, prdDOp) == NULL) {
	err = prtAsspMsg(NULL);
	prdDOp = freeDObj(prdDOp);
	asspFClose(smpDOp, AFC_FREE);
	smpDOp = NULL;
	asspFClose(f0DOp, AFC_FREE);
	f0DOp = NULL;
	remove(f0File);
	break;
      }
    }
    if(computeKSV(smpDOp, &anaOpts, f0DOp, prdDOp) == NULL) {
      err = prtAsspMsg(NULL);
    }
    asspFClose(smpDOp, AFC_FREE);
    smpDOp = NULL;
    asspFClose(f0DOp, AFC_FREE);
    f0DOp = NULL;
    if(prdDOp != NULL) {
      asspFClose(prdDOp, AFC_FREE);
      prdDOp = NULL;
    }
    if(err < 0) {
      remove(f0File); /* contents invalid */
      if(strlen(prdFile) > 0)
	remove(prdFile);
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
  strcpy(outExt, KSV_DEF_SUFFIX);
  strcpy(prdExt, KSV_DEF_PRDEXT);
  X_OPTS = FALSE;
  batFP = NULL;
  setKSVdefaults(&anaOpts);
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
	    printKSVrefs();
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
	    if(anaOpts.channel < 1 || anaOpts.channel > KSV_I_CHANS)
	      return(optError("-C"));
	  }
	  else return(optError("-C"));
	  break;
	case 'g':                                           /* gender */
	  cPtr++;
	  if(*cPtr == '=' && cPtr[1] != EOS) {
	    cPtr++;
	    if(setKSVgenderDefaults(&anaOpts, *cPtr) < 0)
	      return(optError("-g"));
	    cPtr += strlen(cPtr);
	  }
	  else return(optError("-g"));
	  break;
	case 'M':                                       /* maximum F0 */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit(cPtr[1]) || (cPtr[1] == '.' && isdigit(cPtr[2]))) ) {
	    cPtr++;
	    anaOpts.maxF = strtod(cPtr, &cPtr);
	  }
	  else return(optError("-M"));
	  break;
	case 'm':                                       /* minimum F0 */
	  cPtr++;
	  if(*cPtr == '=' &&
	     (isdigit(cPtr[1]) || (cPtr[1] == '.' && isdigit(cPtr[2]))) ) {
	    cPtr++;
	    anaOpts.minF = strtod(cPtr, &cPtr);
	    if(anaOpts.minF < KSV_ABSMIN_F0)
	      return(optError("-m"));
	  }
	  else return(optError("-m"));
	  break;
	case 'o':           /* output modifiers */
	  cPtr++;
	  switch(*cPtr) {
	  case 'A':                                          /* ASCII */
	    strcpy(anaOpts.format, KSV_ASC_FORMAT);
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
	case 'p':           /* period markers */
	  X_OPTS = TRUE;                           /* extended option */
	  anaOpts.options |= KSV_OPT_PRD_OUT;
	  cPtr++;
	  if(*cPtr != EOS) {
	    switch(*cPtr) {
	    case 'M':                                   /* MIX format */
	    case 'm':
	      anaOpts.options |= KSV_OPT_PRD_MIX;
	      cPtr++;
	      break;
	    case 'x':                                    /* extension */
	      cPtr++;
	      if(*cPtr == '=') {
		cPtr++;
		if(*cPtr != EOS) {
		  if(*cPtr != '.') {
		    strcpy(prdExt, ".");
		    strcat(prdExt, cPtr);
		  }
		  else
		    strcpy(prdExt, cPtr);
		  cPtr += strlen(cPtr);
		}
		else return(optError("-px"));
	      }
	      else return(optError("-px"));
	      break;
	    default:
	      return(optError(--cPtr));
	    }
	  }
	  break;            /* END period markers */
	case 'v':           /* voicing parameters */
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
    fprintf(stderr, "\nERROR: maximum F0 (%.1f Hz) too low\n",\
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
  printf("Release : %d.%d (%s)\n", KSV_MAJOR, KSV_MINOR, __DATE__);
  printf("Function: F0 analysis of the speech signal in <file> using the\n");
  printf("          K. Schaefer-Vincent periodicity detection algorithm.\n");
  printf("          Analysis results will be written to a file with the\n");
  printf("          base name of the input file and extension '%s'.\n",\
	 KSV_DEF_SUFFIX);
  printf("          Default output is in SSFF binary format (track '%s').\n",\
         KSV_SSFF_ID);
  printf("          Optionally, location and type of the signal extrema on\n");
  printf("          which the F0 data are based, may be stored in a label\n");
  printf("          file. The name of this file will consist of the base\n");
  printf("          name of the F0 file and the extension '%s'.\n",\
	 KSV_DEF_PRDEXT);
  printf("Options :\n");
  printf(" -h/--help  print this text\n");
  printf(" --refs     print references for the algorithm\n");
  printf(" -b=<time>  set begin of analysis interval to <time> seconds\n");
  printf("            (default: begin of file)\n");
  printf(" -e=<time>  set end of analysis interval to <time> seconds\n");
  printf("            (default: end of file)\n");
  printf(" -s=<dur>   set analysis window shift to <dur> ms\n");
  printf("            (default: %.1f)\n", KSV_DEF_SHIFT);
  printf(" -g=<code>  set gender-specific F0 ranges; <code> may be:\n");
  printf("            \"f[emale]\" (%.1f - %.1f Hz)\n",\
	 KSV_MINF0_f, KSV_MAXF0_f);
  printf("            \"m[ale]\" (%.1f - %.1f Hz)\n",\
	 KSV_MINF0_m, KSV_MAXF0_m);
  printf("            \"u[nknown]\" (default; %.1f - %.1f Hz)\n",\
	 KSV_MINF0_u, KSV_MAXF0_u);
  printf(" -M=<freq>  set maximum F0 value to <freq> Hz\n");
  printf(" -m=<freq>  set minimum F0 value to <freq> Hz\n");
  printf("            (absolute minimum: %.1f)\n", KSV_ABSMIN_F0);
  printf(" -oA        output in XASSP ASCII format\n");
  printf(" -od        store output files in directory of input files\n");
  printf("            (default: current working directory)\n");
  printf(" -od=<dir>  store output files in directory <dir>\n");
  printf(" -ox=<ext>  set extension of output file names to <ext>\n");
  printf(" -va=<amp>  set amplitude threshold for voiced samples to <amp>\n");
  printf("            (default: %.0f)\n", KSV_DEF_VOIMAG); 
  printf(" -vz=<freq> set maximum zero crossing rate to <freq> Hz\n");
  printf("            (default: %.0f; set to zero to switch off)\n",\
	 KSV_DEF_VOIZCR);
  printf(" -z=<list>  batch mode: take input file names from file <list>;\n");
  printf("            overrules arguments <file>\n");
  printf(" -X         show extended options\n");
if(TRACE[0] || X_OPTS) {
  printf("         extended options\n");
  printf(" -C=<num>   for multi-channel input files: analyse channel <num>\n");
  printf("            with 1 <= <num> <= %d  (default: channel %d)\n",\
	 KSV_I_CHANS, KSV_DEF_CHANNEL);
  printf(" -of=<file> set name of output file to <file>\n");
  printf("            (only in single-file mode; overrules -od option)\n");
  printf(" -op=<num>  set precision of ASCII output to <num> digits\n");
  printf("            (default: %d)\n", KSV_DEF_DIGITS);
  printf(" -p         output period markers in ESPS xlabel format\n");
  printf(" -pM        as above but output in IPdS MIX format\n");
  printf(" -px=<ext>  set extension of marker file name(s) to <ext>\n");
}
if(TRACE[0]) {
  printf("         debugging options\n");
  printf(" -Y<abc>    output trace info to 'stderr' for options <abc>\n");
  printf("            the following options are available:\n");
  printf("   'F'   append trace output to file '%s'\n", traceFile);
  printf("   'f'   write trace output to file '%s'\n", traceFile);
  printf("   'A'   analysis parameters\n");
  printf("   'c'   current alive chain\n");
  printf("   'R'   periods rejected by ZCR test\n");
  printf("   'r'   twins rejected by AMV test\n");
  printf("   'x'   extension of twin buffer\n");
}
  printf("NOTE:\n");
  printf(" o   Existing output files will be overwritten without notice.\n");
  printf(" o   The analysis interval will be rounded to the nearest frame\n");
  printf("     boundary (integral multiple of the window shift).\n");
  printf(" o   The -g option overrules earlier specification of the F0\n");
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
   * construct name of output file(s)
   */
  if(fixOutFile != NULL) {
    strcpy(f0File, fixOutFile);
    fixOutFile = NULL;
  }
  else {
    parsepath(smpFile, &dPtr, &bPtr, NULL);
    if(anaOpts.options & AOPT_IN_DIR)
      strcpy(f0File, dPtr);        /* store output in input directory */
    else if(outDir != NULL) {
      strcpy(f0File, outDir);
      n = strlen(outDir);
      if(n > 0 && outDir[n-1] != DIR_SEP_CHR)
	strcat(f0File, DIR_SEP_STR);
    }
    else strcpy(f0File, "");            /* store in current directory */
    strcat(f0File, bPtr);   /* append base name */
    strcat(f0File, outExt); /* append extension */
  }
  if(strcmp(f0File, smpFile) == 0) {
    setAsspMsg(AEC_IO_CLASH, f0File);
    return(prtAsspMsg(NULL));
  }
  strcpy(prdFile, "");
  if(anaOpts.options & KSV_OPT_PRD_OUT) {
    parsepath(f0File, &dPtr, &bPtr, NULL);
    strcpy(prdFile, dPtr);
    strcat(prdFile, bPtr);
    strcat(prdFile, prdExt);
    if(strcmp(prdFile, smpFile) == 0 || strcmp(prdFile, f0File) == 0) {
      setAsspMsg(AEC_IO_CLASH, prdFile);
      return(prtAsspMsg(NULL));
    }
  }
  if(TRACE['F'] || TRACE['f']) {
    if(strlen(prdFile) > 0)
      fprintf(traceFP, "Input file: %s\nOutput files: %s and %s\n",\
	      smpFile, f0File, prdFile);
    else
      fprintf(traceFP, "Input file: %s\nOutput file: %s\n",\
	      smpFile, f0File);
  }
  return(1); /* we got one */
}
