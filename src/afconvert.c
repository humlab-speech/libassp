/***********************************************************************
*                                                                      *
* This file is part of the Advanced Speech Signal Processor program    *
* package.                                                             *
*                                                                      *
* Copyright (C) 1989 - 2010  Michel Scheffers                          *
*                            IPdS, CAU Kiel                            *
*                            Leibnizstr. 10                            *
*                            24118 Kiel, Germany                       *
*                            ms@ipds.uni-kiel.de                       *
*                                                                      *
* This package is free software: you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* This package is distributed in the hope that it will be useful,      *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with this package. If not, see <http://www.gnu.org/licenses/>. *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* File:     afconvert.c                                                *
* Contents: Program for converting the format and/or data contents of  *
*           audio files. 8-bit and unsigned data are automatically     *
*           expanded/converted to signed 16-bit linear PCM.            *
*           Optionally, some signal conditoning may be performed.      *
*           Output formats are restricted to the ones most common for  *
*           speech processing.                                         *
* Author:   Michel T.M. Scheffers                                      *
*                                                                      *
*-- revision history --------------------------------------------------*
*  0.0   derivate of 'hh'                                    MS 140103 *
*  0.1   pre-release                                         MS 240103 *
*  0.2   tested; minor modifications                         MS 050303 *
*  0.3   bug fix: -c option didn't extract channel 1 but only copied   *
*        initial numRecords samples                          MS 080403 *
*  0.4   added warning for non-zero start time and option for setting  *
*        output directory                                    MS 150503 *
*  0.5   moved data conditoning options to extended group; added       *
*        tapering; restructured file name generation for split-channel *
*        mode; included bug fix in generating AIFC header    MS 210703 *
*  0.6   'begRecNr' in DDESC replaced by 'startTime'; included fix for *
*        bug in reading/writing CSRE-ADF header              MS 300903 *
*  0.7   default output format fixed to WAV ;-(              MS 150404 *
*  0.8   added option to invert signal                       MS 180504 *
*  0.9   unified range definition: endTime/SmpNr 0 also EOF  MS 290904 *
*  1.0   reworked: use limited buffer size; added conversions of       *
*        8-bit and unsigned data; added re-muxing option     MS 211204 *
*  1.1   used "pathlims.h" for unified file lengths etc.     MS 250205 *
*  1.2   new winfuncs.c and asspdsp.h                        MS 190805 *
*  1.3   using fixed-size integers; adapted to new headers.[ch],       *
*        dataobj.[ch] and asspfio.[ch]                       MS 160606 *
*  1.4   adapted to changes in auconv.[ch]                   MS 091006 *
*  1.5   bug fix: asspFOpen() in READ mode, always opened file in      *
*        APPEND mode                                         MS 160207 *
*  1.6   adapted to new headers.[ch] and dataobj.[ch]        MS 170407 *
*  1.7   also use extension of 'fixOutFile' to guess format  MS 220507 *
*  2.0   adapted for linking with libmisc and libassp        MS 171007 *
*  2.1   optError() now local; extension via -ox option      MS 311007 *
*  2.2   included new makeWF() syntax                        MS 180408 *
*  2.3   adjusted to changes in aucheck.h                    MS 190109 *
*  2.4   added (hidden) output format 'KTH'                  MS 220610 *
*  2.5   unified extension of and info line in trace file    MS 090710 *
*                                                                      *
***********************************************************************/
/* $Id: afconvert.c,v 1.8 2010/07/09 07:38:43 mtms Exp $ */

#define PROG_MAJOR 2
#define PROG_MINOR 4

#include <stdio.h>      /* stdout sprintf() FILE NULL */
#include <stddef.h>     /* size_t */
#include <stdlib.h>     /* calloc() labs() */
#include <string.h>     /* str... */
#include <inttypes.h>   /* fixed size integer types */
#include <ctype.h>      /* isspace() isdigit() */
#include <math.h>       /* fabs() */

#include <misc.h>       /* system definitions and other handy defines */
#include <mylimits.h>   /* PATH_MAX NAME_MAX INT16_MIN/MAX*/
#include <assp.h>       /* message and trace handler */
#include <assptime.h>   /* standard conversion macros */
#include <asspfio.h>    /* file handler */
#include <asspdsp.h>    /* window functions */
#include <asspendian.h> /* byte order stuff */
#include <dataobj.h>    /* DOBJ DDESC & data object handling */
#include <headers.h>    /* file header handling */
#include <aucheck.h>    /* audio verification stuff */
#include <auconv.h>     /* alaw_t ulaw_t audio conversion functions */

/*
 * constants
 */
#define MAXFILES 2000      /* multiple input files */
#define MAXCHANS 8         /* maximum number of input channels */
#define MAXBUFBYTES 53760  /* use fixed buffer size (3*5*7*512) */

/*
 * option flags
 */
#define OPT_NONE       (0L)
#define OPT_ERROR      (-1L)
#define OPT_TREAT_RAW  0x00000001
#define OPT_IN_DIR     0x00000002
#define OPT_APPEND     0x00000004 /* not yet implemented */
#define OPT_OVERLAP    0x00000008 /* not yet implemented */
#define OPT_RE_MUX     0x00000010
#define OPT_DE_MUX     0x00000020
#define OPT_TO_MONO    0x00000040
#define OPT_MERGE      0x00000080
#define OPT_TAPER      0x00000100
#define OPT_REMOVE_DC  0x00000200
#define OPT_INVERT     0x00000400
#define OPT_AUTOGAIN   0x00001000
#define OPT_BALANCED   0x00002000
#define OPT_DITHER     0x00010000 /* not yet implemented */

/*
 * default option values
 */
#define DEF_FORM FF_WAVE   /* output file format */
#define DEF_EXT  ".wav"    /* output file extension */
#define DEF_TDUR 20.0      /* tapering duration in ms */
#define DEF_GAIN 100.0     /* % of full range for autogain */

/*
 * global arrays and variables
 */
char     progName[NAME_MAX+1];
char    *inpFile[MAXFILES], outDir[PATH_MAX+1];
char     ext[MAXCHANS][SUFF_MAX+1];
char    *fixOutFile=NULL;
void    *readBuf=NULL;
int16_t *workBuf=NULL;
int16_t *outBuf=NULL;
fform_e  outFormat;
int      ASC_OUT, X_OPTS, VERBOSE, ABORT;
int      mux[MAXCHANS], numMux, maxMuxNr;
int      numInpFiles, numOutFiles, numExt;
int      numInpChans, numWorkChans, numOutChans;
long     opts, auCaps, auProps;
long     inpBegRecNr, outBegRecNr, totRecords, maxBufRecs;
long     taperLen;
long     begSmpNr, endSmpNr;
double   begTime, endTime;
double   relMax, taperDur;
double  *tf=NULL;
ENDIAN   sysEndian, outEndian;
DOBJ    *idop, rawDO, outDO[MAXCHANS];

/*
 * prototypes of local functions
 */
LOCAL void    initGlobals(void);
LOCAL int     evalArgs(int argc, char *argv[]);
LOCAL int     optError(char *opt);
LOCAL void    usage(void);
LOCAL int     openFiles(int num);
LOCAL void    closeFiles(void);
LOCAL int     setGlobals(void);
LOCAL int     allocBufs(void);
LOCAL void    freeBufs(void);
LOCAL long    loadWorkBuf(long numRecords);
LOCAL long    saveWorkBuf(long numRecords);
LOCAL int     getExts(char *list);
LOCAL int     getMuxs(char *list);
LOCAL fform_e getFormat(char *str, char *ext);
LOCAL fform_e ext2format(char *ext);
LOCAL char   *format2ext(fform_e format);
LOCAL long    writeASC(void *buffer, long numRecords, DOBJ *dop);


/**********************************************************************/

int main(int argc, char *argv[])
{
  register int16_t *sPtr, *dPtr;
  register int    cn, n;
  register long   rn, orn, numRecs;
  char    info[128];
  int     SUM_CHANNELS, OVRFLW;
  int     err, fn;
  int     pass, numPasses;
  long    block, numBlocks, restRecs;
  long    endTaper;
  double  temp, sumMin, sumMax, ampMul[MAXCHANS];
  double  minAmp[MAXCHANS], maxAmp[MAXCHANS], mean[MAXCHANS];

  strcpy(progName, mybasename(argv[0]));
  strcpy(traceFile, progName);
  strcat(traceFile, TRACE_SUFFIX);
  initGlobals();
  if((err=evalArgs(argc, argv)) < 0)
    exit(1);
  if(err)
    exit(0);
  if(TRACE['F'] || TRACE['f']) {
    sprintf(info, "%s version %d.%d", progName, PROG_MAJOR, PROG_MINOR);
    openTrace(info);
  }
  else openTrace(NULL);
  if(TRACE['C']) {
    for(n = 0; n < argc; n++)
      fprintf(traceFP, "argv[%d] = %s\n", n, argv[n]);
  }
  ABORT = FALSE;
/*
 * Loop over input files.
 */
  for(fn = 0; fn < numInpFiles && !ABORT; fn++) {
    if((err=openFiles(fn)) != 0) {
      prtAsspMsg(NULL);
      if(err < 0 || ABORT)
	break;
    }
  /*
   * Read and process the data block-wise so we can handle 
   * very long files with very little memory load.
   */
    numBlocks = (totRecords / maxBufRecs) + 1;
    restRecs = totRecords % maxBufRecs; /* # of records in last block */
    if(restRecs == 0) {
      restRecs = maxBufRecs;
      numBlocks--;
    }
  /*
   * Because we read the data block-wise we may have to go
   * through the input file in more than one pass.
   */
    numPasses = 1;
    if(opts & (OPT_REMOVE_DC | OPT_AUTOGAIN | OPT_INVERT))
      numPasses++;     /* to gather statistics and check for overflow */
    SUM_CHANNELS = ((opts & OPT_TO_MONO) && numWorkChans > 1);
    if(SUM_CHANNELS)
      numPasses++;                              /* may cause overflow */
    OVRFLW = FALSE;
  /*
   * Initialize statistical variables (even when we do not need them).
   */
    sumMin = INT16_MAX;
    sumMax = INT16_MIN;
    for(cn = 0; cn < MAXCHANS; cn++) {
      minAmp[cn] = INT16_MAX;
      maxAmp[cn] = INT16_MIN;
      mean[cn] = 0.0;
      ampMul[cn] = 1.0;
    }
  /*
   * Position file pointer(s) of output files
   */
    for(cn = 0; cn < numOutFiles && !ABORT; cn++) {
      if(asspFSeek(&(outDO[cn]), outBegRecNr) < 0) {
	err = prtAsspMsg(NULL);
	ABORT = TRUE;
	break;
      }
    }
    if(ABORT) {
      closeFiles();
      freeBufs();
      break; /* leave loop over files */
    }
  /*
   * Loop over passes.
   */
    err = 0;                                    /* clear warning flag */
    for(pass = 1; pass <= numPasses && !ABORT; pass++) {
      if(pass == 1 || numBlocks > 1) {  /* (re-)position file pointer */
	if(asspFSeek(idop, inpBegRecNr) < 0) {
	  err = prtAsspMsg(NULL);
	  ABORT = TRUE;
	  break;
	}
      }
    /*
     * Loop over data blocks.
     */
      for(block = 1; block <= numBlocks && !ABORT; block++) {
	if(block >= numBlocks)
	  numRecs = restRecs;
	else
	  numRecs = maxBufRecs;
	if(pass == 1 || numBlocks > 1) {            /* (re-)load data */
	  if(loadWorkBuf(numRecs) < 0) {
	    err = prtAsspMsg(NULL);
	    ABORT = TRUE;
	    break;
	  }
	}
	if(pass < numPasses) {
        /*
         * Evaluate the sample values.
         */
	  if(pass == 1) {
	    if(opts & (OPT_REMOVE_DC | OPT_INVERT | OPT_AUTOGAIN)) {
	      for(sPtr = workBuf, rn = 0; rn < numRecs; rn++) {
		for(cn = 0; cn < numWorkChans; cn++, sPtr++) {
		  temp = (double)(*sPtr);
		  if(temp < minAmp[cn])
		    minAmp[cn] = temp;
		  if(temp > maxAmp[cn])
		    maxAmp[cn] = temp;
		  if(opts & OPT_REMOVE_DC)
		    mean[cn] += temp;
		}
	      }
	    }
	    else if(SUM_CHANNELS) {
	      for(sPtr = workBuf, rn = 0; rn < numRecs; rn++) {
		temp = 0.0;                  /* evaluate summed value */
		for(cn = 0; cn < numWorkChans; cn++, sPtr++)
		  temp += (double)(*sPtr);
		if(temp < sumMin)
		  sumMin = temp;
		if(temp > sumMax)
		  sumMax = temp;
	      }
	    }
	  }
	  else if(pass == 2) {
	    if(SUM_CHANNELS) {
	      for(sPtr = workBuf, rn = 0; rn < numRecs; rn++) {
		temp = 0.0; /* evaluate summed value after processing */
		for(cn = 0; cn < numWorkChans; cn++, sPtr++)
		  temp += (ampMul[cn] * ((double)(*sPtr) - mean[cn]));
		if(temp < sumMin)
		  sumMin = temp;
		if(temp > sumMax)
		  sumMax = temp;
	      }
	    } /* else no other options yet */
	  } /* else no more passes yet */
	}
	else {
        /*
	 * pass == numPasses : Do the actual processing if required.
	 */
	  orn = (block - 1) * maxBufRecs; /* relative output record # */
	  if(pass == 1) {
	    if(opts & OPT_TAPER) {        /* tapering only (in-place) */
	      endTaper = totRecords - taperLen - 1;
	      if(orn < taperLen || (orn + numRecs) > endTaper) {
		for(sPtr = workBuf, rn = 0; rn < numRecs; rn++, orn++) {
		  for(cn = 0; cn < numWorkChans; cn++, sPtr++) {
		    if(orn < taperLen && orn <= endTaper)
		      *sPtr = (int16_t)myrint((double)(*sPtr) * tf[orn]);
		    else if(orn >= taperLen && orn > endTaper)
		      *sPtr = (int16_t)myrint((double)(*sPtr) *\
					    tf[totRecords - orn - 1]);
		    /* else just count (most of the time) */
		  }
		}
	      }
	    } /* else reformatting only */
	  }
	  else if(pass == 2) {
	    if(opts & (OPT_REMOVE_DC | OPT_INVERT | OPT_AUTOGAIN)) {
	      for(sPtr = workBuf, rn = 0; rn < numRecs; rn++, orn++) {
		for(cn = 0; cn < numWorkChans; cn++, sPtr++) {
		  temp = ampMul[cn] * ((double)(*sPtr) - mean[cn]);
		  if(opts & OPT_TAPER) {
		    endTaper = totRecords - taperLen - 1;
		    if(orn < taperLen && orn <= endTaper)
		      temp *= tf[orn];
		    else if(orn >= taperLen && orn > endTaper)
		      temp *= tf[totRecords - orn -1];
		  }
		  *sPtr = (int16_t)myrint(temp);
		}
	      }
	    }
	    else if(SUM_CHANNELS) {
	      sPtr = dPtr = workBuf;
	      for(rn = 0; rn < numRecs; rn++, orn++, dPtr++) {
		temp = 0.0;
		for(cn = 0; cn < numWorkChans; cn++, sPtr++)
		  temp += (ampMul[cn] * (double)(*sPtr));
		if(opts & OPT_TAPER) {
		  endTaper = totRecords - taperLen - 1;
		  if(orn < taperLen && orn <= endTaper)
		    temp *= tf[orn];
		  else if(orn >= taperLen && orn > endTaper)
		    temp *= tf[totRecords - orn -1];
		}
		*dPtr = (int16_t)myrint(temp);
	      }
	    }
	  }
	  else { /* pass > 2 */
	    if(SUM_CHANNELS) {
	      sPtr = dPtr = workBuf;
	      for(rn = 0; rn < numRecs; rn++, orn++, dPtr++) {
		temp = 0.0;
		for(cn = 0; cn < numWorkChans; cn++, sPtr++)
		  temp += (ampMul[cn] * ((double)(*sPtr) - mean[cn]));
		if(opts & OPT_TAPER) {
		  endTaper = totRecords - taperLen - 1;
		  if(orn < taperLen && orn <= endTaper)
		    temp *= tf[orn];
		  else if(orn >= taperLen && orn > endTaper)
		    temp *= tf[totRecords - orn -1];
		}
		*dPtr = (int16_t)myrint(temp);
	      }
	    }
	  }
        /*
	 * Do the final processing like de-multiplexing
	 * and/or swapping and write the data.
	 */
	  if(ABORT)           /* catch signal before writing anything */
	    break;
	  if(saveWorkBuf(numRecs) < 0) {
	    err = prtAsspMsg(NULL);
	    ABORT = TRUE;
	    break;
	  }
	}
      }
    /*
     * END loop over blocks.
     */
      if(ABORT) /* leave loop over passes */
	break;
    /*
     * Evaluate statistical variables, check for overflow, etc.
     */
      if(pass == 1) {
	if(opts & (OPT_REMOVE_DC | OPT_INVERT | OPT_AUTOGAIN)) {
	  if(TRACE['s']) {
	    if(numWorkChans > 1) {
	      for(cn = 0; cn < numWorkChans; cn++)
		fprintf(traceFP, "   ch%d: min = %.0f  max = %.0f\n",\
			cn+1, minAmp[cn], maxAmp[cn]);
	    }
	    else
	      fprintf(traceFP, "   min = %.0f  max = %.0f\n",\
		      minAmp[0], maxAmp[0]);
	  }
	  if(opts & OPT_REMOVE_DC) {
	    /* compute mean and adjust amplitude range */
	    for(cn = 0; cn < numWorkChans; cn++) {
	      mean[cn] /= (double)totRecords;
	      minAmp[cn] -= mean[cn];     /* minimum after DC removal */
	      maxAmp[cn] -= mean[cn];     /* maximum after DC removal */
	    }
	    if(TRACE['s']) {
	      if(numWorkChans > 1) {
		for(cn = 0; cn < numWorkChans; cn++)
		  fprintf(traceFP, "   ch%d: mean = %.1f  min = %.1f  "\
			  "max = %.1f\n",\
			  cn+1, mean[cn], minAmp[cn], maxAmp[cn]);
	      }
	      else
		fprintf(traceFP, "   mean = %.0f  min = %.0f  "\
			"max = %.0f\n",	mean[0], minAmp[0], maxAmp[0]);
	    }
	  }
	  for(cn = 0; cn < numWorkChans; cn++) {
	    minAmp[cn] = fabs(minAmp[cn]);
	    maxAmp[cn] = fabs(maxAmp[cn]);          /* you never know */
	    if(minAmp[cn] > maxAmp[cn])
	      maxAmp[cn] = minAmp[cn];  /* now contains max magnitude */
	    if(maxAmp[cn] > INT16_MAX)
	      OVRFLW = TRUE;                    /* due to removing DC */
	  }
	  if((opts & OPT_BALANCED) ||
	     (OVRFLW && !(opts & OPT_AUTOGAIN))) {
	    /* preserve balance between channels */
	    temp = maxAmp[0];
	    for(cn = 1; cn < numWorkChans; cn++) { /* get largest mag */
	      if(maxAmp[cn] > temp)
		temp = maxAmp[cn];
	    }
	    for(cn = 0; cn < numWorkChans; cn++)/* set all to max mag */
	      maxAmp[cn] = temp;
	  }
	  if((opts & OPT_AUTOGAIN) || OVRFLW) {
	    /* compute multiplication factors */
	    for(cn = 0; cn < numWorkChans; cn++) {
	      ampMul[cn] = (double)(INT16_MAX-1) / maxAmp[cn];
	      if(opts & OPT_AUTOGAIN)
		ampMul[cn] *= (relMax/100.0);
	    }
	  }
	  if(opts & OPT_INVERT) {
	    /* change sign of multiplication factors */
	    for(cn = 0; cn < numWorkChans; cn++)
	      ampMul[cn] = -ampMul[cn];
	  }
	}
	else if(SUM_CHANNELS) {                 /* check for overflow */
	  if(TRACE['s']) {
	    fprintf(traceFP, "   sum: min = %.0f  max = %.0f\n",\
		    sumMin, sumMax);
	  }
	  sumMin = fabs(sumMin);
	  sumMax = fabs(sumMax);
	  if(sumMin > sumMax)
	    sumMax = sumMin;       /* max. magnitude of summed signal */
	  if(sumMax > INT16_MAX) {       /* set multiplication factors */
	    for(cn = 0; cn < numWorkChans; cn++)
	      ampMul[cn] = (double)(INT16_MAX-1) / sumMax;
	  }
	}
      }
      else if(pass == 2) {
	if(SUM_CHANNELS) {     /* check for overflow after processing */
	  if(TRACE['s']) {
	    fprintf(traceFP, "   sum: min = %.1f  max = %.1f\n",\
		    sumMin, sumMax);
	  }
	  sumMin = fabs(sumMin);
	  sumMax = fabs(sumMax);
	  if(sumMin > sumMax)
	    sumMax = sumMin;       /* max. magnitude of summed signal */
	  if(opts & OPT_AUTOGAIN)
	    sumMax *= (100.0/relMax);
	  if(opts & OPT_AUTOGAIN || sumMax > INT16_MAX) {
	    /* adjust multiplication factors */
	    for(cn = 0; cn < numWorkChans; cn++)
	      ampMul[cn] *= ((double)(INT16_MAX-1) / sumMax);
	  }
	}
      }
    }
  /*
   * END of loop over passes.
   */
    closeFiles();
    freeBufs();
    if(err < 0 || ABORT)
      break;                       /* leave for loop over input files */
  }
/*
 * END loop over input files.
 */
  clearDObj(&rawDO);
  closeTrace();
  if(err < 0 || ABORT)
    exit(1);
  exit(0);
}
/***********************************************************************
* initialize global parameters, pointers and file descriptors          *
***********************************************************************/
LOCAL void initGlobals(void)
{
  int n;

  opts = OPT_NONE;
  ASC_OUT = X_OPTS = VERBOSE = FALSE;
  SETENDIAN(sysEndian);
  SETENDIAN(outEndian);
  strcpy(outDir, "");
  strcpy(ext[0], "");
  numExt = 0;
  mux[0] = 0;
  numMux = 0;
  maxMuxNr = -1;
  numInpFiles = 0;
  fixOutFile = NULL;
  begTime = endTime = -1.0;
  begSmpNr = endSmpNr = -1;
  taperDur = DEF_TDUR;
  relMax = DEF_GAIN;
  readBuf = NULL;
  workBuf = outBuf = NULL;
  tf = NULL;
  outBegRecNr = 0;
  /* set audio capabilities of program */
  auCaps = MAXCHANS;
  auCaps |= (AUC_ALAW | AUC_uLAW | AUC_U8 | AUC_I8 | AUC_U16 | AUC_I16);
  auCaps |= AUC_MSB_X;

  idop = NULL;
  setRawSMP(&rawDO);
  for(n = 0; n < MAXCHANS; n++)
    initDObj(&(outDO[n]));
  clrAsspMsg();

  return;
}
/***********************************************************************
* check program options and arguments                                  *
***********************************************************************/
LOCAL int evalArgs(int argc, char *argv[])
{
  char *cPtr, *list, tmpExt[8];
  int   i, j, l;
  
  list = NULL;
  strcpy(tmpExt, "");
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
	  default:
	    return(optError(cPtr-2));
	  }
	  break;                /* END long option */
	case 'h':                                             /* help */
	  return(optError(NULL));
	case 'B':                                     /* begin sample */
	  cPtr++;
	  if(*cPtr == '=' &&  isdigit((int)cPtr[1])) {
	    cPtr++;
	    begSmpNr = strtol(cPtr, &cPtr, 10);
	  }
	  else return(optError("-B"));
	  break;
	case 'b':                                       /* begin time */
	  cPtr++;
	  if(*cPtr == '=' && (isdigit((int)cPtr[1]) || cPtr[1] == '.' )) {
	    cPtr++;
	    begTime = strtod(cPtr, &cPtr);
	  }
	  else return(optError("-b"));
	  break;
	case 'E':                                       /* end sample */
	  cPtr++;
	  if(*cPtr == '=' &&  isdigit((int)cPtr[1])) {
	    cPtr++;
	    endSmpNr = strtol(cPtr, &cPtr, 10);
	  }
	  else return(optError("-E"));
	  break;
	case 'e':                                         /* end time */
	  cPtr++;
	  if(*cPtr == '=' && (isdigit((int)cPtr[1]) || cPtr[1] == '.' )) {
	    cPtr++;
	    endTime = strtod(cPtr, &cPtr);
	  }
	  else return(optError("-e"));
	  break;
	case 'F':                                    /* output format */
	  cPtr++;
	  if(*cPtr == '=' &&  isalpha((int)cPtr[1])) {
	    cPtr++;
	    outFormat = getFormat(cPtr, tmpExt);
	    if(outFormat > 0)
	      cPtr += strlen(cPtr);
	    else return(optError("-F"));
	  }
	  else return(optError("-F"));
	  break;
	case 'c':                    /* extract/re-multiplex channels */
	  opts |= OPT_RE_MUX;
	  cPtr++;
	  if(*cPtr == '=' &&  isdigit((int)cPtr[1])) {
	    cPtr++;
	    list = cPtr;
	    cPtr += strlen(cPtr);
	    numMux = getMuxs(list);
	    if(numMux < 0)
	      return(prtAsspMsg(NULL));	
	    if(numMux == 0)
	      return(optError("-c"));
	  }
	  else return(optError("-c"));
	  break;
	case 'D':                              /* remove DC component */
	  X_OPTS = TRUE;                           /* extended option */
	  cPtr++;
	  if(*cPtr == 'C') {
	    opts |= OPT_REMOVE_DC;
	    cPtr++;
	  }
	  else return(optError(--cPtr));
	  break;
/* 	case 'd':                         dither signal if gain > 1 */
/* 	  X_OPTS = TRUE;                            extended option */
/* 	  opts |= OPT_DITHER; */
/* 	  cPtr++; */
/* 	  break; */
	case 'G':                                        /* auto gain */
	case 'g':
	  X_OPTS = TRUE;                           /* extended option */
	  opts |= OPT_AUTOGAIN;
	  if(*cPtr == 'G')
	    opts |= OPT_BALANCED;  /* retain balance between channels */
	  else
	    opts &= ~OPT_BALANCED;                      /* clear flag */
	  relMax = DEF_GAIN;
	  cPtr++;
	  if(*cPtr == '=') {
	    cPtr++;
	    if(isdigit((int)*cPtr)) {
	      relMax = strtod(cPtr, &cPtr);
	      if(relMax > 100.0) {
		if(opts & OPT_BALANCED)
		  return(optError("-G"));
		else
		  return(optError("-g"));
	      }
	    }
	    else {
	      if(opts & OPT_BALANCED)
		return(optError("-G"));
	      else
		return(optError("-g"));
	    }
	  }
	  break;
	case 'i':                   /* invert signal (multiply by -1) */
	  X_OPTS = TRUE;                           /* extended option */
	  opts |= OPT_INVERT;
	  cPtr++;
	  break;
	case 'm':             /* convert to mono by channel summation */
	  if(opts & OPT_DE_MUX) {
	    asspMsgNum = AEG_ERR_APPL;
	    sprintf(applMessage, "Can't combine split channels with "\
		    "sum channels option.");
	    return(prtAsspMsg(NULL));
	  }
	  opts |= OPT_TO_MONO;
	  cPtr++;
	  break;
	case 's':                                   /* split channels */
	  if(opts & OPT_TO_MONO) {
	    asspMsgNum = AEG_ERR_APPL;
	    sprintf(applMessage, "Can't combine sum channels with "\
		    "split channels option.");
	    return(prtAsspMsg(NULL));
	  }
	  opts |= OPT_DE_MUX;
	  cPtr++;
	  break;
	case 't':                                         /* tapering */
	  X_OPTS = TRUE;                           /* extended option */
	  opts |= OPT_TAPER;
	  cPtr++;
	  if(*cPtr == '=') {
	    cPtr++;
	    if(isdigit((int)*cPtr)) {
	      taperDur = strtod(cPtr, &cPtr);
	      if(taperDur == 0.0)
		opts &= ~OPT_TAPER;
	    }
	    else
	      return(optError("-t"));
	  }
	  break;
	case 'o':               /* output modifiers */
	  cPtr++;
	  switch(*cPtr) {
	  case 'd':                                      /* directory */
	    cPtr++;
	    if(*cPtr == '=') {
	      cPtr++;
	      if(*cPtr) {
		strcpy(outDir, cPtr);
		l = strlen(outDir);
		if(outDir[l-1] != DIR_SEP_CHR)
		  strcat(outDir, DIR_SEP_STR);
		cPtr += l;
		opts &= ~OPT_IN_DIR;                     /* overruled */
	      }
	      else return(optError("-od"));
	    }
	    else
	      opts |= OPT_IN_DIR;          /* directory of input file */
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
	  case 'm':                                            /* MSB */
	    X_OPTS = TRUE;                         /* extended option */
	    cPtr++;
	    if(*cPtr == '=' && isalpha((int)cPtr[1])) {
	      cPtr++;
	      switch(*cPtr) {
	      case 'B':                                 /* big endian */
	      case 'b':
	      case 'F':                                      /* first */
	      case 'f':
	      case 'M':                                   /* Motorola */
	      case 'm':
		SETMSBFIRST(outEndian);
		cPtr += strlen(cPtr);
		break;
	      case 'L':                              /* little endian */
	      case 'l':                                       /* last */
	      case 'I':                                      /* Intel */
	      case 'i':
		SETMSBLAST(outEndian);
		cPtr += strlen(cPtr);
		break;
	      default:
		return(optError("-om"));
	      }
	    }
	    else return(optError("-om"));
	    break;
	  case 'x':                                   /* extension(s) */
	    cPtr++;
	    if(*cPtr == '=' && cPtr[1]) {
	      cPtr++;
	      list = cPtr;
	      cPtr += strlen(cPtr);
	      numExt = getExts(list);
	      if(numExt < 0)
		return(prtAsspMsg(NULL));
	      if(numExt == 0)
		return(optError("-ox"));
	    }
	    else return(optError("-ox"));
	    break;
	  default:
	    return(optError(--cPtr));
	  }
	  break;                /* END output modifiers */
	case 'R':               /* raw format */
	  opts |= OPT_TREAT_RAW;
	  cPtr++;
	  if(*cPtr) {
	    switch(*cPtr) {
	    case 'c':                                  /* data coding */
	      cPtr++;
	      if(*cPtr == '=' && isalpha((int)cPtr[1])) {
		cPtr++;
		switch(*cPtr) {
		case 'A':                                    /* A-law */
		case 'a':
		  rawDO.ddl.format = DF_UINT8;
		  rawDO.ddl.coding = DC_ALAW;
		  rawDO.ddl.numBits = 8;
		  rawDO.ddl.zeroValue = 0;
		  cPtr += strlen(cPtr);
		  break;
		case 'I':                           /* signed integer */
		case 'i':
		  cPtr++;
		  l = strtol(cPtr, &cPtr, 10);          /* get length */
		  if(l < 8 || l > 16)
		    return(optError("-Rc"));
		  if(l == 8)
		    rawDO.ddl.format = DF_INT8;
		  else
		    rawDO.ddl.format = DF_INT16;
		  rawDO.ddl.coding = DC_PCM;
		  rawDO.ddl.numBits = (uint16_t)l;
		  rawDO.ddl.zeroValue = 0;
		  break;
		case 'U':                 /* u-law / unsigned integer */
		case 'u':
		  cPtr++;
		  if(*cPtr == 'L' || *cPtr == 'l') {
		    rawDO.ddl.format = DF_UINT8;
		    rawDO.ddl.coding = DC_uLAW;
		    rawDO.ddl.numBits = 8;
		    rawDO.ddl.zeroValue = 0;
		    cPtr += strlen(cPtr);
		  }
		  else {
		    l = strtol(cPtr, &cPtr, 10);        /* get length */
		    if(l < 8 || l > 16)
		      return(optError("-Rc"));
		    if(l == 8)
		      rawDO.ddl.format = DF_UINT8;
		    else
		      rawDO.ddl.format = DF_UINT16;
		    rawDO.ddl.coding = DC_PCM;
		    rawDO.ddl.numBits = (uint16_t)l;
		    rawDO.ddl.zeroValue = 1L << (l - 1);
		    break;
		  }
		  break;
		default:
		  return(optError("-Rc"));
		}
	      }
	      else return(optError("-Rc"));
	      break;
	    case 'h':                                  /* header size */
	      cPtr++;
	      if(*cPtr == '=' && isdigit((int)cPtr[1])) {
		cPtr++;
		rawDO.headerSize = strtol(cPtr, &cPtr, 10);
	      }
	      else return(optError("-Rh"));
	      break;
	    case 'm':                                          /* MSB */
	      cPtr++;
	      if(*cPtr == '=' && isalpha((int)cPtr[1])) {
		cPtr++;
		switch(*cPtr) {
		case 'B':       /* big endian */
		case 'b':
		case 'F':       /* first */
		case 'f':
		case 'M':       /* Motorola */
		case 'm':
		  SETMSBFIRST(rawDO.fileEndian);
		  cPtr += strlen(cPtr);
		  break;
		case 'L':       /* last/little endian */
		case 'l':
		case 'I':       /* Intel */
		case 'i':
		  SETMSBLAST(rawDO.fileEndian);
		  cPtr += strlen(cPtr);
		  break;
		default:
		  return(optError("-Rm"));
		}
	      }
	      else return(optError("-Rm"));
	      break;
	    case 'n':                           /* number of channels */
	      cPtr++;
	      if(*cPtr == '=' && isdigit((int)cPtr[1])) {
		cPtr++;
		rawDO.ddl.numFields = strtol(cPtr, &cPtr, 10);
	      }
	      else return(optError("-Rn"));
	      break;
	    case 's':                           /* sampling frequency */
	      cPtr++;
	      if(*cPtr == '=' && isdigit((int)cPtr[1])) {
		cPtr++;
		rawDO.sampFreq = strtod(cPtr, &cPtr);
	      }
	      else return(optError("-Rs"));
	      break;
	    default:
	      return(optError(--cPtr));
	    }
	  }
	  break;                /* END raw format */
        case 'X':                            /* show extended options */
          X_OPTS = TRUE;
          return(optError(NULL));
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
    }                           /* END OPTION */
    else {                      /* ARGUMENT */
      if(numInpFiles < MAXFILES) {
	inpFile[numInpFiles] = cPtr;
	numInpFiles++;
      }
      else {
	fprintf(stderr, "\nERROR: more than %d input files\n", MAXFILES);
	return(-1);
      }
    }                           /* END ARGUMENT */
  }

  if(numInpFiles == 0)          /* naked call */
    return(optError(NULL));
  if(fixOutFile != NULL) {
    if(numInpFiles > 1) {
      asspMsgNum = AEG_ERR_APPL;
      sprintf(applMessage, "Specification of output file only "\
	      "possible for a single input file.");
      return(prtAsspMsg(NULL));
    }
    parsepath(fixOutFile, NULL, NULL, &cPtr);
    strcpy(ext[0], cPtr);
    strcpy(outDir, "");                                  /* overruled */
    opts &= ~OPT_IN_DIR;
  }

  if(opts & OPT_DE_MUX) {
    if(fixOutFile) {
      if(strcmp(fixOutFile, "stdout") == 0) {
	asspMsgNum = AEG_ERR_APPL;
	sprintf(applMessage, "Can't combine split channels option "\
		"with screen output.");
	return(prtAsspMsg(NULL));
      }
    }
    if(numExt == 0) {                /* use channel code as extension */
      for(i = 0; i < MAXCHANS; i++)            /* fill out completely */
	sprintf(ext[i], ".ch%d", i+1);
      numExt = MAXCHANS;  /* number of different extensions available */
      if(outFormat < 0 || outFormat == FF_UNDEF)
	outFormat = DEF_FORM;                  /* take default format */
    }
    else if(numExt == 1) {/* signal: append channel code to base name */
      for(i = 1; i < MAXCHANS; i++)
	strcpy(ext[i], ext[0]);                          /* duplicate */
    }
    else {
      i = numExt;
      while(i < MAXCHANS) {       /* add extensions with channel code */
	sprintf(ext[i], ".ch%d", i+1);
	for(j = 0; j < i; j++) {                            /* verify */
	  if(strcmp(ext[j], ext[i]) == 0)
	    break;                                   /* stop at clash */
	  i++;
	}
	numExt = i;
      }
    }
  }

  if(fixOutFile == NULL && numExt <= 0)
    strcpy(ext[0], tmpExt);
  if(outFormat <= FF_UNDEF) {                   /* evaluate extension */
    if(numExt <= 0 && fixOutFile == NULL) {    /* take default format */
      outFormat = DEF_FORM;
      strcpy(ext[0], DEF_EXT);
    }
    else {
      outFormat = ext2format(ext[0]);
      if(outFormat <= FF_ERROR) {
	asspMsgNum = AEG_ERR_APPL;
	sprintf(applMessage, "Can't determine output format");
	return(prtAsspMsg(NULL));
      }
    }
  }
  if(strlen(ext[0]) == 0) {                /* take default for format */
    cPtr = format2ext(outFormat);
    if(cPtr == NULL) {
      asspMsgNum = AEG_ERR_BUG;
      sprintf(applMessage, "Invalid output format");
      return(prtAsspMsg(NULL));
    }
    strcpy(ext[0], cPtr);
    numExt = 1;
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
  printf("Release : %d.%d (%s)\n", PROG_MAJOR, PROG_MINOR, __DATE__);
  printf("Function: Converts format and/or data of the audio file(s) <file>.\n");
  printf("          A-law and u-Law compressed data will automatically be\n");
  printf("          expanded to 16-bit linear PCM; unsigned data converted\n");
  printf("          to signed.\n");
  printf("          Output will be written to a file with the base name of\n");
  printf("          the input file and the extension corresponding to the\n");
  printf("          output format (see -F option) or the one defined with\n");
  printf("          the -x option. If neither of these options is used, the\n");
  printf("          output format defaults to %s and the extension to '%s'.\n",\
	 (DEF_FORM == FF_WAVE) ? "WAV" : "AIFF", DEF_EXT);
  printf("          If only the extension has been defined the program will\n");
  printf("          try to guess the desired output format on that basis.\n");
  printf("          When splitting multi-channel files, the construction of\n");
  printf("          the names of the output files is as follows:\n");
  printf("          When no extension has been defined, '.ch1', '.ch2', etc.\n");
  printf("          will be used. If one extension has been specified, the\n");
  printf("          channel code ('_l' c.q. '_r' for stereo files and '_1',\n");
  printf("          '_2', etc. for multi-channel files) will be appended to\n");
  printf("          the base name of the output files. Alternatively, you may\n");
  printf("          specify a series of different extensions to be used.\n");
  printf("          Apart from the output formats listed below, the program\n");
  printf("          can handle files in CSL, CSRE, KTH/SWELL and SSFF format.\n");
/*   printf("          can handle input files in CSL, KTH/SWELL and SSFF format.\n"); */
  printf("          The -R options can be used to convert headerless files or\n");
  printf("          those with a faulty header or in an unsupported format.\n");
  printf("Options:\n");
  printf(" -h/--help  print this text\n");
  printf(" -F=<type>  set format of output file to <type> which may be:\n");
  printf("            AIF[F]        (default extension '.aif')\n");
  printf("            AIFC          (default extension '.afc')\n");
  printf("            AU   / SND    (default extension '.au'  / '.snd')\n");
/*   printf("            CSRE          (default extension '.adf')\n"); */
/*   printf("            KTH           (default extension '.smp')\n"); */
  printf("            NIST / SPHERE (default extension '.nst' / '.sph')\n");
  printf("            WAV[E]        (default extension '.wav')\n");
  printf("            RAW  / PCM    (default extension '.raw' / '.pcm')\n");
  printf("            ASC           (default extension '.asc')\n");
  printf(" -b=<time>  extract data starting from <time> seconds\n");
  printf("            (default: begin of data)\n");
  printf(" -B=<num>   as above, but begin at sample number <num>\n");
  printf(" -e=<time>  extract data up to <time> seconds\n");
  printf("            (default: end of data)\n");
  printf(" -E=<num>   as above, but end at sample number <num>\n");
  printf(" -c=<num>   extract channel number <num> (1 <= <num> <= %d)\n",\
	 MAXCHANS);
  printf(" -m         convert multi-channel file to mono by summation\n");
  printf(" -s         split multi-channel file into separate mono files\n");
  printf(" -od        store output file in directory of input file\n");
  printf("            (default: current working directory)\n");
  printf(" -od=<dir>  store output file(s) in directory <dir>\n");
  printf(" -ox=<ext>  set extension of output file name(s) to <ext>\n");
  printf(" -ox=<ext1>,<ext2>,...   use extensions <ext1>, <ext2>, etc. for\n");
  printf("            the output files when splitting channels\n");
  printf(" -R         ignore header of input file (if present) and use the\n");
  printf("            values below for interpreting the file/data format\n");
  printf(" -Rh=<num>  header size in byte (default: %ld)\n",\
	 rawDO.headerSize);
  printf(" -Rc=<type> coding of the data; <type> may be: \n");
  printf("            Al[aw]  8-bit A-law compression\n");
  printf("            ul[aw]  8-bit u-law compression\n");
  printf("            i<num>  <num>-bit signed linear PCM (8 <= <num> <= 16)\n");
  printf("            u<num>  <num>-bit unsigned linear PCM (<num> as above)\n");
  printf("            (default: i16)\n");
  printf(" -Rm=<type> MSB (byte order/endian) of data; <type> may be: \n");
  printf("            f[irst]/b[ig]/M[otorola] or l[ast]/l[ittle]/I[ntel]\n");
  printf("            default is the system's byte order (here MSB %s)\n",\
	 MSBFIRST(outEndian) ? "first" : "last");
  printf(" -Rn=<num>  number of channels (default: %i)\n",\
	 (int)(rawDO.ddl.numFields));
  printf(" -Rs=<freq> sampling frequency in Hz (default: %.1f)\n",\
	 rawDO.sampFreq);
  printf(" -X         show extended options (mainly signal conditioning)\n");
if(TRACE[0] || X_OPTS) {
  printf("         extended options\n");
  printf(" -c=<num1>,<num2>,...   convert multi-channel file to multi-channel\n");
  printf("            file containing channels <num1>, <num2>, etc.\n");
  printf(" -DC        remove DC component (mean) of signal\n");
  printf(" -g         adjust signal amplitude to span the full 16 bit range\n");
  printf(" -g=<num>   adjust signal amplitude to span <num> percent of the range\n");
  printf(" -G[=<num>] as above but preserve loudness balance between channels\n");
  printf(" -i         invert signal (change sign)\n");
  printf(" -t[=<dur>] taper first and last <dur> ms of the signal with a raised\n");
  printf("            cosine; if <dur> is not specified, a value of %.0f ms will\n", DEF_TDUR);
  printf("            be taken\n");
  printf(" -of=<file> set name of output file to <file>\n");
  printf("            (only in single-file mode; overrules -od option)\n");
  printf(" -om=<type> set MSB of the output data to <type> (only for NIST-SPHERE\n");
  printf("            and RAW format; see -Rm option for <type> and default)\n");
}
if(TRACE[0]) {
  printf("         debugging options\n");
  printf(" -Y<abc>    output trace info to 'stderr' for options <abc>\n");
  printf("            the following options are available:\n");
  printf("   'F'   append trace output to file '%s'\n", traceFile);
  printf("   'f'   write trace output to file '%s'\n", traceFile);
  printf("   'C'   items of parsed command line\n");
  printf("   's'   signal statistics (min, max, mean)\n");
}
  printf("NOTE:\n");
  printf(" o   An existing output file will be overwritten without notice.\n");
  printf(" o   The -B/E options overrule the -b/e options.\n");
  printf(" o   The -m and -DC options may result in signal magnitudes\n");
  printf("     exceeding the 16 bit range. If that is the case, overall\n");
  printf("     signal level will be reduced to avoid numerical overflow.\n");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf("This program is free software under the GNU General Public License.\n");
  printf("It comes with ABSOLUTELY NO WARRANTY.\n");
  printf("For details see the file COPYING provided with this program or visit\n");
  printf("http://www.gnu.org/licenses/\n");
  printf("\n");
  return;
}
/***********************************************************************
* open input file and get header; set global variables; allocate memory*
* for buffers; create output file(s) and write header(s)               *
***********************************************************************/
LOCAL int openFiles(int num)
{
  static char outFile[MAXCHANS][PATH_MAX+1];
  char  *bPtr, *cPtr, *dPtr, *xPtr, app[MAXCHANS][6];
  int    err, n, on, TO_SCREEN;
  DOBJ  *odop;
  DDESC *dd;

  clrAsspMsg();
  err = 0;
  idop = NULL;
  if(opts & OPT_TREAT_RAW) {
    idop = &rawDO;
    setRecordSize(idop);
  }

  if(TRACE['F'] || TRACE['f']) {
    fprintf(traceFP, "Input file: %s\n", inpFile[num]);
  }
  
  idop = asspFOpen(inpFile[num], AFO_READ, idop);
  if(idop == NULL) {
    return(-1);
  }
    /* channel count in mux[] starts at zero ! vv ! */
  if((auProps=checkSound(idop, auCaps, maxMuxNr+1)) <= 0) {
    closeFiles();
    return(-1);
  }
  if(idop->Start_Time != 0.0) {
    asspMsgNum = AWG_WARN_APPL;
    sprintf(applMessage, "File %s\n         has non-zero start time "\
	    "(this information will be lost).", idop->filePath);
    prtAsspMsg(NULL);
  }
  if(setGlobals() < 0) {
    closeFiles();
    return(-1);
  }
  if(allocBufs() < 0) {
    closeFiles();
    return(-1);
  }
/*
 * Construct output file name(s).
 */
  TO_SCREEN = FALSE;
  if(num == 0 && fixOutFile) {
    if(strcmp(fixOutFile, "stdout") == 0) {
      TO_SCREEN = TRUE;
      cPtr = fixOutFile;
      strcpy(ext[0], "");
    }
    else {      /* split file path so that we can insert channel code */
      parsepath(fixOutFile, &cPtr, &bPtr, &xPtr);
      strcat(cPtr, bPtr);               /* directory path + base name */
      strcpy(ext[0], xPtr);           /* to make life easier later on */
    }
  }
  else {                          /* get directory path and base name */
    parsepath(inpFile[num], &dPtr, &cPtr, NULL);
    if(opts & OPT_IN_DIR)
      strcpy(outDir, dPtr);        /* store output in input directory */
  }
  if(opts & OPT_DE_MUX) {
    if(numExt < 2) {                   /* insert channel code in name */
      if(numOutFiles < 3) {                            /* stereo file */
	strcpy(app[0], "_l");
	strcpy(app[1], "_r");
      }
      else {
	for(on = 0; on < numOutFiles; on++)
	  sprintf(app[on], "_%d", on+1);
      }
    }
    for(on = 0; on < numOutFiles; on++) {
      strcpy(outFile[on], outDir);
      strcat(outFile[on], cPtr);
      if(numExt < 2)
	strcat(outFile[on], app[on]);
      strcat(outFile[on], ext[on]);
      for(n = 0; n < on; n++) {                      /* verify unique */
	if(strcmp(outFile[n], outFile[on]) == 0) {
	  asspMsgNum = AEG_ERR_APPL;
	  sprintf(applMessage, "Identical file names (%s)\n"\
		  "       for different channels", outFile[on]);
	  closeFiles();
	  freeBufs();
	  return(-1);
	}
      }
    }
  }
  else {
    strcpy(outFile[0], outDir);
    strcat(outFile[0], cPtr);
    strcat(outFile[0], ext[0]);
  }
  
  if(TRACE['F'] || TRACE['f']) {
    fprintf(traceFP, "Output file: %s\n", outFile[0]);
    for(on = 1; on < numOutFiles; on++)
      fprintf(traceFP, "             %s\n", outFile[on]);
  }

  for(on = 0; on < numOutFiles; on++) {
    if(strcmp(inpFile[num], outFile[on]) == 0) {
      setAsspMsg(AEC_IO_CLASH, outFile[on]);
      closeFiles();
      freeBufs();
      return(-1);
    }
  }
/*
 * Fill out header of first file and write to output. We hereby abuse
 * the fact that putHeader() will correct/expand the file descriptor.
 */
  odop = &(outDO[0]);
  copyDObj(odop, idop);           /* will clear 'odop' before copying */
  odop->filePath = outFile[0];
  odop->fileFormat = outFormat;
  CPYENDIAN(odop->fileEndian, outEndian);  /* adjusted in putHeader() */
  odop->version = 0;
  odop->headerSize = 0;
  odop->startRecord = 0;
  odop->numRecords = totRecords;
  dd = &(odop->ddl);
  dd->format = DF_INT16;/* set standard output data format and coding */
  dd->coding = DC_PCM;
  dd->numBits = 16;
  dd->zeroValue = 0;
  dd->numFields = numOutChans;
  setStart_Time(odop);
  setRecordSize(odop);

  if(TO_SCREEN) {            /* have to write header (if any) by hand */
    odop->fp = stdout;
    if((err=putHeader(odop)) != 0) {
      if(err < 0) {
	closeFiles();
	freeBufs();
	return(err);
      }
      prtAsspMsg(NULL);
    }
  }
  else { /* with odop defined, asspFOpen() will also write the header */
    if(asspFOpen(outFile[0], AFO_WRITE, odop) == NULL) {
      closeFiles();
      freeBufs();
      return(-1);
    }
  }
/*
 * Now create the remaining files in split-channel mode.
 * Note that the output file can't be stdout in this case.
 */
  for(on = 1; on < numOutFiles; on++) {
    odop = &(outDO[on]);
    copyDObj(odop, &(outDO[0]));    /* file formats will be identical */
    if(asspFOpen(outFile[on], AFO_WRITE, odop) == NULL) {
      closeFiles();
      freeBufs();
      return(-1);
    }
  }

  return(0);
}
/***********************************************************************
* close all open files; clear/free data objects                        *
***********************************************************************/
LOCAL void closeFiles()
{
  int n, action;

  action = AFC_NONE;
  if(idop) {
    if(opts & OPT_TREAT_RAW)           /* refers to fixed data object */
      action = AFC_KEEP;
    else
      action = AFC_FREE;           /* refers to allocated data object */
    asspFClose(idop, action);
    idop = NULL;
  }
  action = AFC_CLEAR;              /* refers to internal data objects */
  for(n = 0; n < MAXCHANS; n++)
    asspFClose(&(outDO[n]), action);
  return;
}
/***********************************************************************
* set file-dependent global variables.                                 *
***********************************************************************/
LOCAL int setGlobals(void)
{
  long begRecord, endRecord;
/*
 * Convert copy interval to records and clip.
 * NOTE: beg/endSmpNr and beg/endTime refer to
 *       absolute time, NOT to begin of data
 */
  if(idop->startRecord < 0)
    return(-1);

  if(begSmpNr <= 0) {
    if(begTime > 0.0)
      begRecord = TIMEtoSMPNR(begTime, idop->sampFreq);
    else
      begRecord = idop->startRecord;
  }
  else
    begRecord = begSmpNr;
  if(begRecord < idop->startRecord)
    begRecord = idop->startRecord;
  if(endSmpNr <= 0) {
    if(endTime > 0.0)
      endRecord = TIMEtoSMPNR(endTime, idop->sampFreq);
    else
      endRecord = idop->startRecord + idop->numRecords;
  }
  else
    endRecord = endSmpNr;
  if(endRecord > (idop->startRecord + idop->numRecords))
    endRecord = idop->startRecord + idop->numRecords;
  /* set the global variables */
  inpBegRecNr = begRecord;
  totRecords = endRecord - begRecord;
  if(totRecords <= 0) {
    setAsspMsg(AED_NO_DATA, idop->filePath);
    return(-1);
  }
  outBegRecNr = 0; /* future enhancement */
/*
 * Set number of channels at various stages of processing 
 * and the number of output files per input file
 */
  numInpChans = idop->ddl.numFields;
  if(opts & OPT_RE_MUX)
    numWorkChans = numMux;
  else
    numWorkChans = numInpChans;
  if(opts & (OPT_DE_MUX | OPT_TO_MONO))
    numOutChans = 1;
  else
    numOutChans = numWorkChans;
  numOutFiles = 1;
  if(opts & OPT_DE_MUX) {
    numOutFiles = numWorkChans;
    if(numExt > 1 && numExt < numOutFiles) {
      asspMsgNum = AEG_ERR_APPL;
      sprintf(applMessage, "Insufficient extensions specified for "\
	    "splitting\n       file %s.", idop->filePath);
      return(-1);
    }
  }
  return(0);
}
/***********************************************************************
* allocate memory for buffers; always allocate memory for the output   *
* data format (int16_t) so we can map buffers                          *
***********************************************************************/
LOCAL int allocBufs(void)
{
  size_t recSize;

  readBuf = NULL;
  workBuf = outBuf = NULL;
  tf = NULL;

  recSize = (size_t)numInpChans * sizeof(int16_t);
  maxBufRecs = (long)(MAXBUFBYTES / recSize);
  if(((double)totRecords / (double)maxBufRecs) <= 1.2)
    maxBufRecs = totRecords;            /* don't waste memory or time */
  readBuf = calloc((size_t)maxBufRecs, recSize);
  if(auProps & AUC_BYTE_MASK || opts & OPT_RE_MUX) {
    recSize = (size_t)numWorkChans * sizeof(int16_t);
    workBuf = (int16_t *)calloc((size_t)maxBufRecs, recSize);
    if(opts & OPT_DE_MUX)
      outBuf = (int16_t *)readBuf;
    else
      outBuf = workBuf;
  }
  else {
    workBuf = (int16_t *)readBuf;
    if(opts & OPT_DE_MUX)
      outBuf = (int16_t *)calloc((size_t)maxBufRecs, sizeof(int16_t));
    else
      outBuf = workBuf;
  }
  if(readBuf == NULL || workBuf == NULL || outBuf == NULL) {
    setAsspMsg(AEG_ERR_MEM, NULL);
    freeBufs();
    return(-1);
  }
  if(opts & OPT_TAPER) {
    taperLen = TIMEtoSMPNR(taperDur/1000.0, idop->sampFreq);
    if(taperLen < 2)
      taperLen = 2;                    /* at least one sample reduced */
    tf = makeWF(WF_COS_2, (2*taperLen)+1, 0);    /* centre value == 1 */
    if(tf == NULL) {
      setAsspMsg(AEG_ERR_MEM, NULL);
      freeBufs();
      return(-1);
    }
  }
  return(0);
}
/***********************************************************************
* Return memory allocated for buffers. Watch for mappings!             *
***********************************************************************/
LOCAL void freeBufs(void)
{
  if(readBuf) {
    if((void *)workBuf == readBuf)
      workBuf = NULL;
    if((void *)outBuf == readBuf)
      outBuf = NULL;
    free(readBuf);
    readBuf = NULL;
  }
  if(workBuf) {
    if(outBuf == workBuf)
      outBuf = NULL;
    free(workBuf);
    workBuf = NULL;
  }
  if(outBuf) {
    free(outBuf);
    outBuf = NULL;
  }
  if(tf) {
    freeWF(tf);
    tf = NULL;
  }
  return;
}
/***********************************************************************
* Sequentially read 'numRecords' from input file; expand 8-bit codes;  *
* re-mux if desired. Result is 'numRead' * 'numWorkChans' multiplexed  *
* 16-bit integers in 'workBuf'.                                        *
* BEWARE: This funtion needs serious reworking when we allow anything  *
* else than 8- or 16-bit data in the input file.                       *
***********************************************************************/
LOCAL long loadWorkBuf(long numRecords)
{
  register alaw_t     *aPtr;
  register ulaw_t     *uPtr;
  register binoff8_t  *bPtr;
  register binoff16_t *wPtr;
  register int8_t     *cPtr;
  register int16_t    *sPtr, *dPtr, numBits;
  register long        i, n, numValues;
  long format, numRead=0;

  numRead = asspFRead(readBuf, numRecords, idop);
  if(numRead > 0) {
    numValues = numRead * numInpChans;
    if((auProps & AUC_SWAP_MASK) && \
       DIFFENDIAN(idop->fileEndian, sysEndian))
      memswab(readBuf, readBuf, sizeof(int16_t), (size_t)numValues);
    format = auProps & AUC_FORM_MASK;
    numBits = idop->ddl.numBits;
    dPtr = workBuf;                   /* destination always 'workBuf' */
    if(opts & OPT_RE_MUX)
      numValues = numRead;
    switch(format) {
    case AUC_ALAW:
      aPtr = (alaw_t *)readBuf;
      if(opts & OPT_RE_MUX) {
	for(n = 0; n < numValues; n++) {
	  for(i = 0; i < numMux; i++)
	    *(dPtr++) = alaw_to_int16(aPtr[mux[i]]);
	  aPtr += numInpChans;
	}
      }
      else {
	for(n = 0; n < numValues; n++)
	  *(dPtr++) = alaw_to_int16(*(aPtr++));
      }
      break;
    case AUC_uLAW:
      uPtr = (ulaw_t *)readBuf;
      if(opts & OPT_RE_MUX) {
	for(n = 0; n < numValues; n++) {
	  for(i = 0; i < numMux; i++)
	    *(dPtr++) = ulaw_to_int16(uPtr[mux[i]]);
	  uPtr += numInpChans;
	}
      }
      else {
	for(n = 0; n < numValues; n++)
	  *(dPtr++) = ulaw_to_int16(*(uPtr++));
      }
      break;
    case AUC_U8:
      bPtr = (binoff8_t *)readBuf;
      if(opts & OPT_RE_MUX) {
	for(n = 0; n < numValues; n++) {
	  for(i = 0; i < numMux; i++)
	    *(dPtr++) = binoff8_to_int16(bPtr[mux[i]]);
	  bPtr += numInpChans;
	}
      }
      else {
	for(n = 0; n < numValues; n++)
	  *(dPtr++) = binoff8_to_int16(*(bPtr++));
      }
      break;
    case AUC_I8:
      cPtr = (int8_t *)readBuf;
      if(opts & OPT_RE_MUX) {
	for(n = 0; n < numValues; n++) {
	  for(i = 0; i < numMux; i++)
	    *(dPtr++) = (int16_t)cPtr[mux[i]];
	  cPtr += numInpChans;
	}
      }
      else {
	for(n = 0; n < numValues; n++)
	  *(dPtr++) = (int16_t)(*(cPtr++));
      }
      break;
    case AUC_U16:
      wPtr = (binoff16_t *)readBuf;
      if(opts & OPT_RE_MUX) {
	for(n = 0; n < numValues; n++) {
	  for(i = 0; i < numMux; i++)
	    *(dPtr++) = binoff16_to_int16(wPtr[mux[i]], numBits);
	  wPtr += numInpChans;
	}
      }
      else {
	for(n = 0; n < numValues; n++)
	  *(dPtr++) = binoff16_to_int16(*(wPtr++), numBits);
      }
      break;
    case AUC_I16:
      sPtr = (int16_t *)readBuf;
      if(opts & OPT_RE_MUX) {
	for(n = 0; n < numValues; n++) {
	  for(i = 0; i < numMux; i++)
	    *(dPtr++) = sPtr[mux[i]];
	  sPtr += numInpChans;
	}
      }
      else {
	if(dPtr != sPtr) /* shouldn't be the case, but ... */
	  memcpy(dPtr, sPtr, (size_t)numValues * sizeof(int16_t));
      }
      break;
    default:
      asspMsgNum = AEG_ERR_BUG;
      sprintf(applMessage, "Can't handle data format in file %s",\
	      idop->filePath);
      return(-1);
    }
  }
  return(numRead);
}
/***********************************************************************
* sequentially write 'numRecords' to output file(s). Function includes *
* swapping and de-muxing if necessary.                                 *
***********************************************************************/
LOCAL long saveWorkBuf(long numRecords)
{
  register int16_t *sPtr, *dPtr;
  register int      i, SWAP;
  register long     n;
  size_t numShorts;
  long   numWrite=0;
  DOBJ  *odop;

  odop = &(outDO[0]);
  SWAP = DIFFENDIAN(sysEndian, odop->fileEndian);
  if(!(opts & OPT_DE_MUX)) {
    numShorts = (size_t)(numRecords * numOutChans);
    if(SWAP)
      memswab(outBuf, workBuf, sizeof(int16_t), numShorts);
    else if(outBuf != workBuf)
      memcpy(outBuf, workBuf, numShorts * sizeof(int16_t));
    if(ASC_OUT)
      numWrite = writeASC(outBuf, numRecords, odop);
    else
      numWrite = asspFWrite(outBuf, numRecords, odop);
  }
  else {       /* de-multiplex; numOutFiles should equal numWorkChans */
    for(i = 0; i < numOutFiles; i++) {
      odop = &(outDO[i]);
      sPtr = &(workBuf[i]);
      for(dPtr = outBuf, n = 0; n < numRecords; dPtr++, n++) {
	if(SWAP)
	  memswab(dPtr, sPtr, sizeof(int16_t), 1);
	else
	  *dPtr = *sPtr;
	sPtr += numOutFiles;
      }
      if(ASC_OUT)
	numWrite = writeASC(outBuf, numRecords, odop);
      else
	numWrite = asspFWrite(outBuf, numRecords, odop);
      if(numWrite < 0)
	break;
    }
  }
  return(numWrite);
}
/***********************************************************************
* convert string with comma separated extensions to table              *
***********************************************************************/
LOCAL int getExts(char *list)
{
  char *field[MAXCHANS+1];
  int   num, i, j;
  
  num = 0;
  if(list) {
    num = strsplit(list, ',', field, MAXCHANS+1);
    /* NOTE:  + 1 to handle overspecification */
    if(num < 0 || num > MAXCHANS)
      num = MAXCHANS;
    if(num > 0) {
      for(i = 0; i < num; i++) {
	if(strlen(field[i]) == 0) {
	  asspMsgNum = AEG_ERR_APPL;
	  sprintf(applMessage, "Extension missing for channel %d.\n",\
		  i + 1);
	  return(-1);
	}
	if(strlen(field[i]) > SUFF_MAX - 1) {
	  asspMsgNum = AEG_ERR_APPL;
	  sprintf(applMessage, "Extension too long (%s)\n"\
		  "maximally %d characters.", field[i], SUFF_MAX - 1);
	  return(-1);
	}
	for(i = 0; i < num; i++) {                   /* copy to array */
	  if(field[i][0] != '.') {
	    strcpy(ext[i], ".");
	    strcat(ext[i], field[i]);
	  }
	  else
	    strcpy(ext[i], field[i]);
	  for(j = 0; j < i; j++) {                        /* verify */
	    if(strcmp(ext[j], ext[i]) == 0) {
	      asspMsgNum = AEG_ERR_APPL;
	      sprintf(applMessage, "Identical extensions (%s) "\
		      "for different channels.", ext[i]);
	      return(-1);
	    }
	  }
	}
      }
    }
  }
  return(num);
}
/***********************************************************************
* convert string with comma separated channel numbers to table         *
***********************************************************************/
LOCAL int getMuxs(char *list)
{
  char *field[MAXCHANS], *cPtr;
  int   num, i;
  
  num = 0;
  if(list) {
    num = strsplit(list, ',', field, MAXCHANS);
    if(num < 0) {
      asspMsgNum = AEG_ERR_APPL;
      sprintf(applMessage, "-c option: too many channels specified "\
	      "(maximally %d)", MAXCHANS);
      return(-1);
    }
    if(num > 0) {
      for(i = 0; i < num; i++) {
	if(strlen(field[i]) == 0) {
	  asspMsgNum = AEG_ERR_APPL;
	  sprintf(applMessage, "-c option: channel number missing");
	  return(-1);
	}
	mux[i] = strtol(field[i], &cPtr, 10);
	if(*cPtr != EOS) {
	  asspMsgNum = AEG_ERR_APPL;
	  sprintf(applMessage, "-c option: invalid specification");
	  return(-1);
	}
	if(mux[i] < 1 || mux[i] > MAXCHANS) {
	  asspMsgNum = AEG_ERR_APPL;
	  sprintf(applMessage, "-c option: invalid channel number\n"\
		  "       (must be between 1 and %d)", MAXCHANS);
	  return(-1);
	}
	mux[i] -= 1;             /* internal count starts at zero !!! */
      }
      maxMuxNr = -1;     /* channel count in mux[] starts at zero !!! */
      for(i = 0; i < num; i++) {
	if(mux[i] > maxMuxNr)
	  maxMuxNr = mux[i];
      }
    }
  }
  return(num);
}
/***********************************************************************
* convert string to format code and set default extension              *
***********************************************************************/
LOCAL fform_e getFormat(char *str, char *ext)
{
  fform_e format;

  format = FF_ERROR;
  ASC_OUT = FALSE;
  if(strxcmp(str, "AIFC") == 0) {
    format = FF_AIFC;
    if(ext != NULL) strcpy(ext, ".afc");
  }
  else if(strnxcmp(str, "AIF", 3) == 0) {
    format = FF_AIFF;
    if(ext != NULL) strcpy(ext, ".aif");
  }
  else if(strxcmp(str, "AU") == 0) {
    format = FF_AU;
    if(ext != NULL) strcpy(ext, ".au");
  }
  else if(strnxcmp(str, "SND", 2) == 0) {
    format = FF_SND;
    if(ext != NULL) strcpy(ext, ".snd");
  }
  else if(strnxcmp(str, "CSRE", 3) == 0 || strxcmp(str, "ADF") == 0) {
    format = FF_CSRE;
    if(ext != NULL) strcpy(ext, ".adf");
  }
  else if(strnxcmp(str, "KTH", 2) == 0) {
    format = FF_KTH;
    if(ext != NULL) strcpy(ext, ".smp");
  }
  else if(strnxcmp(str, "NIST", 2) == 0) {
    format = FF_NIST;
    if(ext != NULL) strcpy(ext, ".nst");
  }
  else if(strnxcmp(str, "SPH", 2) == 0) {
    format = FF_SPHERE;
    if(ext != NULL) strcpy(ext, ".sph");
  }
  else if(strnxcmp(str, "WAV", 1) == 0) {
    format = FF_WAVE;
    if(ext != NULL) strcpy(ext, ".wav");
  }
  else if(strxcmp(str, "PCM") == 0) {
    format = FF_RAW;
    if(ext != NULL) strcpy(ext, ".pcm");
  }
  else if(strxcmp(str, "RAW") == 0) {
    format = FF_RAW;
    if(ext != NULL) strcpy(ext, ".raw");
  }
  else if(strxcmp(str, "ASC") == 0) {
    format = FF_RAW;
    ASC_OUT = TRUE;
    if(ext != NULL) strcpy(ext, ".asc");
  }
  return(format);
}
/***********************************************************************
* guess file format given extension                                    *
***********************************************************************/
LOCAL fform_e ext2format(char *ext)
{
  fform_e format;

  format = FF_ERROR;
  ASC_OUT = FALSE;
  while(*ext == '.')
    ext++;
  if(strxcmp(ext, "afc") == 0 || strxcmp(ext, "aifc") == 0) {
    format = FF_AIFC;
  }
  else if(strnxcmp(ext, "aif", 3) == 0) {
    format = FF_AIFF;
  }
  else if(strxcmp(ext, "au") == 0) {
    format = FF_AU;
  }
  else if(strxcmp(ext, "snd") == 0) {
    format = FF_SND;
  }
/*   else if(strxcmp(ext, "adf") == 0) { */
/*     format = FF_CSRE; */
/*   } */
/*   else if(strxcmp(ext, "smp") == 0) { */
/*     format = FF_KTH; */
/*   } */
  else if(strxcmp(ext, "nst") == 0 || strxcmp(ext, "nist") == 0) {
    format = FF_NIST;
  }
  else if(strnxcmp(ext, "sph", 3) == 0) {
    format = FF_SPHERE;
  }
  else if(strnxcmp(ext, "wav", 3) == 0) {
    format = FF_WAVE;
  }
  else if(strxcmp(ext, "pcm") == 0 || strxcmp(ext, "raw") == 0) {
    format = FF_RAW;
  }
  else if(strxcmp(ext, "asc") == 0) {
    format = FF_RAW;
    ASC_OUT = TRUE;
  }
  return(format);
}
/***********************************************************************
* return pointer to default extension for file format                  *
***********************************************************************/
LOCAL char *format2ext(fform_e format)
{
  switch(format) {
  case FF_AIFC:
    return(".afc");
  case FF_AIFF:
    return(".aif");
/*   case FF_CSRE: */
/*     return(".adf"); */
/*   case FF_KTH: */
/*     return(".smp"); */
   case FF_SPHERE:
    return(".sph");
  case FF_AU:
    return(".au");
  case FF_WAVE:
    return(".wav");
  case FF_RAW:
    if(!ASC_OUT)
      return(".raw");
    else
      return(".asc");
  default:
    return(NULL);
  }
  return(NULL);
}

/*DOC

Function 'writeASC'

Preliminary ! print in ASCII

NOTE: This function does not change file descriptor items. 

DOC*/

LOCAL long writeASC(void *buffer, long numRecords, DOBJ *dop)
{
  register int16_t *sPtr;
  register long   n, N, i;

  if(TRACE[0]) {
    if(dop == NULL || buffer == NULL || numRecords < 0) {
      setAsspMsg(AEB_BAD_ARGS, "writeAsc");
      return(-1);
    }
    if(dop->fp == NULL || dop->recordSize < 1) {
      setAsspMsg(AEB_BAD_CALL, "writeAsc");
      return(-1);
    }
  }
  if(numRecords > 0) {
    clearerr(dop->fp);    /* because we'll have to test on these later */
    sPtr = (int16_t *)buffer;
    N = dop->ddl.numFields;
    for(n = 0; n < numRecords; n++) {
      fprintf(dop->fp, "%6d", *(sPtr++));
      for(i = 1; i < N; i++)
	fprintf(dop->fp, " %6d", *(sPtr++));
      fprintf(dop->fp, "\n");
    }
    if(feof(dop->fp) || ferror(dop->fp)) {
      setAsspMsg(AEF_ERR_WRIT, dop->filePath);
      return(-1);
    }
  }
/*   clrAsspMsg(); */
  return(numRecords);
}
