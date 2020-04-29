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
* File:     ffinfo.c                                                   *
* Contents: Program to determine format and contents of data files.    *
* Author:   Michel T.M. Scheffers                                      *
*                                                                      *
*-- revision history --------------------------------------------------*
* 0.0   first implementation                                 MS 221098 *
* 0.1   data types etc. in lower case                        MS 231098 *
* 0.2   adapted to new headers.[ch]                          MS 031298 *
* 0.3   adapted to new headers.[ch]                          MS 101298 *
* 0.4   adapted to new headers.[ch]                          MS 250299 *
* 0.5   added NIST-SPHERE format                             MS 150800 *
* 0.6   extended data types and formats                      MS 240101 *
* 0.7   allowed renaming of program                          MS 280301 *
* 1.0   renamed in ffinfo; added non-audio formats           MS 080501 *
* 1.1   installed message handler; adapted to new headers    MS 140801 *
* 1.2   used function setAsspMsg() where appropriate         MS 040901 *
* 1.3   removed copying of label header from getHeader()     MS 161101 *
* 1.4   added SSFF header                                    MS 220302 *
* 1.5   extended data formats                                MS 120402 *
* 2.0   adjusted to new file/data descriptor                 MS 290702 *
* 2.1   added start time for audio and spectrum parameters   MS 211102 *
* 2.2   fixed to compile under MINGW                         MS 240203 *
* 2.3   adapted to data descriptor item 'startTime'          MS 250903 *
* 2.4   minor modifications in argument evaluation           MS 190804 *
* 2.5   added signal duration in output                      MS 181004 *
* 2.6   adapted to new pathlims.h                            MS 280205 *
* 2.7   adapted to new dataobj.[ch] and headers.[ch]         MS 120406 *
* 2.8   adapted to new dataobj.[ch] and headers.[ch]         MS 120307 *
* 3.0   adapted to be linked with libmisc and libassp        MS 170807 *
* 3.1   added support for extended RIFF-Wave format          MS 260208 *
* 3.2   added description for AC1, LP1 and PROB              MS 290509 *
* 3.3   prepared for option to get audio amplitude range     MS 091009 *
* 3.4   adapted to new dataobj.[ch] and headers.[ch]         MS 050110 *
* 3.5   added data type DT_FTSQR; DT_CEP => DT_FTCEP         MS 140610 *
*                                                                      *
***********************************************************************/
/* $Id: ffinfo.c,v 1.12 2010/07/01 14:19:32 mtms Exp $ */

#define PROG_MAJOR 3
#define PROG_MINOR 5

#include <stdio.h>      /* printf() FILE NULL */
#include <stdlib.h>     /* exit */
#include <string.h>     /* str... */

#include <misc.h>       /* TRUE FALSE LOCAL mybasename() */
#include <assp.h>       /* message and trace handler */
#include <mylimits.h>   /* PATH_MAX NAME_MAX */
#include <dataobj.h>    /* DOBJ DDESC */
#include <headers.h>    /* getHeader() and codes */
#include <ipds_lbl.h>   /* MIX_SFR */

/*
 * constants
 */
#define MAXFILES 2000 /* multiple input files */

/*
 * global arrays and variables
 */
char  progName[NAME_MAX+1];
char *argList[MAXFILES];
int   numFiles, X_OPTS;

/*
 * prototypes of local functions
 */
LOCAL int  evalArgs(int argc, char *argv[]);
LOCAL int  optError(char *opt);
LOCAL void usage(void);
LOCAL void smpInfo(DOBJ *dop);
LOCAL void tagInfo(DOBJ *dop);
LOCAL void asspInfo(DOBJ *dop);
LOCAL void ssffInfo(DOBJ *dop);
LOCAL void uwmInfo(DOBJ *dop);

/**********************************************************************/
int main(int argc, char *argv[])
{
  int    n, err;
  DOBJ   dobj, *dop;
  DDESC *dd;

  strcpy(progName, mybasename(argv[0]));
  if((err=evalArgs(argc, argv)) < 0)
    exit(1);
  if(err)
    exit(0);

  dop = &dobj;
  initDObj(dop);
  for(n = 0; n < numFiles; n++) {
    printf("\n");
    clrAsspMsg();
    dop->filePath = argList[n];
    dop->fp = fopen(dop->filePath, "rb");
    if(dop->fp == NULL) {
      setAsspMsg(AEF_ERR_OPEN, dop->filePath);
      prtAsspMsg(NULL);
      continue;
    }
    err = getHeader(dop);
    fclose(dop->fp);
    if(err && asspMsgNum != AEF_BAD_FORM) {
      prtAsspMsg(NULL);
      if(err < 0) {
	clearDObj(dop);
	continue;
      }
    }

    printf("file_name\t%s\n", dop->filePath);
    dd = &(dop->ddl);
    switch(dd->type) {
    case DT_SMP:
      smpInfo(dop);
      break;
    case DT_TAG:
    case DT_LBL:
    case DT_MRK:
    case DT_EPO:
    case DT_PRD:
      tagInfo(dop);
      break;
    default:                                       /* parametric data */
      switch(dop->fileFormat) {
      case FF_XASSP:
	asspInfo(dop);
	break;
      case FF_SSFF:
	ssffInfo(dop);
	break;
      case FF_UWM:
	uwmInfo(dop);
	break;
      default:
	printf("file_format\tUNKNOWN\n");
	printf("data_format\t");
	fflush(stdout);
	if(dop->fileData == FDF_ASC)
	  printf("ASCII\n");
	else if(dop->fileData == FDF_BIN)
	  printf("BINARY\n");
	else
	  printf("UNKNOWN\n");
	err = TRUE;
      }
      break;
    }
    clearDObj(dop);
  }
  printf("\n");
  exit(0);
}

/***********************************************************************
* check program options and arguments                                  *
***********************************************************************/
LOCAL int evalArgs(int argc, char *argv[])
{
  char   *cPtr;
  int     i;

  X_OPTS = FALSE;
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
	case 'X':                            /* show extended options */
	  X_OPTS = TRUE;
	  return(optError(NULL));
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

  if(numFiles == 0) /* naked call */
    return(optError(NULL));
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
/**********************************************************************/
LOCAL void usage(void)
{
  printf("\n");
  printf("Syntax  : %s [<opts>] <file> {<file>}\n", progName);
  printf("Release : %d.%d (%s)\n", PROG_MAJOR, PROG_MINOR, __DATE__);
  printf("Function: Attempts to determine the format of the file(s) <file>.\n");
  printf("          If the file header can be identified, information on\n");
  printf("          the file format and contents will be printed to the\n");
  printf("          screen in the form: <keyword><Tab><value>.\n");
  printf("Options:\n");
  printf(" -h/--help  print this text\n");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf("This program is free software under the GNU General Public License.\n");
  printf("It comes with ABSOLUTELY NO WARRANTY.\n");
  printf("For details see the file COPYING provided with this program or visit\n");
  printf("http://www.gnu.org/licenses/\n");
  printf("\n");
  return;
}
/**********************************************************************/
LOCAL void smpInfo(DOBJ *dop)
{
  int    nd, HANDLE;
  DDESC *dd=&(dop->ddl);

  HANDLE = TRUE;
  printf("file_format\t");
  fflush(stdout);
  switch(dop->fileFormat) {
  case FF_RAW:
    printf("RAW\n");
    break;
  case FF_KTH:
    printf("KTH\n");
    break;
  case FF_AIFF:
    printf("AIFF\n");
    break;
  case FF_AIFC:
    printf("AIFC\n");
    break;
  case FF_AU:
    printf("AU/SND\n");
    break;
  case FF_CSL:
    printf("CSL\n");
    break;
  case FF_CSRE:
    printf("CSRE\n");
    break;
  case FF_NIST:
    printf("NIST-SPHERE\n");
    break;
  case FF_SSFF:
    printf("SSFF\n");
    break;
  case FF_WAVE:
    printf("RIFF-WAVE\n");
    break;
  case FF_WAVE_X:
    printf("RIFF-WAVE (extended)\n");
    break;
  default:
    printf("UNKNOWN\n");
    HANDLE = FALSE;
  }
  printf("header_size\t%ld\n", dop->headerSize);
  if(dop->fileFormat == FF_SSFF)
    printf("ssff_track\t%s\n", dd->ident);
  printf("data_type\taudio\n");
  printf("data_format\t");
  fflush(stdout);
  switch(dd->format) {
    case DF_STR:
      printf("ASCII\n");
      break;
    case DF_BIT:
      printf("BIT\n");
      break;
    case DF_CHAR:
      printf("CHAR\n");
      break;
    case DF_UINT8:
      printf("INT8_unsigned\n");
      break;
    case DF_INT8:
      printf("INT8\n");
      break;
    case DF_UINT16:
      printf("INT16_unsigned\n");
      break;
    case DF_INT16:
      printf("INT16\n");
      break;
    case DF_UINT24:
      printf("INT24_unsigned\n");
      break;
    case DF_INT24:
      printf("INT24\n");
      break;
    case DF_UINT32:
      printf("INT32_unsigned\n");
      break;
    case DF_INT32:
      printf("INT32\n");
      break;
    case DF_UINT64:
      printf("INT64_unsigned\n");
      break;
    case DF_INT64:
      printf("INT64\n");
      break;
    case DF_REAL32:
      printf("FLOAT32\n");
      break;
    case DF_REAL64:
      printf("FLOAT64\n");
      break;
/*     case DF_REAL80: */
/*       printf("FLOAT80\n"); */
/*       break; */
/*     case DF_REAL128: */
/*       printf("FLOAT128\n"); */
/*       break; */
    default:
      printf("UNKNOWN\n");
      HANDLE = FALSE;
  }
  printf("data_coding\t");
  fflush(stdout);
  switch(dd->coding) {
    case DC_PCM:
      printf("PCM\n");
      break;
    case DC_BINOFF:
      printf("binary_offset\n");
      break;
    case DC_ALAW:
      printf("A-LAW\n");
      break;
    case DC_uLAW:
      printf("u-LAW\n");
      break;
    case DC_DELTA:
      printf("DELTA\n");
      HANDLE = FALSE;
      break;
    case DC_ACE2:
    case DC_ACE8:
      printf("Apple_compressed\n");
      HANDLE = FALSE;
      break;
    case DC_MAC3:
    case DC_MAC6:
      printf("Mac_compressed\n");
      HANDLE = FALSE;
      break;
    case DC_ADPCM:
    case DC_G721:
    case DC_G722:
    case DC_G723_3:
    case DC_G723_5:
    case DC_MS_ADPCM:
    case DC_CL_ADPCM:
    case DC_IDVI_ADPCM:
    case DC_OKI_ADPCM:
    case DC_IBM_ADPCM:
      printf("ADPCM\n");
      HANDLE = FALSE;
      break;
    case DC_MPEG3:
      printf("MP3\n");
      HANDLE = FALSE;
      break;
    default:
      printf("UNKNOWN\n");
      HANDLE = FALSE;
      break;
  }
  if(dd->numBits > 8) {
    printf("byte_order\t");
    fflush(stdout);
    if(MSBFIRST(dop->fileEndian))
      printf("MSB_first\n");
    else if(MSBLAST(dop->fileEndian))
      printf("MSB_last\n");
    else
      printf("UNDEFINED\n");
  }
  printf("bits/sample\t%d\n", dd->numBits);
  nd = numDecim(dop->sampFreq, 11) + 1;
  printf("sample_rate\t%.*f\n", nd, dop->sampFreq);
  printf("number_channels\t%ld\n", (long)(dd->numFields));
  if(HANDLE)
    printf("samples/channel\t%ld\n", dop->numRecords);
  if(dop->fileFormat == FF_SSFF) {
    nd = numDecim(dop->Start_Time, 12);/* sampled data need no adjust */
    if(nd <= 0) nd = 1;
    printf("start_time\t%.*f\n", nd, dop->Start_Time);
  }
  if(HANDLE)
    printf("duration\t%s\n", smp2dur(dop->numRecords, dop->sampFreq));
  return;
}
/**********************************************************************/
LOCAL void tagInfo(DOBJ *dop)
{
  int    nd;
  DDESC *dd=&(dop->ddl);

  printf("file_format\t");
  fflush(stdout);
  switch(dop->fileFormat) {
  case FF_RAW:
    printf("RAW\n");
    break;
  case FF_XASSP:
    printf("XASSP\n");
    break;
  case FF_IPDS_M:
    printf("IPdS_MIX\n");
    if(dop->sampFreq <= 0.0)                /* not specified in header */
      dop->sampFreq = MIX_SFR;
    break;
  case FF_IPDS_S:
    printf("IPdS_SAMPA\n");
    if(dop->sampFreq <= 0.0)                /* not specified in header */
      dop->sampFreq = MIX_SFR;
    break;
  case FF_XLABEL:
    printf("ESPS_XLABEL\n");
    break;
  default:
    printf("UNKNOWN\n");
    break;
  }
  printf("header_size\t%ld\n", dop->headerSize);
  printf("data_type\t");
  fflush(stdout);
  switch(dd->type) {
  case DT_TAG:
    if(dd->ident)
      printf("%s\n", dd->ident);
    else
      printf("tags\n");
    break;
  case DT_MRK:
    if(dd->ident)
      printf("%s\n", dd->ident);
    else
      printf("marks\n");
    break;
  case DT_LBL:
    if(dd->ident)
      printf("%s\n", dd->ident);
    else
      printf("labels\n");
    break;
  case DT_EPO:
    if(dd->ident)
      printf("%s\n", dd->ident);
    else
      printf("epochs\n");
    break;
  case DT_PRD:
    if(dd->ident)
      printf("%s\n", dd->ident);
    else
      printf("period_markers\n");
    break;
  default:
    if(dd->ident)
      printf("%s\n", dd->ident);
    else
      printf("UNKNOWN\n");
    break;
  }
  if(dop->fileData == FDF_ASC)
    printf("data_format\tASCII\n");
  else
    printf("data_format\tNON-ASCII\n");
  if(dop->sampFreq > 0.0) {
    nd = numDecim(dop->sampFreq, 11) + 1;
    printf("sample_rate\t%.*f\n", nd, dop->sampFreq);
  }
  else
    printf("sample_rate\tUNDEFINED\n");
  printf("number_tiers\t%ld\n", (long)(dd->numFields));
  return;
}
/**********************************************************************/
LOCAL void asspInfo(DOBJ *dop)
{
  int    nd;
  DDESC *dd=&(dop->ddl);

  printf("file_format\tXASSP\n");
  printf("header_size\t%ld\n", dop->headerSize);
  printf("data_type\t");
  fflush(stdout);
  switch(dd->type) {
  case DT_RMS:
    printf("energy");
    break;
  case DT_PIT:
    printf("pitch");
    break;
  case DT_ZCR:
    printf("zero-crossing_rate");
    break;
  case DT_EPG:
    printf("palatogram");
    break;
  case DT_AMP:
    printf("amplification");
    break;
  case DT_DUR:
    printf("lengthening");
    break;
  case DT_DATA_LOG:
    printf("datalog");
    break;
  default:
    printf("UNKNOWN");
  }
  if(strlen(dd->unit)) {
    if(strlen(dd->factor))
      printf(" (%s%s)\n", dd->factor, dd->unit);
    else
      printf(" (%s)\n", dd->unit);
  }
  else
    printf("\n");
  printf("data_format\tASCII\n");
  if(dop->dataRate > 0.0) {
    nd = numDecim(dop->dataRate, 12);
    if(nd <= 0) nd = 1;
    printf("frame_rate\t%.*f\n", nd, dop->dataRate);
  }
  if(dop->sampFreq > 0.0) {
    nd = numDecim(dop->sampFreq, 12);
    if(nd <= 0) nd = 1;
    printf("sample_rate\t%.*f\n", nd, dop->sampFreq);
  }
  if(dop->frameDur > 1) {
    printf("samples/frame\t%ld\n", dop->frameDur);
  }
  printf("start_frame\t%ld\n", dop->startRecord);
  if(dop->numRecords > 0) {
    printf("number_frames\t%ld\n", dop->numRecords);
    if(dop->sampFreq > 0.0 && dop->frameDur > 0) {
      printf("duration\t%s\n",\
	     smp2dur(dop->numRecords, dop->sampFreq / (double)dop->frameDur));
    }
  }
  return;
}
/**********************************************************************/
LOCAL void ssffInfo(DOBJ *dop)
{
  int    nd;
  DDESC *dd=&(dop->ddl);

  printf("file_format\tSSFF\n");
  printf("header_size\t%ld\n", dop->headerSize);
  printf("byte_order\t");
  if(MSBFIRST(dop->fileEndian))
    printf("MSB_first\n");
  else if(MSBLAST(dop->fileEndian))
    printf("MSB_last\n");
  else
    printf("UNDEFINED\n");
  if(dop->dataRate > 0.0) {
    nd = numDecim(dop->dataRate, 12);
    if(nd <= 0) nd = 1;
    printf("data_rate\t%.*f\n", nd, dop->dataRate);
  }
  else
    printf("data_rate\tUNDEFINED\n");
  if(dop->sampFreq > 0.0) {
    nd = numDecim(dop->sampFreq, 12);
    if(nd <= 0) nd = 1;
    printf("sample_rate\t%.*f\n", nd, dop->sampFreq);
  }
  else if(dop->frameDur < 0)
    printf("sample_rate\tUNDEFINED\n");
  if(dop->numRecords > 0)
    printf("number_frames\t%ld\n", dop->numRecords);
  nd = numDecim(dop->Start_Time, 12);
  if(nd <= 0) nd = 1;
  printf("start_time\t%.*f\n", nd, dop->Start_Time);
  if(dop->numRecords > 0) {
    if(dop->dataRate > 0.0)
      printf("duration\t%s\n", smp2dur(dop->numRecords, dop->dataRate));
    else if(dop->sampFreq > 0.0)
      printf("duration\t%s\n", smp2dur(dop->numRecords, dop->sampFreq));
  }
  while(dd != NULL) {
    printf("ssff_track\t%s\n", dd->ident);
    printf("data_type\t");
    fflush(stdout);
    switch(dd->type) {
    case DT_RMS:
      printf("RMS");
      break;
    case DT_GAIN:
      printf("gain");
      break;
    case DT_PIT:
      printf("pitch");
      break;
    case DT_AC1:
      printf("AC1_coefficient");
      break;
    case DT_LP1:
      printf("LP1_coefficient");
      break;
    case DT_PROB:
      printf("probability");
      break;
    case DT_ZCR:
      printf("zero-crossing_rate");
      break;
    case DT_ACF:
      printf("autocorrelation_function");
      break;
    case DT_LPC:
      printf("LP_coefficients");
      break;
    case DT_RFC:
      printf("reflection_coefficients");
      break;
    case DT_ARF:
      printf("area_function");
      break;
    case DT_LAR:
      printf("log_area_ratios");
      break;
    case DT_LPCEP:
      printf("cepstral_coefficients");
      break;
    case DT_FFR:
      printf("formant_frequency");
      break;
    case DT_FBW:
      printf("formant_bandwidth");
      break;
    case DT_FAM:
      printf("formant_amplitude");
      break;
    case DT_FTAMP:
      printf("amplitude_spectrum");
      break;
    case DT_FTSQR:
    case DT_FTPOW:
      printf("power_spectrum");
      break;
    case DT_FTFTS:
      printf("smoothed_spectrum");
      break;
    case DT_FTLPS:
      printf("LP_smoothed_spectrum");
      break;
    case DT_FTCSS:
      printf("cepstral_smoothed_spectrum");
      break;
    case DT_FTCEP:
      printf("cepstrum");
      break;
    case DT_EPG:
      printf("palatogram");
      break;
    default:
      printf("UNKNOWN");
    }
    if(strlen(dd->unit)) {
      if(strlen(dd->factor))
	printf(" (%s%s)\n", dd->factor, dd->unit);
      else
	printf(" (%s)\n", dd->unit);
    }
    else
      printf("\n");
    printf("data_format\t");
    fflush(stdout);
    switch(dd->format) {
    case DF_CHAR:
      printf("CHAR\n");
      break;
    case DF_UINT8:
      printf("INT8_unsigned\n");
      break;
    case DF_INT16:
      printf("INT16\n");
      break;
    case DF_INT32:
      printf("INT32\n");
      break;
    case DF_REAL32:
      printf("FLOAT\n");
      break;
    case DF_REAL64:
      printf("DOUBLE\n");
      break;
    default:
      printf("UNKNOWN\n");
    }
    printf("number_fields\t%ld\n", (long)(dd->numFields));
    dd = dd->next;
  }   /* end loop over data descriptors */
  return;
}
/**********************************************************************/
LOCAL void uwmInfo(DOBJ *dop)
{
  int    nd;
  DDESC *dd=&(dop->ddl);

  printf("file_format\tWISCONSIN\n");
  printf("header_size\t%ld\n", dop->headerSize);
  printf("data_type\tX-ray_microbeam\n");
  printf("data_format\tASCII\n");
  if(dd->coding == DC_TXY)
    printf("data_coding\tTXY\n");
  else
    printf("data_coding\tXYD\n");
  printf("number_fields\t%ld\n", (long)(dd->numFields));
  nd = numDecim(dop->sampFreq, 11) + 1;
  printf("sample_rate\t%.*f\n", nd, dop->sampFreq);
  printf("samples/frame\t%ld\n", dop->frameDur);
  if(dop->numRecords > 0) {
    printf("number_frames\t%ld\n", dop->numRecords);
    if(dop->sampFreq > 0.0 && dop->frameDur > 0) {
      printf("duration\t%s\n",\
	     smp2dur(dop->numRecords, dop->sampFreq / (double)dop->frameDur));
    }
  }
  return;
}
