
# set the include path found by configure
INCLUDES= $(all_includes)
pkginclude_HEADERS = acf.h assp.h asspana.h asspdsp.h asspendian.h asspfio.h \
				asspmess.h assptime.h aucheck.h auconv.h dataobj.h dlldef.h headers.h ieee.h \
				ipds_lbl.h labelobj.h mhs.h misc.h miscdefs.h mylimits.h rms.h trace.h uwm_xrmb.h \
				zcr.h rfc.h fmt.h filter.h esps_lbl.h diff.h ksv.h spectra.h

# the library search path.
lib_LTLIBRARIES = libmisc.la libassp.la 
libmisc_la_LDFLAGS = -no-undefined
libmisc_la_CFLAGS =  -DBUILDING_MISC
libmisc_la_SOURCES = fgetl.c myrint.c numdecim.c parsepath.c query.c trace.c \
			miscstring.c bitarray.c isgerman.c myrand.c statistics.c chain.c
libassp_la_LDFLAGS = -no-undefined
libassp_la_CFLAGS =  -DBUILDING_ASSP
libassp_la_LIBADD = $(top_builddir)/src/libmisc.la -lm
libassp_la_SOURCES = asspana.c asspfio.c asspmess.c aucheck.c auconv.c \
	dataobj.c dsputils.c filters.c headers.c ieee.c labelobj.c memswab.c rms.c \
	smp2dur.c winfuncs.c zcr.c acf.c fft.c freqconv.c mhs.c math.c rfc.c lpc.c fmt.c \
	filter.c diff.c ksv.c spectra.c
bin_PROGRAMS = acfana afconvert affilter ffinfo forest rfcana rmsana spectrum \
	zcrana f0_mhs msb afdiff f0_ksv
acfana_SOURCES = acfana.c
affilter_SOURCES = affilter.c
forest_SOURCES = forest.c
rfcana_SOURCES = rfcana.c
rmsana_SOURCES = rmsana.c
spectrum_SOURCES = spectrum.c
zcrana_SOURCES = zcrana.c
acfana_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la \
	-lm
affilter_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la \
	-lm
rfcana_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la \
	-lm
rmsana_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la \
	-lm
spectrum_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la \
	-lm
zcrana_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la \
	-lm
forest_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la \
	-lm


libmisc_la_LIBADD = -lm
ffinfo_SOURCES = ffinfo.c
ffinfo_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la \
	-lm
afconvert_SOURCES = afconvert.c
afconvert_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la\
	-lm
f0_mhs_SOURCES = f0_mhs.c
f0_mhs_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la \
	-lm
msb_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la -lm
msb_SOURCES = msb.c
afdiff_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la -lm
afdiff_SOURCES = afdiff.c


f0_ksv_SOURCES = f0_ksv.c
f0_ksv_LDADD = $(top_builddir)/src/libassp.la $(top_builddir)/src/libmisc.la -lm
