SUBDIRS = MyFlavKtPlugin FlavKtPlugin

F77=gfortran

ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS= test run
test_SOURCES = test.cpp
run_SOURCES = run.cpp


lib_LTLIBRARIES=libMINLOlib.la 

libMINLOlib_la_SOURCES = \
	MINLOfunctions.cpp MINLOfunctions.h \
	MinloKeithDeclarations.h \
	InterpolatedFunction.cpp InterpolatedFunction.h\
	indie-minlo.f NLLsudakov.f sudakov.f \
	debug.h debug.cpp \
	pdf.h pdf.cpp

AM_LDFLAGS =  $(FASTJET_LDFLAGS) $(LHAPDF_LDFLAGS) 
AM_LDFLAGS += $(ROOT_LIBS)

libMINLOlib_la_LDFLAGS = $(NTUPLEREADER_LDFLAGS)
#libMINLOlib_la_LIBADD = MyFlavKtPlugin/libMyFlavKtPlugin.la FlavKtPlugin/libFlavKtPlugin.la -lgfortran $(GSL_LDFLAGS)
libMINLOlib_la_LIBADD = MyFlavKtPlugin/libMyCluster.la FlavKtPlugin/libFlavKtPlugin.la $(GSL_LDFLAGS)


AM_CXXFLAGS =  $(NTUPLEREADER_CXXFLAGS) $(FASTJET_CXXFLAGS) $(LHAPDF_CXXFLAGS)
AM_CXXFLAGS += $(ROOT_CXXFLAGS)
AM_CXXFLAGS += $(BOOST_CXXFLAGS)

test_LDFLAGS = $(NTUPLEREADER_LDFLAGS) $(FASTJET_CXXFLAGS) $(LHAPDF_CXXFLAGS)
test_LDFLAGS += $(ROOT_LIBS)
test_LDADD = libMINLOlib.la MyFlavKtPlugin/libMyCluster.la  $(FASTJET_LDFLAGS)

run_LDFLAGS = $(NTUPLEREADER_LDFLAGS) $(FASTJET_CXXFLAGS) $(LHAPDF_CXXFLAGS) $(BOOST_LDFLAGS)
run_LDFLAGS += $(ROOT_LIBS)
run_LDADD = libMINLOlib.la MyFlavKtPlugin/libMyCluster.la  $(FASTJET_LDFLAGS)

LIBS += $(BOOST_PROGRAM_OPTIONS_LIB)

bin_SCRIPTS=MINLOlib-config

pkginclude_HEADERS = MINLOfunctions.h



