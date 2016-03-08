SUBDIRS = MyFlavKtPlugin FlavKtPlugin



bin_PROGRAMS= test
test_SOURCES = test.cpp

lib_LTLIBRARIES=libMINLOlib.la 

libMINLOlib_la_SOURCES = \
	MINLOfunctions.cpp MINLOfunctions.h \
	MinloKeithDeclarations.h \
	InterpolatedFunction.cpp InterpolatedFunction.h\
	indie-minlo.f NLLsudakov.f sudakov.f \
	debug.h debug.cpp \
	pdf.h pdf.cpp

AM_LDFLAGS =  $(FASTJET_LDFLAGS) $(LHAPDF_LDFLAGS) 

libMINLOlib_la_LDFLAGS = $(NTUPLEREADER_LDFLAGS)
libMINLOlib_la_LIBADD = MyFlavKtPlugin/libMyFlavKtPlugin.la FlavKtPlugin/libFlavKtPlugin.la -lgfortran $(GSL_LDFLAGS)


AM_CXXFLAGS =  $(NTUPLEREADER_CXXFLAGS) $(FASTJET_CXXFLAGS) $(LHAPDF_CXXFLAGS)

test_LDFLAGS = $(NTUPLEREADER_LDFLAGS) $(FASTJET_CXXFLAGS) $(LHAPDF_CXXFLAGS)
test_LDADD = libMINLOlib.la MyFlavKtPlugin/libMyFlavKtPlugin.la  $(FASTJET_LDFLAGS)


bin_SCRIPTS=MINLOlib-config

pkginclude_HEADERS = MINLOfunctions.h

