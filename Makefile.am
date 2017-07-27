SUBDIRS = MyFlavKtPlugin FlavKtPlugin

F77=gfortran

ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS= run ExpIntegralTest 
run_SOURCES = run.cpp
ExpIntegralTest_SOURCES = ExpIntegralTest.cpp

lib_LTLIBRARIES=libMINLOlib.la 

libMINLOlib_la_SOURCES = \
	MINLOfunctions.cpp MINLOfunctions.h \
	MinloKeithDeclarations.h \
	InterpolatedFunction.cpp InterpolatedFunction.h\
	indie-minlo.f NLLsudakov.f sudakov.f \
	debug.h debug.cpp \
	pdf.h pdf.cpp \
	utils.h utils.cpp \
	MinloReader.h MinloReader.cpp \
	sudakovs.cpp sudakovs.h \
	MinloInfo.h MinloInfo.cpp \
	nll.h nll.cpp \
	coreScales.cpp coreScales.h\
	ExpIntegral.cpp ExpIntegral.h

AM_LDFLAGS =  $(FASTJET_LDFLAGS) $(LHAPDF_LDFLAGS) 
AM_LDFLAGS += $(ROOT_LIBS)

libMINLOlib_la_LDFLAGS = $(NTUPLEREADER_LDFLAGS)
#libMINLOlib_la_LIBADD = MyFlavKtPlugin/libMyFlavKtPlugin.la FlavKtPlugin/libFlavKtPlugin.la -lgfortran $(GSL_LDFLAGS)
libMINLOlib_la_LIBADD = MyFlavKtPlugin/libMyCluster.la FlavKtPlugin/libFlavKtPlugin.la $(GSL_LDFLAGS)


AM_CXXFLAGS =  $(NTUPLEREADER_CXXFLAGS) $(FASTJET_CXXFLAGS) $(LHAPDF_CXXFLAGS)
AM_CXXFLAGS += $(ROOT_CXXFLAGS)
AM_CXXFLAGS += $(BOOST_CXXFLAGS)


run_LDFLAGS = $(NTUPLEREADER_LDFLAGS) $(FASTJET_CXXFLAGS) $(LHAPDF_CXXFLAGS) $(BOOST_LDFLAGS)
run_LDFLAGS += $(ROOT_LIBS)
run_LDADD = libMINLOlib.la MyFlavKtPlugin/libMyCluster.la  $(FASTJET_LDFLAGS)

ExpIntegralTest_LDFLAGS = $(NTUPLEREADER_LDFLAGS) $(FASTJET_CXXFLAGS) $(LHAPDF_CXXFLAGS) $(BOOST_LDFLAGS)
ExpIntegralTest_LDFLAGS += $(ROOT_LIBS)
ExpIntegralTest_LDADD = libMINLOlib.la 


AM_CPPFLAGS += $(PYTHON_INCLUDE) -I$(PYTHON_INCLUDE_DIR) $(NTUPLEREADER_CPPFLAGS) $(ROOT_CPPFLAGS)

lib_LTLIBRARIES+=libPyMinlo.la

libPyMinlo_la_SOURCES =  MinloInterface.cpp
libPyMinlo_la_LDFLAGS = $(NTUPLEREADER_LDFLAGS) $(FASTJET_CXXFLAGS) $(LHAPDF_CXXFLAGS) $(BOOST_LDFLAGS)
libPyMinlo_la_LDFLAGS += $(ROOT_LIBS)
libPyMinlo_la_LIBADD = libMINLOlib.la MyFlavKtPlugin/libMyCluster.la  $(FASTJET_LDFLAGS)

MinloInterface.d:	Minlo.i
	swig $(AM_CPPFLAGS) -M -Wall -python -c++ -o MinloInterface.cpp $(srcdir)/Minlo.i

@am__include@ ./MinloInterface.d


MinloInterface.cpp:	Minlo.i
	swig $(AM_CPPFLAGS) -Wall -python -c++ -o MinloInterface.cpp $(srcdir)/Minlo.i

CLEANFILES = MinloInterface.cpp




LIBS += $(BOOST_PROGRAM_OPTIONS_LIB)




bin_SCRIPTS=MINLOlib-config

pkginclude_HEADERS = MINLOfunctions.h MinloInfo.h pdf.h nll.h

uninstall-hook:
	rm -f  $(DESTDIR)$(libdir)/_MinloReader.so
	if test -f  $(DESTDIR)$(libdir)/MinloReader.py; then rm -f $(DESTDIR)$(libdir)/MinloReader.py; fi



install-exec-hook:
	ls -l ${prefix}/lib
	ln -fs $(DESTDIR)$(libdir)/libPyMinlo.so $(DESTDIR)$(libdir)/_MINLOReader.so
	if test -f $(top_builddir)/MINLOReader.py; then \
	cp $(top_builddir)/MINLOReader.py $(DESTDIR)$(libdir); \
	else \
	cp $(srcdir)/MinloReader.py $(DESTDIR)$(libdir)	;	\
	fi


