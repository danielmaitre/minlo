
# hepmc has a version of fastjet that is not new enough
AM_LDFLAGS =  $(FASTJET_LDFLAGS) 
# hepmc has a version of fastjet that is not new enough
AM_CXXFLAGS =  $(FASTJET_CXXFLAGS) 



lib_LTLIBRARIES = libFlavKtPlugin.la


libFlavKtPlugin_la_SOURCES = \
	alg_declarations.f90 types.f90 alt_flavclus.f90 assertions.f90 \
	warnings_and_errors.f90 dpsorter.f90 \
	 io_utils.f90 dpktclus_gen.f90 \
	 ktflav.f90 ktflavf77.f90 lcl_dec.f90 \
	NameSelect.f90  \
	FlavKtPlugin.cc

pkginclude_HEADERS = FlavKtPlugin.hh 

libFlavKtPlugin_la_LIBADD = $(FCLIBS)

.NOTPARALLEL:

