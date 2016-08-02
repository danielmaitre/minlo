
# hepmc has a version of fastjet that is not new enough
AM_LDFLAGS =  $(FASTJET_LDFLAGS) $(ROOT_LDFLAGS)
# hepmc has a version of fastjet that is not new enough
AM_CXXFLAGS =  $(NTUPLEREADER_CXXFLAGS) $(FASTJET_CXXFLAGS)


lib_LTLIBRARIES = libMyCluster.la


libMyCluster_la_SOURCES = \
		MyFlavKtPlugin.cc boosts.h boosts.cpp



libMyCluster_la_LIBS = fastjettool
libMyCluster_la_LIBADD = $(FCLIBS) 
