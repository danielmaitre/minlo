AC_INIT([MINLOlib], [0.0], [daniel.maitre@durham.ac.uk])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])
AM_DISABLE_STATIC
AC_PROG_CXX
AC_PROG_F77
AC_PROG_FC
AC_FC_LIBRARY_LDFLAGS
AM_PATH_PYTHON
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
AC_PROG_LIBTOOL
AC_CONFIG_FILES([
    Makefile
    MyFlavKtPlugin/Makefile
    FlavKtPlugin/Makefile
    MINLOlib-config
])

AC_ARG_WITH([fastjet-path], [AS_HELP_STRING([--with-fastjet-path],
          [location of the fastjet installation])],
          [],
          [with_fastjet_path=yes])

AS_IF([test "x$with_fastjet_path" != xyes],
          [AC_SUBST([FASTJET_CXXFLAGS], ["`$with_fastjet_path/bin/fastjet-config --plugins=yes --cxxflags`"])
	   AC_SUBST([FASTJET_LDFLAGS], ["`$with_fastjet_path/bin/fastjet-config --plugins=yes --libs`"])
	  ],
          [AC_SUBST([FASTJET_CXXFLAGS], [`fastjet-config  --plugins=yes --cxxflags`])
	   AC_SUBST([FASTJET_LDFLAGS], ["`fastjet-config --plugins=yes --libs`"])
])

oldCPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $FASTJET_CXXFLAGS"
AC_LANG([C++])
AC_CHECK_HEADERS([fastjet/PseudoJet.hh],
     [],
     [AC_ERROR(
	[Need fastjet headers, please specifiy the path to the fastjet installation with the --with-fastjet-path option.]
	)
	])
CPPFLAGS=$oldCPPFLAGS



AC_ARG_WITH([gsl-path], [AS_HELP_STRING([--with-gsl-path],
          [location of the gsl installation])],
          [],
          [with_gsl_path=yes])

AS_IF([test "x$with_gsl_path" != xyes],
          [AC_SUBST([GSL_CXXFLAGS], ["`$with_gsl_path/bin/gsl-config --cflags`"])
			AC_SUBST([GSL_LDFLAGS], ["`$with_gsl_path/bin/gsl-config --libs`"])
	  ],
        [AC_SUBST([GSL_CXXFLAGS], [`gsl-config --cflags`])
	   AC_SUBST([GSL_LDFLAGS], ["`gsl-config --libs`"])
])



AC_ARG_WITH([nTupleReader-path], [AS_HELP_STRING([--with-nTupleReader-path],
          [location of the nTupleReader installation])],
          [use_ntuplereader_path=yes])

AS_IF([test "x$use_ntuplereader_path" == xyes],
          [AC_SUBST([NTUPLEREADER_CXXFLAGS], ["`$with_nTupleReader_path/bin/nTupleReader-config --cxxflags`"])
	   AC_SUBST([NTUPLEREADER_LDFLAGS], ["`$with_nTupleReader_path/bin/nTupleReader-config --libs`"])
	   AC_SUBST([NTUPLEREADER_LIBS], ["`$with_nTupleReader_path/bin/nTupleReader-config --libs`"])
	  ],
          [AC_SUBST([NTUPLEREADER_CXXFLAGS], [`nTupleReader-config --cxxflags`])
	   AC_SUBST([NTUPLEREADER_LDFLAGS], ["`nTupleReader-config --libs`"])
	   AC_SUBST([NTUPLEREADER_LIBS], ["`nTupleReader-config --libs`"])
])




oldCPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $NTUPLEREADER_CXXFLAGS"
echo $CPPFLAGS
AC_LANG([C++])
AC_CHECK_HEADERS([ntuplereader/nTupleReader.h],
     [],
     [AC_ERROR(
	[Need nTupleReader headers, please specifiy the path to the nTupleReader installation with the --with-nTupleReader-path option.]
	)
	])
CPPFLAGS=$oldCPPFLAGS

AC_ARG_WITH([lhapdf-path], [AS_HELP_STRING([--with-lhapdf-path],
          [location of the lhapdf installation])],
          [],
          [with_lhapdf_path=yes])

AS_IF([test "x$with_lhapdf_path" != xyes],
          [ 
			AC_SUBST([LHAPDF_CXXFLAGS], ["`$with_lhapdf_path/bin/lhapdf-config --cppflags`"])
			AC_SUBST([LHAPDF_LDFLAGS], ["`$with_lhapdf_path/bin/lhapdf-config --ldflags`"])
	  ],
        [
		AC_SUBST([LHAPDF_CXXFLAGS], [`lhapdf-config --cppflags`])
	   AC_SUBST([LHAPDF_LDFLAGS], ["`lhapdf-config --ldflags`"])
])


oldCPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $LHAPDF_CXXFLAGS"
echo $CPPFLAGS
AC_LANG([C++])
AC_CHECK_HEADERS([LHAPDF/LHAPDF.h],
     [],
     [AC_ERROR(
	[Need LHAPDF headers, please specifiy the path to the LHAPDF installation with the --with-lhapdf-path option.]
	)
	])
CPPFLAGS=$oldCPPFLAGS



AC_ARG_WITH([root-path], [AS_HELP_STRING([--with-root-path],
          [location of the root installation])],
          [],
          [with_root_path=yes])

AS_IF([test "x$with_root_path" != xyes],
          [AC_SUBST([ROOT_CXXFLAGS], ["`$with_root_path/bin/root-config --cflags`"])
	   AC_SUBST([ROOT_LDFLAGS], ["`$with_root_path/bin/root-config --libs`"])
	   AC_SUBST([ROOTCINT], ["`$with_root_path/bin/root-config --prefix`/bin/rootcint"])
	   ROOTINCDIR=`$with_root_path/bin/root-config --incdir`
	  ],
          [AC_SUBST([ROOT_CXXFLAGS], [`root-config --cflags`])
	   AC_SUBST([ROOT_LDFLAGS], ["`root-config --ldflags`"])
	   AC_SUBST([ROOT_LIBS], ["`root-config --libs`"])
	   AC_SUBST([ROOTCINT], ["`root-config --prefix`/bin/rootcint"])
	   ROOTINCDIR=`root-config --incdir`
])


oldCPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS -I$ROOTINCDIR"
AC_LANG([C++])
AC_CHECK_HEADERS([TBuffer.h TChain.h TFile.h TH1F.h TKey.h TList.h TObject.h TString.h TSystem.h TTree.h
],[],[AC_ERROR([Need ROOT headers, please specifiy the path to the ROOT installation with the --with-root-path option. This can also fail if you are using ROOT >= 6.0 and do not the c++11 standards are not supported by default, consider using CXXFLAGS=-std=c++11])])
CPPFLAGS=$oldCPPFLAGS


AC_ARG_ENABLE([debug],
     [  --enable-debug    Turn on debugging],
     [case "${enableval}" in
       yes) debug=true ;;
       no)  debug=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
     esac],[debug=false])
AM_CONDITIONAL([DEBUG], [test x$debug = xtrue])

if test x$debug = xtrue ; then
    AM_CPPFLAGS="$AM_CPPFLAGS -DDEBUG_ON "
    echo "Using debugging"
else
    AM_CPPFLAGS="$AM_CPPFLAGS  "
    echo "Not using debugging"
fi

AC_SUBST([AM_CPPFLAGS])



AC_OUTPUT
