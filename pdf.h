#ifndef _H_PDF_H
#define _H_PDF_H

#include "LHAPDF/LHAPDF.h"

#ifdef LHAPDF_MAJOR_VERSION
#if LHAPDF_MAJOR_VERSION == 6
#define LHAPDF_NEW_VERSION
#endif
#endif




#ifdef LHAPDF_NEW_VERSION
class currentPDF {
public:
	static LHAPDF::PDF* s_PDF;
	static LHAPDF::PDFSet* s_PDFSet;
	static void init(const std::string& name, int set);

};

#endif



#endif /* _H_PDF_H */
