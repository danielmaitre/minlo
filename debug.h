/*
 * degug.h
 *
 *  Created on: 31 Jul 2010
 *      Author: daniel
 */

#ifndef DEBUG_H_
#define DEBUG_H_

#include <iostream>

//! Macro that prints its argument and its value (useful to remember what we intended to print)
#define _PRINT(X) std::cout  <<(#X) << ": "  <<  (X) <<  std::endl
//! Macro that prints a message. The difference with _PRINT() is that it will only print the message, _PRit() will print it twice.
#define _MESSAGE(X) std::cout <<  (X) << std::endl
#define _MESSAGE2(X,Y) std::cout << (X)  << (Y) << std::endl
#define _MESSAGE3(X,Y,Z) std::cout << (X) << (Y) <<  (Z) << std::endl
#define _MESSAGE4(W,X,Y,Z) std::cout << (W) << (X) << (Y) <<  (Z) << std::endl
#define _MESSAGE5(V,W,X,Y,Z) std::cout << (V) << (W) << (X) << (Y) <<  (Z) << std::endl
//! Macro that a warning (with cerr instead of cout.)
#define _WARNING(X) std::cerr << (X) << std::endl
#define _WARNING2(X,Y) std::cerr << (X)  << (Y) << std::endl
#define _WARNING3(X,Y,Z) std::cerr << (X) << (Y) <<  (Z) << std::endl
#define _WARNING4(X,Y,Z,A) std::cerr << (X) << (Y) <<  (Z) << (A) << std::endl
#define _WARNING5(X,Y,Z,A,B) std::cerr << (X) << (Y) <<  (Z) << (A) << (B) << std::endl
#define _MESSAGE6(V,W,X,Y,Z,Z1) std::cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << std::endl
#define _MESSAGE7(V,W,X,Y,Z,Z1,Z2) std::cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << std::endl
#define _MESSAGE8(V,W,X,Y,Z,Z1,Z2,Z3) std::cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << std::endl
#define _MESSAGE9(V,W,X,Y,Z,Z1,Z2,Z3,Z4) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << std::endl
#define _MESSAGE10(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << std::endl
#define _MESSAGE11(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << std::endl
#define _MESSAGE12(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << std::endl
#define _MESSAGE13(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << (Z8) << std::endl
#define _MESSAGE14(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << (Z8) << (Z9) << std::endl
#define _MESSAGE15(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,U1) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << (Z8) << (Z9) << (Z10) << std::endl
#define _MESSAGE16(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,U1,U2) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << (Z8) << (Z9) << (U1) << (U2) << std::endl
#define _MESSAGE17(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,U1,U2,U3) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << (Z8) << (Z9) << (U1) << (U2) << (U3) << std::endl



#ifdef DEBUG_ON
#define UNIQUE_ID(X)  unique_ID_ ## X
#define UNIQUE_STR_ID(X)  unique_ID_str ## X

#define DEBUG_MESSAGEIMPL(X,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE2(UNIQUE_STR_ID(UNIQUE),X);}
#define DEBUG_MESSAGE(X) DEBUG_MESSAGEIMPL(X,__LINE__)
#define DEBUG_MESSAGEIMPL2(X,Y,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE3(UNIQUE_STR_ID(UNIQUE),X,Y);}
#define DEBUG_MESSAGE2(X,Y) DEBUG_MESSAGEIMPL2(X,Y,__LINE__)
#define DEBUG_MESSAGEIMPL3(X,Y,Z,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE4(UNIQUE_STR_ID(UNIQUE),X,Y,Z);}
#define DEBUG_MESSAGE3(X,Y,Z) DEBUG_MESSAGEIMPL3(X,Y,Z,__LINE__)
#define DEBUG_MESSAGEIMPL4(X,Y,Z,A,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE5(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A);}
#define DEBUG_MESSAGE4(X,Y,Z,A) DEBUG_MESSAGEIMPL4(X,Y,Z,A,__LINE__)
#define DEBUG_MESSAGEIMPL5(X,Y,Z,A,B,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE6(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B);}
#define DEBUG_MESSAGE5(X,Y,Z,A,B) DEBUG_MESSAGEIMPL5(X,Y,Z,A,B,__LINE__)
#define DEBUG_MESSAGEIMPL6(X,Y,Z,A,B,C,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE7(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B,C);}
#define DEBUG_MESSAGE6(X,Y,Z,A,B,C) DEBUG_MESSAGEIMPL6(X,Y,Z,A,B,C,__LINE__)
#define DEBUG_MESSAGEIMPL7(X,Y,Z,A,B,C,D,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE8(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B,C,D);}
#define DEBUG_MESSAGE7(X,Y,Z,A,B,C,D) DEBUG_MESSAGEIMPL7(X,Y,Z,A,B,C,D,__LINE__)
#define DEBUG_MESSAGEIMPL8(X,Y,Z,A,B,C,D,E,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE9(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B,C,D,E);}
#define DEBUG_MESSAGE8(X,Y,Z,A,B,C,D,E) DEBUG_MESSAGEIMPL8(X,Y,Z,A,B,C,D,E,__LINE__)
#define DEBUG_MESSAGEIMPL9(X,Y,Z,A,B,C,D,E,F,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE10(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B,C,D,E,F);}
#define DEBUG_MESSAGE9(X,Y,Z,A,B,C,D,E,F) DEBUG_MESSAGEIMPL9(X,Y,Z,A,B,C,D,E,F,__LINE__)
#define DEBUG_MESSAGEIMPL10(X,Y,Z,A,B,C,D,E,F,G,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE11(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B,C,D,E,F,G);}
#define DEBUG_MESSAGE10(X,Y,Z,A,B,C,D,E,F,G) DEBUG_MESSAGEIMPL10(X,Y,Z,A,B,C,D,E,F,G,__LINE__)

#define DEBUG_PRINTIMPL(X,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){ std::cout << UNIQUE_STR_ID(UNIQUE); _PRINT(X);}
#define DEBUG_PRINT(X) DEBUG_PRINTIMPL(X,__LINE__)
#define DEBUG_IMPL(X,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){ if (UNIQUE_STR_ID(UNIQUE).size()!=0){ _MESSAGE2("::begin of ",UNIQUE_STR_ID(UNIQUE)); }; X; if (UNIQUE_STR_ID(UNIQUE).size()!=0){ _MESSAGE2("::end of ",UNIQUE_STR_ID(UNIQUE)); }; }
#define DEBUG(X) DEBUG_IMPL(X,__LINE__)

#define NAMED_DEBUG_IMPL(NAME,X,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,NAME) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,NAME,UNIQUE); if (UNIQUE_ID(UNIQUE)){ if (UNIQUE_STR_ID(UNIQUE).size()!=0){ _MESSAGE2("::begin of ",UNIQUE_STR_ID(UNIQUE)); }; X; if (UNIQUE_STR_ID(UNIQUE).size()!=0){ _MESSAGE2("::end of ",UNIQUE_STR_ID(UNIQUE)); }; }
#define NAMED_DEBUG(NAME,X) NAMED_DEBUG_IMPL(NAME,X,__LINE__)


#else
#define DEBUG_MESSAGE(X)
#define DEBUG_MESSAGE2(X,Y)
#define DEBUG_MESSAGE3(X,Y,Z)
#define DEBUG_MESSAGE4(X,Y,Z,A)
#define DEBUG_MESSAGE5(X,Y,Z,A,B)
#define DEBUG_MESSAGE6(X,Y,Z,A,B,C)
#define DEBUG_MESSAGE7(X,Y,Z,A,B,C,D)
#define DEBUG_MESSAGE8(X,Y,Z,A,B,C,D,E)
#define DEBUG_MESSAGE9(X,Y,Z,A,B,C,D,E,F)
#define DEBUG_MESSAGE10(X,Y,Z,A,B,C,D,E,F,G)

#define DEBUG_PRINT(X)
#define DEBUG(X)
#define NAMED_DEBUG(NAME,X)
#endif

#include <string>

namespace Debug {

bool need_debug(const char* FileName,const char* FuncName);
std::string get_info_str(const char* FileName,const char* FuncName,int line);

} /* Debug */

using namespace Debug;

#endif /* DEBUG_H_ */
