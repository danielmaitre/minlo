#ifndef __KTFLAVF77_H__
#define __KTFLAVF77_H__ 1

#include "JetIModes.h"

#define ktflavf77  ktflavf77_
extern "C" void ktflavf77(
                   const int & imode,
		   const double & rparam,
                   const int & np,
                   const double * pp,    // p(4,np)
                   const int & nflav,  
                   const int * flav,     // flav(nflav,np)
                   const int * beamflav, // beamflav(nflav,-2:-1)
                   int    & njet,
                   double * jetpp,       // jetpp(4,np) (filled up to njet)
                   int    * jetflav,     // jetflav(nflav,np)
                   int    * incomingflav,// incomingflav(nflav,-2:-1)
                   double & d3           // d3 parameter (not normalised)
                   );

#define filterflavour filterflavour_
extern "C" void filterflavour(const int & np, 
			      const int & nflav, 
                              const int * flav,  // flav(nflav,np) 
                              const int & id  
                             );

#define ktf77MakeNodeList ktf77makenodelist_
extern "C" void ktf77MakeNodeList(
                   const int    & imode,   // algorithm
                   const double & rparam,  // R
                   const int    & np,      // number of particles
                   const double * pp       // p(4,np)
                  );

#define ktf77MakeNodeListFlav ktf77makenodelistflav_
extern "C" void ktf77MakeNodeListFlav(
                   const int    & imode,    // algorithm
                   const double & rparam,   // R
                   const int    & np,       // number of particles
                   const double * pp,       // p(4,np)
		   const int    & nflav,    // number of flavours we consider
		   const int    * flav,     // flav(nflav,np)
		   const int    * beamflav  // beamflav(nflav,2)
                  );

#define ktf77DeleteNodeList ktf77deletenodelist_
extern "C" void ktf77DeleteNodeList();

#define ktf77ReadNodeList  ktf77readnodelist_
extern "C" void ktf77ReadNodeList (
                        const int  & i, // the index of the node list to access 
                        double     * p, // p(4) the momentum at the nodelist
                        int  & parent1, // index of the first parent
                        int  & parent2, // index of the second parent
                        double & y      // dij???
                        );

#endif //  __KTFLAVF77_H__ 
