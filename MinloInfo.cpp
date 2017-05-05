/*
 * MinloInfo.cpp
 *
 *  Created on: 13 Apr 2017
 *      Author: daniel
 */




#include "MinloInfo.h"
#include "boost/program_options.hpp"
#include <iostream>

namespace po = boost::program_options;

void MinloInfo::readFromStream(std::istream& is){

	po::options_description desc;

	desc.add_options()
	("minlo.njetsOrig", po::value<int>(), "number of jets to start with")
	("minlo.njetsClus", po::value<int>(), "number of jets to cluster to")
	("minlo.beamEnergy", po::value<double>(), "beam energy")
	("minlo.type",po::value<std::string>(),"type of contribution")
	("minlo.radius",po::value<double>()->default_value(1.0),"radius for the kt clustering")
	("minlo.alltheway", po::value<int>(), " non-zero for the raising-all-the-way policy")
	("minlo.coreScale",po::value<std::string>()->default_value("shat"),"core scale choice")
	("minlo.useModifiedR",po::value<bool>()->default_value(false),"whether to use the modified R definition")
	("minlo.scaleMode",po::value<std::string>()->default_value("geometric"),"mode for the nlo scale, either geometric or inverseAlpha")
	("minlo.stopAfterFirstDrop",po::value<bool>()->default_value(false),"whether to stop clustering when a nodal scale becomes smaller than its predecessor.")
	("minlo.usePDFalphas",po::value<bool>()->default_value(false),"whether to use alphas from the pdf rather than the build in one.")
	("minlo.useSherpa",po::value<bool>()->default_value(false),"whether to use sherpa sudakovs.")
	("minlo.useLOalphas",po::value<bool>()->default_value(false),"whether to use LO alphas running.")
	("minlo.lambda",po::value<double>()->default_value(0.226),"lambda_QCD for five light flavours")
	("minlo.sherpaMode", po::value<int>()->default_value(3), "mode for sherpa sudakov")
	("minlo.useRapidityInClustering", po::value<bool>()->default_value(true), "use rapidity in clustering")
	;


	po::variables_map vm;

	po::store(po::parse_config_file(is,desc,true), vm);
	po::notify(vm);

	// number of jets at the beginning
	d_njetsOrig=vm["minlo.njetsOrig"].as<int>();
	// number of jets to cluster to
	d_njetsClus=vm["minlo.njetsClus"].as<int>();
	d_energy=vm["minlo.beamEnergy"].as<double>();
	if ( vm["minlo.type"].as<std::string>() == "born"){
		d_type=MinloInfo::born;
	}
	if ( vm["minlo.type"].as<std::string>() == "real"){
		d_type=MinloInfo::real;
	}
	if ( vm["minlo.type"].as<std::string>() == "nlo"){
		d_type=MinloInfo::nlo;
	}
	if ( vm["minlo.type"].as<std::string>() == "bornLO"){
		d_type=MinloInfo::bornLO;
	}
	d_R=vm["minlo.radius"].as<double>();
	d_alltheway=vm["minlo.alltheway"].as<int>();
	d_useModifiedR=vm["minlo.useModifiedR"].as<bool>();
	d_stopAfterFirstDrop=vm["minlo.stopAfterFirstDrop"].as<bool>();
	d_usePDFalphas=vm["minlo.usePDFalphas"].as<bool>();
	d_useLOalphas=vm["minlo.useLOalphas"].as<bool>();
	d_lambda=vm["minlo.lambda"].as<double>();
	d_useSherpa=vm["minlo.useSherpa"].as<bool>();
	d_sherpaMode=vm["minlo.sherpaMode"].as<int>();
	d_useRapidityInClustering=vm["minlo.useRapidityInClustering"].as<bool>();
	if ( vm["minlo.scaleMode"].as<std::string>() == "geometric"){
		d_scaleMode=MinloInfo::geometric;
	} else if ( vm["minlo.scaleMode"].as<std::string>() == "inverseAlpha"){
		d_scaleMode=MinloInfo::inverseAlpha;
	} else {
		std::cout << "wrong setting for scaleMode!" << std::endl;
	}
	if ( vm["minlo.coreScale"].as<std::string>() == "shat"){
		d_coreScaleType=MinloInfo::shat;
	} else if ( vm["minlo.coreScale"].as<std::string>() == "hthalf"){
		d_coreScaleType=MinloInfo::hthalf;
	} else 	if ( vm["minlo.coreScale"].as<std::string>() == "stefan"){
		d_coreScaleType=MinloInfo::stefan;
	}


		else{
		std::cout << "wrong setting for scaleMode!" << std::endl;
	}


}
