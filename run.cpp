#include "ntuplereader/nTupleReader.h"
#include "ntuplereader/NtupleInfo.h"
#include "ntuplereader/EventReaderBase.h"
#include "MINLOfunctions.h"
#include "pdf.h"
#include "boost/program_options.hpp"
#include "MinloReader.h"
using namespace std;



namespace po = boost::program_options;

int main(int ac,char** av){

	po::options_description desc;

	desc.add_options()
    	("help", "produce help")
    	("run.pdf", po::value<std::string>()->default_value(""), "pdf")
    	("run.ntupleFile", po::value<std::string>()->default_value(""), "nTuple file to be read")
    	("config-file", po::value<std::string>()->default_value("run.cfg"), "config file to be read")
    	("minlo.njetsOrig", po::value<int>(), "number of jets to start with")
    	("minlo.njetsClus", po::value<int>(), "number of jets to cluster to")
    	("minlo.beamEnergy", po::value<double>(), "beam energy")
		("minlo.type",po::value<std::string>(),"type of contribution")
    	("keith.flg_bornonly", po::value<int>(), " Are we feeding through only Born stuff (1), or NLO (0)?")
    	("keith.imode", po::value<int>(), " imode=1 for Born, imode=2 for all NLO contribs")
    	("keith.isReal", po::value<int>(), " Set isReal=1 for real kinematics, 0 otherwise.")

    	;


	po::variables_map vm;


	po::positional_options_description p;
	p.add("config-file", -1);

	po::store(po::command_line_parser(ac, av).options(desc).positional(p).run(), vm);
	po::notify(vm);


	ifstream is;
	is.open(vm["config-file"].as<std::string>());
	po::store(po::parse_config_file(is,desc), vm);
	po::notify(vm);


	if (vm.count("help")) {
	    cout << desc << "\n";
	    return 1;
	}

	if (vm.count("run.pdf")) {
	    cout << "pdf set to " << vm["run.pdf"].as<std::string>() << ".\n";
	} else {
	    cout << "No pdf set.\n";
	}
	if (vm.count("run.ntupleFile")) {
	    cout << "ntupleFile set to " << vm["run.ntupleFile"].as<std::string>() << ".\n";
	} else {
	    cout << "No nTuple file  set.\n";
	}



	currentPDF::init(vm["run.pdf"].as<std::string>(),0);

	std::vector<std::string> fs;
	fs.push_back(vm["run.ntupleFile"].as<std::string>());

	MINLOreader r;

	r.addFiles(fs);

	int weightType=1;

	MinloInfo MI;
	// number of jets at the beginning
	MI.d_njetsOrig=vm["minlo.njetsOrig"].as<int>();
	// number of jets to cluster to
	MI.d_njetsClus=vm["minlo.njetsClus"].as<int>();
	MI.d_energy=vm["minlo.beamEnergy"].as<double>();
	if ( vm["minlo.type"].as<std::string>() == "born"){
		MI.d_type=MinloInfo::born;
	}
	if ( vm["minlo.type"].as<std::string>() == "real"){
		MI.d_type=MinloInfo::real;
	}
	if ( vm["minlo.type"].as<std::string>() == "nlo"){
		MI.d_type=MinloInfo::nlo;
	}
	if ( vm["minlo.type"].as<std::string>() == "bornLO"){
		MI.d_type=MinloInfo::bornLO;
	}

	MI.print(std::cout);

	while(r.nextEntry()){
		int id=r.d_NI.id;

		double q0,scaleForNLO;
		double alphaFactor=MINLOcomputeSudakov(MI,r.d_NI,weightType,q0,scaleForNLO);



		int flg_bornonly=vm["keith.flg_bornonly"].as<int>();    //! Are we feeding through only Born stuff (1), or NLO (0)?
		int imode=vm["keith.imode"].as<int>();           //! imode=1 for Born, imode=2 for all NLO contribs
		int isReal=vm["keith.isReal"].as<int>();          //! Set isReal=1 for real kinematics, 0 otherwise.
		double keith = MINLO_computeSudakovKeith(r.d_NI,flg_bornonly,imode,isReal);
		cout << "Keith weight: " << keith << endl;
		cout << "(" << r.d_NI.id<<") Keith weight/my weight: " << keith/alphaFactor << endl;


	}

	return 0;

}

