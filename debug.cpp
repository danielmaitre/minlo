/*
 * debug.cpp
 *
 *  Created on: 31 Jul 2010
 *      Author: daniel
 */

#include "debug.h"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>


using namespace std;

namespace Debug {

bool DoPrintFileName=false;

vector<string> ReadDebugFile( ){

	vector<string> all;
	string path="debug.dat";

	std::ifstream file;
	file.open(path.c_str());
	if (file.fail()){
		_WARNING3("Could not open ",path,": done nothing. ");
	} else {
		_MESSAGE3("Reading debug output information from ",path,". ");
		while ( file ){
			char next_line[100];
			file.getline(next_line,100);
//			_MESSAGE3("next line '",next_line,"'");
			stringstream ss(next_line);
			if (ss && next_line[0]!='#' && next_line[0]!=':'){
				string new_path,new_func;
				ss >> new_path;
//				_MESSAGE3("new path '",new_path,"'");
				if (ss){
					ss >> new_func;
				}
				if ( new_func.size() == 0){
					new_func="*";
				}
				if (new_path.size() != 0){   // prevents empty lines to be taken as path
					all.push_back(new_path+string("|")+new_func);
				}
			} else if (ss && next_line[0]==':'){
				string new_opt,value;
				ss >> new_opt;
//				_MESSAGE3("new path '",new_path,"'");
				if (ss){
					ss >> value;
				}
				if ( new_opt == ":PrintFileName" ){
					DoPrintFileName=true;
					_MESSAGE("Filenames will be printed. To switch this option off, remove \":PrintFileName\" from the BH_debug.dat file.");
				}
			}

		}

	}
	return all;
};

// returns the filename from a path (e.g. XXX from /../../../XXX )
string GetFileName(const char* path){
	string filename(path);
	int pos=filename.rfind('/');
	if (pos != string::npos){
		return filename.substr(pos+1,filename.size()-pos-1);
	} else {
		return filename;
	}
}

bool need_debug(const char* FileName,const char* FuncName){

	static vector<string> DebugFileList=ReadDebugFile();

	string filename=GetFileName(FileName);
	string all=filename+string("|*");
	vector<string>::iterator all_pos=find(DebugFileList.begin(),DebugFileList.end(),all);
	if ( all_pos!=DebugFileList.end()){
		return true;
	} else {
		string specific=filename+string("|")+string(FuncName);
		vector<string>::iterator specific_pos=find(DebugFileList.begin(),DebugFileList.end(),specific);
		if ( specific_pos!=DebugFileList.end()){
			return true;
		} else {
//			_MESSAGE2("No need to output message for ",all);
			return false;
		}
	}
};

std::string get_info_str(const char* FileName,const char* FuncName,int line){
	if (DoPrintFileName){
		string str;
		stringstream ss(str);
		ss << GetFileName(FileName) << "|" << FuncName << "(" << line << "): ";
		return ss.str();
	} else return string("");
};


} /* Debug */
