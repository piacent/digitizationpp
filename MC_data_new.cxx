//compiling it seems I need -lstdc++fs
//CALL: ./progname.exe Configfile -I inputdir -O outputdir
#include <iostream>
#include <fstream>
#include <chrono>
#include <map>
#include <random>
#include <string>
#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <cstdlib>
#include <experimental/filesystem>
#include <memory>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

using namespace std;


////Global and constant value;
double GEM1_gain;
double GEM2_gain;
double GEM3_gain;
double extraction_eff_GEM1;
double extraction_eff_GEM2;
double extraction_eff_GEM3;
int x_ini = 0;
int y_ini = 0;
int z_ini = 0;

void ReadConfig(string name, map<string,string>& options);
void SaveValues(map<string,string>& options, shared_ptr<TFile>& outfile);
int NelGEM2(const vector<double>& energyDep,const vector<double>& z_hit, map<string,string>& options);

int main(int argc, char** argv)
{
	string outfolder;
	string infolder;
	if(argc<2) {cerr<<"No Configfile given!!\nSuggested use: ./progname.exe Configfile -I inputdir -O outputdir"; exit(EXIT_FAILURE);}
	string nome=argv[1];
	if(argc<3)
	{
		infolder="./"
		outfolder="OutDir/"
	}
	else
	{
		if(argc!=6) {cerr<<"Wrong parameter inputs given!!\nCorrect use: ./progname.exe Configfile -I inputdir -O outputdir"; exit(EXIT_FAILURE);}
		string parseop=argv[2];
		string parseop2=argv[4];
		//cout<< parseop << "   " << parseop2<< endl;
		if(parseop!="-I" && parseop2!="-I") {cerr<<"Wrong options name!!\nSuggested use: ./progname.exe Configfile -I inputdir -O outputdir"; exit(EXIT_FAILURE);}
		if(parseop!="-O" && parseop2!="-O") {cerr<<"Wrong options name!!\nSuggested use: ./progname.exe Configfile -I inputdir -O outputdir"; exit(EXIT_FAILURE);}
	
		if(parseop=="-I") {infolder=argv[3]; outfolder=argv[5];}
		else {infolder=argv[5]; outfolder=argv[3];}
	}
	
	map<string,string> options;
	ReadConfig(nome,options);						//Function to be checked
	
	GEM1_gain = 0.0347*exp((0.0209)*stod(options["GEM1_HV"]));
    GEM2_gain = 0.0347*exp((0.0209)*stod(options["GEM2_HV"]));
    GEM3_gain = 0.0347*exp((0.0209)*stod(options["GEM3_HV"]));
    cout<< "GEM1_gain = " << GEM1_gain;
    cout<< "\nGEM2_gain = " << GEM2_gain;
    cout<< "\nGEM3_gain = " << GEM3_gain << endl;
    
    //dividing Fernando's to Francesco&Karolina's single GEM gain measurement
    extraction_eff_GEM1 = 0.87319885*exp(-0.0020000000*stod(options["GEM1_HV"]));
    extraction_eff_GEM2 = 0.87319885*exp(-0.0020000000*stod(options["GEM2_HV"]));
    extraction_eff_GEM3 = 0.87319885*exp(-0.0020000000*stod(options["GEM3_HV"]));
    cout<< "extraction eff GEM1 = " << extraction_eff_GEM1;
    cout<< "\nextraction eff GEM2 = " << extraction_eff_GEM2;
    cout<< "\nextraction eff GEM3 = " << extraction_eff_GEM3 << endl;
	
	double y_dim=stod(options["y_dim"]);
	double demag=y_dim/stod(options["sensor_size"]);
	double aperture=stod(options["camera_aperture"]);
	double omega=1./pow(4.*(demag+1)*aperture,2);
	
	//Code execution
	
	int runcount=stoi(options["historunstart"]);		// historunstart does not exist now but make sense to add
    auto t0 = std::chrono::steady_clock::now();	
    
    if(options["fixed_seed"]=="True" || options["fixed_seed"]=="true") gRandom->SetSeed(0);
    
    vector<int> eventnumber;
    vector<int> particle_type;
    vector<float> energy_ini;
    vector<float> theta_ini;
    vector<float> phi_ini;
    
    if(! experimental::filesystem::exists(outfolder)) system(("mkdir " + outfolder).c_str() );
    
    string ending=".root";
    for(const auto& entry : experimental::filesystem::directory_iterator(infolder))
    {
		string filename=entry.path();
		cout<<filename<<endl;
		bool ends= false;
		if (ending.size() <= filename.size())    ends=equal(ending.rbegin(), ending.rend(), filename.rbegin());
		cout<<ends<<endl;
		//if(filename.ends_with(".root"))   //c++20 fix it wtf
		//{
		if(ends)
		{
			z_ini=0.;
			//int zbins= int(stod(options["zcloud"])/stod(options["z_vox_dim"]));   //line disappeared on original python script?
			auto f= unique_ptr<TFile> {TFile::Open(filename.c_str()) };
			auto inputtree = unique_ptr<TTree> {(TTree*)f->Get("nTuple")};
			string pathread=Form("/media/giorgio/DATA/Giorgio/Documenti/Uni/Dottorato/Lavoro/Saturation_Mango_Sim/%s/%s:/",infolder.c_str(),filename.c_str());
			string fileoutname= Form("histogram_Runs%05d_MC.root",runcount);
			string pathwrite=Form("/media/giorgio/DATA/Giorgio/Documenti/Uni/Dottorato/Lavoro/Saturation_Mango_Sim/digitization/%s/%s:/",outfolder.c_str(),fileoutname.c_str());
			
			auto outfile= shared_ptr<TFile> {TFile::Open(Form("%s/%s",outfolder.c_str(),fileoutname.c_str()),"RECREATE") };
			outfile->mkdir("event_info");
			SaveValues(options,outfile);				//To be checked
			
			
			
			
			
			
			
			f->Close();
		}
	} 
	
	
	
	
	
	auto t1 = std::chrono::steady_clock::now();	
	std::chrono::duration<double> dur=t1-t0;
	cout << "Time taken in seconds is: " << dur.count() << endl; 
	return 0;
}

///////////FUNCTIONS DEFINITION
void ReadConfig(string name, map<string,string>& options)
{
	ifstream config(name.c_str());
	
	string line;
	while(getline(config,line))
	{
		line.erase(remove_if(line.begin(),line.end(),[] (char c){return isspace(c);}),line.end());
		line.erase(remove_if(line.begin(),line.end(),[] (char c){return c=='\'';}),line.end());
		if(line[0] == '#' || line.empty() || line[0] == '{' || line[0] == '}') continue;
		
		auto delim1= line.find(":");
		auto delim2= line.find(",");
		if(delim2==string::npos) delim2=line.size();
		auto index= line.substr(0,delim1);
		auto val= line.substr(delim1+1,min(delim2,line.size())-delim1-1 );
		//cout<<index <<"   " << val << endl;
		options[index]=val;
	}
	
}

int NelGEM2(vector<double> energyDep,const vector<double>& z_hit, map<string,string>& options)
{
	vector<double> n_ioniz_el;
	double opt_pot=stod(options["ion_pot"]);
	transform(energyDep.begin(),energyDep.end(),back_inserter(n_ioniz_el), [&] (double a) { return a/opt_pot;});
	vector<double> drift_l;
	int opt_gem=stoi(options["z_gem"]);
	transform(z_hit.begin(),z_hit.end(),back_inserter(drift_l), [&] (double a) { return abs(a-opt_gem);}); 
	vector<double> n_ioniz_el_mean(n_ioniz_el.size(),0);
	double optabsorption_l=stod(options["absorption_l"]);
	for(int i=0;i<n_ioniz_el_mean.size();i++) n_ioniz_el_mean[i]=abs(n_ioniz_el[i]*exp(-drift_l[i]/optabsorption_l));
	
	
	return 0;
}

void SaveValues(map<string,string>& options, shared_ptr<TFile>& outfile)
{
	outfile->cd();
	outfile->mkdir("param_dir");
	gDirectory->cd("param_dir");
	for(auto const& [key, val] : options)
	{
		if(key!="tag")
		{
			TH1F h(string(key).c_str(),"",1,0,1);
			double value=stod(val);
			h.SetBinContent(1,value);
			h.Write();
		}
	}
	outfile->cd();
	return;
}
