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
#include <filesystem>
#include <memory>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

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
		infolder="../input/other/";
		outfolder="../OutDir/";
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
	
    // DEBUG
    cout<<"Input Folder: "<<infolder<<endl;
    cout<<"Output Folder: "<<outfolder<<endl;
    
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
    
    if(! filesystem::exists(outfolder)){
        //DEBUG
        cout<<"Creating oufolder..."<<
        system(("mkdir " + outfolder).c_str() );
    }
    
    string ending=".root";
    for(const auto& entry : filesystem::directory_iterator(infolder))
    {
		string filename=entry.path();
		cout<<filename<<endl;
		bool ends= false;
		if (ending.size() <= filename.size())    ends=equal(ending.rbegin(), ending.rend(), filename.rbegin());
        
        //DEBUG
		if(!ends) cout<<"Not ending with "<<ending<<endl;
        else cout<<filename<<" is ending with "<<ending<<"!"<<endl;
        
		//if(filename.ends_with(".root"))   //c++20 fix it wtf
		//{
        
		if(ends) {
			z_ini=0.;
			
            auto f         = unique_ptr<TFile> {TFile::Open(filename.c_str())};
            auto inputtree = (TTree*)f->Get("nTuple");
			
			string fileoutname= Form("histogram_Runs%05d_MC.root",runcount);
			auto outfile= shared_ptr<TFile> {TFile::Open(Form("%s/%s",
                                                              outfolder.c_str(),
                                                              fileoutname.c_str()),
                                                         "RECREATE") };
			outfile->mkdir("event_info");
			SaveValues(options,outfile);
            
            // Input file branches
            Int_t eventnumber;
            Int_t numhits;
            Double_t energyDep;
            Double_t energyDep_NR;
            Int_t particle_type;
            vector<int>    pdgID_hits;
            vector<double> tracklen_hits;
            vector<double> px_particle;
            vector<double> py_particle;
            vector<double> pz_particle;
            vector<double> energyDep_hits;
            vector<double> x_hits;
            vector<double> y_hits;
            vector<double> z_hits;
            
            //Output file branches
            int eventnumber_out;
            int particle_type_out;
            float energy;
            float theta;
            float phi;
            float track_length_3D;
            float x_vertex;
            float y_vertex;
            float z_vertex;
            float x_vertex_end;
            float y_vertex_end;
            float z_vertex_end;
            float x_min;
            float x_max;
            float y_min;
            float y_max;
            float z_min;
            float z_max;
            float px;
            float py;
            float pz;
            float proj_track_2D;
            float nhits_og;
            float N_photons;
			
            int max_events = inputtree->GetEntries();
            int totev = (stod(options["events"])==-1) ? max_events : stod(options["events"]);
            totev = min(totev, max_events);
            
            // DEBUG
            //cout<<max_events<<" - "<<totev<<endl;
            
            unique_ptr<TH2D> VignMap;
            if(options["Vignetting"]=="True") {
                string vignfilename = Form("../VignettingMap/%s", options["Vig_Map"].c_str());
                cout<<"Opening "<<vignfilename<<"..."<<endl;
                auto VignFile = unique_ptr<TFile> {TFile::Open(vignfilename.c_str())};
                
                VignMap = std::make_unique<TH2D>(*(TH2D*)VignFile->Get("normmap_lime"));
                
                VignMap->Smooth();
                
                VignMap->SetDirectory(0); // To keep object in memory outside scope
                VignFile->Close();
            }
            
            //DEBUG
            //cout<<"DEBUG: "<<VignMap->GetBinContent(0,0)<<endl;
            
            
			f->Close();
            outfile->Close();
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
	for(unsigned int i=0;i<n_ioniz_el_mean.size();i++) n_ioniz_el_mean[i]=abs(n_ioniz_el[i]*exp(-drift_l[i]/optabsorption_l));
	
	
	return 0;
}

void SaveValues(map<string,string>& options, shared_ptr<TFile>& outfile)
{
	outfile->cd();
	outfile->mkdir("param_dir");
	gDirectory->cd("param_dir");
	for(auto const& [key, val] : options)
	{
		if(key!="tag" and key !="Vig_Map")
		{
			TH1F h(string(key).c_str(),"",1,0,1);
            
            double value;
            if(val == "True")       value = 1.;
            else if(val == "False") value = 0.;
            else                   value  = stod(val);
            
			h.SetBinContent(1,value);
			h.Write();
		}
	}
	outfile->cd();
	return;
}
