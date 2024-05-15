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
#include "s3.h"
#include "cygnolib.h"

using namespace std;


// Global and constant value;
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
void ReadG4Isotopes(string name, map<string,string>& dict_isotopes);
void ReadIonlist(string name, vector<vector<string>>& ionlist);
void SaveValues(map<string,string>& options, shared_ptr<TFile>& outfile);
int NelGEM2(const vector<double>& energyDep,const vector<double>& z_hit, map<string,string>& options);

void AddBckg(map<string,string>& options, int entry, TH2F& background);

// Old approach
//string rootlocation(string tag, int run);   // inconsistency for the MC-old tag!!!
//string root_filename(string tag, int run);

int main(int argc, char** argv)
{
    string outfolder;
    string infolder;
    if(argc<2) {cerr<<"No Configfile given!!\nSuggested use: ./progname.exe Configfile -I inputdir -O outputdir"; exit(EXIT_FAILURE);}
    string nome=argv[1];
    if(argc<3)
    {
        infolder="../input/";
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
    // cout<<"Input Folder: "<<infolder<<endl;
    // cout<<"Output Folder: "<<outfolder<<endl;
    
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
    
    map<string, string> dict_isotopes;
    if(options["GEANT4isotopes"]=="True") ReadG4Isotopes("../Z_A_isotopes.csv", dict_isotopes);
    
    
    
    
    
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
            
            
            
            vector<vector<string>> SRIM_events;
            
            if(options["NR"]=="True" and options["NR_list"]!="") {
                auto delim1 = filename.find("part");
                auto delim2 = filename.find(".root");
                if(delim1==string::npos) throw runtime_error("Cannot determine the 'part' of the file.\n");
                if(delim2==string::npos) throw runtime_error("Input file is not a root file.\n");
                auto part= filename.substr(delim1+4, delim2-delim1-4);
                cout<<"Using NR list from "<<options["NR_list"]<<"_part"<<part<<".py"<<endl;
                
                ReadIonlist(Form("%s/%s_part%s.py", infolder.c_str(), options["NR_list"].c_str(), part.c_str()),
                            SRIM_events);
                
            }
            
            
            z_ini=0.;
            
            auto f         = unique_ptr<TFile> {TFile::Open(filename.c_str())};
            auto inputtree = (TTree*)f->Get("nTuple");
            
            string fileoutname= Form("histogram_Runs%05d_MC.root",runcount);
            auto outfile = shared_ptr<TFile> {TFile::Open(Form("%s/%s",
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
            vector<int>    *pdgID_hits = 0;
            vector<double> *tracklen_hits = 0;
            vector<double> *px_particle = 0;
            vector<double> *py_particle = 0;
            vector<double> *pz_particle = 0;
            vector<double> *energyDep_hits = 0;
            vector<double> *x_hits = 0;
            vector<double> *y_hits = 0;
            vector<double> *z_hits = 0;
            
            // Some of the following variables not present in (old?) NR simulations???
            inputtree->SetBranchAddress("eventnumber", &eventnumber);
            inputtree->SetBranchAddress("numhits", &numhits);
            inputtree->SetBranchAddress("energyDep", &energyDep);
            inputtree->SetBranchAddress("energyDep_NR", &energyDep_NR);
            inputtree->SetBranchAddress("pdgID_hits", &pdgID_hits);
            inputtree->SetBranchAddress("tracklen_hits", &tracklen_hits);
            inputtree->SetBranchAddress("px_particle", &px_particle);
            inputtree->SetBranchAddress("py_particle", &py_particle);
            inputtree->SetBranchAddress("pz_particle", &pz_particle);
            inputtree->SetBranchAddress("energyDep_hits", &energyDep_hits);
            inputtree->SetBranchAddress("x_hits", &x_hits);
            inputtree->SetBranchAddress("y_hits", &y_hits);
            inputtree->SetBranchAddress("z_hits", &z_hits);
            
            if(options["NR"]=="True") {
                inputtree->SetBranchAddress("particle_type", &particle_type);
            }
            
            
            
            
            
            //Output file branches
            Int_t eventnumber_out   = -999;
            Int_t particle_type_out = -999;
            Float_t energy = -999;
            Float_t theta  = -999;
            Float_t phi    = -999;
            Float_t track_length_3D = -1;
            Float_t x_vertex = -1;
            Float_t y_vertex = -1;
            Float_t z_vertex = -1;
            Float_t x_vertex_end = -1;
            Float_t y_vertex_end = -1;
            Float_t z_vertex_end = -1;
            Float_t x_min = -1;
            Float_t x_max = -1;
            Float_t y_min = -1;
            Float_t y_max = -1;
            Float_t z_min = -1;
            Float_t z_max = -1;
            Float_t px = -1;
            Float_t py = -1;
            Float_t pz = -1;
            Float_t proj_track_2D = -1;
            Int_t nhits_og = -1;
            Int_t N_photons = -1;
            
            auto outtree = new TTree("event_info", "event_info");
            
            outtree->Branch("eventnumber", &eventnumber_out, "eventnumber/I");
            outtree->Branch("particle_type", &particle_type_out, "particle_type/I");
            outtree->Branch("energy", &energy, "energy/F");
            outtree->Branch("theta", &theta, "theta/F");
            outtree->Branch("phi", &phi, "phi/F");
            outtree->Branch("track_length_3D", &track_length_3D, "track_length_3D/F");
            outtree->Branch("proj_track_2D", &proj_track_2D, "proj_track_2D/F");
            outtree->Branch("x_vertex", &x_vertex, "x_vertex/F");
            outtree->Branch("y_vertex", &y_vertex, "y_vertex/F");
            outtree->Branch("z_vertex", &z_vertex, "z_vertex/F");
            outtree->Branch("x_vertex_end", &x_vertex_end, "x_vertex_end/F");
            outtree->Branch("y_vertex_end", &y_vertex_end, "y_vertex_end/F");
            outtree->Branch("z_vertex_end", &z_vertex_end, "z_vertex_end/F");
            outtree->Branch("x_min", &x_min, "x_min/F");
            outtree->Branch("x_max", &x_max, "x_max/F");
            outtree->Branch("y_min", &y_min, "y_min/F");
            outtree->Branch("y_max", &y_max, "y_max/F");
            outtree->Branch("z_min", &z_min, "z_min/F");
            outtree->Branch("z_max", &z_max, "z_max/F");
            outtree->Branch("N_photons", &N_photons, "N_photons/I");
            outtree->Branch("px", &px, "px/F");
            outtree->Branch("py", &py, "py/F");
            outtree->Branch("pz", &pz, "pz/F");
            outtree->Branch("nhits_og", &nhits_og, "nhits_og/I");
            
            
            int max_events = inputtree->GetEntries();
            int totev = (stod(options["events"])==-1) ? max_events : stod(options["events"]);
            totev = min(totev, max_events);
            
            
            unique_ptr<TH2F> VignMap;
            if(options["Vignetting"]=="True") {
                string vignfilename = Form("../VignettingMap/%s", options["Vig_Map"].c_str());
                cout<<"Opening "<<vignfilename<<"..."<<endl;
                auto VignFile = unique_ptr<TFile> {TFile::Open(vignfilename.c_str())};
                
                VignMap = std::make_unique<TH2F>(*(TH2F*)VignFile->Get("normmap_lime"));
                
                VignMap->Smooth();
                
                VignMap->SetDirectory(0); // To keep object in memory outside scope
                VignFile->Close();
            }
            
            //DEBUG
            //cout<<"DEBUG: "<<VignMap->GetBinContent(0,0)<<endl;
            
            for(int entry=0; entry<totev; entry++) {  // RUNNING ON ENTRIES
                
                inputtree->GetEntry(entry);
                
                //DEBUG
                //if(options["NR"]=="True") cout<<particle_type<<endl;
                //if (entry==0) {
                //    cout<<numhits<<" - "<<pdgID_hits->size()<<endl;
                //    for(unsigned int i=0; i<pdgID_hits->size(); i++) {
                //        cout<<"---"<< (*pdgID_hits)[i] <<endl;
                //    }
                //}
                
                //cout<<"Entry "<<entry<<" of "<<totev<<endl;
                //cout<<"Energy "<<energyDep<<" keV"<<endl;
                
                if (energyDep>200) {
                    continue;
                }
                
                TH2F background;
                AddBckg(options, entry, background);
                
                
                
                
                eventnumber_out = eventnumber;
                
                
                outtree->Fill();
                gDirectory->cd();
                
                // DEBUG
                //background.Write();
                
            }
            gDirectory->cd("event_info");
            outtree->Write();
            
            
            f->Close();
            outfile->Close();
        }
        
    }
	
    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> dur=t1-t0;
    cout << "Time taken in seconds is: " << dur.count() << endl;
    return 0;
    
}

// ====== FUNCTIONS DEFINITION ======

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

void ReadG4Isotopes(string name, map<string,string>& dict_isotopes) {
    ifstream config(name.c_str());
    
    string line;
    while(getline(config, line))
    {
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return isspace(c);}),line.end());
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return c=='\'';}),line.end());
        if(line[0] == '#' || line.empty() || line[0] == '{' || line[0] == '}') continue;
        
        stringstream to_split(line.c_str());
        string element;
        
        int counter = 0;
        string index;
        string val;
        while(getline(to_split, element, ',')) {
            if(counter == 0) index = element;
            else if (counter == 1) val = element;
            else if (counter == 2) val += Form("%03d", stoi(element));
            counter++;
        }
        dict_isotopes[index]=val;
        
        // DEBUG
        // cout<<"---"<<index<<": "<<val<<endl;
    }
    
    return;
}


void ReadIonlist(string name, vector<vector<string>>& ionlist) {
    
    cout<<"Opening and parsing "<<name<<" ..."<<endl;
    ifstream config(name.c_str());
    
    string line;
    while(getline(config, line))
    {
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return isspace(c);}),line.end());
        if(line.find("ionlist")!=string::npos) continue;
        stringstream to_split(line.c_str());
        
        string element;
        while(getline(to_split, element, '[')) {
            if(!element.empty()){
                element.erase(remove_if(element.begin(),element.end(),[=] (char c){return c==']';}), element.end());
                stringstream content(element.c_str());

                string var;
                vector<string> row;
                while(getline(content, var, ',')) {
                    if(!var.empty()){
                        row.push_back(var);
                    }
                }
                ionlist.push_back(row);
            }
        }
        break;
    }
    return;
}

int NelGEM2(vector<double> energyDep,const vector<double>& z_hit, map<string,string>& options)
{
    vector<double> n_ioniz_el;
    double opt_pot=stod(options["ion_pot"]);
    transform(energyDep.begin(),energyDep.end(),back_inserter(n_ioniz_el), [&] (double a) { return a/opt_pot;});
    
    vector<double> drift_l;
    int opt_gem=stod(options["z_gem"]);
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
        // DEBUG
        // cout<<key<<": "<<val<<endl;
        
        if(key!="tag"       and key !="Vig_Map" and
           key!="bckg_path" and key !="ped_cloud_dir" and
           key!="NR_list"
           )
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

void AddBckg(map<string,string>& options, int entry, TH2F& background) {
    
    string tmpfolder = options["bckg_path"];
    int x_pix = stoi(options["x_pix"]);
    int y_pix = stoi(options["y_pix"]);
    
    
    if(! filesystem::exists(tmpfolder)){
        //DEBUG
        cout<<"Creating tmpfolder..."<<
        system(("mkdir " + tmpfolder).c_str() );
    }
    
    
    if(options["bckg"]=="True") {
        
        string tmpname = Form( "%s/run%05d.mid.gz", tmpfolder.c_str(), stoi(options["noiserun"]));
        
        
        if(! filesystem::exists(tmpname)) {
            // DEBUG
            cout << "Looking for file: "<<tmpname<<endl;
            bool verbose = false;
            bool cloud   = true;
            
            int run = stoi(options["noiserun"]);
            
            //Download or find midas file
            string filename = s3::cache_file(s3::mid_file(run, options["ped_cloud_dir"], cloud, verbose),
                                             options["bckg_path"],
                                             cloud,
                                             options["ped_cloud_dir"],
                                             verbose);
        }
        
        int pic_index = gRandom->Integer(100);
        
        // DEBUG
        cout<<"Using pic # "<<pic_index<<" as a pedestal..."<<endl;
        
        //reading data from midas file
        cout<<"Opening midas file "<<tmpname<<" ..."<<endl;
        TMReaderInterface* reader = cygnolib::OpenMidasFile(tmpname);
        bool reading = true;
        
        int counter = 0;
        bool found = false;
        while (reading) {
            
            TMidasEvent event = TMidasEvent();
            reading = TMReadEvent(reader, &event);
            if (!reading) {
                // DEBUG
                // std::cout<<"EOF reached."<<std::endl;
                break;
            }
            
            bool cam_found = cygnolib::FindBankByName(event, "CAM0");
            if(cam_found) {
                if(counter == pic_index) {
                    cygnolib::Picture pic=cygnolib::daq_cam2pic(event, "fusion");
                    
                    TH2F rootpic(Form("pic_run1_ev%d", entry) , Form("pic_run1_ev%d", entry) , x_pix, -0.5, x_pix-0.5, y_pix, -0.5, y_pix-0.5);
                    
                    vector<vector<uint16_t>> vecpic = pic.GetFrame();
                    for(unsigned int i = 0; i<vecpic.size(); i++) {
                        for (unsigned int j =0; j<vecpic[0].size(); j++) {
                            rootpic.SetBinContent(j, i, vecpic[i][j]);
                        }
                    }
                    background=rootpic;
                    
                    found = true;
                    break;
                }
                
                counter++;
            }
            
        }
        if(!found) {
            cerr<<"AddBckg: Error: Cannot find pic # "<<pic_index<<" in pedestal run."<<endl;
            exit(EXIT_FAILURE);
        }
        
    }
    
    return;
}

// Old approach:
//string rootlocation(string tag, int run){
//    string sel;
//    if (tag == "Data") {
//        if ((run>=936) && (run<=1601)) {
//            sel = "Data/LTD/Data_Camera/ROOT";
//        } else if ((run>=1632) && (run<4505)) {
//            sel = "Data/LAB";
//        } else if ((run>=4470) && (run<10000)) {
//            sel = "LAB";
//        } else {
//            cerr<<"rootlocation: Error: Data taken with another DAQ or not yet uploaded to the cloud"<<endl;
//            exit(EXIT_FAILURE);
//        }
//    } else if (tag == "DataMango") {
//        sel = (run<3242) ? "Data/Man" : "MAN";
//    } else if (tag ==  "MC") {
//        sel = "Simulation";
//        cerr<<"rootlocation: Error: automatic download for Simulated data not implemented yet"<<endl;
//        exit(EXIT_FAILURE);
//    }
//
//    return sel;
//}

//string root_filename(string tag, int run) {
//
//    string sel = rootlocation(tag,run);
//
//    string BASE_URL = "https://s3.cloud.infn.it/v1/AUTH_2ebf769785574195bde2ff418deac08a/";
//    string bucket;
//    if(tag.find("MC") != std::string::npos) {
//        bucket = (tag=="MC-old") ? "cygnus" : "cygno-sim";
//    } else if (tag == "Data") {
//        bucket = (run<4505) ? "cygnus" : "cygno-data";
//    } else if (tag == "DataMango") {
//        bucket = (run<3242) ? "cygnus" : "cygno-data";
//    }
//
//    BASE_URL = Form("%s%s/", BASE_URL.c_str(), bucket.c_str());
//    string file_root = Form("%s/histogr_Run%05d.root", sel.c_str(), run);
//
//    return Form("%s%s", BASE_URL.c_str(), file_root.c_str());
//}
//
