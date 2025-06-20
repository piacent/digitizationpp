#include "DigitizationRunner.h"
#include "Globals.h"
#include "Utils.h"
#include <filesystem>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;

DigitizationRunner::DigitizationRunner(ConfigManager& configMgr,
                                   const std::string& inputDir_,
                                   const std::string& outputDir_)
    : config(configMgr),
      inputDir(inputDir_),
      outputDir(outputDir_),
      runCount(config.getInt("start_run_number"))
{
    initializeGlobals();
    prepareCameraSettings();
    checkDimensionConsistency();
    
    num_threads = config.getInt("Parallel_threads");

    // DEBUG
    // cout<<"Initialized globals:"<<endl;
    // cout<<"GEM1_gain = "<< GEM1_gain <<" - GEM2_gain = " << GEM2_gain<<" - GEM3_gain = " << GEM3_gain << endl;
    // cout<<"extraction_eff_GEM1 = "<< extraction_eff_GEM1 << " - extraction_eff_GEM2 = "
    //     << extraction_eff_GEM2 << " - extraction_eff_GEM3 = " << extraction_eff_GEM3 << endl;
    // cout<<"omega = "<< omega <<endl;
    // cout<<"x_pix, y_pix = "<<x_pix<<", "<<y_pix<<endl;
    // cout<<"optcounts_per_photon = "<<optcounts_per_photon<<endl;
    // cout<<"y_sensor_size = "<<y_sensor_size<<endl;
    // cout<<"readout_time = "<<readout_time<<endl;
    // cout<<"num_threads = "<<num_threads<<endl;
    
}

void DigitizationRunner::run() {
    auto t0 = std::chrono::steady_clock::now();

    if (config.getBool("fixed_seed")) {
        gRandom->SetSeed(10);
    }

    processRootFiles();

    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> dur = t1 - t0;
    cout << "Total digitization time: " << dur.count() << " seconds" << endl;
}

void DigitizationRunner::initializeGlobals() {
    GEM1_gain = 0.0347 * exp(0.0209 * config.getDouble("GEM1_HV"));
    GEM2_gain = 0.0347 * exp(0.0209 * config.getDouble("GEM2_HV"));
    GEM3_gain = 0.0347 * exp(0.0209 * config.getDouble("GEM3_HV"));

    extraction_eff_GEM1 = 0.87319885 * exp(-0.002 * config.getDouble("GEM1_HV"));
    extraction_eff_GEM2 = 0.87319885 * exp(-0.002 * config.getDouble("GEM2_HV"));
    extraction_eff_GEM3 = 0.87319885 * exp(-0.002 * config.getDouble("GEM3_HV"));
}

void DigitizationRunner::prepareCameraSettings() {
    std::string camera = config.get("Camera_type");

    y_pix = 2304;
    if (camera == "Fusion") {
        x_pix = 2304;
        optcounts_per_photon = 2.; // apparently calibrated with LEMOn in the past
        y_sensor_size = 14.976;    // mm
        readout_time  = 184.4;     // ms in in ultra quiet scan (UQS)
    } else if (camera == "Quest1" || camera == "Quest2") {
        x_pix = 4096;
        optcounts_per_photon = 4.49; // equal to 2.*2.245 which is the ratio of e/count of the two cameras
        y_sensor_size = 10.598;      // mm
        readout_time = (camera == "Quest1") ? 199.0 : 39.0; // ms in in ultra quiet scan (UQS)
    }

    double y_dim = config.getDouble("y_dim");
    double demag = y_dim / y_sensor_size;
    double aperture = config.getDouble("camera_aperture");
    omega = 1. / pow(4. * (demag + 1) * aperture, 2);

}

void DigitizationRunner::checkDimensionConsistency() {
    double x_dim = config.getDouble("x_dim");
    double y_dim = config.getDouble("y_dim");
    double testvalue = (x_dim / y_dim) - (double(x_pix) / y_pix);

    if (testvalue > 1e-1) {
        cerr << "Using Quest image dimension with Fusion simulation! Digitization FAILED!" << endl;
        exit(EXIT_FAILURE);
    } else if (testvalue < -1e-1) {
        cerr << "Using Fusion image dimension with Quest simulation! Digitization FAILED!" << endl;
        exit(EXIT_FAILURE);
    }
}

bool DigitizationRunner::isValidInputFile(const std::string& filename) const {
    bool fn_size = (filename.size() >= 5);
    if (!fn_size) {
        return false;
    } else if (fn_size && filename.substr(filename.size() - 5) == ".root") {
        return true;
    }
    return false;
}

void DigitizationRunner::setSeed(int seed) {
    if(config.getBool("fixed_seed")) gRandom->SetSeed(seed);
    return;
}

void DigitizationRunner::processRootFiles() {

    DigitizationRunner::setSeed(10);
    
    string infolder  = Utils::resolvePath(inputDir);
    string outfolder = Utils::resolvePath(outputDir);

    // Probably not needed anymore
    vector<int> eventnumber;
    vector<int> particle_type;
    vector<float> energy_ini;
    vector<float> theta_ini;
    vector<float> phi_ini;

    //std::string inDir = ConfigManager::resolvePath(inputDir);
    for (const auto& fentry : filesystem::directory_iterator(infolder)) {
        const string filename = fentry.path().string();
        if (!isValidInputFile(filename)) continue;

        cout << "Processing: " << filename << endl;

        ///// aaaaaaaaaaa

        //DEBUG
        if(filename.find("HeCF4gas_AmBe_part") != string::npos) {
        //if(filename.find("LIME_CADshield") != string::npos) {
            continue;
        }

        // For SRIM: TO BE TESTED ON SRIM SIMS
        vector<vector<string>> SRIM_events;
        if(config.getBool("NR") && config.get("NR_list")!="") {
            cout<<"debug"<<endl;
            auto delim1 = filename.find("part");
            auto delim2 = filename.find(".root");
            if(delim1==string::npos) throw runtime_error("Cannot determine the 'part' of the file.\n");
            if(delim2==string::npos) throw runtime_error("Input file is not a root file.\n");
            auto part= filename.substr(delim1+4, delim2-delim1-4);
            cout<<"Using NR list from "<<config.get("NR_list")<<"_part"<<part<<".py"<<endl;


            string nrlist = config.get("NR_list").c_str();
            config.loadIonList(Form("%s/%s_part%s.py",
                                    infolder.c_str(), nrlist.c_str(), part.c_str()),
                               SRIM_events);
        }
        
        auto f         = unique_ptr<TFile> {TFile::Open(filename.c_str())};
        auto inputtree = (TTree*)f->Get("nTuple");
            
        int max_events = inputtree->GetEntries();
        int totev = (config.getInt("events") == -1) ? max_events : config.getInt("events");
        totev = min(totev, max_events);
            
        int firstentry = config.getInt("start_event");
        if (firstentry>totev) throw runtime_error("Error: First entry is larger than last entry, exiting!");
        cout << "Processing entries from "<<firstentry<<" to "<<totev<<"."<<endl;
            
            
        // output saved in outfolder/filename/
        stringstream ssfilename(filename.c_str());
        string ssfilename_el, filename_tmp;
        while(getline(ssfilename, ssfilename_el, '/')) {
            filename_tmp = ssfilename_el;
        }
        auto delimFN = filename_tmp.find(".root");
        string basefilename   = filename_tmp.substr(0, delimFN);
        string fnameoutfolder = outfolder + "/" + basefilename;
        if(! filesystem::exists(fnameoutfolder) && !config.getBool("queue")){
            system(("mkdir -p " + fnameoutfolder).c_str() );
        }
            
        // standard: name of output file = histograms_RunRRRRR.root (R run number)
        if(config.getBool("queue")) fnameoutfolder="./";
        string fileoutname= Form("%s/histograms_Run%05d.root",
                                 fnameoutfolder.c_str(),
                                 runCount);
            
        // for radioisotope simulation: histograms_RunZZAAANN.root (Z=atomic number, A=mass numbe$
        // NOTE: this is a 7 digit run number, while reconstruction currently looks for 5
        string isot_numb = "0000000";
        if(config.getBool("GEANT4isotopes")) {
            cout<<"GEANT4isotopes option is active."<<endl;
            
            stringstream ssinfile(basefilename);
            string tmpstr;
            int counter = 0;
            while(getline(ssinfile, tmpstr, '_')) {
                if(counter == 1) {
                    isot_numb = config.getIsotope(tmpstr);
                    // DEBUG
                    //isot_numb = "00000";
                }
                counter++;
            }
                
            auto delimBFN = basefilename.find("part");
            if(delimBFN==string::npos) throw runtime_error("Cannot determine the 'part' of the file.\n");
            auto part = basefilename.substr(delimBFN+4, basefilename.size() - delimBFN - 4);
                    
            if(filename.find("part")!= string::npos) {
                fileoutname = Form("%s/histograms_Run%05d%02d.root",
                                   fnameoutfolder.c_str(),
                                   stoi(isot_numb),
                                   stoi(part)
                                   );
                isot_numb = Form("%05d%02d", stoi(isot_numb), stoi(part));
            } else {
                fileoutname = Form("%s/histograms_Run%05d00.root",
                                   fnameoutfolder.c_str(),
                                   stoi(isot_numb));
                
                isot_numb = Form("%05d00", stoi(isot_numb));
            }
            //DEBUG
            //cout<<"DEBUG: fileoutname = "<<fileoutname<<endl;
            //cout<<"DEBUG: isot_numb = "<<isot_numb<<endl;
        }

        if(config.get("start_event")!="0"){
            cout<<"out folder "<<fnameoutfolder<<endl;
            int newpart = (int)(config.getInt("start_event")/500);
            int oldpart = stoi(isot_numb);
            int partnum = oldpart + newpart;
            fileoutname = Form("%s/histograms_Run%05d.root",
                               fnameoutfolder.c_str(),
                               partnum);
        }
                        
            
        auto outfile = shared_ptr<TFile> {TFile::Open(fileoutname.c_str(),
                                                          "RECREATE") };
        //outfile->mkdir("event_info");
        
        /*
        SaveValues(options,outfile);
            
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
            
        Int_t row_cut = -1;
        Int_t N_photons_cut = -1;
        Float_t cut_energy = -1;
        Float_t x_min_cut = -1;
        Float_t x_max_cut = -1;
        Float_t y_min_cut = -1;
        Float_t y_max_cut = -1;
        Float_t z_min_cut = -1;
        Float_t z_max_cut = -1;
        Float_t proj_track_2D_cut = -1;
            
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
        if (config.getBool("exposure_time_effect")) {
            outtree->Branch("N_photons_cut", &N_photons_cut, "N_photons_cut/I");
            outtree->Branch("row_cut", &row_cut, "row_cut/I");
            outtree->Branch("cut_energy", &cut_energy, "cut_energy/F");
            outtree->Branch("x_min_cut", &x_min_cut, "x_min_cut/F");
            outtree->Branch("x_max_cut", &x_max_cut, "x_max_cut/F");
            outtree->Branch("y_min_cut", &y_min_cut, "y_min_cut/F");
            outtree->Branch("y_max_cut", &y_max_cut, "y_max_cut/F");
            outtree->Branch("z_min_cut", &z_min_cut, "z_min_cut/F");
            outtree->Branch("z_max_cut", &z_max_cut, "z_max_cut/F");
            outtree->Branch("proj_track_2D_cut", &proj_track_2D_cut, "proj_track_2D_cut/F");
        }
            
        // Input file branches
        Int_t eventnumber;
        Int_t numhits;
        Double_t energyDep;
        Double_t energyDep_NR;
        Float_t  ekin_particle;
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
            
            
        inputtree->SetBranchAddress("eventnumber", &eventnumber);
        inputtree->SetBranchAddress("numhits", &numhits);
        
        if(!config.getBool("SRIM")) { // TO BE CHECK ON REAL SRIM SIMULATIONS
            inputtree->SetBranchAddress("energyDep",    &energyDep);
            inputtree->SetBranchAddress("energyDep_NR", &energyDep_NR);
            inputtree->SetBranchAddress("pdgID_hits",    &pdgID_hits);
            inputtree->SetBranchAddress("tracklen_hits", &tracklen_hits);
            inputtree->SetBranchAddress("px_particle", &px_particle);
            inputtree->SetBranchAddress("py_particle", &py_particle);
            inputtree->SetBranchAddress("pz_particle", &pz_particle);
        }
            
        inputtree->SetBranchAddress("energyDep_hits", &energyDep_hits);
        inputtree->SetBranchAddress("x_hits", &x_hits);
        inputtree->SetBranchAddress("y_hits", &y_hits);
        inputtree->SetBranchAddress("z_hits", &z_hits);
            
        if(config.getBool("NR")) { // TO BE CHECK ON REAL SRIM SIMULATIONS
            inputtree->SetBranchAddress("particle_type", &particle_type);
            inputtree->SetBranchAddress("ekin_particle", &ekin_particle);
        }

        TH2F VignMap;
        /*if(config.getBool("Vignetting")) {
            string vignfilename = Form("%sVignettingMap/%s",SOURCE_DIR.c_str(), options["Vig_Map"].c_str());
                cout<<"Opening "<<vignfilename<<"..."<<endl;
                auto VignFile = unique_ptr<TFile> {TFile::Open(vignfilename.c_str())};
                
                VignMap = (*(TH2F*)VignFile->Get("normmap"));
                
                VignMap.Smooth();
                
                VignFile->Close();

                int x_vign=VignMap.GetNbinsX();
                int y_vign=VignMap.GetNbinsY();

                if(x_vign!=y_vign)
                {
                    if(x_pix==y_pix)
                    {
                        cerr<<"You are using a Quest vignette map with a Fusion simulation! Digitization FAILED!\n";
                        exit(0);
                    }
                }
                else
                {
                    if(x_pix!=y_pix)
                    {
                        cerr<<"You are using a Fusion vignette map with a Quest simulation! Digitization FAILED!\n";
                        exit(0);
                    }

                }
            }
            
            //DEBUG
            //cout<<"DEBUG: "<<VignMap.GetBinContent(100,100)<<endl;
            
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
                
                cout<<"\nEntry "<<entry<<" of "<<totev<<endl;
                
                if (options["SRIM"]=="True"){
                    cout<<"Energy "<<ekin_particle<<" keV"<<endl;
                } else {
                    cout<<"Energy "<<energyDep    <<" keV"<<endl;
                }

                bool NR_flag=false;
                if(options["SRIM"] == "True") {
                    energy            = ekin_particle;
                    particle_type_out = particle_type;
                    //if(particle_type_out==??) NR_flag=true;       Not known output from SRIM
                } else {
                    // this would be the energy of the primary particle - equal to
                    // deposited energy only if it is completely contained in the sensitive volume
                    // energy = ekin_particle * 1000;
                    energy = energyDep;
                    if (energyDep_NR>0){
                        particle_type_out = is_NR(*pdgID_hits, int(1.e9));
                        NR_flag = true;
                    } 
                    else particle_type_out = (*pdgID_hits)[0];  // this will tell us if the deposit was
                    // started by a photon or an electron
                }

                if(options["NR"]=="False" && NR_flag==true ) continue;
                if(options["SRIM"]=="True"  && ekin_particle>900) continue;     //not corrected for SRIM
                
                //initialize array values - to save info also if the track is skipped (background only)
                row_cut         = -1;
                eventnumber_out = eventnumber;
                
                
                // DEBUG
                //cout<<"energyDep_NR = "<<energyDep_NR<<endl;
                //cout<<"particle_type_out = "<<particle_type_out<<endl;
                
                //particle_type_out = -999;
                //energy = -1;
                cut_energy      = -1;
                theta           =  0;
                phi             =  0;
                track_length_3D = -1;
                proj_track_2D   = -1;
                x_vertex        = -1;
                y_vertex        = -1;
                z_vertex        = -1;
                x_vertex_end    = -1;
                y_vertex_end    = -1;
                z_vertex_end    = -1;
                x_min           = -1;
                x_max           = -1;
                y_min           = -1;
                y_max           = -1;
                z_min           = -1;
                z_max           = -1;
                N_photons       =  0;
                x_min_cut       = -1;
                x_max_cut       = -1;
                y_min_cut       = -1;
                y_max_cut       = -1;
                z_min_cut       = -1;
                z_max_cut       = -1;
                N_photons_cut   =  0;
                proj_track_2D_cut = -1;
                px    = 0;
                py    = 0;
                pz    = 0;
                nhits_og = numhits;


                vector<vector<int>> background(x_pix,
                                               vector<int>(y_pix, 0));
                                               
                AddBckg(options, background);
                
                if (energy < stod(options["ion_pot"])){
                    energy = 0;
                    TH2I final_image(Form("pic_run%d_ev%d", run_count, entry), "",
                                     x_pix, -0.5, x_pix -0.5,
                                     y_pix, -0.5, y_pix -0.5);
                    
                    for(unsigned int xx =0; xx < background.size(); xx++) {
                        for(unsigned int yy =0; yy < background[0].size(); yy++) {
                            final_image.SetBinContent(xx+1, yy+1, background[xx][yy]);
                        }
                    }
                    
                    outtree->Fill();
                    outfile->cd();
                    final_image.Write();
                    
                    continue;
                }
                
                vector<double> x_hits_tr;
                vector<double> y_hits_tr;
                vector<double> z_hits_tr;
                

                if (options["SRIM"]=="True") {
                    // x_hits_tr = np.array(tree.x_hits) + opt.x_offset
                    // y_hits_tr = np.array(tree.y_hits) + opt.y_offset
                    // z_hits_tr = np.array(tree.z_hits) + opt.z_offset
                    vector<double> v1 = {1.,0.,0.};
                    vector<double> v2 = {stod(SRIM_events[entry][3])-stod(SRIM_events[entry][2]),
                                         stod(SRIM_events[entry][5])-stod(SRIM_events[entry][4]),
                                         stod(SRIM_events[entry][7])-stod(SRIM_events[entry][6]),
                                        };
                    
                    double        angle = angle_between(v1, v2);
                    vector<double> axis =  crossProduct(v1, v2);
                    
                    
                    double norm = sqrt(inner_product(axis.begin(), axis.end(), axis.begin(), 0.0));
                    vector<double> uaxis = {axis[0]/norm, axis[1]/norm, axis[2]/norm};
                    // DEBUG
                    //std::cout<<angle<<endl;
                    //std::cout<<"-"<<axis[0]<<","<<axis[1]<<","<<axis[2]<<endl;
                    
                    
                    for(int ihit=0; ihit < numhits; ihit++) {
                        vector<double> tmpvec = {(*x_hits)[ihit], (*y_hits)[ihit], (*z_hits)[ihit]};
                        vector<double> rotvec = rotateByAngleAndAxis(tmpvec, angle, uaxis);
                        
                        x_hits_tr.push_back(rotvec[0]+stod(SRIM_events[entry][2])+stod(options["x_offset"]));
                        y_hits_tr.push_back(rotvec[1]+stod(SRIM_events[entry][4])+stod(options["y_offset"]));
                        z_hits_tr.push_back(rotvec[2]+stod(SRIM_events[entry][6])+stod(options["z_offset"]));
                    }
                    
                } else {
                    transform(z_hits->begin(),
                              z_hits->end(),
                              back_inserter(x_hits_tr),
                              [&] (double a) {return a + stod(options["x_offset"]);});
                    transform(y_hits->begin(),
                              y_hits->end(),
                              back_inserter(y_hits_tr),
                              [&] (double a) {return a + stod(options["y_offset"]);});
                    transform(x_hits->begin(),
                              x_hits->end(),
                              back_inserter(z_hits_tr),
                              [&] (double a) {return a + stod(options["z_offset"]);});
                    
                    // FIXME: [Check which is the z axis orientation]
                    
                }
                
                // DEBUG
                //for(int ihit=0; ihit < numhits; ihit++) {
                //    cout<<x_hits_tr[ihit]<<",";
                //    cout<<y_hits_tr[ihit]<<",";
                //    cout<<z_hits_tr[ihit]<<"\n";
                //}
                
                vector<double> energy_hits = (*energyDep_hits);
                
                // add random Z to tracks
                if (stod(options["randZ_range"]) != 0) {
                    double rand = (gRandom->Uniform() - 0.5) * stod(options["randZ_range"]);
                    //DEBUG
                    //cout<<"rand = "<<rand<<endl;
                    transform(z_hits_tr.begin(), z_hits_tr.end(), z_hits_tr.begin(),
                              [&] (double a) {return a + rand;}
                              );
                    
                }
                
                //Compute length and extremes of the track before the cut
                proj_track_2D = 0;
                for(int ihit=0; ihit < numhits-1; ihit++){
                    proj_track_2D += sqrt((x_hits_tr[ihit+1]-x_hits_tr[ihit])*(x_hits_tr[ihit+1]-x_hits_tr[ihit])+
                                          (y_hits_tr[ihit+1]-y_hits_tr[ihit])*(y_hits_tr[ihit+1]-y_hits_tr[ihit])
                                          );
                }
                // DEBUG
                //cout<<"proj_track_2D = "<<Form("%.10f", proj_track_2D)<<endl;
                
                
                x_vertex = (x_hits_tr[0] + 0.5 * stod(options["x_dim"]) )*static_cast<double>(x_pix)/stod(options["x_dim"]); //in pixels
                y_vertex = (y_hits_tr[0] + 0.5 * stod(options["y_dim"]) )*static_cast<double>(y_pix)/stod(options["y_dim"]); //in pixels
                z_vertex = (z_hits_tr[0]+stod(options["z_extra"])); //distance from GEMs in mm
                // DEBUG
                //cout<<"x_vertex = "<<x_vertex<<" ### y_vertex = "<<y_vertex<<" ### z_vertex = "<<z_vertex<<endl;
                
                x_vertex_end = (x_hits_tr[numhits-1] + 0.5 * stod(options["x_dim"])) * static_cast<double>(x_pix) / stod(options["x_dim"]); //in pixels
                y_vertex_end = (y_hits_tr[numhits-1] + 0.5 * stod(options["y_dim"])) * static_cast<double>(y_pix) / stod(options["y_dim"]); //in pixels
                z_vertex_end = (z_hits_tr[numhits-1]+stod(options["z_extra"])); //distance from GEMs in mm
                //DEBUG
                //cout<<"x_vertex_end = "<<x_vertex_end<<" ### y_vertex_end = "<<y_vertex_end<<" ### z_vertex_end = "<<z_vertex_end<<endl;
                
                x_min = (*min_element(x_hits_tr.begin(), x_hits_tr.end()) + 0.5 * stod(options["x_dim"])) * static_cast<double>(x_pix) / stod(options["x_dim"]);
                x_max = (*max_element(x_hits_tr.begin(), x_hits_tr.end()) + 0.5 * stod(options["x_dim"])) * static_cast<double>(x_pix) / stod(options["x_dim"]);
                y_min = (*min_element(y_hits_tr.begin(), y_hits_tr.end()) + 0.5 * stod(options["y_dim"])) * static_cast<double>(y_pix) / stod(options["y_dim"]);
                y_max = (*max_element(y_hits_tr.begin(), y_hits_tr.end()) + 0.5 * stod(options["y_dim"])) * static_cast<double>(y_pix) / stod(options["y_dim"]);
                z_min = min((*max_element(z_hits_tr.begin(),
                                             z_hits_tr.end()) + stod(options["z_extra"])),
                            (*min_element(z_hits_tr.begin(),
                                             z_hits_tr.end()) + stod(options["z_extra"])));
                z_max = max((*max_element(z_hits_tr.begin(),
                                             z_hits_tr.end()) + stod(options["z_extra"])),
                            (*min_element(z_hits_tr.begin(),
                                             z_hits_tr.end()) + stod(options["z_extra"])));
                //DEBUG
                //cout<<" x_min = "<<x_min<<" x_max = "<<x_max<<" y_min = "<<y_min<<" y_max = "<<y_max<<" z_min = "<<z_min<<" z_max = "<<z_max<<endl;
                
                
                //CUT TRACKS due to exposure of camera
                double randcut = gRandom->Uniform(stod(options["exposure_time"])+readout_time);
                //randcut = 390.0;
                
                if (options["exposure_time_effect"] == "True") {
                    if (randcut<readout_time) {
                        
                        double y_cut_tmp = stod(options["y_dim"]) * (0.5 - randcut/readout_time)-3.;
                        
                        // Removing elements from x_hits_tr
                        x_hits_tr.erase(std::remove_if(x_hits_tr.begin(), x_hits_tr.end(), [&](const double& x) {
                            return y_hits_tr[&x-&*x_hits_tr.begin()] < y_cut_tmp;
                        }), x_hits_tr.end());
                        // Removing elements from z_hits_tr
                        z_hits_tr.erase(std::remove_if(z_hits_tr.begin(), z_hits_tr.end(), [&](const double& z) {
                            return y_hits_tr[&z-&*z_hits_tr.begin()] < y_cut_tmp;
                        }), z_hits_tr.end());
                        // Removing elements from energy_hits
                        energy_hits.erase(std::remove_if(energy_hits.begin(), energy_hits.end(), [&](const double& e) {
                            return y_hits_tr[&e-&*energy_hits.begin()] < y_cut_tmp;
                        }), energy_hits.end());
                        
                        // Removing elements from y_hits_tr [must be done after the previous ones]
                        y_hits_tr.erase(std::remove_if(y_hits_tr.begin(), y_hits_tr.end(),[&](const double& y) {
                            return y < y_cut_tmp;
                        }), y_hits_tr.end());
                        
                        row_cut = y_pix - (int)(randcut * static_cast<double>(y_pix) / readout_time);
                        
                        //DEBUG
                        //cout<<"y_cut_tmp = "<<y_cut_tmp<<endl;
                        //cout<<"row_cut = "<<row_cut<<endl;
                        //cout<<"sizes = ["<< x_hits_tr.size()<<","<<y_hits_tr.size()<<","
                        //    <<z_hits_tr.size()<<","<<energy_hits.size()<<"]"<<endl;
                        
                    } else if (randcut>stod(options["exposure_time"])) {
                        double y_cut_tmp = stod(options["y_dim"]) * (0.5 - (randcut - stod(options["exposure_time"])) / readout_time)+3.;
                        
                        // Removing elements from x_hits_tr
                        x_hits_tr.erase(std::remove_if(x_hits_tr.begin(), x_hits_tr.end(), [&](const double& x) {
                            return y_hits_tr[&x-&*x_hits_tr.begin()] > y_cut_tmp;
                        }), x_hits_tr.end());
                        // Removing elements from z_hits_tr
                        z_hits_tr.erase(std::remove_if(z_hits_tr.begin(), z_hits_tr.end(), [&](const double& z) {
                            return y_hits_tr[&z-&*z_hits_tr.begin()] > y_cut_tmp;
                        }), z_hits_tr.end());
                        // Removing elements from energy_hits
                        energy_hits.erase(std::remove_if(energy_hits.begin(), energy_hits.end(), [&](const double& e) {
                            return y_hits_tr[&e-&*energy_hits.begin()] > y_cut_tmp;
                        }), energy_hits.end());
                        
                        // Removing elements from y_hits_tr [must be done after the previous ones]
                        y_hits_tr.erase(std::remove_if(y_hits_tr.begin(), y_hits_tr.end(),[&](const double& y) {
                            return y > y_cut_tmp;
                        }), y_hits_tr.end());
                        
                        row_cut = y_pix - (int)((randcut-stod(options["exposure_time"])) * static_cast<double>(y_pix) / readout_time);
                        
                        //DEBUG
                        //cout<<"y_cut_tmp = "<<y_cut_tmp<<endl;
                        //cout<<"row_cut = "<<row_cut<<endl;
                        //cout<<"sizes = ["<< x_hits_tr.size()<<","<<y_hits_tr.size()<<","
                        //    <<z_hits_tr.size()<<","<<energy_hits.size()<<"]"<<endl;
                        
                    }
                    if(x_hits_tr.size()==0){
                        cut_energy = 0;
                        cout<<"The track was completely cut"<<endl;
                        TH2I final_image(Form("pic_run%d_ev%d", run_count, entry), "",
                                         x_pix, -0.5, x_pix -0.5,
                                         y_pix, -0.5, y_pix -0.5);
                        
                        for(unsigned int xx =0; xx < background.size(); xx++) {
                            for(unsigned int yy =0; yy < background[0].size(); yy++) {
                                final_image.SetBinContent(xx+1, yy+1, background[xx][yy]);
                            }
                        }
                        
                        outtree->Fill();
                        outfile->cd();
                        final_image.Write();
                        
                        continue;
                    }
                }
                
                vector<vector<double>> array2d_Nph(x_pix,
                                                   vector<double>(y_pix, 0.0));
                
                auto ta = std::chrono::steady_clock::now();
                // with saturation
                if(options["saturation"]=="True") {
                    cout<<"Starting compute_cmos_with_saturation with size = "<<x_hits_tr.size()<<"..."<<endl;
                    compute_cmos_with_saturation(x_hits_tr,
                                                 y_hits_tr,
                                                 z_hits_tr,
                                                 energy_hits,
                                                 options,
                                                 array2d_Nph,
                                                 energy,
                                                 NR_flag
                                                 );
                } else {// no saturation
                    compute_cmos_without_saturation(x_hits_tr,
                                                    y_hits_tr,
                                                    z_hits_tr,
                                                    energy_hits,
                                                    options,
                                                    array2d_Nph
                                                    );
                }
                
                // Integral of the track - if opt.exposure_effect, it's computed anyway after the cut on the original hits (to save time we digitize only the part that will be visible)
                N_photons = accumulate(array2d_Nph.cbegin(), array2d_Nph.cend(), 0, [](auto sum, const auto& row) {
                    return accumulate(row.cbegin(), row.cend(), sum);
                });
                // DEBUG
                cout<<"N_photons = "<<N_photons<<endl;
                
                if(options["Vignetting"]=="True") {
                    TrackVignetting(array2d_Nph,
                                    x_pix,
                                    y_pix,
                                    VignMap);
                }
                
                
                if (options["exposure_time_effect"]=="True") { //cut the track post-smearing
                    if (randcut<readout_time) {
                        for(unsigned int xx=0; xx < array2d_Nph.size(); xx++) {
                            for(int yy=0; yy < row_cut; yy++) {
                                array2d_Nph[xx][yy] = 0.0;
                            }
                        }
                    } else if(randcut> stod(options["exposure_time"]) ) {
                        for(unsigned int xx=0; xx < array2d_Nph.size(); xx++) {
                            for(int yy=row_cut; yy < (int)array2d_Nph[0].size(); yy++) {
                                array2d_Nph[xx][yy] = 0.0;
                            }
                        }
                    }
                    
                    //integral after camera exposure cut
                    N_photons_cut = accumulate(array2d_Nph.cbegin(),
                                               array2d_Nph.cend(), 0, [](auto sum, const auto& row) {
                        return accumulate(row.cbegin(), row.cend(), sum);
                    });
                }
                
                TH2I final_image(Form("pic_run%d_ev%d", run_count, entry), "",
                                 x_pix, -0.5, x_pix -0.5,
                                 y_pix, -0.5, y_pix -0.5);
                
                for(unsigned int xx =0; xx < array2d_Nph.size(); xx++) {
                    for(unsigned int yy =0; yy < array2d_Nph[0].size(); yy++) {
                        
                        int binc = background[xx][yy]+(int)array2d_Nph[xx][array2d_Nph[0].size()-1-yy];
                        final_image.SetBinContent(xx+1, yy+1, binc);
                    }
                }
                
                //Cut again the hits to save the effective length and energy which is visible in the final image,
                //and compute the number of photons post-cut
                if (options["exposure_time_effect"]=="True") {
                    if (randcut<readout_time) {
                        double y_cut_tmp = stod(options["y_dim"]) * (0.5 - randcut/readout_time);

                        // Removing elements from x_hits_tr
                        x_hits_tr.erase(std::remove_if(x_hits_tr.begin(), x_hits_tr.end(), [&](const double& x) {
                            return y_hits_tr[&x-&*x_hits_tr.begin()] < y_cut_tmp;
                        }), x_hits_tr.end());
                        // Removing elements from z_hits_tr
                        z_hits_tr.erase(std::remove_if(z_hits_tr.begin(), z_hits_tr.end(), [&](const double& z) {
                            return y_hits_tr[&z-&*z_hits_tr.begin()] < y_cut_tmp;
                        }), z_hits_tr.end());
                        // Removing elements from energy_hits
                        energy_hits.erase(std::remove_if(energy_hits.begin(), energy_hits.end(), [&](const double& e) {
                            return y_hits_tr[&e-&*energy_hits.begin()] < y_cut_tmp;
                        }), energy_hits.end());
                        
                        // Removing elements from y_hits_tr [must be done after the previous ones]
                        y_hits_tr.erase(std::remove_if(y_hits_tr.begin(), y_hits_tr.end(),[&](const double& y) {
                            return y < y_cut_tmp;
                        }), y_hits_tr.end());
                        
                    } else if (randcut> stod(options["exposure_time"]) ) {
                        double y_cut_tmp = stod(options["y_dim"]) * (0.5 - (randcut - stod(options["exposure_time"])) / readout_time);

                        // Removing elements from x_hits_tr
                        x_hits_tr.erase(std::remove_if(x_hits_tr.begin(), x_hits_tr.end(), [&](const double& x) {
                            return y_hits_tr[&x-&*x_hits_tr.begin()] > y_cut_tmp;
                        }), x_hits_tr.end());
                        // Removing elements from z_hits_tr
                        z_hits_tr.erase(std::remove_if(z_hits_tr.begin(), z_hits_tr.end(), [&](const double& z) {
                            return y_hits_tr[&z-&*z_hits_tr.begin()] > y_cut_tmp;
                        }), z_hits_tr.end());
                        // Removing elements from energy_hits
                        energy_hits.erase(std::remove_if(energy_hits.begin(), energy_hits.end(), [&](const double& e) {
                            return y_hits_tr[&e-&*energy_hits.begin()] > y_cut_tmp;
                        }), energy_hits.end());
                        
                        // Removing elements from y_hits_tr [must be done after the previous ones]
                        y_hits_tr.erase(std::remove_if(y_hits_tr.begin(), y_hits_tr.end(),[&](const double& y) {
                            return y > y_cut_tmp;
                        }), y_hits_tr.end());
                        
                    }
                    
                    cut_energy = accumulate(energy_hits.begin(), energy_hits.end(), 0.0);
                    
                    if(x_hits_tr.size()==0){
                        cut_energy = 0;
                        cout<<"The track was completely cut after smearing"<<endl;
                        
                        TH2I final_image_cut(Form("pic_run%d_ev%d", run_count, entry), "",
                                         x_pix, -0.5, x_pix -0.5,
                                         y_pix, -0.5, y_pix -0.5);
                        
                        for(unsigned int xx =0; xx < background.size(); xx++) {
                            for(unsigned int yy =0; yy < background[0].size(); yy++) {
                                final_image_cut.SetBinContent(xx+1, yy+1, background[xx][yy]);
                            }
                        }
                        
                        outtree->Fill();
                        outfile->cd();
                        final_image_cut.Write();
                        
                        continue;
                    }
                    
                }
                
                
                // Compute variables to be saved in the tree
                proj_track_2D_cut = 0;
                for(int ihit=0; ihit < numhits-1; ihit++){
                    proj_track_2D_cut += sqrt((x_hits_tr[ihit+1]-x_hits_tr[ihit])*(x_hits_tr[ihit+1]-x_hits_tr[ihit])+
                                              (y_hits_tr[ihit+1]-y_hits_tr[ihit])*(y_hits_tr[ihit+1]-y_hits_tr[ihit])
                                             );
                }
                // DEBUG
                //cout<<"proj_track_2D_cut = "<<Form("%.10f", proj_track_2D_cut)<<endl;
                
                
                x_min_cut = (*min_element(x_hits_tr.begin(), x_hits_tr.end()) + 0.5 * stod(options["x_dim"])) * static_cast<double>(x_pix) / stod(options["x_dim"]);
                x_max_cut = (*max_element(x_hits_tr.begin(), x_hits_tr.end()) + 0.5 * stod(options["x_dim"])) * static_cast<double>(x_pix) / stod(options["x_dim"]);
                y_min_cut = (*min_element(y_hits_tr.begin(), y_hits_tr.end()) + 0.5 * stod(options["y_dim"])) * static_cast<double>(y_pix) / stod(options["y_dim"]);
                y_max_cut = (*max_element(y_hits_tr.begin(), y_hits_tr.end()) + 0.5 * stod(options["y_dim"])) * static_cast<double>(y_pix) / stod(options["y_dim"]);
                z_min_cut = min((*max_element(z_hits_tr.begin(),
                                                 z_hits_tr.end()) + stod(options["z_extra"])),
                                (*min_element(z_hits_tr.begin(),
                                                 z_hits_tr.end()) + stod(options["z_extra"])));
                z_max_cut = max((*max_element(z_hits_tr.begin(),
                                                 z_hits_tr.end()) + stod(options["z_extra"])),
                                (*min_element(z_hits_tr.begin(),
                                                 z_hits_tr.end()) + stod(options["z_extra"])));
                //DEBUG
                //cout<<" x_min_cut = "<<x_min_cut<<" x_max_cut = "<<x_max_cut<<" y_min_cut = "<<y_minv<<" y_max_cut = "<<y_max_cut<<" z_min_cut = "<<z_min_cut<<" z_max_cut = "<<z_max_cut<<endl;
                
                
                // if there are at least 2 hits compute theta and phi
                if(x_hits_tr.size() > 1) {
                    phi   = atan2( y_hits_tr[1] - y_hits_tr[0],
                                   x_hits_tr[1] - x_hits_tr[0]
                                  );
                    theta = acos((z_hits_tr[1] - z_hits_tr[0]) /
                                 sqrt(
                                     (x_hits_tr[1] - x_hits_tr[0])*(x_hits_tr[1] - x_hits_tr[0])+
                                     (y_hits_tr[1] - y_hits_tr[0])*(y_hits_tr[1] - y_hits_tr[0])+
                                     (z_hits_tr[1] - z_hits_tr[0])*(z_hits_tr[1] - z_hits_tr[0]))
                                 );
                } else {
                    phi   = -999;
                    theta = -999;
                }
                
                if(options["SRIM"]=="False") {
                    track_length_3D = accumulate(tracklen_hits->begin(), tracklen_hits->end(), 0.0);
                    px              = (*px_particle)[0];
                    py              = (*py_particle)[0];
                    pz              = (*pz_particle)[0];
                }
                
                auto tb = std::chrono::steady_clock::now();
                std::chrono::duration<double> durtmp=tb-ta;
                cout << "Time taken in seconds to compute_cmos_with_saturation is: " << durtmp.count() << endl;
                
                
                outtree->Fill();
                outfile->cd();
                
                // DEBUG
                final_image.Write();
                
            }
            //outfile->cd("event_info");
            outtree->Write();
            
            cout<<Form("COMPLETED RUN %d",run_count)<<endl;
            
            f->Close();
            outfile->Close();
            
            run_count ++;*/


        
        ///// bbbbbbbbbbb



        

        
    }
}
