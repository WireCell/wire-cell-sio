#include "WireCellSio/XinFileSink.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TFile.h"
#include "TTree.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(XinFileSink, WireCell::Sio::XinFileSink, WireCell::IFrameSink, WireCell::IConfigurable);

using namespace WireCell;

Sio::XinFileSink::XinFileSink()
{
}

Sio::XinFileSink::~XinFileSink()
{
}

void Sio::XinFileSink::configure(const WireCell::Configuration& config)
{
    m_cfg = config;
}

WireCell::Configuration Sio::XinFileSink::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = "";
    cfg["rebin"] = 1;
    return cfg;
}

void save_into_file(const char* filename,IFrame::pointer frame_decon, int nrebin);

bool Sio::XinFileSink::operator()(const IFrame::pointer& frame)
{
    std::string fname = m_cfg["filename"].asString();
    if (fname.empty()) { return false; }
    int nrebin = m_cfg["rebin"].asInt();
    save_into_file(fname.c_str(), frame, nrebin);
    return true;
}




void save_into_file(const char* filename,IFrame::pointer frame_decon, int nrebin){
    TFile *file1 = new TFile(filename);

    TFile *file = new TFile("temp.root","RECREATE");
    TTree *Trun = ((TTree*)file1->Get("Trun"))->CloneTree();
    Trun->SetDirectory(file);

    TH2I *hu_orig = (TH2I*)file1->Get("hu_orig");
    TH2I *hv_orig = (TH2I*)file1->Get("hv_orig");
    TH2I *hw_orig = (TH2I*)file1->Get("hw_orig");

    hu_orig->SetDirectory(file);
    hv_orig->SetDirectory(file);
    hw_orig->SetDirectory(file);
  
    int nwire_u = hu_orig->GetNbinsX();
    int nwire_v = hv_orig->GetNbinsX();
    int nwire_w = hw_orig->GetNbinsX();
    int nticks = hu_orig->GetNbinsY();

  
    TH2F *hu_raw = (TH2F*)file1->Get("hu_raw");
    TH2F *hv_raw = (TH2F*)file1->Get("hv_raw");
    TH2F *hw_raw = (TH2F*)file1->Get("hw_raw");

    hu_raw->SetDirectory(file);
    hv_raw->SetDirectory(file);
    hw_raw->SetDirectory(file);

  
    TH1F *hu_baseline = (TH1F*)file1->Get("hu_baseline");
    TH1F *hv_baseline = (TH1F*)file1->Get("hv_baseline");
    TH1F *hw_baseline = (TH1F*)file1->Get("hw_baseline");
  
    hu_baseline->SetDirectory(file);
    hv_baseline->SetDirectory(file);
    hw_baseline->SetDirectory(file);

  
    // temporary ...  need a data structure to load the threshold ... 
    TH1F *hu_threshold = (TH1F*)file1->Get("hu_threshold");
    TH1F *hv_threshold = (TH1F*)file1->Get("hv_threshold");
    TH1F *hw_threshold = (TH1F*)file1->Get("hw_threshold");
  
    hu_threshold->SetDirectory(file);
    hv_threshold->SetDirectory(file);
    hw_threshold->SetDirectory(file);
  


  
    // temporary ...

  
    TH2F *hu_decon = new TH2F("hu_decon","hu_decon",nwire_u,-0.5,nwire_u-0.5,int(nticks/nrebin),0,nticks);
    TH2F *hv_decon = new TH2F("hv_decon","hv_decon",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,int(nticks/nrebin),0,nticks);
    TH2F *hw_decon = new TH2F("hw_decon","hw_decon",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,int(nticks/nrebin),0,nticks);

  
    auto traces = frame_decon->traces();
    for (auto trace : *traces.get()) {
        int tbin = trace->tbin();
        int ch = trace->channel();
        auto charges = trace->charge();
        if (ch < nwire_u){
            int counter = 0;
            int rebin_counter = 0;
            float acc_charge = 0;
      
            for (auto q : charges) {
                if (rebin_counter < nrebin){
                    acc_charge += q;
                    rebin_counter ++;
                }
                if (rebin_counter == nrebin){
                    counter ++;
                    hu_decon->SetBinContent(ch+1,tbin+counter,acc_charge); 
                    //reset ... 
                    rebin_counter = 0;
                    acc_charge = 0;
                }
            }
        }else if (ch < nwire_v + nwire_u){

            int counter = 0;
            int rebin_counter = 0;
            float acc_charge = 0;
      
            for (auto q : charges) {
	
                if (rebin_counter < nrebin){
                    acc_charge += q;
                    rebin_counter ++;
                }
                if (rebin_counter == nrebin){
                    counter ++;
                    hv_decon->SetBinContent(ch+1-nwire_u,tbin+counter,acc_charge); 
                    //reset ... 
                    rebin_counter = 0;
                    acc_charge = 0;
                }
            }

      
            // int counter = 0;
            // for (auto q : charges) {
            // 	counter ++;
            // 	hv_decon->SetBinContent(ch+1-nwire_u,tbin+counter,q); 
	
            // }
        }else{

            int counter = 0;
            int rebin_counter = 0;
            float acc_charge = 0;
      
            for (auto q : charges) {
	
                if (rebin_counter < nrebin){
                    acc_charge += q;
                    rebin_counter ++;
                }
                if (rebin_counter == nrebin){
                    counter ++;
                    hw_decon->SetBinContent(ch+1-nwire_u-nwire_v,tbin+counter,acc_charge); 
                    //reset ... 
                    rebin_counter = 0;
                    acc_charge = 0;
                }
            }
      
            // int counter = 0;
            // for (auto q : charges) {
            // 	counter ++;
            // 	hw_decon->SetBinContent(ch+1-nwire_u-nwire_v,tbin+counter,q); 
	
            // }
        }
    }


  
    // save bad channels 
    TTree *T_bad = new TTree("T_bad","T_bad");
    int chid, plane, start_time,end_time;
    T_bad->Branch("chid",&chid,"chid/I");
    T_bad->Branch("plane",&plane,"plane/I");
    T_bad->Branch("start_time",&start_time,"start_time/I");
    T_bad->Branch("end_time",&end_time,"end_time/I");
    T_bad->SetDirectory(file);

    TTree *T_lf = new TTree("T_lf","T_lf");
    int channel;
    T_lf->Branch("channel",&channel,"channel/I");
  

    Waveform::ChannelMaskMap input_cmm = frame_decon->masks();
    for (auto const& it: input_cmm) {

        if (it.first == "bad"){ // save bad ... 
            //std::cout << "Xin1: " << it.first << " " << it.second.size() << std::endl;
            for (auto const &it1 : it.second){
                chid = it1.first;
                if (chid < nwire_u){
                    plane = 0;
                }else if (chid < nwire_v + nwire_u){
                    plane = 1;
                }else{
                    plane = 2;
                }
                //std::cout << "Xin1: " << chid << " " << plane << " " << it1.second.size() << std::endl;
                for (size_t ind = 0; ind < it1.second.size(); ++ind){
                    start_time = it1.second[ind].first;
                    end_time = it1.second[ind].second;
                    T_bad->Fill();
                }
            }
        }else if (it.first =="lf_noisy"){
            for (auto const &it1 : it.second){
                channel = it1.first;
                T_lf->Fill();
            }
      
        }else if (it.first=="threshold"){
            for (auto const &it1 : it.second){
                chid = it1.first;
                float threshold = it1.second[0].first/it1.second[0].second;
                if (chid < nwire_u){
                    hu_threshold->SetBinContent(chid+1,threshold*nrebin*3.0);
                }else if (chid < nwire_u+nwire_v){
                    hv_threshold->SetBinContent(chid+1-nwire_u,threshold*nrebin*3.0);
                }else{
                    hw_threshold->SetBinContent(chid+1-nwire_u-nwire_v,threshold*nrebin*3.0);
                }
	 
            }
        }

    
    }

    file->Write();
    file->Close();
  
}
