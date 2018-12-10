#include "WireCellSio/CelltreeSource.h"
#include "WireCellIface/SimpleTrace.h"
#include "WireCellIface/SimpleFrame.h"

#include "WireCellUtil/NamedFactory.h"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1F.h"


WIRECELL_FACTORY(CelltreeSource, WireCell::Sio::CelltreeSource,
                 WireCell::IFrameSource, WireCell::IConfigurable)

using namespace WireCell;

Sio::CelltreeSource::CelltreeSource()
    : m_calls(0)
{
}

Sio::CelltreeSource::~CelltreeSource()
{
}

void Sio::CelltreeSource::configure(const WireCell::Configuration& cfg)
{
    if (cfg["filename"].empty()) {
        THROW(ValueError() << errmsg{"CelltreeSource: must supply input \"filename\" configuration item."});
    }
    m_cfg = cfg;
}

WireCell::Configuration Sio::CelltreeSource::default_configuration() const
{
    Configuration cfg;
    // Give a URL for the input file.
    cfg["filename"] = "";

    // which event in the celltree file to be processed 
    cfg["EventNo"] = "0";
    
    // Give a list of frame/tree tags. 

    // just raw waveform and no other choice at present
    // Tree: Sim, Wf: raw_wf
    cfg["frames"][0] = "orig";


    return cfg;
}

bool Sio::CelltreeSource::operator()(IFrame::pointer& out)
{
    out = nullptr;
    ++m_calls;
    std::cerr << "CelltreeSource: called " << m_calls << " times\n";
    if (m_calls > 2) {
        std::cerr << "CelltreeSource: past EOS\n";
        return false;
    }
    if (m_calls > 1) {
        std::cerr << "CelltreeSource: EOS\n";
        return true;            // this is to send out EOS
    }
    

  std::string url = m_cfg["filename"].asString();

  TFile* tfile = TFile::Open(url.c_str());


  
  TTree *tree = (TTree*)tfile->Get("/Event/Sim");
  if (!tree) {
    THROW(IOError() << errmsg{"No tree: /Event/Sim in input file"});
  }

  tree->SetBranchStatus("*",0);

  int run_no, subrun_no, event_no;
  tree->SetBranchStatus("eventNo",1);
  tree->SetBranchAddress("eventNo" , &event_no);
  tree->SetBranchStatus("runNo",1);
  tree->SetBranchAddress("runNo"   , &run_no);
  tree->SetBranchStatus("subRunNo",1);
  tree->SetBranchAddress("subRunNo", &subrun_no);
  

  std::vector<int> *channelid = new std::vector<int>;
  TClonesArray* esignal = new TClonesArray;
  
  tree->SetBranchStatus("raw_channelId",1);
  tree->SetBranchAddress("raw_channelId", &channelid);
  tree->SetBranchStatus("raw_wf",1);
  tree->SetBranchAddress("raw_wf", &esignal);

  int frame_number = std::stoi(m_cfg["EventNo"].asString());

  // loop entry 
  int siz = 0;
  unsigned int entries = tree->GetEntries();
  for(unsigned int ent = 0; ent<entries; ent++)
  {
      siz = tree->GetEntry(ent);
      if(event_no == frame_number)
      {
          break;
      }
  }
    

  // need run number and subrunnumber
  int frame_ident = event_no;
  double frame_time=0;

  // some output using eventNo, runNo, subRunNO, ...
  std::cerr<<"CelltreeSource: frame "<<frame_number<<"\n";
  std::cerr << "CelltreeSource: runNo "<<run_no<<", subrunNo "<<subrun_no<<", eventNo "<<event_no<<"\n";
  
  ITrace::vector all_traces;
  std::unordered_map<IFrame::tag_t, IFrame::trace_list_t> tagged_traces;
  // celltree input now is raw data, no information about any noisy or bad channels
  // leave cmm empty.
  WireCell::Waveform::ChannelMaskMap cmm;

  
  if (siz>0 && frame_number == frame_ident){
    int nchannels = channelid->size();

    auto frametag = m_cfg["frames"][0].asString();
    int channel_number = 0;
    
	std::cerr<<"CelltreeSource: loading "<<frametag<<" "<<nchannels<<" channels \n";

    // fill waveform ... 
     for (int ind=0; ind < nchannels; ++ind) {
	TH1F* signal = dynamic_cast<TH1F*>(esignal->At(ind));
	if (!signal) continue;
    channel_number = channelid->at(ind);

	 ITrace::ChargeSequence charges;
     int nticks = signal->GetNbinsX();
     //std::cerr<<"CelltreeSource: tick "<<nticks<<"\n";
     //nticks = 9600,  this could be an issue cause just 9594 have non-ZERO value around baseline
     for (int itickbin = 0; itickbin != nticks; itickbin++){
        if(signal->GetBinContent(itickbin+1)!=0) {
            charges.push_back(signal->GetBinContent(itickbin+1));
        }
	 }

	 const size_t index = all_traces.size();
	 tagged_traces[frametag].push_back(index);
	 
	 all_traces.push_back(std::make_shared<SimpleTrace>(channel_number, 0, charges));
     }
     auto sframe = new SimpleFrame(frame_ident, frame_time,
				   all_traces, 0.5*units::microsecond, cmm);
     for (auto const& it : tagged_traces) {
       sframe->tag_traces(it.first, it.second);
       std::cerr << "CelltreeSource: tag " << it.second.size() << " traces as: \"" << it.first << "\"\n";
     }
     
     out = IFrame::pointer(sframe);
     
     return true;
  }else{
    std::cerr << "CelltreeSource: Event Number is out of boundary! " << std::endl;
    return false;
  }
  

  

  
  
  // // runNo, subRunNo??
  // trun->SetBranchAddress("eventNo", &frame_ident);
  // trun->SetBranchAddress("total_time_bin", &nticks);
  // //trun->SetBranchAddress("time_offset", &frame_time); use this??
  // trun->GetEntry(0);
  
  //  std::cerr << "CelltreeSource: frame ident="<<frame_ident<<", time=0, nticks="<<nticks<<"\n";
  
    
}
