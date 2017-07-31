#include "WireCellSio/CelltreeSource.h"
#include "WireCellIface/SimpleTrace.h"
#include "WireCellIface/SimpleFrame.h"

#include "WireCellUtil/NamedFactory.h"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1F.h"


WIRECELL_FACTORY(CelltreeSource, WireCell::Sio::CelltreeSource, WireCell::IFrameSource, WireCell::IConfigurable);

using namespace WireCell;

Sio::CelltreeSource::CelltreeSource()
{
}

Sio::CelltreeSource::~CelltreeSource()
{
}

void Sio::CelltreeSource::configure(const WireCell::Configuration& cfg)
{
    if (cfg["filename"].empty()) {
        THROW(ValueError() << errmsg{"MagnifySource: must supply input \"filename\" configuration item."});
    }
    m_cfg = cfg;
}

WireCell::Configuration Sio::CelltreeSource::default_configuration() const
{
    Configuration cfg;
    // Give a URL for the input file.
    cfg["filename"] = "";

    cfg["EventNo"] = 0;
    
    // Give a list of frame tags.  These are translated to histogram
    cfg["frames"][0] = "orig";

    return cfg;
}

bool Sio::CelltreeSource::operator()(IFrame::pointer& out)
{
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

  int frame_number = m_cfg["EventNo"].asDouble();
  
  int siz = tree->GetEntry(frame_number);

  // meed run number and subrunnumber
  int frame_ident = event_no;
  double frame_time=0;

  ITrace::vector all_traces;
  std::unordered_map<IFrame::tag_t, IFrame::trace_list_t> tagged_traces;
  WireCell::Waveform::ChannelMaskMap cmm;

  int nticks=0;
  TH1F* signal = dynamic_cast<TH1F*>(esignal->At(0));
  nticks = signal->GetNbinsX();
  
  
  if (siz>0 && frame_number < siz){
    int nchannels = channelid->size();

    auto frametag = m_cfg["frames"][0].asString();
    int channel_number = 0;
    
    // fill waveform ... 
     for (int ind=0; ind < nchannels; ++ind) {
	TH1F* signal = dynamic_cast<TH1F*>(esignal->At(ind));
	if (!signal) continue;
	channel_number = ind;

	 ITrace::ChargeSequence charges;
	 for (int itickbin = 0; itickbin != signal->GetNbinsX(); itickbin++){
	   charges.push_back(signal->GetBinContent(itickbin+1));
	 }
	 
	 const size_t index = all_traces.size();
	 tagged_traces[frametag].push_back(index);
	 
	 all_traces.push_back(std::make_shared<SimpleTrace>(channel_number, 0, charges));
     }
     auto sframe = new SimpleFrame(frame_ident, frame_time,
				   all_traces, 0.5*units::microsecond, cmm);
     for (auto const& it : tagged_traces) {
       sframe->tag_traces(it.first, it.second);
       std::cerr << "MagnifySource: tag " << it.second.size() << " traces as: \"" << it.first << "\"\n";
     }
     
     out = IFrame::pointer(sframe);
     
    return true;
  }else{
    std::cerr << "Frame Number is out of boundary! " << std::endl;
    return false;
  }
  

  

  
  
  // // runNo, subRunNo??
  // trun->SetBranchAddress("eventNo", &frame_ident);
  // trun->SetBranchAddress("total_time_bin", &nticks);
  // //trun->SetBranchAddress("time_offset", &frame_time); use this??
  // trun->GetEntry(0);
  
  //  std::cerr << "CelltreeSource: frame ident="<<frame_ident<<", time=0, nticks="<<nticks<<"\n";
  
    
}
