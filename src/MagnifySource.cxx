#include "WireCellSio/MagnifySource.h"
#include "WireCellIface/SimpleTrace.h"
#include "WireCellIface/SimpleFrame.h"

#include "WireCellUtil/NamedFactory.h"

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"


WIRECELL_FACTORY(MagnifySource, WireCell::Sio::MagnifySource, WireCell::IFrameSource, WireCell::IConfigurable);

using namespace WireCell;

Sio::MagnifySource::MagnifySource()
{
}

Sio::MagnifySource::~MagnifySource()
{
}

void Sio::MagnifySource::configure(const WireCell::Configuration& cfg)
{
    if (cfg["filename"].empty()) {
        THROW(ValueError() << errmsg{"MagnifySource: must supply input \"filename\" configuration item."});
    }
    m_cfg = cfg;
}

WireCell::Configuration Sio::MagnifySource::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = "";       // actually any URL
    cfg["histtype"] = "raw";
    return cfg;
}


bool Sio::MagnifySource::operator()(IFrame::pointer& out)
{
    std::string url = m_cfg["filename"].asString();
    std::string histtype = m_cfg["histtype"].asString();

    TFile* tfile = TFile::Open(url.c_str());
    TH2* hists[3];
    for (int iplane=0; iplane<3; ++iplane) {
	std::string name = Form("h%c_%s", 'u'+iplane, histtype.c_str());
        std::cerr << "MagnifySource: loading " << name << std::endl;
	hists[iplane] = (TH2*)tfile->Get(name.c_str());
    }

    int frame_ident=0;
    double frame_time=0;
    TTree *trun = (TTree*)tfile->Get("Trun");
    if (trun) {
        trun->SetBranchAddress("eventNo", &frame_ident);
        //trun->SetBranchAddress("time_offset", &frame_time);??
        trun->GetEntry(0);
    }
    else {
        std::cerr << "MagnifySource: \"run\" tree not in input, using frame ident=0, time=0\n";
    }


    WireCell::Waveform::ChannelMaskMap cmm;
    TTree *T_bad = (TTree*)tfile->Get("T_bad");
    if (T_bad) { 
        int chid=0, plane=0, start_time=0, end_time=0;
        T_bad->SetBranchAddress("chid",&chid);
        T_bad->SetBranchAddress("plane",&plane);
        T_bad->SetBranchAddress("start_time",&start_time);
        T_bad->SetBranchAddress("end_time",&end_time);
      
        for (int i=0;i!=T_bad->GetEntries();i++){
            T_bad->GetEntry(i);
            WireCell::Waveform::BinRange chirped_bins;
            chirped_bins.first = start_time;
            chirped_bins.second = end_time;
            cmm["bad"][chid].push_back(chirped_bins);
        }
    }
    else {
        std::cerr << "MagnifySource: \"bad\" tree not in input, not setting \"bad\" channel mask map\n";
    }
      
    TTree *T_lf = (TTree*)tfile->Get("T_lf");
    if (T_lf) {
        int channel=0;
        T_lf->SetBranchAddress("channel",&channel);
        for (int i=0;i!=T_lf->GetEntries();i++){
            T_lf->GetEntry(i);
            WireCell::Waveform::BinRange chirped_bins;
            chirped_bins.first = 0;
            chirped_bins.second = hists[0]->GetNbinsY();
            cmm["lf_noisy"][channel].push_back(chirped_bins);
        }
    }
    else {
        std::cerr << "MagnifySource:  \"lf\" tree not in input, not setting \"lf_noisy\" channel mask map\n";
    }

    ITrace::vector traces;
    int channel_number = 0;     // fixme: this likely breaks for non-microboone
    for (int iplane=0; iplane<3; ++iplane) {
        TH2* h = hists[iplane];

        int nchannels = h->GetNbinsX();
        int nticks = h->GetNbinsY();

        double qtot = 0;
        for (int ichbin = 1; ichbin <= nchannels; ++ichbin) {
            ITrace::ChargeSequence charges;
            for (int itickbin = 1; itickbin <= nticks; ++itickbin) {
                auto q = h->GetBinContent(ichbin, itickbin);
                charges.push_back(q);
                qtot += q;
            }
            traces.push_back(std::make_shared<SimpleTrace>(channel_number, 0, charges));
            ++channel_number;
        }
        std::cerr << "MagnifySource: plane " << iplane
                  << ": " << nchannels << " X " << nticks
                  << " qtot=" << qtot
                  << std::endl;
    }

    out = std::make_shared<SimpleFrame>(frame_ident, frame_time,
                                        traces, 0.5*units::microsecond, cmm);
    return true;
}


