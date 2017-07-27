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
    // Give a URL for the input file.
    cfg["filename"] = "";

    // Give a list of frame tags.  These are translated to histogram
    // names h[uvw]_<tag> to be used for input.  Missing histograms
    // raise an exception.  Each frame found will be loaded as tagged
    // traces into the output frame.  The ensemble of traces will be
    // tagged by both <tag> and the per-plane subset as
    // <planeletter>plane.
    cfg["frames"][0] = "raw";
    return cfg;
}


bool Sio::MagnifySource::operator()(IFrame::pointer& out)
{
    std::string url = m_cfg["filename"].asString();

    TFile* tfile = TFile::Open(url.c_str());

    int frame_ident=0;
    int nticks=0;
    double frame_time=0;
    {
        TTree *trun = (TTree*)tfile->Get("Trun");
        if (!trun) {
            THROW(IOError() << errmsg{"No tree: Trun in input file"});
        }

        // runNo, subRunNo??
        trun->SetBranchAddress("eventNo", &frame_ident);
        trun->SetBranchAddress("total_time_bin", &nticks);
        //trun->SetBranchAddress("time_offset", &frame_time); use this??
        trun->GetEntry(0);
    }
    std::cerr << "MagnifySource: frame ident="<<frame_ident<<", time=0, nticks="<<nticks<<"\n";


    WireCell::Waveform::ChannelMaskMap cmm;
    {
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
    }
      
    { 
        TTree *T_lf = (TTree*)tfile->Get("T_lf");
        if (T_lf) {             
            int channel=0;
            T_lf->SetBranchAddress("channel",&channel);
            const int nentries = T_lf->GetEntries();
            for (int ind = 0; ind != nentries; ++ind){
                T_lf->GetEntry(ind);
                WireCell::Waveform::BinRange all_bins;
                all_bins.first = 0;
                all_bins.second = nticks-1;
                cmm["lf_noisy"][channel].push_back(all_bins);
            }
            std::cerr << "MagnifySource:  "<<nentries<<" \"lf_noisy\" noisy channels\n";
        }
        else {
            std::cerr << "MagnifySource:  \"lf\" tree not in input, not setting \"lf_noisy\" channel mask map\n";
        }
    }


    ITrace::vector all_traces;
    std::unordered_map<IFrame::tag_t, IFrame::trace_list_t> tagged_traces;

    {
        for (auto jtag : m_cfg["frames"]) {
            auto frametag = jtag.asString();
            int channel_number = 0;
            for (int iplane=0; iplane<3; ++iplane) {
                std::string plane_name = Form("%cplane", 'u'+iplane);
                std::string hist_name = Form("h%c_%s", 'u'+iplane, frametag.c_str());
                std::cerr << "MagnifySource: loading " << hist_name << std::endl;
                TH2* hist = (TH2*)tfile->Get(hist_name.c_str());
            
                ITrace::vector plane_traces;
            
                int nchannels = hist->GetNbinsX();
                int nticks = hist->GetNbinsY();
                double qtot = 0;

                // loop over channels in plane
                for (int chip = 0; chip<nchannels; ++chip) {
                    const int ichbin = chip+1;

                    ITrace::ChargeSequence charges;
                    for (int itickbin = 1; itickbin <= nticks; ++itickbin) {
                        auto q = hist->GetBinContent(ichbin, itickbin);
                        charges.push_back(q);
                        qtot += q;
                    }
                    const size_t index = all_traces.size();
                    tagged_traces[frametag].push_back(index);
                    all_traces.push_back(std::make_shared<SimpleTrace>(channel_number, 0, charges));

                    ++channel_number;
                }
                std::cerr << "MagnifySource: plane " << iplane
                          << ": " << nchannels << " X " << nticks
                          << " qtot=" << qtot
                          << std::endl;
            } // plane

        }
    }

    auto sframe = new SimpleFrame(frame_ident, frame_time,
                                  all_traces, 0.5*units::microsecond, cmm);
    for (auto const& it : tagged_traces) {
        sframe->tag_traces(it.first, it.second);
        std::cerr << "MagnifySource: tag " << it.second.size() << " traces as: \"" << it.first << "\"\n";
    }

    out = IFrame::pointer(sframe);
    return true;
}


