#include "WireCellSio/MagnifySink.h"
#include "WireCellIface/ITrace.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TFile.h"
#include "TTree.h"

#include "WireCellUtil/NamedFactory.h"

#include <vector>
#include <string>

WIRECELL_FACTORY(MagnifySink, WireCell::Sio::MagnifySink, WireCell::IFrameSink, WireCell::IConfigurable);

using namespace WireCell;


Sio::MagnifySink::MagnifySink()
{
}

Sio::MagnifySink::~MagnifySink()
{
}

void Sio::MagnifySink::configure(const WireCell::Configuration& cfg)
{
    std::string fn;

    fn = cfg["input_filename"].asString();
    if (fn.empty()) {
        if (cfg["shunk"].empty()) {
            std::cerr << "MagnifySink no objects to copy but given input: " << fn << std::endl;
        }
        else {
            THROW(ValueError() << errmsg{"MagnifySink: must provide input filename to shunt objects to output"});
        }
    }

    fn = cfg["output_filename"].asString();
    if (fn.empty()) {
        THROW(ValueError() << errmsg{"Must provide output filename to MagnifySink"});
    }

    auto anode_tn = get<std::string>(cfg, "anode", "AnodePlane");
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

    m_cfg = cfg;
}

WireCell::Configuration Sio::MagnifySink::default_configuration() const
{
    Configuration cfg;

    cfg["anode"] = "AnodePlane";

    cfg["input_filename"] = ""; // fixme: this TOTALLY violates the design of wire cell DFP

    cfg["shunt"] = Json::arrayValue;

    cfg["output_filename"] = "";
    cfg["frames"] = Json::arrayValue;

    return cfg;
}

typedef std::unordered_set<std::string> string_set_t;
string_set_t getset(const WireCell::Configuration& cfg)
{
    string_set_t ret;
    for (auto jone : cfg) {
        ret.insert(jone.asString());
    }
    return ret;
}

// fixme: this little helper needs to move to FrameUtil
ITrace::vector get_tagged_traces(IFrame::pointer frame, IFrame::tag_t tag)
{
    ITrace::vector ret;
    auto const& all_traces = frame->traces();
    for (size_t index : frame->tagged_traces(tag)) {
        ret.push_back(all_traces->at(index));
    }
    return ret;
}


std::vector<WireCell::Binning> collate_byplane(const ITrace::vector traces, const IAnodePlane::pointer anode,
                                               ITrace::vector byplane[])
{
    std::vector<int> uvwt[4];
    for (auto trace : traces) {
        const int chid = trace->channel();
        auto wpid = anode->resolve(chid);
        const int iplane = wpid.index();
        if (iplane<0 || iplane>=3) {
            THROW(RuntimeError() << errmsg{"Illegal wpid"});
        }
        uvwt[iplane].push_back(chid);
        byplane[iplane].push_back(trace);
        uvwt[3].push_back(trace->tbin());
        uvwt[3].push_back(trace->tbin() + trace->charge().size());
    }

    std::vector<Binning> binnings;
    for (int ind=0; ind<4; ++ind) {
        auto const& one = uvwt[ind];
        auto mme = std::minmax_element(one.begin(), one.end());
        const int vmin = *mme.first;
        const int vmax = *mme.second;
        if (ind == 3) {
            const int n = vmax - vmin;
            binnings.push_back(Binning(n, vmin, vmax));
        }
        else {
            // Channel-centered binning
            const double diff = vmax - vmin;
            binnings.push_back(Binning(diff+1, vmin-0.5, vmax+0.5));
        } 
    }
    return binnings;
}


bool Sio::MagnifySink::operator()(const IFrame::pointer& frame)
{
    if (!frame) {
        // eos 
        return true;
    }

    std::string ofname = m_cfg["output_filename"].asString();
    std::cerr << "MagnifySink: opening for output: " << ofname << std::endl;
    TFile* output_tf = TFile::Open(ofname.c_str(), "RECREATE");

    for (auto tag : getset(m_cfg["frames"])) {

        ITrace::vector traces_byplane[3], traces = get_tagged_traces(frame, tag);
        if (traces.empty()) {
            std::cerr << "MagnifySink: no tagged traces for \"" << tag << "\"\n";
            THROW(ValueError() << errmsg{"MagnifySink: no tagged traces"});
        }
        auto binnings = collate_byplane(traces, m_anode, traces_byplane);

        Binning tbin = binnings[3];
        for (int iplane=0; iplane < 3; ++iplane) {

            const std::string name = Form("h%c_%s", 'u'+iplane, tag.c_str());
            Binning cbin = binnings[iplane];
            std::cerr << "MagnifySink:"
                      << " cbin:"<<cbin.nbins()<<"["<<cbin.min() << "," << cbin.max() << "]"
                      << " tbin:"<<tbin.nbins()<<"["<<tbin.min() << "," << tbin.max() << "]\n";

            TH2F* hist = new TH2F(name.c_str(), name.c_str(),
                                  cbin.nbins(), cbin.min(), cbin.max(),
                                  tbin.nbins(), tbin.min(), tbin.max());

            hist->SetDirectory(output_tf);
            
            for (auto trace : traces_byplane[iplane]) {
                const int tbin = trace->tbin();
                const int ch = trace->channel();
                auto const& charges = trace->charge();
                for (size_t itick=0; itick < charges.size(); ++itick) {
                    hist->Fill(ch, tbin+itick, charges[itick]);
                }
            }
        }
    }


    // Handle any trace summaries
    for (auto tag : getset(m_cfg["summaries"])) {
        auto traces = get_tagged_traces(frame, tag);
        auto summary = frame->trace_summary(tag);

        const size_t nchannels = summary.size();
        const int channel0 = traces.front()->channel();

        // Warning: makes huge assumption that summary is defined over
        // a contiguous and ordered span of channel numbers!  Works
        // for MicroBooNE.
        TH1F* hist = new TH1F(("h"+tag).c_str(),("h"+tag).c_str(),
                              nchannels, channel0, channel0 + nchannels);
        for (size_t ch=channel0; ch < channel0+nchannels; ++ch) {
            hist->SetBinContent(1+ch, summary[ch]);
        }
        hist->SetDirectory(output_tf);
    }


    // Now deal with "shunting" input Magnify data to output.
    std::string ifname = m_cfg["input_filename"].asString();
    if (ifname.empty()) {
        // good, we shouldn't be peeking into the input file anyways.
        return true;
    }
    auto toshunt = getset(m_cfg["shunt"]);
    if (toshunt.empty()) {
        std::cerr << "MagnifySink no objects to copy but given input: " << ifname << std::endl;
        return true;
    }
    std::cerr << "MagnifySink: sneaking peaks into input file: " << ifname << std::endl;
    TFile *input_tf = TFile::Open(ifname.c_str());

    for (auto name : toshunt) {
        TObject* obj = input_tf->Get(name.c_str());

        TTree* tree = dynamic_cast<TTree*>(obj);
        if (tree) {
            tree = tree->CloneTree();
            tree->SetDirectory(output_tf);
            continue;
        }

        TH1* hist = dynamic_cast<TH1*>(obj);
        if (hist) {
            hist->SetDirectory(output_tf);
            continue;
        }

    }

    
    {
        Waveform::ChannelMaskMap input_cmm = frame->masks();
        for (auto const& it: input_cmm) {
            
            if (it.first == "bad"){

                // save "bad" channels

                TTree *T_bad = new TTree("T_bad","T_bad");
                int chid, plane, start_time,end_time;
                T_bad->Branch("chid",&chid,"chid/I");
                T_bad->Branch("plane",&plane,"plane/I");
                T_bad->Branch("start_time",&start_time,"start_time/I");
                T_bad->Branch("end_time",&end_time,"end_time/I");
                T_bad->SetDirectory(output_tf);

                for (auto const &it1 : it.second){
                    chid = it1.first;
                    plane = m_anode->resolve(chid).index();
                    for (size_t ind = 0; ind < it1.second.size(); ++ind){
                        start_time = it1.second[ind].first;
                        end_time = it1.second[ind].second;
                        T_bad->Fill();
                    }
                }
                continue;
            }

            if (it.first =="lf_noisy"){

                // save "noisy" channels

                TTree *T_lf = new TTree("T_lf","T_lf");
                int channel;
                T_lf->Branch("channel",&channel,"channel/I");
                for (auto const &it1 : it.second){
                    channel = it1.first;
                    T_lf->Fill();
                }
                continue;
            }
        }
    }


    std::cerr << "MagnifySink: closing output file " << ofname << std::endl;
    output_tf->Write();
    output_tf->Close();

    return true;
}
