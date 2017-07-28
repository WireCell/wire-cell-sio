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

    // fixme: this TOTALLY violates the design of wire cell DFP.
    cfg["input_filename"] = ""; 

    // List of TObjects to copy from input file to output file.
    cfg["shunt"] = Json::arrayValue;

    // Name of ROOT file to write.
    cfg["output_filename"] = "";

    // A list of trace tags defining which waveforms are saved to Magnify histograms.
    cfg["frames"] = Json::arrayValue;

    // A list of pairs mapping a cmm key name to a ttree name.
    cfg["cmmtree"] = Json::arrayValue;

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


std::vector<WireCell::Binning> collate_byplane(const ITrace::vector& traces, const IAnodePlane::pointer anode,
                                               ITrace::vector byplane[])
{
    std::vector<int> uvwt[4];
    for (auto trace : traces) {
        const int chid = trace->channel();
        auto wpid = anode->resolve(chid);
        const int iplane = wpid.index();
        //std::cerr << "\tchid="<<chid<<" iplane="<<iplane<<" wpid="<<wpid<<std::endl;
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
        if (one.empty()) {
            std::cerr << "MagnifySink: bogus bounds " << ind << "\n";
            THROW(ValueError() << errmsg{"MagnifySink: bogus bounds"});
        }
            
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

        std::cerr << "MagnifySink: tag: \"" << tag << "\" with " << traces.size() << " traces\n";
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
                    // Using what should be an identical call to
                    // Fill() ends up with a file that is more than
                    // two times bigger:
                    // 826M Jul 27 18:22 orig-bl-nf-fill-tbin.root
                    // 342M Jul 27 18:28 orig-bl-nf-setbincontent-tplus1.root
                    //hist->Fill(ch, tbin+itick+0.5, charges[itick]);

                    hist->SetBinContent(cbin.bin(ch)+1, tbin+itick+1, charges[itick]);
                }
            }
        }
    }


    // Handle any trace summaries
    for (auto tag : getset(m_cfg["summaries"])) {
        auto traces = get_tagged_traces(frame, tag);
        if (traces.empty()) {
            std::cerr << "MagnifySink: warning: no traces tagged with \"" << tag << "\", skipping summary\n";
            continue;
        }
        auto const& summary = frame->trace_summary(tag);
        if (summary.empty()) {
            std::cerr << "MagnifySink: warning: empty summary tagged with \"" << tag << "\", skipping summary\n";
            continue;
        }
            
        std::cerr << "MagnifySink: saving summaries tagged with \"" << tag << "\" into per-plane hists\n";

        // warning: this is going to get ugly.  we jump through these
        // hoops to avoid hard coding microboone channel ranges and
        // allow for sparse summaries on the domain of channels;
        const int ntot = traces.size();
        std::vector<int> perplane_channels[3];
        std::vector<double> perplane_values[3];
        for (int ind=0; ind<ntot; ++ind) {
            const int chid = traces[ind]->channel();
            const int iplane = m_anode->resolve(chid).index();
            perplane_channels[iplane].push_back(chid);
            perplane_values[iplane].push_back(summary[ind]);
        }
        for (int iplane=0; iplane<3; ++iplane) {
            std::vector<int>& chans = perplane_channels[iplane];
            std::vector<double>& vals = perplane_values[iplane];
            auto mme = std::minmax_element(chans.begin(), chans.end());
            const int ch0 = *mme.first;
            const int chf = *mme.second;
            const std::string hname = Form("h%c_%s", 'u'+iplane, tag.c_str());
            TH1F* hist = new TH1F(hname.c_str(), hname.c_str(),
                                  chf-ch0+1, ch0, chf);
            for (size_t ind=0; ind<chans.size(); ++ind) {
                const int ch = chans[ind];
                const double val = vals[ind];
                hist->Fill(ch+0.5, val);
            }
            hist->SetDirectory(output_tf);
        }
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

        if (!obj) {
            std::cerr << "MagnifySink: warning: failed to find input object: \"" << name << "\" for copying to output\n";
        }

        TTree* tree = dynamic_cast<TTree*>(obj);
        if (tree) {
            std::cerr << "MagnifySink: copying tree: \"" << name << "\"\n";
            tree = tree->CloneTree();
            tree->SetDirectory(output_tf);
            continue;
        }

        TH1* hist = dynamic_cast<TH1*>(obj);
        if (hist) {
            std::cerr << "MagnifySink: copying hist: \"" << name << "\"\n";
            hist->SetDirectory(output_tf);
            continue;
        }

        std::cerr << "MagnifySink: warning: object of unknown type: \"" << name << "\", will not copy\n";

    }

    
    {
        std::unordered_map<std::string, std::string> cmmkey2treename;
        for (auto jcmmtree : m_cfg["cmmtree"]) {
            cmmkey2treename[jcmmtree[0].asString()] = jcmmtree[1].asString();
        }

        Waveform::ChannelMaskMap input_cmm = frame->masks();
        for (auto const& it: input_cmm) {
            auto cmmkey = it.first;
            auto ct = cmmkey2treename.find(cmmkey);
            if (ct == cmmkey2treename.end()) {
                std::cerr << "MagnifySink: warning: no tree configured to save channel mask \"" << cmmkey << "\"\n";
                continue;
            }
            
            auto treename = ct->second;
            
            std::cerr << "MagnifySink: saving channel mask \"" << cmmkey << "\" to tree \"" << treename << "\"\n";

            TTree *tree = new TTree(treename.c_str(), treename.c_str());
            int chid=0, plane=0, start_time=0, end_time=0;
            tree->Branch("chid",&chid,"chid/I");
            tree->Branch("plane",&plane,"plane/I");
            tree->Branch("start_time",&start_time,"start_time/I");
            tree->Branch("end_time",&end_time,"end_time/I");
            tree->SetDirectory(output_tf);

            for (auto const &chmask : it.second){
                chid = chmask.first;
                plane = m_anode->resolve(chid).index();
                auto mask = chmask.second;
                for (size_t ind = 0; ind < mask.size(); ++ind){
                    start_time = mask[ind].first;
                    end_time = mask[ind].second;
                    tree->Fill();
                }
            }
        }
    }


    std::cerr << "MagnifySink: closing output file " << ofname << std::endl;
    output_tf->Write();
    output_tf->Close();

    return true;
}
