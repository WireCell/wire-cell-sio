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

static std::vector<std::string> known_cateogries{"threshold","baseline","orig","raw","decon"};

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
        THROW(ValueError() << errmsg{"Must provide input filename to MagnifySink"});
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

    cfg["rebin"] = 1;           // fixme: so does this.  Rebinning should be done in an IFrameFilter
    cfg["input_filename"] = ""; // fixme: this TOTALLY violates the design of wire cell DFP
    // More evilness: these are the histogram types that I know about and will shunt by default.
    for (size_t ind=0; ind<known_cateogries.size(); ++ind) {
        cfg["shunt"][(int)ind] = known_cateogries[ind];
    }

    cfg["output_filename"] = "";
    cfg["histtype"] = "decon";

    return cfg;
}


bool Sio::MagnifySink::operator()(const IFrame::pointer& frame)
{
    std::string ifname = m_cfg["input_filename"].asString();
    std::string ofname = m_cfg["output_filename"].asString();
    const int nrebin = m_cfg["rebin"].asInt();
    std::string histtype = m_cfg["histtype"].asString();

    // figure out what histograms the user wants to "shunt" to the
    // output.
    std::vector<std::string> toshunt = known_cateogries;
    auto jshunt = m_cfg["shunt"];
    if (!jshunt.empty()) {
        toshunt.clear();
        for (auto jht : jshunt) {
            toshunt.push_back(jht.asString());
        }
    }
    std::vector<std::string> tmp;
    for (auto ht : toshunt) {
        if (ht != histtype) {   // exclude the one we are actually
            tmp.push_back(ht);  // writing from the frame
        }
    }
    toshunt = tmp;
        
    // sus out channel and tick bins
    std::vector<int> uvwt[4];
    ITrace::vector traces[3];
    for (auto trace : *frame->traces()) {
        const int chid = trace->channel();
        auto wpid = m_anode->resolve(chid);
        const int iplane = wpid.index();
        if (iplane<0 || iplane>=3) {
            THROW(RuntimeError() << errmsg{"Illegal wpid"});
        }
        uvwt[iplane].push_back(chid);
        traces[iplane].push_back(trace);
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
            // Keep histogram bounds still in original tick number but rebin
            const int n = vmax - vmin;
            binnings.push_back(Binning(n/nrebin, vmin, vmax));
        }
        else {
            // Channel-centered binning
            const double diff = vmax - vmin;
            binnings.push_back(Binning(diff+1, vmin-0.5, vmax+0.5));
        } 
    }

    TFile *input_tf = TFile::Open(ifname.c_str());
    TFile* output_tf = TFile::Open(ofname.c_str(), "RECREATE");

    // do evil shunting from input file
    {
        TTree *Trun = ((TTree*)input_tf->Get("Trun"))->CloneTree();
        Trun->SetDirectory(output_tf);
    }

    // more evilness.  thresholds get rewritten, if they exist, with
    // stuff stashed in the channel mask maps, themselves an unding
    // source of evilness.
    std::vector<TH1F*> thresholds(3, nullptr);

    for (auto ht : toshunt) {
        for (int iplane=0; iplane<3; ++iplane) {
            const std::string name = Form("h%c_%s", 'u'+iplane, ht.c_str());
            TH1* obj = (TH1*)input_tf->Get(name.c_str());
            if (!obj) {
                std::cerr <<"MagnifySink: warning \"" << name << "\" not found in " << ifname << std::endl;
                continue;
            }
            std::cerr <<"MagnifySink: evilly shunting \"" << name << "\"\n";
            obj->SetDirectory(output_tf);
            if (ht == "threshold") {
                std::cerr <<iplane<<": " << obj->IsA()->ClassName() << std::endl;
                thresholds[iplane] = (TH1F*)obj;
            }
        }
    }
                
    // make actual output histogram
    std::vector<TH2F*> hists;
    for (int ind=0; ind<3; ++ind) {
        const std::string name = Form("h%c_%s", 'u'+ind, histtype.c_str());
        Binning cbin = binnings[ind];
        Binning tbin = binnings[3];
        std::cerr << "MagnifySink:"
                  << " cbin:"<<cbin.nbins()<<"["<<cbin.min() << "," << cbin.max() << "]"
                  << " tbin:"<<tbin.nbins()<<"["<<tbin.min() << "," << tbin.max() << "]\n";

        TH2F* hist = new TH2F(name.c_str(), name.c_str(),
                              cbin.nbins(), cbin.min(), cbin.max(),
                              tbin.nbins(), tbin.min(), tbin.max());
        hist->SetDirectory(output_tf);
        for (auto trace : traces[ind]) {
            const int tbin = trace->tbin();
            const int ch = trace->channel();
            auto const& charges = trace->charge();
            for (size_t itick=0; itick < charges.size(); ++itick) {
                hist->Fill(ch, tbin+itick, charges[itick]);
                // note: Xin's original jumped through hoops to use
                // SetBinContent().  Do I miss something?
            }
        }
        if (thresholds.size() < 3) { // make these fresh if we don't already have them.
            const std::string name = Form("h%c_threshold", 'u'+ind);
            TH1F* thresh = new TH1F(name.c_str(), name.c_str(),
                                    cbin.nbins(), cbin.min(), cbin.max());
            thresholds.push_back(thresh);
            thresh->SetDirectory(output_tf);
        }        
    }
    
    {
        // save "bad" channels 
        TTree *T_bad = new TTree("T_bad","T_bad");
        int chid, plane, start_time,end_time;
        T_bad->Branch("chid",&chid,"chid/I");
        T_bad->Branch("plane",&plane,"plane/I");
        T_bad->Branch("start_time",&start_time,"start_time/I");
        T_bad->Branch("end_time",&end_time,"end_time/I");
        T_bad->SetDirectory(output_tf);

        // save "noisy" channels
        TTree *T_lf = new TTree("T_lf","T_lf");
        int channel;
        T_lf->Branch("channel",&channel,"channel/I");

        Waveform::ChannelMaskMap input_cmm = frame->masks();
        for (auto const& it: input_cmm) {
            
            if (it.first == "bad"){ // save bad ... 
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
                for (auto const &it1 : it.second){
                    channel = it1.first;
                    T_lf->Fill();
                }
                continue;
            }
            if (it.first=="threshold"){ 
                for (auto const &it1 : it.second){
                    chid = it1.first;
                    TH1F* hthresh = thresholds[plane]; // evil
                    const float tval = it1.second[0].first/it1.second[0].second; // evilevil
                    hthresh->SetBinContent(chid+1, tval*nrebin*3.0); // evilevilevil
                }
                continue;
            }
        }
    }


    std::cerr << "MagnifySink: closing output file\n";
    output_tf->Write();
    output_tf->Close();

    return true;
}
