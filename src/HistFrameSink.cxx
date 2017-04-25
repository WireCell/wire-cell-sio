#include "WireCellSio/HistFrameSink.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/BoundingBox.h"
#include "WireCellUtil/Waveform.h"

#include "TFile.h"
#include "TH2F.h"

#include <iostream>
#include <algorithm>
#include <unordered_map>

WIRECELL_FACTORY(HistFrameSink, WireCell::Sio::HistFrameSink, WireCell::IFrameSink, WireCell::IConfigurable);


using namespace std;
using namespace WireCell;

Sio::HistFrameSink::HistFrameSink()
    : m_filepat("histframe-%02d.root")
    , m_anode_tn("AnodePlane")
    , m_units(1.0)
{
}

Sio::HistFrameSink::~HistFrameSink()
{
}


WireCell::Configuration Sio::HistFrameSink::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_filepat;
    cfg["anode"] = m_anode_tn;
    cfg["units"] = units::mV;
    return cfg;
}

void Sio::HistFrameSink::configure(const WireCell::Configuration& cfg)
{
    m_filepat = get<std::string>(cfg, "filename", m_filepat);
    m_units = get<double>(cfg, "units", m_units);
    m_anode_tn = get<std::string>(cfg, "anode", m_anode_tn);
    m_anode = Factory::lookup_tn<IAnodePlane>(m_anode_tn);
    if (!m_anode) {
        cerr << "Sio::HistFrameSink: failed to get anode: \"" << m_anode_tn << "\"\n";
        return;
    }
    cerr << "Sio::HistFrameSink: configured with: "
         << "file:" << m_filepat << ", "
         << "units:" << m_units << ", "
         << "anode:" << m_anode_tn << endl;
}



bool Sio::HistFrameSink::operator()(const IFrame::pointer& frame)
{
    std::string fname = Form(m_filepat.c_str(), frame->ident());
    TFile* file = TFile::Open(fname.c_str(), "recreate");
    

    typedef std::tuple< ITrace::vector, std::vector<int>, std::vector<int> > tct_tuple;
    std::unordered_map<int, tct_tuple> perplane; 

    // collate traces into per plane and also calculate bounds 
    ITrace::shared_vector traces = frame->traces();
    for (ITrace::pointer trace : *traces) {
        int ch = trace->channel();
        auto wpident = m_anode->resolve(ch).ident(); 
        double tmin = trace->tbin();
        double tlen = trace->charge().size();

        auto& tct = perplane[wpident];
        get<0>(tct).push_back(trace);
        get<1>(tct).push_back(ch);
        get<2>(tct).push_back(tmin);
        get<2>(tct).push_back(tmin+tlen);
    }

    const double t0 = frame->time();
    const double tick = frame->tick();


    for (auto& thisplane : perplane) {
        int wpident = thisplane.first;
        auto& tct = thisplane.second;
        auto& traces = get<0>(tct);
        auto& chans = get<1>(tct);
        auto chmm = std::minmax_element(chans.begin(), chans.end());
        auto& tbins = get<2>(tct);
        auto tbmm = std::minmax_element(tbins.begin(), tbins.end());


        const double tmin = t0 + tick*(*tbmm.first);
        const double tmax = t0 + tick*(*tbmm.second);
        const int ntbins = (*tbmm.second)-(*tbmm.first);

        const float chmin = *chmm.first;
        const float chmax = *chmm.second + 1;
        const int nchbins = (*chmm.second) - (*chmm.first) + 1;
        
        TH2F* hist = new TH2F(Form("plane%d", wpident),
                              Form("Plane %d", wpident),
                              ntbins, tmin/units::us, tmax/units::us,
                              nchbins, chmin, chmax);
        hist->SetDirectory(file); // needed?
        hist->SetXTitle("time [us]");
        hist->SetYTitle("channel");

        for (auto& trace : traces) {
            double fch = trace->channel() + 0.5; // make sure we land in bin-center.
            int tbin = trace->tbin();
            auto& charge = trace->charge();
            int nbins = charge.size();
            for (int ibin=0; ibin<nbins; ++ibin) {
                const double t = t0 + (tick)*(tbin+ibin+0.5); // 0.5 to land in bin-center
                hist->Fill(t/units::us, fch, charge[ibin]/m_units);
            }
        }

        cerr << wpident
             << " qunit:" << m_units << " "
             << " integ:" << hist->Integral()
             << " min:" << hist->GetMinimum()
             << " max:" << hist->GetMaximum()
             << " chan:[" << chmin << "," << chmax << "] "
             << " time:[" << tmin/units::us << "," << tmax/units::us <<"]us\n";


        hist->Write();
    }
    file->Close();
    delete file;
    file = nullptr;
    return true;
}
