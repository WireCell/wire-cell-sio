#include "WireCellSio/HistFrameSink.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/BoundingBox.h"

#include "TFile.h"
#include "TH2F.h"

#include <iostream>
#include <unordered_map>

using namespace std;
using namespace WireCell;

Sio::HistFrameSink::HistFrameSink()
    : m_filepat("histframe-%02d.root")
    , m_anode_tn("AnodePlane")
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
    return cfg;
}

void Sio::HistFrameSink::configure(const WireCell::Configuration& cfg)
{
    m_filepat = get<std::string>(cfg, "filename", m_filepat);

    m_anode_tn = get<std::string>(cfg, "anode", m_anode_tn);
    m_anode = Factory::lookup_tn<IAnodePlane>(m_anode_tn);
    if (!m_anode) {
        cerr << "Sio::HistFrameSink: failed to get anode: \"" << m_anode_tn << "\"\n";
        return;
    }
}


bool Sio::HistFrameSink::operator()(const IFrame::pointer& frame)
{
    std::string fname = Form(m_filepat.c_str(), frame->ident());
    TFile* file = TFile::Open(fname.c_str(), "recreate");
    

    // fixme: index by packed integer as WirePlaneId can not be hashed.

    // plane ident to its traces
    std::unordered_map<int, ITrace::vector> perplane; 

    // fixme: this is an abuse of bounding box.  A general,
    // N-dimensional std::tuple based BB should be developed!
    std::unordered_map<int, BoundingBox> limits;      

    ITrace::shared_vector traces = frame->traces();
    for (ITrace::pointer trace : *traces) {
        int ch = trace->channel();
        auto wpident = m_anode->resolve(ch).ident(); 
        perplane[wpident].push_back(trace);
        double tmin = trace->tbin();
        double tlen = trace->charge().size();
        limits[wpident](Point(ch, tmin, tmin+tlen));
    }

    const double t0 = frame->time();
    const double tick = frame->tick();


    for (auto wpid_bb : limits) {
        int wpident = wpid_bb.first;
        auto wpid = WirePlaneId(wpident);
        Ray bb = wpid_bb.second.bounds();

        // fixme: this is rife for off-by-one bugs
        const double tmin = t0 + tick*bb.first[1];
        const double tmax = t0 + tick*bb.second[2];
        const int ntbins = (tmax-tmin)/tick;
        const int chmin = round(bb.first[0]);
        const int chmax = round(bb.first[1]);
        const int nchbins = chmax - chmin + 1;

        
        TH2F* hist = new TH2F(Form("plane%d", wpident),
                              Form("Plane %d", wpident),
                              ntbins, tmin/units::us, tmax/units::us,
                              nchbins, chmin, chmax);
        hist->SetDirectory(file); // needed?

        auto& traces = perplane[wpident];
        for (auto& trace : traces) {
            double fch = trace->channel() + 0.5; // make sure we land in bin-center.
            int tbin = trace->tbin();
            auto& charge = trace->charge();
            int nbins = charge.size();
            for (int ibin=0; ibin<nbins; ++ibin) {
                const double t = t0 + (tick+0.5)*(tbin+ibin); // 0.5 to land in bin-center
                hist->Fill(t/units::us, fch, charge[ibin]);
            }
        }
        hist->Write();
    }
    file->Close();
    delete file;
    file = nullptr;
}
