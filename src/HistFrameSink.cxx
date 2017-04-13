#include "WireCellSio/HistFrameSink.h"
#include "WireCellUtil/NamedFactory.h"

#include <iostream>

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
    std::string fname = Form(m_filepath.c_str(), frame->ident());
    TFile* file = TFile::Open(fname.c_str, "recreate");
    
    // need way to go from channel to frame.

    auto traces = frame->traces();
    const int ntraces = traces->size();
    
    cerr << "Frame: #" << frame->ident()
         << " @" << frame->time()/units::ms
         << " with " << ntraces << " traces" << endl;

    for (int iplane=0; iplane<3; ++iplane) {

    TH2F* hist = new TH2F(Form("plane%d", iplane),
                          Form("Plane %d", iplane),
                          ntbins, tmin/units::us, tmax/units::us,
                          nwbins, minwire, maxwire);
    hist->SetDirectory(file)
        for (auto it : frame) {
            const int iwire = it.first;
            const auto& wave = it.second;
            for (int itick = mintick; itick < maxtick; ++itick) {
                hist.Fill(tbins.center(itick)/units::us, iwire, wave[itick]);
            }
        }
        hist.Write();
    }

}
