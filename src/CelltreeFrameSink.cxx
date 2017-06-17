#include "WireCellSio/CelltreeFrameSink.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/BoundingBox.h"
#include "WireCellUtil/Waveform.h"

#include "TFile.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TTree.h"
#include "TH1F.h"
#include "TClonesArray.h"


#include <iostream>
#include <algorithm>
#include <unordered_map>

WIRECELL_FACTORY(CelltreeFrameSink, WireCell::Sio::CelltreeFrameSink, WireCell::IFrameSink, WireCell::IConfigurable);


using namespace std;
using namespace WireCell;

Sio::CelltreeFrameSink::CelltreeFrameSink()
    : m_filepat("celltreeframe-%02d.root")
    , m_anode_tn("AnodePlane")
    , m_units(1.0)
    , m_readout(5.0*units::ms)
{
}

Sio::CelltreeFrameSink::~CelltreeFrameSink()
{
}


WireCell::Configuration Sio::CelltreeFrameSink::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_filepat;
    cfg["anode"] = m_anode_tn;
    cfg["units"] = units::uV;
    cfg["readout_time"] = m_readout;
    return cfg;
}

void Sio::CelltreeFrameSink::configure(const WireCell::Configuration& cfg)
{
    m_filepat = get<std::string>(cfg, "filename", m_filepat);
    m_units = get<double>(cfg, "units", m_units);
    m_readout = get<double>(cfg, "readout_time", m_readout);
    m_anode_tn = get<std::string>(cfg, "anode", m_anode_tn);
    m_anode = Factory::lookup_tn<IAnodePlane>(m_anode_tn);
    if (!m_anode) {
        cerr << "Sio::CelltreeFrameSink: failed to get anode: \"" << m_anode_tn << "\"\n";
        return;
    }
    cerr << "Sio::CelltreeFrameSink: configured with: "
         << "file:" << m_filepat << ", "
         << "units:" << m_units << ", "
         << "anode:" << m_anode_tn << ", "
         << "readout: "<< m_readout/units::ms <<endl;
}



bool Sio::CelltreeFrameSink::operator()(const IFrame::pointer& frame)
{
    std::string fname = Form(m_filepat.c_str(), frame->ident());
    
    // [HY] Simulation output celltree
    // Sim Tree
    TTree *Sim = new TTree("Sim","Wire-cell toolkit simulation output");
    // Branches init
    Int_t runNo = 0; 
    Int_t subRunNo = 0;
    Int_t eventNo = 0;
    Int_t raw_nChannel = 8256;
    std::vector<int> *raw_channelId = new std::vector<int>;
    TClonesArray *sim_wf = new TClonesArray("TH1I");
    TH1::AddDirectory(kFALSE);
    // Brances set
    Sim->Branch("runNo", &runNo, "runNo/I");
    Sim->Branch("subRunNo", &subRunNo, "subRunNo/I");
    Sim->Branch("eventNo", &eventNo, "eventNo/I");
    Sim->Branch("raw_nChannel", &raw_nChannel, "raw_nChannel/I");

    Sim->Branch("raw_channelId", &raw_channelId);
    Sim->Branch("raw_wf", &sim_wf, 256000, 0);
    // for sim_wf fill one by one
    int sim_wf_ind = 0;  

    /* // fill channel: if all channels need to be filled */
    /* auto& chanchan = m_anode->channels(); */
    /* raw_nChannel = chanchan.size(); */
    /* std::cout<<"nChannel: "<<raw_nChannel<<std::endl; */
    /* for(int i=0; i<raw_nChannel; i++){ */
    /*     raw_channelId->push_back(chanchan.At(i)); */   
    /*     std::cout<<"channelId: "<<chanchan.At(i)<<std::endl; */
    /* } */


    typedef std::tuple< ITrace::vector, std::vector<int>, std::vector<int> > tct_tuple;
    std::unordered_map<int, tct_tuple> perplane; 

    const double tick = frame->tick();
    //std::cout<<"tick: "<<tick/units::us<<std::endl;
    const int nsamples = m_readout/tick;
    //std::cout<<"nsamples: "<<nsamples<<std::endl;


    /** [DEBUG] storage for searching min/max channel and tbin
    typedef std::tuple< std::vector<int>, std::vector<int> > tc_tuple;
    std::unordered_map<int, tc_tuple> planes;
     */


    // collate traces into per plane and also calculate bounds 
    ITrace::shared_vector traces = frame->traces();
    for (ITrace::pointer trace : *traces) {
        int ch = trace->channel();
        //std::cout<<"channel number: "<<ch<<std::endl;
        raw_channelId->push_back(ch);
        // fill raw_wf
        TH1I *htemp = new ( (*sim_wf)[sim_wf_ind] ) TH1I("", "",  nsamples, 0,  nsamples);
        auto& wave = trace->charge();
        int nbins = wave.size();
        //std::cout<<"waveform size: "<<nbins<<std::endl;
        int tmin = trace->tbin();
        //std::cout<<"tmin: "<<tmin<<std::endl;
        for(Int_t i=0; i<nbins; i++){
            if(tmin+i+1<=nsamples){ 
            htemp->SetBinContent(tmin+i+1, (int)wave[i]/m_units);
            }
        } 
        sim_wf_ind ++;

        /** [DEBUG]
        auto& tc = planes[m_anode->resolve(ch).ident()];
        get<0>(tc).push_back(ch);
        get<1>(tc).push_back(tmin);
        get<1>(tc).push_back(tmin+nbins);
          */


    }
    raw_nChannel = sim_wf_ind;
    //cout<<"nChannel: "<<raw_nChannel<<std::endl;

    Sim->Fill();

    /** [DEBUG]
     * Histogram copy of CelltreeFrameSink waveforms (TClonesArray)
     * easy used for validation against HistoFrameSink output
    const double t0 = frame->time();
    TH2F* hist[3]; // not auto, need to know wpident list first of all
    int tbng[3];
    int tlen[3];
    for (auto& thisplane : planes) {
        int wpident = thisplane.first;
        auto& tc = thisplane.second;
        auto& chans = get<0>(tc);
        auto chmm = std::minmax_element(chans.begin(), chans.end());
        auto& tbins = get<1>(tc);
        auto tbmm = std::minmax_element(tbins.begin(), tbins.end());


        const double tmin = t0 + tick*(*tbmm.first);
        const double tmax = t0 + tick*(*tbmm.second);
        const int ntbins = (*tbmm.second)-(*tbmm.first);

        const int chmin = round(*chmm.first);
        const int chmax = round(*chmm.second + 1);
        const int nchbins = chmax - chmin;

        int hind = wpident/2;
        tbng[hind] = (*tbmm.first);
        tlen[hind] = ntbins;
        hist[hind] = new TH2F(Form("plane%d", wpident),
                Form("Plane %d", wpident),
                ntbins, tmin/units::us, tmax/units::us,
                nchbins, chmin, chmax);
        hist[hind]->SetXTitle("time [us]");
        hist[hind]->SetYTitle("channel");
    }

    for(int ich=0; ich<raw_nChannel; ich++)
    {
        int chh = raw_channelId->at(ich);
        TH1F* h = dynamic_cast<TH1F*>(sim_wf->At(ich));
        int pind = m_anode->resolve(chh).ident()/2;
        //cout<<"Plane ID: "<<pind<<endl;

        for(int it = 0; it<tlen[pind]; it++){
            int abstind = tbng[pind]+it;
            hist[pind]->Fill((t0+(tick)*(abstind+0.5))/units::us, chh+0.5, h->GetBinContent(abstind+1)/m_units);
        }
    }

    

    hist[0]->Write();
    hist[1]->Write();
    hist[2]->Write();
    /// end
    */
   

    TFile* file = TFile::Open(fname.c_str(), "recreate");
    TDirectory *event = file->mkdir("Event");
    event->cd();
    Sim->Write();
    
    file->Close();
    sim_wf->Delete();
    delete raw_channelId;
    delete file;
    file = nullptr;
    return true;
}
