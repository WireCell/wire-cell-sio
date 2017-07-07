#include "WireCellSio/MagnifyFileIterator.h"
#include "WireCellIface/SimpleTrace.h"
#include "WireCellIface/SimpleFrame.h"

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"

using namespace WireCell;
using namespace std;

Sio::MagnifyFileIterator::MagnifyFileIterator(const char* filename, const char* histtype)
{
    file = TFile::Open(filename);
    string uvw = "uvw";
    for (int ind=0; ind<3; ++ind) {
	auto c = uvw[ind];
	std::string name = Form("h%c_%s", c, histtype);
	cerr << "Loading " << name << endl;
	hist[ind] = (TH2*)file->Get(name.c_str());
    }
      
    TTree *T_bad = (TTree*)file->Get("T_bad");
    int chid, plane, start_time,end_time;
    T_bad->SetBranchAddress("chid",&chid);
    T_bad->SetBranchAddress("plane",&plane);
    T_bad->SetBranchAddress("start_time",&start_time);
    T_bad->SetBranchAddress("end_time",&end_time);
      
    for (int i=0;i!=T_bad->GetEntries();i++){
	T_bad->GetEntry(i);
	WireCell::Waveform::BinRange chirped_bins;
	chirped_bins.first = start_time;
	chirped_bins.second = end_time;
	ret["bad"][chid].push_back(chirped_bins);
    }
      
      
    TTree *T_lf = (TTree*)file->Get("T_lf");
    int channel;
    T_lf->SetBranchAddress("channel",&channel);
    for (int i=0;i!=T_lf->GetEntries();i++){
	T_lf->GetEntry(i);
	WireCell::Waveform::BinRange chirped_bins;
	chirped_bins.first = 0;
	chirped_bins.second = hist[0]->GetNbinsY();
	ret["lf_noisy"][channel].push_back(chirped_bins);
    }
    delete T_lf;
    delete T_bad;
      
    //file->Close();
    //delete file;
}

int Sio::MagnifyFileIterator::plane(int ch)
{
    if (ch < 2400) return 0;
    if (ch < 2400+2400) return 1;
    return 2;
}
int Sio::MagnifyFileIterator::index(int ch)
{
    if (ch < 2400) return ch;
    if (ch < 2400+2400) return ch-2400;
    return ch-2400-2400;
}

vector<float> Sio::MagnifyFileIterator::at(int ch) 
{
    TH2* h = hist[plane(ch)];
    int ind = index(ch);
    vector<float> ret(9600);
    for (int itick=0; itick<9600; ++itick) {
        ret[itick] = h->GetBinContent(ind+1, itick+1);
    }
    return ret;
}

void Sio::MagnifyFileIterator::clear()
{
    delete hist[0];
    delete hist[1];
    delete hist[2];
    
    file->Close();
    delete file;
}
  
/// Return a frame, the one and only in the file.
IFrame::pointer Sio::MagnifyFileIterator::frame()
{
    ITrace::vector traces;

    int chindex=0;
    for (int iplane=0; iplane<3; ++iplane) {
        TH2* h = hist[iplane];

        int nchannels = h->GetNbinsX();
        int nticks = h->GetNbinsY();

        cerr << "plane " << iplane << ": " << nchannels << " X " << nticks << endl;

        double qtot = 0.0;
        for (int ich = 1; ich <= nchannels; ++ich) {
            ITrace::ChargeSequence charges;
            for (int itick = 1; itick <= nticks; ++itick) {
                auto q = h->GetBinContent(ich, itick);
                charges.push_back(q);
                qtot += q;
            }
            SimpleTrace* st = new SimpleTrace(chindex, 0.0, charges);
            traces.push_back(ITrace::pointer(st));
            ++chindex;
            //cerr << "qtot in plane/ch/index "
            //     << iplane << "/" << ich << "/" << chindex << " = " << qtot << endl;
        }
    }
    SimpleFrame* sf = new SimpleFrame(0, 0, traces, 0.5*units::microsecond, ret);
    return IFrame::pointer(sf);
}
    


