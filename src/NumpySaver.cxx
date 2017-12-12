#include "WireCellSio/NumpySaver.h"

#include "WireCellIface/FrameTools.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/cnpy.h"

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <tuple>

WIRECELL_FACTORY(NumpySaver, WireCell::Sio::NumpySaver,
                 WireCell::IFrameFilter, WireCell::IDepoFilter, WireCell::IConfigurable)

using namespace WireCell;

Sio::NumpySaver::NumpySaver()
    : m_save_count(0)
{
}

Sio::NumpySaver::~NumpySaver()
{
}


WireCell::Configuration Sio::NumpySaver::default_configuration() const
{
    Configuration cfg;

    // If digitize is true, then samples as 16 bit ints.  Otherwise
    // save as 32 bit floats.
    cfg["digitize"] = false;

    // This number is set to the waveform sample array before any
    // charge is added.
    cfg["baseline"] = 0.0;

    // This number will be multiplied to each waveform sample before
    // casting to dtype.
    cfg["scale"] = 1.0;

    // This number will be added to each scaled waveform sample before
    // casting to dtype.
    cfg["offset"] = 0.0;

    // The frame tags to consider for saving
    cfg["frame_tags"] = Json::arrayValue;
    // The summary tags to consider for saving
    //cfg["summary_tags"] = Json::arrayValue;    
    // The channel mask maps to consider for saving
    //cfg["chanmaskmaps"] = Json::arrayValue;

    // The output file name to write.  Only compressed (zipped) Numpy
    // files are supported.  Writing is always in "append" mode.  It's
    // up to the user to delete a previous instance of the file if
    // it's old contents are not wanted.
    cfg["filename"] = "wct-frame.npz";
        
    return cfg;
}

void Sio::NumpySaver::configure(const WireCell::Configuration& config)
{
    m_cfg = config;
}


typedef std::tuple<IDepo::pointer, size_t, size_t> depo_gen_child;
typedef std::vector< depo_gen_child > depos_with_prior;

static void push_depo(depos_with_prior& dp, WireCell::IDepo::pointer depo, size_t gen=0, size_t childid=0)
{
    dp.push_back(depo_gen_child(depo, gen, childid));
    auto prior = depo->prior();
    if (!prior) {
        return;
    }
    push_depo(dp, prior, gen+1, dp.size());
}
static depos_with_prior flatten_depos(std::vector<WireCell::IDepo::pointer> depos)
{
    depos_with_prior ret;
    for (auto depo : depos) {
        push_depo(ret, depo);
    }
    return ret;
}



bool Sio::NumpySaver::operator()(const IFrame::pointer& inframe,
                                 IFrame::pointer& outframe)
{
    outframe = inframe;         // pass through actual frame

    const std::string mode = "a";

    const float baseline = m_cfg["baseline"].asFloat();
    const float scale = m_cfg["scale"].asFloat();
    const float offset = m_cfg["offset"].asFloat();
    const bool digitize = m_cfg["digitize"].asBool();

    const std::string fname = m_cfg["filename"].asString();

    // Eigen3 array is indexed as (irow, icol) or (ichan, itick)
    // one row is one channel, one column is a tick.
    // Numpy saves reversed dimensions: {ncols, nrows} aka {ntick, nchan} dimensions.

    if (! m_cfg["frame_tags"].isNull()) {
        for (auto jtag : m_cfg["frame_tags"]) {
            const std::string tag = jtag.asString();
            auto traces = FrameTools::tagged_traces(inframe, tag);
            if (traces.empty()) {
                continue;
            }
            auto channels = FrameTools::channels(traces);
            std::sort(channels.begin(), channels.end());
            auto chbeg = channels.begin();
            auto chend = std::unique(chbeg, channels.end());
            auto tbinmm = FrameTools::tbin_range(traces);

            // fixme: may want to give user some config over tbin range to save.
            const size_t ncols = tbinmm.second-tbinmm.first;
            const size_t nrows = std::distance(chbeg, chend);
            Array::array_xxf arr = Array::array_xxf::Zero(nrows, ncols) + baseline;
            FrameTools::fill(arr, traces, channels.begin(), chend, tbinmm.first);
            arr = arr * scale + offset;

            {                   // the 2D frame array
                const std::string aname = String::format("frame_%s_%d", tag.c_str(), m_save_count);
                if (digitize) {
                    Array::array_xxs sarr = arr.cast<short>();
                    const short* sdata = sarr.data();
                    cnpy::npz_save(fname, aname, sdata, {ncols, nrows}, mode);
                }
                else {
                    cnpy::npz_save(fname, aname, arr.data(), {ncols, nrows}, mode);
                }
                std::cerr << "Saved " << aname << " with " << nrows << " channels "
                          << ncols << " ticks @t=" << inframe->time() / units::ms << "ms\n";

            }

            {                   // the channel array
                const std::string aname = String::format("channels_%s_%d", tag.c_str(), m_save_count);
                cnpy::npz_save(fname, aname, channels.data(), {nrows}, mode);
            }

            {                   // the tick array
                const std::string aname = String::format("tickinfo_%s_%d", tag.c_str(), m_save_count);
                const std::vector<double> tickinfo{inframe->time(), inframe->tick(), (double)tbinmm.first};
                cnpy::npz_save(fname, aname, tickinfo.data(), {3}, mode);
            }

        }

    } // if any frame tags to save

    const size_t ndepos = m_depos.size();
    if (ndepos) {
        
        auto fdepos = flatten_depos(m_depos);
        const size_t nfdepos = fdepos.size();

        // time, charge, x, y, z, dlong, dtran
        const size_t ndata=7;
        Array::array_xxf data(nfdepos, ndata);
        // ID, pdg, gen, child
        const size_t ninfo = 4;
        Array::array_xxi info(nfdepos, ninfo);
        for (size_t idepo=0; idepo != nfdepos; ++idepo) {
            auto depogc = fdepos[idepo];
            auto depo = std::get<0>(depogc);
            auto gen = std::get<1>(depogc);
            auto child = std::get<2>(depogc);
            data(idepo, 0) = depo->time();
            data(idepo, 1) = depo->charge();
            data(idepo, 2) = depo->pos().x();
            data(idepo, 3) = depo->pos().y();
            data(idepo, 4) = depo->pos().z();
            data(idepo, 5) = depo->extent_long();
            data(idepo, 6) = depo->extent_tran();
            info(idepo, 0) = depo->id();
            info(idepo, 1) = depo->pdg();
            info(idepo, 2) = gen;
            info(idepo, 3) = child;
        }
        const std::string data_name = String::format("depo_data_%d", m_save_count);
        const std::string info_name = String::format("depo_info_%d", m_save_count);

        cnpy::npz_save(fname, data_name, data.data(), {ndata, ndepos}, mode);
        cnpy::npz_save(fname, info_name, info.data(), {ninfo, ndepos}, mode);
        m_depos.clear();
    }


    ++m_save_count;
    return true;
}

bool Sio::NumpySaver::operator()(const WireCell::IDepo::pointer& indepo,
                                 WireCell::IDepo::pointer& outdepo)
{
    outdepo = indepo;
    if (indepo) {
        m_depos.push_back(indepo);
    }
    return true;
}
