#include "WireCellSio/NumpySaver.h"

#include "WireCellIface/FrameTools.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/cnpy.h"

#include <string>
#include <vector>
#include <algorithm>

WIRECELL_FACTORY(NumpySaver, WireCell::Sio::NumpySaver,
                 WireCell::IFrameFilter, WireCell::IConfigurable)

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
    cfg["summary_tags"] = Json::arrayValue;    
    // The channel mask maps to consider for saving
    cfg["chanmaskmaps"] = Json::arrayValue;

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

    if (! m_cfg["frame_tags"].isNull()) {
        for (auto jtag : m_cfg["frame_tags"]) {
            const std::string tag = jtag.asString();
            auto traces = FrameTools::tagged_traces(inframe, tag);
            if (traces.empty()) {
                continue;
            }
            auto channels = FrameTools::channels(traces);
            std::sort(channels.begin(), channels.end());
            auto chend = std::unique(channels.begin(), channels.end());
            auto tbinmm = FrameTools::tbin_range(traces);

            // fixme: may want to give user some config over tbin range to save.
            const size_t ncols = tbinmm.second-tbinmm.first;
            const size_t nrows = std::distance(channels.begin(), chend);
            Array::array_xxf arr = Array::array_xxf::Zero(ncols, nrows) + baseline;
            FrameTools::fill(arr, traces, channels.begin(), chend, tbinmm.first);
            arr = arr * scale + offset;

            const std::string arrname = String::format("%s_%d", tag.c_str(), m_save_count);
            if (digitize) {
                Array::array_xxs sarr = arr.cast<short>();
                const short* sdata = sarr.data();
                cnpy::npz_save(fname, arrname, sdata, {nrows, ncols}, mode);
            }
            else {
                cnpy::npz_save(fname, arrname, arr.data(), {nrows, ncols}, mode);
            }
        }

    }
    ++m_save_count;
    return true;
}

