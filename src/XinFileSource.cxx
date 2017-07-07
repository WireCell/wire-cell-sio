#include "WireCellSio/XinFileSource.h"
#include "WireCellSio/XinFileIterator.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(XinFileSource, WireCell::Sio::XinFileSource, WireCell::IFrameSource, WireCell::IConfigurable);

using namespace WireCell;

Sio::XinFileSource::XinFileSource()
{
}

Sio::XinFileSource::~XinFileSource()
{
}

void Sio::XinFileSource::configure(const WireCell::Configuration& config)
{
    m_cfg = config;
}

WireCell::Configuration Sio::XinFileSource::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = "";       // actually any URL
    cfg["histtype"] = "raw";
    return cfg;
}


bool Sio::XinFileSource::operator()(IFrame::pointer& out)
{
    std::string url = m_cfg["filename"].asString();
    if (url.empty()) { return false; }
    std::string histtype = m_cfg["histtype"].asString();
    XinFileIterator fs(url.c_str(), histtype.c_str());
    out = fs.frame();
    fs.clear();
    return true;
}


