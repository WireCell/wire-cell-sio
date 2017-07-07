#include "WireCellSio/MagnifySource.h"
#include "WireCellSio/MagnifyFileIterator.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(MagnifySource, WireCell::Sio::MagnifySource, WireCell::IFrameSource, WireCell::IConfigurable);

using namespace WireCell;

Sio::MagnifySource::MagnifySource()
{
}

Sio::MagnifySource::~MagnifySource()
{
}

void Sio::MagnifySource::configure(const WireCell::Configuration& config)
{
    m_cfg = config;
}

WireCell::Configuration Sio::MagnifySource::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = "";       // actually any URL
    cfg["histtype"] = "raw";
    return cfg;
}


bool Sio::MagnifySource::operator()(IFrame::pointer& out)
{
    std::string url = m_cfg["filename"].asString();
    if (url.empty()) { return false; }
    std::string histtype = m_cfg["histtype"].asString();
    MagnifyFileIterator fs(url.c_str(), histtype.c_str());
    out = fs.frame();
    fs.clear();
    return true;
}


