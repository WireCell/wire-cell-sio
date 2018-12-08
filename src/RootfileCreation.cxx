#include "WireCellSio/RootfileCreation.h"

#include "TFile.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(RootfileCreation_depos, WireCell::Sio::RootfileCreation_depos,
                 WireCell::IDepoFilter,  WireCell::IConfigurable)

WIRECELL_FACTORY(RootfileCreation_frames, WireCell::Sio::RootfileCreation_frames,
                 WireCell::IDepoFilter,  WireCell::IConfigurable)


using namespace WireCell;

Sio::RootfileCreation_depos::RootfileCreation_depos(){
}

Sio::RootfileCreation_depos::~RootfileCreation_depos(){
}


void Sio::RootfileCreation_depos::configure(const WireCell::Configuration& cfg)
{
  m_cfg = cfg;
}

WireCell::Configuration Sio::RootfileCreation_depos::default_configuration() const
{
  Configuration cfg;
  cfg["output_filename"] = "";
  cfg["root_file_mode"] = "RECREATE";
  return cfg;
}

bool Sio::RootfileCreation_depos::operator()(const WireCell::IDepo::pointer& indepo,
				       WireCell::IDepo::pointer& outdepo)
{
  outdepo = indepo;
  create_file();
  return true;
}

void Sio::RootfileCreation_depos::create_file(){
  const std::string ofname = m_cfg["output_filename"].asString();
  const std::string mode = m_cfg["root_file_mode"].asString();
  TFile* output_tf = TFile::Open(ofname.c_str(), mode.c_str());
  output_tf->Close("R");
  delete output_tf;
  output_tf = nullptr;
}


Sio::RootfileCreation_frames::RootfileCreation_frames(){
}

Sio::RootfileCreation_frames::~RootfileCreation_frames(){
}


void Sio::RootfileCreation_frames::configure(const WireCell::Configuration& cfg)
{
  m_cfg = cfg;
}

WireCell::Configuration Sio::RootfileCreation_frames::default_configuration() const
{
  Configuration cfg;
  cfg["output_filename"] = "";
  cfg["root_file_mode"] = "RECREATE";
  return cfg;
}

void Sio::RootfileCreation_frames::create_file(){
  const std::string ofname = m_cfg["output_filename"].asString();
  const std::string mode = m_cfg["root_file_mode"].asString();

  TFile* output_tf = TFile::Open(ofname.c_str(), mode.c_str()); 
  output_tf->Close("R");
  delete output_tf;
  output_tf = nullptr;
}


bool Sio::RootfileCreation_frames::operator()(const WireCell::IFrame::pointer& in, WireCell::IFrame::pointer& out){
  out = in;
  create_file();
  return true;
}
 
