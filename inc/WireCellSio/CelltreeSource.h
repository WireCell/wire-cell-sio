/* Celltree input for ( noise filter + signal processing ) */
#ifndef WIRECELLSIO_CELLTREEFILESOURCE
#define WIRECELLSIO_CELLTREEFILESOURCE

#include "WireCellIface/IFrameSource.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell {
  namespace Sio {
    class CelltreeSource : public IFrameSource, public IConfigurable {
    public:
      CelltreeSource();
      virtual ~CelltreeSource();

      virtual bool operator()(IFrame::pointer& out);

      virtual WireCell::Configuration default_configuration() const;
      virtual void configure(const WireCell::Configuration& config);
      
    private:
      Configuration m_cfg;
      
    };
  }
}


#endif
