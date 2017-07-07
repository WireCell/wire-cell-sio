/** A "Magnify" file is one used for viewing in the Magnify display. */

#ifndef WIRECELLSIO_MAGNIFYFILESOURCE
#define WIRECELLSIO_MAGNIFYFILESOURCE

#include "WireCellIface/IFrameSource.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell {
    namespace Sio {

        class MagnifySource : public IFrameSource, public IConfigurable {
        public:

            MagnifySource();
            virtual ~MagnifySource();

            /// IFrameSource
            virtual bool operator()(IFrame::pointer& out);

            /// IConfigurable
            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);
        private:
            Configuration m_cfg;

        };
    }
}

#endif
