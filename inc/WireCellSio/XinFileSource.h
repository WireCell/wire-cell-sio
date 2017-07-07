/** A "Xin" file is one used in the WC prototype and some large tests in the toolkit. */

#ifndef WIRECELLSIO_XINFILESOURCE
#define WIRECELLSIO_XINFILESOURCE

#include "WireCellIface/IFrameSource.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell {
    namespace Sio {

        class XinFileSource : public IFrameSource, public IConfigurable {
        public:

            XinFileSource();
            virtual ~XinFileSource();

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
