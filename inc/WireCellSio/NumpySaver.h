/** Save frames to Numpy array files. */

#ifndef WIRECELLSIO_NUMPYSAVER
#define WIRECELLSIO_NUMPYSAVER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell {
    namespace Sio {
        class NumpySaver : public WireCell::IFrameFilter,
                           public WireCell::IConfigurable {
        public:
            NumpySaver();
            virtual ~NumpySaver();

            /// IFrameFilter
            virtual bool operator()(const WireCell::IFrame::pointer& inframe,
                                    WireCell::IFrame::pointer& outframe);

            /// IConfigurable
            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);
        private:

            Configuration m_cfg;
            int m_save_count;
        };
    }
}
#endif
