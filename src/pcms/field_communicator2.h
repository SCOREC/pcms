#ifndef FIELD_COMMUNICATOR2_H_
#define FIELD_COMMUNICATOR2_H_

#include "field_layout_communicator.h"
#include "pcms/field_layout.h"
#include "pcms/field.h"
#include "pcms/profile.h"
#include "pcms/assert.h"
#include "pcms/inclusive_scan.h"
#include "pcms/arrays.h"
#include <redev.h>

namespace pcms
{

template <typename T>
class FieldCommunicator2
{
public:
  FieldCommunicator2(FieldLayoutCommunicator<T>& layout_comm, FieldT<T>& field)
    : comm_buffer_{}, layout_comm_(layout_comm), field_(field)
  {
    PCMS_ALWAYS_ASSERT(&layout_comm.GetLayout() == &field.GetLayout());
    comm_buffer_.resize(layout_comm.GetMsgSize());
  }

  void Send(redev::Mode mode = redev::Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(layout_comm_.GetChannel().InSendCommunicationPhase());
    auto n = field_.Serialize({}, {});
    auto buffer = make_array_view(comm_buffer_);
    field_.Serialize(buffer, layout_comm_.GetPermutationArray());
    layout_comm_.Send(buffer.data_handle(), mode);
  }

  void Receive()
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(layout_comm_.GetChannel().InReceiveCommunicationPhase());
    // Current implementation requires that Receive is always called in Sync
    // mode because we make an immediate call to deserialize after a call to
    // receive.
    auto data = layout_comm_.Recv(redev::Mode::Synchronous);
    field_.Deserialize(make_const_array_view(data),
                       layout_comm_.GetPermutationArray());
  }

private:
  std::vector<T> comm_buffer_;
  FieldLayoutCommunicator<T>& layout_comm_;
  FieldT<T>& field_;
};
} // namespace pcms

#endif // FIELD_COMMUNICATOR2_H_
