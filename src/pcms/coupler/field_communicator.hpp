#ifndef FIELD_COMMUNICATOR2_H_
#define FIELD_COMMUNICATOR2_H_

#include "pcms/coupler/field_layout_communicator.h"
#include "pcms/field/field.h"
#include "pcms/field/field_layout.h"
#include "pcms/coupler/field_serializer.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/inclusive_scan.h"
#include "pcms/utility/arrays.h"
#include <redev.h>
#include <memory>

namespace pcms
{

template <typename T>
class FieldCommunicator
{
public:
  FieldCommunicator(const std::string& name,
                    FieldLayoutCommunicator& layout_comm, Field<T>& field)
    : FieldCommunicator(name, layout_comm, field,
                        std::make_unique<FieldSerializer<T>>())
  {
  }

  FieldCommunicator(const std::string& name,
                    FieldLayoutCommunicator& layout_comm, Field<T>& field,
                    std::unique_ptr<FieldSerializer<T>> serializer)
    : comm_buffer_{},
      layout_comm_(layout_comm),
      field_(field),
      serializer_(std::move(serializer))
  {
    PCMS_ALWAYS_ASSERT(&layout_comm.GetLayout() == &field.GetLayout());
    // The exchange plan is per DOF holder; the field payload carries
    // num_components values per holder, so the buffer is scaled accordingly.
    comm_buffer_.resize(
      layout_comm.GetMsgSize() *
      static_cast<size_t>(field.GetLayout().GetNumComponents()));
    comm_ =
      layout_comm_.GetChannel().CreateComm<T>(name, layout_comm.GetMPIComm());
    layout_comm_.SetOutMessageLayout<T>(comm_);
  }

  void Send(redev::Mode mode = redev::Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(layout_comm_.GetChannel().InSendCommunicationPhase());
    auto buffer = make_array_view(comm_buffer_);
    serializer_->Serialize(field_.GetData(), field_.GetLayout(), buffer,
                           layout_comm_.GetPermutationArray());
    comm_.Send(buffer.data_handle(), mode);
  }

  void Receive()
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(layout_comm_.GetChannel().InReceiveCommunicationPhase());
    // Current implementation requires that Receive is always called in Sync
    // mode because we make an immediate call to deserialize after a call to
    // receive.
    auto data = comm_.Recv(redev::Mode::Synchronous);
    serializer_->Deserialize(field_.GetData(), field_.GetLayout(),
                             make_const_array_view(data),
                             layout_comm_.GetPermutationArray());
  }

private:
  std::vector<T> comm_buffer_;
  redev::BidirectionalComm<T> comm_;
  FieldLayoutCommunicator& layout_comm_;
  Field<T>& field_;
  std::unique_ptr<FieldSerializer<T>> serializer_;
};

template <typename T>
using FieldCommunicator2PtrT = std::unique_ptr<FieldCommunicator<T>>;

using FieldCommunicator2Ptr =
  std::variant<FieldCommunicator2PtrT<int8_t>, FieldCommunicator2PtrT<int32_t>,
               FieldCommunicator2PtrT<int64_t>, FieldCommunicator2PtrT<float>,
               FieldCommunicator2PtrT<double>>;
} // namespace pcms

#endif // FIELD_COMMUNICATOR2_H_
