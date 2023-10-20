#ifndef PCMS_COUPLING_CLIENT_H
#define PCMS_COUPLING_CLIENT_H
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pcms/common.h"
#include "pcms/field_communicator.h"
#include "pcms/profile.h"
#include "../client.h"


namespace py = pybind11;

PYBIND11_MODULE(pcms_client, m) {
  m.doc() = "CMS Coupling Client C++ Pybind11 Bindings";

  py::enum_<Mode>(m, "Mode")
        .value("Deferred", Mode::Deferred);
        .value("Synchronous", Mode::Synchronous);


  py::class_<CoupledField>(m, "CoupledField")
        .def(py::init<const std::string&, FieldAdapterT, MPI_Comm, redev::Redev&, redev::Channel&, bool>());
        .def("Send", &CoupledField::Send, py::arg("mode") = Mode::Synchronous);
        .def("Receive", &CoupledField::Receive);

  py::class_<CouplerClient>(m, "CouplerClient")
        .def(py::init<std::string, MPI_Comm, redev::TransportType, adios2::Params, std::string>());
        .def("GetPartition", &CouplerClient::GetPartition, py::return_value_policy::reference); //redev::Partition&
        .def("AddField", &CouplerClient::AddField, py::arg("name"), py::arg("field_adapter"), py::arg("participates") = true, 
        py::return_value_policy::reference); //CoupledField*
        .def("SendField", &CouplerClient::SendField, const std::string&, py::arg("mode") = Mode::Synchronous); //void
        .def("ReceiveField", &CouplerClient::ReceiveField, const std::string&); //void 
        .def("InSendPhase", &CouplerClient::InSendPhase); //bool
        .def("InReceivePhase", &CouplerClient::InReceivePhase); //bool
        .def("BeginSendPhase", &CouplerClient::BeginSendPhase); //void
        .def("EndSendPhase", &CouplerClient::EndSendPhase); //void
        .def("BeginReceivePhase", &CouplerClient::BeginReceivePhase); //void
        .def("EndReceivePhase", &CouplerClient::EndReceivePhase); //void
}
 // namespace pcms

#endif // PCMS_COUPLING_CLIENT_H
