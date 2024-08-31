#include <pcms.h>
#include <fusion_io.h>

#pragma once

namespace fusion_io
{
  class Library
  {
    public:
    MPI_Comm comm;
    pcms::CouplerClient* coupler;
    pcms::Application* application;

    Library(int argc, char** argv, std::string name) {
      MPI_Init(&argc, &argv);
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);
      coupler = new pcms::CouplerClient(name, comm);
    }

    ~Library() {
      MPI_Finalize();
      delete coupler;
    }
  };

  int fio_open_source(fio_source** src, const int type, const char* filename, Library& library) {
    int ierr;
    *src = 0;

    switch(type) {
    case(FIO_GATO_SOURCE):
      *src = new gato_source();
      ierr = (*src)->open(filename);
      break;

    case(FIO_GEQDSK_SOURCE):
      *src = new geqdsk_source();
      ierr = (*src)->open(filename);
      break;

    case(FIO_GPEC_SOURCE):
      *src = new gpec_source();
      ierr = (*src)->open(filename);
      break;

    case(FIO_M3DC1_SOURCE):
      *src = new m3dc1_source();
      ierr = (*src)->open(filename);
      break;

    case(FIO_MARS_SOURCE):
      *src = new mars_source();
      ierr = (*src)->open(filename);
      break;

    default:
      std::cerr << "Source type " << type << " unsupported." << std::endl;
      return FIO_UNSUPPORTED;
    }

    if(ierr != FIO_SUCCESS) {
      delete(*src);
      return ierr;
    }
    
    return FIO_SUCCESS;
  };

  // int gato_source::get_field(const field_type t,fio_field** f, const fio_option_list* opt, Library library);

  // int m3dc1_scalar_field::eval(const double* x, double* v, fio_hint s, Library library);
};
