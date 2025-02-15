#include <iostream>

#include "Omega_h_adapt.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_shape.hpp"
#include "Omega_h_timer.hpp"
#include "Omega_h_recover.hpp" //recover_hessians
#include <Omega_h_file.hpp> //Omega_h::binary
#include <Omega_h_atomics.hpp> //Omega_h::atomic_fetch_add
#include <sstream> //ostringstream
#include <iomanip> //precision
#include <Omega_h_dbg.hpp>
#include <Omega_h_for.hpp>

//detect floating point exceptions
#include <fenv.h>

using namespace Omega_h;

const Real ICE_DENSITY = 910; //kg/m^3

template <Int dim>
Reals get_eigen_values_dim(Reals metrics) {
  auto n = divide_no_remainder(metrics.size(), symm_ncomps(dim));
  auto isDegen = Write<LO>(n);
  Write<Real> eigenVals(dim*n);
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto m = get_symm<dim>(metrics, i);
    auto m_dc = decompose_eigen(m);
    Few<Int, dim> m_ew_is_degen;
    for (Int j = 0; j < dim; ++j) {
      eigenVals[i*2+j] = m_dc.l[j];
    }
  };
  parallel_for(n, f, "get_eigen_vals");
  return read(eigenVals);
}

Reals get_eigen_values(Mesh* mesh, Reals metrics) {
  if( mesh->dim() == 1)
    return get_eigen_values_dim<1>(metrics); // is this supported?
  else if( mesh->dim() == 2)
    return get_eigen_values_dim<2>(metrics);
  else if( mesh->dim() == 3)
    return get_eigen_values_dim<3>(metrics);
  else
    return Reals(0);
}

template <typename T>
void printTagInfo(Omega_h::Mesh& mesh, std::ostringstream& oss, int dim, int tag, std::string type) {
    auto tagbase = mesh.get_tag(dim, tag);
    auto array = Omega_h::as<T>(tagbase)->array();

    Omega_h::Real min = get_min(array);
    Omega_h::Real max = get_max(array);

    oss << std::setw(18) << std::left << tagbase->name().c_str()
        << std::setw(5) << std::left << dim
        << std::setw(7) << std::left << type
        << std::setw(5) << std::left << tagbase->ncomps()
        << std::setw(10) << std::left << min
        << std::setw(10) << std::left << max
        << "\n";
}

Reals getVolumeFromIceThickness(Omega_h::Mesh& mesh) {
  auto coords = mesh.coords();
  auto triArea = mesh.ask_sizes();
  auto iceThickness = mesh.get_array<Real>(OMEGA_H_VERT, "ice_thickness");
  auto faces2verts = mesh.ask_elem_verts();
  Write<Real> vol(mesh.nfaces());
  auto getVol = OMEGA_H_LAMBDA(LO faceIdx) {
    auto triVerts = Omega_h::gather_verts<3>(faces2verts, faceIdx);
    auto triThickness = Omega_h::gather_scalars<3>(iceThickness, triVerts);
    Real avgThickness = 0;
    for(LO i=0; i<3; i++)
       avgThickness += triThickness[i];
    avgThickness /= 3;
    vol[faceIdx] = triArea[faceIdx]*avgThickness;
  };
  parallel_for(mesh.nfaces(), getVol);
  return vol;
}

Reals getIceMass(Mesh& mesh, Reals vol) {
  return multiply_each_by(vol, ICE_DENSITY);
}

void printTags(Mesh& mesh) {
    std::ostringstream oss;
    // always print two places to the right of the decimal
    // for floating point types (i.e., imbalance)
    oss.precision(2);
    oss << std::fixed;

    if (!mesh.comm()->rank()) {
        oss << "\nTag Properties by Dimension: (Name, Dim, Type, Number of Components, Min. Value, Max. Value)\n";
        for (int dim=0; dim <= mesh.dim(); dim++)
        for (int tag=0; tag < mesh.ntags(dim); tag++) {
            auto tagbase = mesh.get_tag(dim, tag);
            if (tagbase->type() == OMEGA_H_I8)
                printTagInfo<Omega_h::I8>(mesh, oss, dim, tag, "I8");
            if (tagbase->type() == OMEGA_H_I32)
                printTagInfo<Omega_h::I32>(mesh, oss, dim, tag, "I32");
            if (tagbase->type() == OMEGA_H_I64)
                printTagInfo<Omega_h::I64>(mesh, oss, dim, tag, "I64");
            if (tagbase->type() == OMEGA_H_F64)
                printTagInfo<Omega_h::Real>(mesh, oss, dim, tag, "F64");
        }

        std::cout << oss.str();
    }

}

static void check_total_mass(Mesh& mesh) {
  auto vol = getVolumeFromIceThickness(mesh);
  auto masses = getIceMass(mesh, vol);
  auto owned_masses = mesh.owned_array(mesh.dim(), masses, 1);
  auto mass = get_sum(mesh.comm(), owned_masses);
  if (!mesh.comm()->rank()) {
    std::cout << "mass " << mass << '\n';
  }
}

void printTriCount(Mesh* mesh) {
  const auto nTri = mesh->nglobal_ents(2);
  if (!mesh->comm()->rank())
    std::cout << "nTri: " << nTri << "\n";
}

void setupFieldTransfer(AdaptOpts& opts) {
  opts.xfer_opts.type_map["velocity"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["ice_thickness"] = OMEGA_H_LINEAR_INTERP;
  const int numLayers = 11;
  for(int i=1; i<=numLayers; i++) {
    std::stringstream ss;
    ss << "temperature_" << std::setfill('0') << std::setw(2) << i;
    opts.xfer_opts.type_map[ss.str()] = OMEGA_H_LINEAR_INTERP;
  }
}

Reals exponentialVelocity(Mesh& mesh) {
  auto coords = mesh.coords();
  auto coords_h = HostRead<Real>(coords);
  auto vel_h = HostWrite<Real>(mesh.nverts());
  for(int i=0; i<mesh.nverts(); i++) {
    auto x = get_vector<2>(coords_h, i);
    if(std::fabs(x[0]) < 100) {
      vel_h[i] = 20 * std::exp(-0.02 * std::fabs(x[0]));
    } else {
      vel_h[i] = 0;
    }
  }
  return Reals(vel_h);
}

Reals squareVelocity(Mesh& mesh) {
  auto coords = mesh.coords();
  auto coords_h = HostRead<Real>(coords);
  auto vel_h = HostWrite<Real>(mesh.nverts());
  for(int i=0; i<mesh.nverts(); i++) {
    auto x = get_vector<2>(coords_h, i);
    vel_h[i] = 20 * (-1.0e-6 * std::pow(x[0],2));
  }
  return Reals(vel_h);
}

Reals normSquaredVelocity(Mesh& mesh) {
  auto coords = mesh.coords();
  auto coords_h = HostRead<Real>(coords);
  auto vel_h = HostWrite<Real>(mesh.nverts());
  for(int i=0; i<mesh.nverts(); i++) {
    auto x = get_vector<2>(coords_h, i);
    vel_h[i] = norm_squared(x);
  }
  return Reals(vel_h);
}

Reals magnitudeOfVelocity(Mesh& mesh) {
  auto vel_x = mesh.get_array<Real>(VERT, "solution_1");
  auto vel_y = mesh.get_array<Real>(VERT, "solution_2");
  auto vel_mag = Write<Real>(mesh.nverts());
  auto f = OMEGA_H_LAMBDA(LO i) {
    Vector<2> v{vel_x[i], vel_y[i]};
    vel_mag[i] = norm(v);
  };
  parallel_for(mesh.nverts(), f, "magnitude_of_velocity");
  return Read(vel_mag);
}

int main(int argc, char** argv) {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);  // Enable all floating point exceptions but FE_INEXACT
  auto lib = Library(&argc, &argv);
  if( argc != 3 ) {
    fprintf(stderr, "Usage: %s inputMesh.osh outputMeshPrefix\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  auto world = lib.world();
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[1], world, &mesh);

  Omega_h::vtk::write_parallel("beforeClassFix_edges.vtk", &mesh, 1);
  classify_by_angles(&mesh,60*(Omega_h::PI/180));

  //create analytic velocity field
  auto velocity = magnitudeOfVelocity(mesh);
  mesh.add_tag(VERT, "velocity", 1, Reals(velocity));
  const auto minSizeThreshold = 0.5;
  auto max_velocity = get_max(velocity)*minSizeThreshold;
  auto norm_velocity = divide_each_by(velocity, max_velocity);
  mesh.add_tag(VERT, "norm_velocity", 1, Reals(norm_velocity));
  const auto minScaleFactor = 0.25;
  const auto maxScaleFactor = 4.0;
  std::cout << "minScaleFacor " << minScaleFactor << " " <<
               "maxScaleFacor " << maxScaleFactor << "\n";
  auto scale = Write<Real>(mesh.nverts());
  auto set_metric_scaling = OMEGA_H_LAMBDA(LO i) {
    const auto a = 3.85939018034835;
    const auto b = 20.0;
    const auto d = 0.25;
    const auto x = norm_velocity[i];
    scale[i] = a * std::log(b*x+1) + d;
    if( x > minSizeThreshold ) scale[i] = 12;
  };
  parallel_for(mesh.nverts(), set_metric_scaling, "set_metric_scaling");
  auto scale_r = Read(scale);
  mesh.add_tag(VERT, "scale", 1, scale_r);

  Omega_h::add_implied_metric_tag(&mesh);
  auto metric = mesh.get_array<Real>(VERT, "metric");
  auto tgt_metric = multiply_each(metric, scale_r);
  mesh.add_tag(VERT, "tgt_metric", 3, Reals(tgt_metric));

  auto genopts = Omega_h::MetricInput();
  auto metric_scaling = 1.0; //1.0: no scaling
  genopts.sources.push_back(
      Omega_h::MetricSource{OMEGA_H_GIVEN, metric_scaling, "tgt_metric"});
  genopts.nsmoothing_steps = 2;
  Omega_h::generate_target_metric_tag(&mesh, genopts);

  mesh.set_parting(OMEGA_H_GHOSTED);

  //DEBUG { 
  auto edge_lengths_source = measure_edges_metric(&mesh, mesh.get_array<Real>(VERT, "metric"));
  mesh.add_tag(EDGE, "lengths_source", 1, edge_lengths_source);
  auto edge_lengths_target = measure_edges_metric(&mesh, mesh.get_array<Real>(VERT, "target_metric"));
  mesh.add_tag(EDGE, "lengths_target", 1, edge_lengths_target);
  auto elen = measure_edges_real(&mesh);
  mesh.add_tag(EDGE, "lengths_real", 1, elen);
  Omega_h::vtk::write_parallel("beforeAdapt_edges.vtk", &mesh, 1);
  // }

  Omega_h::AdaptOpts opts(&mesh);
  setupFieldTransfer(opts);

  printTriCount(&mesh);
  printTags(mesh);
  Omega_h::vtk::write_parallel("beforeAdapt.vtk", &mesh, 2);
  check_total_mass(mesh);

  while (Omega_h::approach_metric(&mesh, opts)) {
    adapt(&mesh, opts);
  }

  auto post_elen = measure_edges_real(&mesh);
  mesh.add_tag(EDGE, "post_lengths_real", 1, post_elen);

  printTriCount(&mesh);
  printTags(mesh);
  check_total_mass(mesh);
  const std::string vtkFileName = std::string(argv[2]) + ".vtk";
  Omega_h::vtk::write_parallel(vtkFileName, &mesh, 2);
  const std::string vtkFileName_edges = std::string(argv[2]) + "_edges.vtk";
  Omega_h::vtk::write_parallel(vtkFileName_edges, &mesh, 1);
  return 0;
}
