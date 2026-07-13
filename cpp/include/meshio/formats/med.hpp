#pragma once
//
// MED/Salome (.med) reader/writer, mirroring meshio's h5py dialect:
// ENS_MAA/<mesh> (with a -..1-..1 time-step group), NOE/COO + FAM,
// MAI/<MED type>/NOD + FAM (Fortran-order flattened, 1-based), FAS families
// (point/cell tags), and CHA fields (NOEU nodal, ELEM cell, ELNO
// nodes-per-cell; ELGA and non-default profiles are deferred to Python).
//
// point_tags / cell_tags are custom attributes on the Python Mesh and
// field_data["med:nom"] is a list of string-lists, so they travel through
// the MedInfo side-channel rather than the Mesh conversion layer.
// Available only when built with MESHIO_HAS_HDF5.

#ifdef MESHIO_HAS_HDF5

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "meshio/mesh.hpp"

namespace meshio {

struct MedInfo {
    std::map<std::int64_t, std::vector<std::string>> point_tags;
    std::map<std::int64_t, std::vector<std::string>> cell_tags;
    std::vector<std::vector<std::string>> med_nom;  // field_data["med:nom"]
};

Mesh read_med(const std::string& path, MedInfo& info);
void write_med(const std::string& path, const Mesh& mesh, const MedInfo& info);

}  // namespace meshio

#endif  // MESHIO_HAS_HDF5
