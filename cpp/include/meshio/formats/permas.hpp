#pragma once
//
// PERMAS dat (.post/.dato) reader/writer. Plain text: $COOR node block and
// $ELEMENT TYPE=<permas type> connectivity blocks with the permas<->meshio
// type maps and the write-side node reorders for the second-order elements.
// $NSET/$ESET and other keywords are ignored (meshio drops them too). The
// gzip `.gz` containers are left to the Python fallback.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_permas(const std::string& path, const Mesh& mesh);
Mesh read_permas(const std::string& path);

}  // namespace meshio
