#include "dgt_field.hpp"

namespace dgt {

static void verify_info(FieldInfo& info) {
  if (info.name == "") {
    throw std::runtime_error("Field - unset name");
  }
  if (info.ent_dim < 0) {
    throw std::runtime_error("Field - invalid ent dim");
  }
  if (info.ncomps <= 0) {
    throw std::runtime_error("Field - invalid ncomps");
  }
}

static void verify_ent_dim(int ent_dim, int mesh_dim) {
  if (ent_dim != mesh_dim-1) {
    throw std::runtime_error("Field - only side fields supported");
  }
}

void Field::set_info(FieldInfo const& info) {
  m_info = info;
}

std::string Field::name() const {
  return m_info.name;
}

int Field::ncomps() const {
  return m_info.ncomps;
}

int Field::ent_dim() const {
  return m_info.ent_dim;
}

View<double**> Field::data(int axis) const {
  return m_data[axis];
}

static std::string view_name(FieldInfo const& info, int axis) {
  return "dgt::Field::m_data[" +
    std::to_string(axis) +
    "].(" + info.name + ")";
}

void Field::allocate(p3a::grid3 const& cell_grid) {
  verify_info(m_info);
  int const mesh_dim = get_dim(cell_grid);
  int const ent_dim = m_info.ent_dim;
  int const ncomps = m_info.ncomps;
  verify_ent_dim(ent_dim, mesh_dim);
  p3a::grid3 const cgrid = generalize(cell_grid);
  for (int axis = 0; axis < mesh_dim; ++axis) {
    int const nsides = get_side_grid(cgrid, axis).size();
    m_data[axis] = View<double**>(view_name(m_info, axis), nsides, ncomps);
  }
}

void Field::deallocate() {
  for (int axis = 0; axis < DIMS; ++axis) {
    m_data[axis] = View<double**>();
  }
}

}
