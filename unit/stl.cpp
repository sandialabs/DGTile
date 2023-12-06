#include <dgt_stl.hpp>

#include <gtest/gtest.h>

#include <dgt_print.hpp> // debug

using namespace dgt;
using namespace dgt::stl;

TEST(stl, read)
{
  std::filesystem::path const data_dir(std::getenv("UNIT_DATA_DIR"));
  std::filesystem::path const stl_path = data_dir / "ex.stl";
  HostView<Triangle*> triangles = read(stl_path);
  EXPECT_EQ(triangles.size(), 2);
  // normal of first triangle
  EXPECT_EQ(triangles[0][0].x(), 1.);
  EXPECT_EQ(triangles[0][0].y(), 0.);
  EXPECT_EQ(triangles[0][0].z(), 0.);
  // first vertex
  EXPECT_EQ(triangles[0][1].x(), 0.);
  EXPECT_EQ(triangles[0][1].y(), 0.);
  EXPECT_EQ(triangles[0][1].z(), 0.);
  // second vertex
  EXPECT_EQ(triangles[0][2].x(), 0.);
  // stl stores floats in its binary format
  EXPECT_NEAR(triangles[0][2].y(), 1.23, 1.e-7);
  EXPECT_EQ(triangles[0][2].z(), 0.);
  // third vertex
  EXPECT_EQ(triangles[0][3].x(), 0.);
  EXPECT_EQ(triangles[0][3].y(), 0.);
  EXPECT_NEAR(triangles[0][3].z(), 1.23, 1.e-7);
  // normal of second triangle
  EXPECT_EQ(triangles[1][0].x(), -1.);
  EXPECT_EQ(triangles[1][0].y(), 0.);
  EXPECT_EQ(triangles[1][0].z(), 0.);
  // first vertex
  EXPECT_EQ(triangles[1][1].x(), 0.);
  EXPECT_EQ(triangles[1][1].y(), 0.);
  EXPECT_EQ(triangles[1][1].z(), 0.);
  // second vertex
  EXPECT_EQ(triangles[1][2].x(), 0.);
  EXPECT_EQ(triangles[1][2].y(), 0.);
  EXPECT_NEAR(triangles[1][2].z(), 1.23, 1.e-7);
  // third vertex
  EXPECT_EQ(triangles[1][3].x(), 0.);
  EXPECT_NEAR(triangles[1][3].y(), 1.23, 1.e-7);
  EXPECT_EQ(triangles[1][3].z(), 0.);
}

TEST(stl, bounding_box)
{
  std::filesystem::path const data_dir(std::getenv("UNIT_DATA_DIR"));
  std::filesystem::path const stl_path = data_dir / "ex.stl";
  HostView<Triangle*> const triangles = read(stl_path);
  Box3<real> const bounding_box = compute_bounding_box(triangles);

  std::cout << bounding_box << "\n";

}

TEST(stl, to_device)
{
  std::filesystem::path const data_dir(std::getenv("UNIT_DATA_DIR"));
  std::filesystem::path const stl_path = data_dir / "ex.stl";
  HostView<Triangle*> host_triangles = read(stl_path);
  View<Triangle*> device_triangles = to_device(host_triangles);
  EXPECT_EQ(device_triangles.size(), 2);
}
