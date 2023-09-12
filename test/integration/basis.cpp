#include <dgt_dg.hpp>
#include <dgt_lua.hpp>
#include <dgt_lua_interface.hpp>

#include <gtest/gtest.h>

using namespace dgt;

static void check_basis(
    Basis<View> const& B,
    int const dim,
    int const p,
    int const q,
    bool const tensor)
{
  EXPECT_EQ(B.dim, dim);
  EXPECT_EQ(B.p, p);
  EXPECT_EQ(B.q, q);
  EXPECT_EQ(B.tensor, tensor);
}

extern "C" {

static int check_basis(lua_State* L)
{
  return lua::function_wrapper(L,
    [] (lua::stack stack) {
      auto input = lua::table(stack.function_argument(1, "dgt_check_basis"));
      Basis<View> basis = make_basis(input);
      int const dim = input.get_integer("dimension");
      int const p = input.get_integer("polynomial_order");
      int const q = input.get_integer("quadrature_rule");
      bool const tensor = input.get_boolean("tensor_product");
      check_basis(basis, dim, p, q, tensor);
      return std::vector<lua::stack_object>({});
    });
}

}

TEST(basis, from_lua)
{
  auto stack = lua::stack::newstate();
  stack.setglobal("dgt_check_basis", check_basis);
  std::filesystem::path const data_dir(std::getenv("INTEGRATION_DATA_DIR"));
  std::filesystem::path const basis_path = data_dir / "basis.lua";
  stack.dofile(basis_path);
}
