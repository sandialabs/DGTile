#pragma once

// This code was originally developed by:
// Dan Ibanez, Sandia National Laboratories

#include <filesystem>
#include <iterator>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include <fmt/core.h>

#include <mpi.h>

#include <lua.hpp>

namespace dgt {
namespace lua {

class stack;

class stack_object {
  lua::stack* m_stack_pointer{nullptr};
  int m_index{std::numeric_limits<int>::max()};
  std::string m_name;
 public:
  stack_object() = default;
  stack_object(lua::stack* stack_pointer_arg, int index_arg, std::string&& name_arg);
  stack_object(stack_object&&);
  stack_object& operator=(stack_object&&);
  stack_object(stack_object const&);
  stack_object& operator=(stack_object const&);
  ~stack_object();
  lua::stack*& stack() { return m_stack_pointer; }
  lua::stack* stack() const { return m_stack_pointer; }
  lua_State* state() const;
  int index() const { return m_index; }
  int& index() { return m_index; }
  std::string const& name() const { return m_name; }
  int type() const;
  char const* type_name() const;
  void check_type(int expected_type) const;
};

class boolean : public stack_object {
 public:
  boolean(stack_object&& so_arg);
  bool value() const;
};

class function : public stack_object {
 public:
  function() = default;
  function(stack_object const& so_arg);
  template <typename ... Args>
  std::vector<stack_object> operator()(Args&&... args) const
  {
    int reference_point = push_self();
    push_args(std::forward<Args>(args)...);
    return raw_call(reference_point);
  }
  template <typename ... Args>
  stack_object first_result(Args&&... args) const
  {
    auto results = operator()(std::forward<Args>(args)...);
    if (results.size() < 1) {
      throw std::runtime_error(fmt::format("Lua function {} called from C++ should have returned at least one result but returned none", this->name()));
    }
    return std::move(results[0]);
  }
 private:
  int push_self() const;
  void push_args() const {};
  template <typename ... Args>
  void push_args(stack_object const& arg, Args&&... args) const
  {
    push_stack_object(arg);
    push_args(std::forward<Args>(args)...);
  }
  template <typename ... Args>
  void push_args(lua_Number arg, Args&&... args) const
  {
    push_number(arg);
    push_args(std::forward<Args>(args)...);
  }
  template <typename ... Args>
  void push_args(lua_Integer arg, Args&&... args) const
  {
    push_integer(arg);
    push_args(std::forward<Args>(args)...);
  }
  template <typename ... Args>
  void push_args(std::string const& arg, Args&&... args) const
  {
    push_string(arg);
    push_args(std::forward<Args>(args)...);
  }
  template <typename ... Args>
  void push_args(char const* arg, Args&&... args) const
  {
    push_string(arg);
    push_args(std::forward<Args>(args)...);
  }
  void push_stack_object(stack_object const& so) const;
  void push_number(lua_Number value) const;
  void push_integer(lua_Integer value) const;
  void push_string(std::string const& value) const;
  std::vector<stack_object> raw_call(int reference_point) const;
};

class integer : public stack_object {
 public:
  integer(stack_object&& so_arg);
  lua_Integer value() const;
};

class number : public stack_object {
 public:
  number();
  number(stack_object const& so_arg);
  lua_Number value() const;
};

class string : public stack_object {
 public:
  string(stack_object const& so_arg);
  std::string value() const;
};

class table_iterator {
  stack_object m_table;
  stack_object m_key;
  bool m_is_end{false};
 public:
  using iterator_category = std::input_iterator_tag;
  using value_type = std::pair<stack_object, stack_object>;
  using difference_type = int;
  using pointer = void*;
  using reference = value_type;
  table_iterator() = default;
  table_iterator(table_iterator&&) = default;
  table_iterator& operator=(table_iterator&&) = default;
  table_iterator(table_iterator const&) = default;
  table_iterator& operator=(table_iterator const&) = default;
  table_iterator(
      stack_object table_arg,
      stack_object key_arg,
      bool is_end_arg);
  bool operator==(table_iterator const& other) const;
  bool operator!=(table_iterator const& other) const;
  reference operator*() const;
  pointer operator->() const { return nullptr; }
  table_iterator& operator++();
  table_iterator operator++(int)
  {
    auto result = *this;
    operator++();
    return result;
  }
};

class table : public stack_object {
 public:
  table() = default;
  table(stack_object const& so_arg);
  std::optional<stack_object> get_optional(char const* key) const;
  std::optional<stack_object> get_optional(lua_Integer key) const;
  bool has(char const* key) const;
  bool has(lua_Integer key) const;
  stack_object get(char const* key) const;
  stack_object get(lua_Integer key) const;
  table get_table(char const* key) const;
  table get_table(lua_Integer key) const;
  lua_Number get_number(char const* key) const;
  lua_Number get_number(lua_Integer key) const;
  lua_Integer get_integer(char const* key) const;
  lua_Integer get_integer(lua_Integer key) const;
  std::string get_string(char const* key) const;
  std::string get_string(lua_Integer key) const;
  bool get_boolean(char const* key) const;
  bool get_boolean(lua_Integer key) const;
  double get_or(char const* key, double value) const;
  double get_or(lua_Integer key, double value) const;
  int get_or(char const* key, int value) const;
  int get_or(lua_Integer key, int value) const;
  std::string get_or(char const* key, std::string const& value) const;
  std::string get_or(lua_Integer key, std::string const& value) const;
  std::string get_or(char const* key, char const* value) const;
  std::string get_or(lua_Integer key, char const* value) const;
  bool get_or(char const* key, bool value) const;
  bool get_or(lua_Integer key, bool value) const;
  template <class T>
  void replace_if(char const* key, T& value) const
  {
    value = get_or(key, value);
  }
  table get_or_table(char const* key) const;
 private:
  void set_stack_object(char const* key, stack_object const& so);
  void set_stack_object(lua_Integer key, stack_object const& so);
 public:
  void set_number(char const* key, lua_Number value);
  void set_number(lua_Integer key, lua_Number value);
  void set_integer(char const* key, lua_Integer value);
  void set_string(char const* key, std::string const& value);
  void set_boolean(char const* key, bool value);
  void set_cfunction(char const* key, lua_CFunction value);
  void set(char const* key, lua_CFunction value);
  void set(char const* key, char const* value);
  void set(char const* key, int value);
  void set(char const* key, double value);
  void set(char const* key, lua::table const& value);
  lua_Integer size() const;
  lua_Integer num_entries() const;
  table set_table(char const* key);
  table set_table(lua_Integer key);
  table_iterator begin() const;
  table_iterator end() const;
};

class stack {
  lua_State* m_state;
  bool m_owns_state{false};
  std::vector<stack_object*> m_references;
 public:
  stack(lua_State* state_arg);
  stack(stack&&);
  stack& operator=(stack&&);
  stack(stack const&) = delete;
  stack& operator=(stack const&) = delete;
  ~stack();
  lua_State* state() const { return m_state; }
  void add_reference(stack_object* reference);
  void remove_reference(stack_object* reference);
  stack_object function_argument(int i, std::string const& function_name);
  lua::table table(std::string&& name_arg);
  lua::stack_object nil();
  lua::number number(lua_Number value, std::string&& name_arg);
  lua::integer integer(lua_Integer value, std::string&& name_arg);
  lua::string string(std::string const& value, std::string&& name_arg);
  lua::boolean boolean(bool value, std::string&& name_arg);
  lua::function cfunction(lua_CFunction value, std::string&& name_arg);
  stack_object push(stack_object const& other);
  // call this before any function that pushes values onto the Lua stack
  void prepush() const;
  // call this on a stack object that references a stack index that is about
  // to be consumed (popped/removed) by Lua C API functions
  void release(stack_object&& so);
  static stack newstate();
  void setglobal(char const* key, stack_object const& so);
  void setglobal(char const* key, lua_CFunction value);
  void dofile(std::filesystem::path const& path);
  stack_object getglobal(char const* key);
  lua::table require(char const* key);
};

template <class Functor>
int function_wrapper(lua_State* L, Functor const& f)
{
  std::vector<lua::stack_object> results;
  try {
    results = f(lua::stack(L));
  } catch (std::exception const& e) {
    luaL_error(L, "\nLua detected an error in wrapped C++ code:\n%s", e.what()); 
  }
  return int(results.size());
}

}
}
