#include <algorithm>

#include "dgt_lua.hpp"

namespace dgt {
namespace lua {

stack_object::stack_object(
      lua::stack* stack_pointer_arg,
      int index_arg,
      std::string&& name_arg)
  :m_stack_pointer(stack_pointer_arg)
  ,m_index(index_arg)
  ,m_name(std::move(name_arg))
{
  if (m_stack_pointer) m_stack_pointer->add_reference(this);
}

stack_object::stack_object(stack_object&& other)
  :m_stack_pointer(other.m_stack_pointer)
  ,m_index(other.m_index)
  ,m_name(std::move(other.m_name))
{
  if (m_stack_pointer) m_stack_pointer->add_reference(this);
  if (other.m_stack_pointer) other.m_stack_pointer->remove_reference(&other);
  other.m_stack_pointer = nullptr;
}

stack_object& stack_object::operator=(stack_object&& other)
{
  if (m_stack_pointer) {
    m_stack_pointer->remove_reference(this);
  }
  m_stack_pointer = other.m_stack_pointer;
  m_index = other.m_index;
  m_name = std::move(other.m_name);
  if (m_stack_pointer) m_stack_pointer->add_reference(this);
  if (other.m_stack_pointer) other.m_stack_pointer->remove_reference(&other);
  other.m_stack_pointer = nullptr;
  return *this;
}

stack_object::stack_object(stack_object const& other)
  :m_stack_pointer(other.m_stack_pointer)
  ,m_index(other.m_index)
  ,m_name(other.m_name)
{
  if (m_stack_pointer) m_stack_pointer->add_reference(this);
}

stack_object& stack_object::operator=(stack_object const& other)
{
  if (m_stack_pointer) {
    m_stack_pointer->remove_reference(this);
  }
  m_stack_pointer = other.m_stack_pointer;
  m_index = other.m_index;
  m_name = other.m_name;
  if (m_stack_pointer) m_stack_pointer->add_reference(this);
  return *this;
}

stack_object::~stack_object()
{
  if (m_stack_pointer) {
    m_stack_pointer->remove_reference(this);
  }
}

lua_State* stack_object::state() const
{
  return m_stack_pointer->state();
}

int stack_object::type() const
{
  return lua_type(state(), index());
}

char const* stack_object::type_name() const
{
  return lua_typename(state(), type());
}

void stack_object::check_type(int expected_type) const
{
  int got_type = type();
  if (got_type == expected_type) return;
  throw std::runtime_error(fmt::format("expected {} to be a {} but it is a {}",
        name(),
        lua_typename(state(), expected_type),
        lua_typename(state(), got_type)));
}

boolean::boolean(stack_object&& so_arg)
  :stack_object(std::move(so_arg))
{
  check_type(LUA_TBOOLEAN);
}

bool boolean::value() const
{
  return lua_toboolean(state(), index()) == 1;
}

function::function(stack_object const& so_arg)
  :stack_object(so_arg)
{
  check_type(LUA_TFUNCTION);
}

int function::push_self() const
{
  stack()->prepush();
  lua_pushvalue(state(), index());
  return lua_gettop(state());
}

void function::push_stack_object(stack_object const& so) const
{
  stack()->prepush();
  lua_pushvalue(state(), so.index());
}

void function::push_number(lua_Number value) const
{
  stack()->prepush();
  lua_pushnumber(state(), value);
}

void function::push_integer(lua_Integer value) const
{
  stack()->prepush();
  lua_pushinteger(state(), value);
}

void function::push_string(std::string const& value) const
{
  stack()->prepush();
  lua_pushstring(state(), value.c_str());
}

std::vector<stack_object> function::raw_call(int reference_point) const
{
  int nargs = lua_gettop(state()) - reference_point;
  int result = lua_pcall(state(), nargs, LUA_MULTRET, 0);
  if (result != LUA_OK) {
    if (lua_gettop(state()) != reference_point) {
      throw std::logic_error(fmt::format(
            "lua_pcall was not LUA_OK but there is not one error object on the stack, instead gettop is {} and reference_point is {}",
            lua_gettop(state()), reference_point));
    }
    std::string error_string = lua_tostring(state(), -1);
    lua_pop(state(), 1);
    throw std::runtime_error(fmt::format(
        "calling Lua function {} produced error: {}",
        name(), error_string));
  }
  std::vector<stack_object> results;
  for (int i = reference_point; i <= lua_gettop(state()); ++i) {
    results.push_back(stack_object(stack(), i,
          fmt::format("(result {} of {})", results.size() + 1, name())));
  }
  return results;
}

integer::integer(stack_object&& so_arg)
  :stack_object(std::move(so_arg))
{
  check_type(LUA_TNUMBER);
  if (!lua_isinteger(state(), index())) {
    throw std::runtime_error(fmt::format("expected number {} to be an integer but it is not", name()));
  }
}

lua_Integer integer::value() const
{
  return lua_tointeger(state(), index());
}

number::number()
  :stack_object()
{
}

number::number(stack_object const& so_arg)
  :stack_object(so_arg)
{
  check_type(LUA_TNUMBER);
}

lua_Number number::value() const
{
  return lua_tonumber(state(), index());
}

string::string(stack_object const& so_arg)
  :stack_object(so_arg)
{
  check_type(LUA_TSTRING);
}

std::string string::value() const
{
  return lua_tostring(state(), index());
}

table_iterator::table_iterator(
    stack_object table_arg,
    stack_object key_arg,
    bool is_end_arg)
  :m_table(table_arg)
  ,m_key(key_arg)
  ,m_is_end(is_end_arg)
{
}

bool table_iterator::operator==(table_iterator const& other) const
{
  if (m_is_end || other.m_is_end) {
    return m_is_end == other.m_is_end;
  } else {
    return 1 == lua_compare(m_key.state(), m_key.index(), other.m_key.index(), LUA_OPEQ);
  }
}

bool table_iterator::operator!=(table_iterator const& other) const
{
  return !(operator==(other));
}

table_iterator::reference table_iterator::operator*() const
{
  auto* state = m_key.state();
  auto* stack = m_key.stack();
  stack->release(stack->push(m_key));
  lua_gettable(m_key.state(), m_table.index());
  auto value = stack_object(stack, lua_gettop(state),
      fmt::format("{}[{}]", m_table.name(), m_key.name()));
  return std::make_pair(m_key, value);
}

table_iterator& table_iterator::operator++()
{
  auto* state = m_key.state();
  auto* stack = m_key.stack();
  stack->release(stack->push(m_key));
  int result = lua_next(state, m_table.index());
  if (result == 0) {
    m_is_end = true;
  } else {
    lua_pop(state, 1);
    m_key = stack_object(stack, lua_gettop(state), "iterator_key");
  }
  return *this;
}

table::table(stack_object const& so_arg)
  :stack_object(std::move(so_arg))
{
  check_type(LUA_TTABLE);
}

static std::string make_value_name(
    std::string const& table_name, std::string const& key)
{
  return fmt::format("{}.{}", table_name, key);
}

static std::string make_value_name(
    std::string const& table_name, lua_Integer index)
{
  return fmt::format("{}[{}]", table_name, index);
}

std::optional<stack_object> table::get_optional(char const* key) const
{
  using namespace std::string_literals;
  stack()->prepush();
  int field_type = lua_getfield(state(), index(), key);
  auto value_name = make_value_name(name(), key);
  if (field_type == LUA_TNIL) {
    lua_pop(state(), 1);
    return std::optional<stack_object>();
  }
  return stack_object(stack(), lua_gettop(state()), std::move(value_name));
}

std::optional<stack_object> table::get_optional(lua_Integer key) const
{
  using namespace std::string_literals;
  stack()->prepush();
  int field_type = lua_geti(state(), index(), key);
  auto value_name = make_value_name(name(), key);
  if (field_type == LUA_TNIL) {
    lua_pop(state(), 1);
    return std::optional<stack_object>();
  }
  return stack_object(stack(), lua_gettop(state()), std::move(value_name));
}

bool table::has(char const* key) const
{
  return bool(get_optional(key));
}

bool table::has(lua_Integer key) const
{
  return bool(get_optional(key));
}

stack_object table::get(char const* key) const
{
  auto optional = get_optional(key);
  if (!optional) {
    throw std::runtime_error(fmt::format("{} doesn't exist", make_value_name(name(), key)));
  }
  return optional.value();
}

stack_object table::get(lua_Integer index_arg) const
{
  using namespace std::string_literals;
  int field_type = lua_geti(state(), index(), index_arg);
  auto value_name = make_value_name(name(), index_arg);
  if (field_type == LUA_TNIL) {
    throw std::runtime_error(fmt::format("{} doesn't exist", value_name));
  }
  return stack_object(stack(), lua_gettop(state()), std::move(value_name));
}

table table::get_table(char const* key) const
{
  return table(get(key));
}

table table::get_table(lua_Integer key) const
{
  return table(get(key));
}

lua_Number table::get_number(char const* key) const
{
  return number(get(key)).value();
}

lua_Number table::get_number(lua_Integer key) const
{
  return number(get(key)).value();
}

lua_Integer table::get_integer(char const* key) const
{
  return integer(get(key)).value();
}

lua_Integer table::get_integer(lua_Integer key) const
{
  return integer(get(key)).value();
}

std::string table::get_string(char const* key) const
{
  return string(get(key)).value();
}

std::string table::get_string(lua_Integer key) const
{
  return string(get(key)).value();
}

bool table::get_boolean(char const* key) const
{
  return boolean(get(key)).value();
}

bool table::get_boolean(lua_Integer key) const
{
  return boolean(get(key)).value();
}

double table::get_or(char const* key, double value) const
{
  auto optional = get_optional(key);
  if (optional) return number(std::move(optional.value())).value();
  return value;
}

double table::get_or(lua_Integer key, double value) const
{
  auto optional = get_optional(key);
  if (optional) return number(std::move(optional.value())).value();
  return value;
}

int table::get_or(char const* key, int value) const
{
  auto optional = get_optional(key);
  if (optional) return number(std::move(optional.value())).value();
  return value;
}

int table::get_or(lua_Integer key, int value) const
{
  auto optional = get_optional(key);
  if (optional) return number(std::move(optional.value())).value();
  return value;
}

std::string table::get_or(char const* key, std::string const& value) const
{
  auto optional = get_optional(key);
  if (optional) return string(std::move(optional.value())).value();
  return value;
}

std::string table::get_or(lua_Integer key, std::string const& value) const
{
  auto optional = get_optional(key);
  if (optional) return string(std::move(optional.value())).value();
  return value;
}

std::string table::get_or(char const* key, char const* value) const
{
  auto optional = get_optional(key);
  if (optional) return string(std::move(optional.value())).value();
  return value;
}

std::string table::get_or(lua_Integer key, char const* value) const
{
  auto optional = get_optional(key);
  if (optional) return string(std::move(optional.value())).value();
  return value;
}

bool table::get_or(char const* key, bool value) const
{
  auto optional = get_optional(key);
  if (optional) return boolean(std::move(optional.value())).value();
  return value;
}

bool table::get_or(lua_Integer key, bool value) const
{
  auto optional = get_optional(key);
  if (optional) return boolean(std::move(optional.value())).value();
  return value;
}

table table::get_or_table(char const* key) const
{
  auto optional = get_optional(key);
  if (optional) return table(std::move(optional.value()));
  return stack()->table(make_value_name(name(), key));
}

void table::set_stack_object(char const* key, stack_object const& so)
{
  stack()->release(stack()->push(so));
  lua_setfield(state(), index(), key);
}

void table::set_stack_object(lua_Integer key, stack_object const& so)
{
  stack()->release(stack()->push(so));
  lua_seti(state(), index(), key);
}

void table::set_number(char const* key, lua_Number value)
{
  set_stack_object(key, stack()->number(value, make_value_name(name(), key)));
}

void table::set_number(lua_Integer key, lua_Number value)
{
  set_stack_object(key, stack()->number(value, make_value_name(name(), key)));
}

void table::set_integer(char const* key, lua_Integer value)
{
  set_stack_object(key, stack()->integer(value, make_value_name(name(), key)));
}

void table::set_string(char const* key, std::string const& value)
{
  set_stack_object(key, stack()->string(value, make_value_name(name(), key)));
}

void table::set_boolean(char const* key, bool value)
{
  set_stack_object(key, stack()->boolean(value, make_value_name(name(), key)));
}

void table::set_cfunction(char const* key, lua_CFunction value)
{
  set_stack_object(key, stack()->cfunction(value, make_value_name(name(), key)));
}

void table::set(char const* key, lua_CFunction value)
{
  set_cfunction(key, value);
}

void table::set(char const* key, char const* value)
{
  set_string(key, value);
}

void table::set(char const* key, int value)
{
  set_integer(key, value);
}

void table::set(char const* key, double value)
{
  set_number(key, value);
}

void table::set(char const* key, lua::table const& value)
{
  set_stack_object(key, value);
}

lua_Integer table::size() const
{
  stack()->prepush();
  lua_len(state(), index());
  return lua::integer(lua::stack_object(stack(), lua_gettop(state()),
        fmt::format("#{}", name()))).value();
}

lua_Integer table::num_entries() const
{
  lua_Integer table_size = 0;
  lua_pushnil(state());
  while (lua_next(state(), index()) != 0) {
    lua_pop(state(), 1);
    table_size++;
  }
  return table_size;
}

table table::set_table(char const* key)
{
  auto result = stack()->table(make_value_name(name(), key));
  set_stack_object(key, result);
  return result;
}

table table::set_table(lua_Integer key)
{
  auto result = stack()->table(make_value_name(name(), key));
  set_stack_object(key, result);
  return result;
}

table_iterator table::begin() const
{
  return ++table_iterator(*this, stack()->nil(), false);
}

table_iterator table::end() const
{
  return table_iterator(*this, stack_object(), true);
}

stack::stack(lua_State* state_arg)
  :m_state(state_arg)
{
}

stack::stack(stack&& other)
  :m_state(other.m_state)
  ,m_references(std::move(other.m_references))
{
  for (auto r : m_references)
  {
    r->stack() = this;
  }
}

stack& stack::operator=(stack&& other)
{
  m_state = other.m_state;
  m_references = std::move(other.m_references);
  for (auto r : m_references)
  {
    r->stack() = this;
  }
  return *this;
}

stack::~stack()
{
  for (auto r : m_references)
  {
    r->stack() = nullptr;
  }
  if (m_owns_state) {
    lua_close(m_state);
  }
}

void stack::add_reference(stack_object* reference)
{
  auto it = std::find(m_references.begin(), m_references.end(), reference);
  if (it != m_references.end()) {
    throw std::logic_error(fmt::format(
        "trying to add stack_object reference {} but it already exists!",
        static_cast<void*>(reference)));
  }
  m_references.push_back(reference);
}

void stack::remove_reference(stack_object* reference)
{
  auto new_end = std::remove(
      m_references.begin(),
      m_references.end(),
      reference);
  m_references.erase(new_end, m_references.end());
  bool is_still_referenced = false;
  for (auto r : m_references) {
    if (r->index() == reference->index()) is_still_referenced = true;
  }
  if (!is_still_referenced) {
    lua_remove(m_state, reference->index());
    for (auto r : m_references) {
      if (r->index() > reference->index()) --(r->index());
    }
  }
}

stack_object stack::function_argument(int i, std::string const& function_name)
{
  if (i > lua_gettop(m_state)) {
    throw std::runtime_error(fmt::format("tried to access argument {} of {} but it has {} arguments",
          i, function_name, lua_gettop(m_state)));
  }
  return stack_object(this, i, fmt::format("(argument {} of {})", i, function_name));
}

lua::table stack::table(std::string&& name_arg)
{
  prepush();
  lua_newtable(m_state);
  return lua::table(stack_object(this, lua_gettop(m_state), std::move(name_arg)));
}

lua::number stack::number(lua_Number value, std::string&& name_arg)
{
  prepush();
  lua_pushnumber(m_state, value);
  return lua::number(stack_object(this, lua_gettop(m_state), std::move(name_arg)));
}

lua::stack_object stack::nil()
{
  prepush();
  lua_pushnil(m_state);
  return stack_object(this, lua_gettop(m_state), "nil");
}

lua::integer stack::integer(lua_Integer value, std::string&& name_arg)
{
  prepush();
  lua_pushinteger(m_state, value);
  return lua::integer(stack_object(this, lua_gettop(m_state), std::move(name_arg)));
}

lua::string stack::string(std::string const& value, std::string&& name_arg)
{
  prepush();
  lua_pushstring(m_state, value.c_str());
  return lua::string(stack_object(this, lua_gettop(m_state), std::move(name_arg)));
}

lua::boolean stack::boolean(bool value, std::string&& name_arg)
{
  prepush();
  lua_pushboolean(m_state, value ? 1 : 0);
  return lua::string(stack_object(this, lua_gettop(m_state), std::move(name_arg)));
}

lua::function stack::cfunction(lua_CFunction value, std::string&& name_arg)
{
  prepush();
  lua_pushcfunction(m_state, value);
  return lua::function(stack_object(this, lua_gettop(m_state), std::move(name_arg)));
}

stack_object stack::push(stack_object const& other)
{
  prepush();
  lua_pushvalue(m_state, other.index());
  std::string name = other.name();
  return stack_object(this, lua_gettop(m_state), std::move(name));
}

void stack::prepush() const
{
  int result = lua_checkstack(m_state, 1);
  if (!result) {
    throw std::runtime_error("reached the maximum stack limit for Lua");
  }
}

void stack::release(stack_object&& reference)
{
  auto new_end = std::remove(
      m_references.begin(),
      m_references.end(),
      &reference);
  m_references.erase(new_end, m_references.end());
  for (auto r : m_references) {
    if (r->index() == reference.index()) {
      throw std::logic_error("stack::release: there was another object referencing this index");
    }
  }
  reference.stack() = nullptr;
}

stack stack::newstate()
{
  lua_State* L = luaL_newstate();
  luaL_openlibs(L);
  auto result = stack(L);
  result.m_owns_state = true;
  return result;
}

void stack::setglobal(char const* key, stack_object const& so)
{
  release(push(so));
  lua_setglobal(m_state, key);
}

void stack::setglobal(char const* key, lua_CFunction value)
{
  setglobal(key, cfunction(value, std::string(key)));
}

void stack::dofile(std::filesystem::path const& path)
{
  if (luaL_dofile(m_state, path.c_str())) {
    throw std::runtime_error(fmt::format(
          "executing Lua file {} failed saying:\n{}", 
          path.string(), lua_tostring(m_state, -1)));
  }
}

stack_object stack::getglobal(char const* key)
{
  prepush();
  lua_getglobal(m_state, key);
  return stack_object(this, lua_gettop(m_state), std::string(key));
}

lua::table stack::require(char const* key)
{
  auto global_require_function = lua::function(getglobal("require"));
  return lua::table(global_require_function.first_result(key));
}

}
}
