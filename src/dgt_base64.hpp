#pragma once

#include <string>

namespace dgt {
namespace base64 {

std::string encode(void const* data, std::size_t size);
void decode(std::string const& text, void* data, std::size_t size);

}
}
