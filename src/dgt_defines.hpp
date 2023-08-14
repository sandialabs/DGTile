#pragma once

namespace dgt {

static constexpr int X = 0;
static constexpr int Y = 1;
static constexpr int Z = 2;
static constexpr int DIMENSIONS = 3;

static constexpr int MIN = 0;
static constexpr int MAX = 1;
static constexpr int LEFT = 0;
static constexpr int RIGHT = 1;
static constexpr int DIRECTIONS = 2;

static constexpr int SEND = 0;
static constexpr int RECV = 1;
static constexpr int OWNED = 0;
static constexpr int GHOST = 1;

extern int errors;

using real = double;

}
