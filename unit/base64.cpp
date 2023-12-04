#include <dgt_base64.hpp>

#include <gtest/gtest.h>

TEST(base64, encode_integers)
{
  int in_data[3];
  in_data[0] = 100;
  in_data[1] = 12;
  in_data[2] = 12;
  std::uint64_t const bytes = sizeof(int) * static_cast<std::uint64_t>(3);
  std::string const encoded = dgt::base64::encode(in_data, bytes);
  EXPECT_EQ(encoded, "ZAAAAAwAAAAMAAAA");
}

TEST(base64, decode_integers)
{
  int out_data[3];
  std::string const encoded = "ZAAAAAwAAAAMAAAA";
  std::uint64_t const bytes = sizeof(int) * static_cast<std::uint64_t>(3);
  dgt::base64::decode(encoded, out_data, bytes);
  EXPECT_EQ(out_data[0], 100);
  EXPECT_EQ(out_data[1], 12);
  EXPECT_EQ(out_data[2], 12);
}

TEST(base64, encode_doubles)
{
  double in_data[3];
  in_data[0] = 100.;
  in_data[1] = 12.;
  in_data[2] = 12.;
  std::uint64_t const bytes = sizeof(double) * static_cast<std::uint64_t>(3);
  std::string const encoded = dgt::base64::encode(in_data, bytes);
  EXPECT_EQ(encoded, "AAAAAAAAWUAAAAAAAAAoQAAAAAAAAChA");
}

TEST(base64, decode_doubles)
{
  double out_data[3];
  std::string const encoded = "AAAAAAAAWUAAAAAAAAAoQAAAAAAAAChA";
  std::uint64_t const bytes = sizeof(double) * static_cast<std::uint64_t>(3);
  dgt::base64::decode(encoded, out_data, bytes);
  EXPECT_EQ(out_data[0], 100.);
  EXPECT_EQ(out_data[1], 12.);
  EXPECT_EQ(out_data[2], 12.);
}
