#include "doctest.h"
#include <cstdio>
#include <cstdlib>
#include "../src/ThrainConfig.h"

struct SetUp {
  SetUp() {
    ThrainConfig::reconfigure("./tests/tests.config");
    printf("test of testing\n");
  }
  ~SetUp() = default;
};

TEST_SUITE("basic") {

TEST_CASE_FIXTURE(SetUp, "sanity") {
  REQUIRE( true );
  REQUIRE( !false );
  REQUIRE( 1 == 1 );
  REQUIRE( 1 != 2 );
}

TEST_CASE("addition") {
  REQUIRE( 1 + 1 > 1 );
  REQUIRE( 1 + 1 == 2 );
  REQUIRE( 1 + 1 != 3 );
}

}