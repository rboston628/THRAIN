#include "doctest.h"
#include <cstdio>
#include <cstdlib>
#include "../src/ThrainConfig.h"

TEST_SUITE("basic [unit]") {

TEST_CASE("basic:sanity") {
  REQUIRE( true );
  REQUIRE( !false );
  REQUIRE( 1 == 1 );
  REQUIRE( 1 != 2 );
}

TEST_CASE("basic:addition") {
  REQUIRE( 1 + 1 > 1 );
  REQUIRE( 1 + 1 == 2 );
  REQUIRE( 1 + 1 != 3 );
}

}