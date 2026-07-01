// Tests of the simple string functions
// which are intended to enable better support
// for simple c-style strings and c-style formatting,
// mixed with the more modern std::string objects

#include "../src/constants.h"
#include "../lib/string.h"
#include "doctest.h"

TEST_SUITE("StringFormatter [unit]") {
TEST_CASE("string: string_addition") {
  // test adding string and c-string
  char x[] = "_c_string_1";
  char y[] = "_c_string_2";
  std::string s = "_std_string";

  CHECK( x[11] == '\0' );
  CHECK( y[11] == '\0' );

  std::string res1 = "_std_string_c_string_1"; // string + cstring
  std::string res2 = "_c_string_1_std_string"; // cstring + string
  std::string res3 = "_c_string_1_c_string_2"; // cstring + cstring

  CHECK( addstring(s,x) == res1 );
  CHECK( addstring(x,s) == res2 );
  CHECK( addstring(x,y) == res3 );

  CHECK( addstring(s,x)[22] == '\0');
  CHECK( addstring(x,s)[22] == '\0');
  CHECK( addstring(x,y)[22] == '\0');
}

TEST_CASE("string: string_format") {
  // test formating integers
  std::string integers = strmakef("x%d%d", 2,7);
  CHECK( integers == "x27" );

  // test formatting doubles
  std::string doubles = strmakef("%0.4lf + %2.6le", m_pi, sqrt(117));
  CHECK( doubles == "3.1416 + 1.081665e+01" );

  // the implementation uses a buffer of size 256
  // try to deliberately make something larger than this
  constexpr std::size_t really_big_size = 400;
  char really_big_char[really_big_size+1];
  std::string exp_result = "";
  for(int i=0; i<really_big_size; i++){
    really_big_char[i] = 'a';
    exp_result.push_back('a');
  }
  // must add null-terminator
  really_big_char[really_big_size] = 0;
  // compare results
  std::string test_big_string = strmakef("%s", really_big_char);
  CHECK( test_big_string == exp_result );
}
}
