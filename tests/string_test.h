#include "../src/constants.h"
#include "../lib/string.h"
#include <cxxtest/TestSuite.h>

class StringFormatterTest : public CxxTest::TestSuite {
public:

    void test_string_addition( void ){

        // test adding string and c-string
        char x[] = "_c_string_1";
        char y[] = "_c_string_2";
        std::string s = "_std_string";

        TS_ASSERT_EQUALS(x[11], '\0');
        TS_ASSERT_EQUALS(y[11], '\0');
        TS_ASSERT_EQUALS(s[11], '\0');

        std::string res1 = "_std_string_c_string_1"; // string + cstring
        std::string res2 = "_c_string_1_std_string"; // cstring + string
        std::string res3 = "_c_string_1_c_string_2"; // cstring + cstring

        TS_ASSERT_EQUALS( addstring(s,x), res1 );
        TS_ASSERT_EQUALS( addstring(x,s), res2 );
        TS_ASSERT_EQUALS( addstring(x,y), res3 );

        TS_ASSERT_EQUALS( addstring(s,x)[22], '\0');
        TS_ASSERT_EQUALS( addstring(x,s)[22], '\0');
        TS_ASSERT_EQUALS( addstring(x,y)[22], '\0');

    }

    void test_string_format( void ) {
        printf("TEST OF STRING FORMAT\n");

        // test formating integers
        std::string integers = strmakef("x%d%d", 2,7);
        TS_ASSERT_EQUALS(integers, "x27");

        // test formatting doubles
        std::string doubles = strmakef("%0.4lf + %2.6le", m_pi, sqrt(117));
        TS_ASSERT_EQUALS(doubles, "3.1416 + 1.081665e+01");
    }
};