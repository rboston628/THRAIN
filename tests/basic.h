// for installing through brew
// https://github.com/CxxTest/cxxtest/issues/133

 #include <cxxtest/TestSuite.h>

class MostBasicTest : public CxxTest::TestSuite {
public:
    void testAddition( void ) {
        TS_ASSERT( 1 + 1 > 1 );
        TS_ASSERT_EQUALS( 1 + 1, 2 );
        TS_ASSERT_DIFFERS( 1 + 1, 3 );
    }
};