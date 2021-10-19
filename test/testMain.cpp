//
// Created by Alec Glisman on 7/30/21
//

/* Tell Catch2 to provide a main() - only do this in one cpp file */
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

/* REVIEW: Introductory material
 * @SOURCE: https://github.com/catchorg/Catch2/blob/master/docs/assertions.md#top

 * NATURAL EXPRESSIONS
 *     REQUIRE( expression ) -> Test expression and aborts if it fails.
 *     CHECK( expression )   -> Same as require but execution continues.

 * FLOATING POINT COMPARISONS
 *     Use Approx( double ) macro to take tolerance into account
 *     REQUIRE( performComputation() == Approx( 2.1 ) );

 * EXCEPTIONS
 *     REQUIRE_NOTHROW( expression ), REQUIRE_THROWS( expression ),
 *     REQUIRE_THROWS_AS( expression, exception type ),
 *     REQUIRE_THROWS_WITH( expression, string or string matcher )
 *     Also CHECK equivalents
 */

/* REVIEW: this is just example code

TEST_CASE( "vectors can be sized and resized", "[vector]" ) {

    std::vector<int> v( 5 );

    REQUIRE( v.size() == 5 );
    REQUIRE( v.capacity() >= 5 );

    SECTION( "resizing bigger changes size and capacity" ) {
        v.resize( 10 );

        REQUIRE( v.size() == 10 );
        REQUIRE( v.capacity() >= 10 );
    }
    SECTION( "resizing smaller changes size but not capacity" ) {
        v.resize( 0 );

        REQUIRE( v.size() == 0 );
        REQUIRE( v.capacity() >= 5 );
    }
    SECTION( "reserving bigger changes capacity but not size" ) {
        v.reserve( 10 );

        REQUIRE( v.size() == 5 );
        REQUIRE( v.capacity() >= 10 );
    }
    SECTION( "reserving smaller does not change size or capacity" ) {
        v.reserve( 0 );

        REQUIRE( v.size() == 5 );
        REQUIRE( v.capacity() >= 5 );
    }
}

 */
