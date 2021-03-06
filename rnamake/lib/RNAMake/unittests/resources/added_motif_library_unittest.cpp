
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "resources/added_motif_library.h"

TEST_CASE( "Test Added Motif Library", "[AddedMotifLibrary]" ) {
    
    auto path = motif_dirs() + "base.motif";
    auto m = file_to_motif(path);
    
    auto mlib = AddedMotifLibrary();
    mlib.add_motif(m);
    
    SECTION("try getting back motif that was stored in library") {
        REQUIRE_NOTHROW(mlib.get(m->name()));
        REQUIRE_NOTHROW(mlib.get("", m->end_ids()[0]));
        REQUIRE_NOTHROW(mlib.get("", "", m->ends()[0]->name()));
        REQUIRE_NOTHROW(mlib.get(m->name(), m->end_ids()[0], m->ends()[0]->name()));
        
        REQUIRE_THROWS_AS(mlib.get(), AddedMotifLibraryException);
        REQUIRE_THROWS_AS(mlib.get("TEST"), AddedMotifLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "TEST"), AddedMotifLibraryException);
        REQUIRE_THROWS_AS(mlib.get("", "", "TEST"), AddedMotifLibraryException);

    }
    
    SECTION("trying to catch adding the same motif twice") {
        REQUIRE_THROWS_AS(mlib.add_motif(m), AddedMotifLibraryException);
    }
    
    SECTION("test detecting whether a motif exists in the library") {
        REQUIRE(mlib.contains(m->name()) == 1);
        REQUIRE(mlib.contains("", m->end_ids()[0]) == 1);

        REQUIRE(mlib.contains("TEST") == 0);
        
    }
    
}