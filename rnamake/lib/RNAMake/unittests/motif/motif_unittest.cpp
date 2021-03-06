
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "structure/is_equal.hpp"
#include "motif/motif.h"

TEST_CASE( "Test Motifs the core of everything!", "[Motif]" ) {
    
    auto path = motif_dirs() + "base.motif";
    auto m = file_to_motif(path);
    path = motif_dirs() + "ref.motif";
    auto ref_m = file_to_motif(path);
    auto rts = ResidueTypeSet();
    
    REQUIRE(m->ends().size() == 2);
    REQUIRE(m->basepairs().size() == 5);
    
    SECTION("test copy constuctor for motif") {
        auto m_copy = Motif(*m);
        
        REQUIRE(m_copy.ends().size() == 2);
        REQUIRE(m_copy.basepairs().size() == 5);

        for(auto const & r : m->residues()) {
            REQUIRE(m_copy.get_residue(r->uuid()) != nullptr);
        }
        
        for(auto const & bp : m->basepairs()) {
            REQUIRE(m_copy.get_basepair(bp->uuid()).size() == 1);
        }
        
        
        REQUIRE(m_copy.atoms().size() == m->atoms().size());
        auto atoms1 = m->atoms();
        auto atoms2 = m_copy.atoms();
        
        REQUIRE(are_atom_vectors_equal(atoms1, atoms2));

    }
    
    SECTION("test stringifying motif") {
        auto s = m->to_str();
        auto m_copy = Motif(s, rts);
        
        REQUIRE(m_copy.ends().size() == 2);
        REQUIRE(m_copy.basepairs().size() == 5);

        for(auto const & r : m->residues()) {
            REQUIRE(m_copy.get_residue(r->num(), r->chain_id(), r->i_code()) != nullptr);
        }
        
        REQUIRE(m_copy.atoms().size() == m->atoms().size());
        auto atoms1 = m->atoms();
        auto atoms2 = m_copy.atoms();
        
        REQUIRE(are_atom_vectors_equal(atoms1, atoms2));
    }
    
    SECTION("test aligning motifs") {
        m->move(Point(10, 10, 10));
        auto m_aligned = get_aligned_motif(ref_m->ends()[0], m->ends()[0], m);
        
        auto dist = ref_m->ends()[0]->d().distance(m_aligned->ends()[0]->d());
        auto r_dist = ref_m->ends()[0]->r().difference(m_aligned->ends()[0]->r());

        REQUIRE(dist < 1.0);
        REQUIRE(r_dist < 0.0001);
        
    }
    
    SECTION("test that repeat aligning does not cause error") {
        m->move(Point(10, 10, 10));
        auto ref_bp = ref_m->ends()[0]->state();
        
        auto dist = 0.0f, r_dist = 0.0f;
        
        for(int i = 0; i < 100; i++) {
            align_motif(ref_bp, m->ends()[0], m);
            
            dist = ref_bp->d().distance(m->ends()[0]->d());
            r_dist = ref_bp->r().difference(m->ends()[0]->r());
            
        }
        
        REQUIRE(dist < 1.0);
        REQUIRE(r_dist < 0.0001);

    }
    
    SECTION("test getting new unique indentifers") {
        auto m_copy = Motif(*m);
        m_copy.new_res_uuids();
        
        SECTION("all residues and basepairs have new indentifers cant search with originals anymore") {
        
            for(auto const & r : m->residues()) {
                REQUIRE(m_copy.get_residue(r->uuid()) == nullptr);
            }
        
            for(auto const & bp : m->basepairs()) {
                REQUIRE(m_copy.get_basepair(bp->uuid()).size() == 0);
            }
        }
        
        SECTION("should still be able to find them using names to find them") {
            for(auto const & r : m->residues()) {
                REQUIRE(m_copy.get_residue(r->num(), r->chain_id(), r->i_code()) != nullptr);
            }
        }

        
    }
    
    SECTION("test copying uuids from one motif to another") {
        auto path = motif_dirs() + "base.motif";
        auto m1 = file_to_motif(path);
        auto m2 = file_to_motif(path);

        REQUIRE(m1->id() != m2->id());
        m1->copy_uuids_from_motif(*m2);
        
        REQUIRE(m1->id() == m2->id());
        
        int i = 0;
        for(auto const & bp : m1->basepairs()) {
            REQUIRE(bp->uuid() == m2->basepairs()[i]->uuid());
            i++;
        }
        
        i = 0;
        for(auto const & r : m1->residues()) {
            REQUIRE(r->uuid() == m2->residues()[i]->uuid());
            i++;
        }
        
    }
    
}







