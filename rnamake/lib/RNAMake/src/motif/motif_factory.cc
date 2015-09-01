//
//  motif_factory.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>

#include "motif/motif_factory.h"
#include "structure/chain.h"
#include "util/file_io.h"
#include "util/x3dna.h"


MotifOP
MotifFactory::motif_from_file(
    String const & path) {
    
    auto fname = filename(path);
    auto pdb_path = path;
    StructureOP structure;
    if(is_dir(path)) {
        structure = sf_.get_structure(path + "/" + fname + ".pdb");
        pdb_path = path + "/" + fname + ".pdb";
    }
    else {
        structure = sf_.get_structure(path);
        fname = fname.substr(0, -4);
    }
    
    auto basepairs = _setup_basepairs(pdb_path, structure);
    auto ends = _setup_basepair_ends(structure, basepairs);
    auto m = std::make_shared<Motif>(structure, basepairs, ends);
    m->name(fname);
    m->path(path);
    _setup_secondary_structure(m);
    
    return m;
}

MotifOP
MotifFactory::motif_from_res(
    ResidueOPs & res,
    BasepairOPs const & bps) {
    
    auto chains = sf_.build_chains(res);
    auto structure = std::make_shared<Structure>(chains);
    auto ends = _setup_basepair_ends(structure, bps);
    auto m = std::make_shared<Motif>(structure, bps, ends);
    _setup_secondary_structure(m);
    return m;
}


void
MotifFactory::standardize_motif(
    MotifOP & m) {
    
    _align_chains(m);
    _align_ends(m);
    _setup_secondary_structure(m);
}


MotifOP
MotifFactory::can_align_motif_to_end(
    MotifOP const & m,
    int ei) {
    
    auto m_end = m->ends()[ei];
    auto m_added = get_aligned_motif(base_motif_->ends()[1], m_end, m);
    if(_steric_clash(base_motif_, m_added)) {
        m_end->flip();
        m_added = get_aligned_motif(base_motif_->ends()[1], m_end, m);
        if(_steric_clash(base_motif_, m_added)) { return nullptr; }
    }
    
    int fail = 0;
    MotifOP m2_added;
    for(auto & end : m_added->ends()) {
        if(end == m_end) { continue; }
        m2_added = get_aligned_motif(end, added_helix_->ends()[0], added_helix_);
        if(_steric_clash(m_added, m2_added) == 0 &&
           _steric_clash(base_motif_, m2_added) == 0) {
            continue;
        }
        
        end->flip();
        m2_added = get_aligned_motif(end, added_helix_->ends()[0], added_helix_);
        if(_steric_clash(m_added, m2_added) || _steric_clash(base_motif_, m2_added)) {
            fail = 1; break;
        }

    }
    
    if (fail) { return nullptr; }
    else      { return m_added; }
    
}

MotifOP
MotifFactory::align_motif_to_common_frame(
    MotifOP const & m,
    int ei) {
    
    auto m_added = get_aligned_motif(ref_motif_->ends()[0], m->ends()[ei], m);
    standardize_motif(m_added);
    return m_added;
}

BasepairOPs
MotifFactory::_setup_basepairs(
    String const & path,
    StructureOP const & structure) {
    
    BasepairOPs basepairs;
    X3dna x3dna_parser;
    X3Basepairs x_basepairs = x3dna_parser.get_basepairs(path);
    ResidueOP res1, res2;
    BasepairOP bp;
    for(auto const & xbp : x_basepairs) {
        res1 = structure->get_residue(xbp.res1.num, xbp.res1.chain_id, xbp.res1.i_code);
        res2 = structure->get_residue(xbp.res2.num, xbp.res2.chain_id, xbp.res2.i_code);
        if (res1 == nullptr || res2 == nullptr) {
            throw "cannot find residues in basepair during setup";
        }
        
        bp = BasepairOP(new Basepair(res1, res2, xbp.r, xbp.bp_type));
        basepairs.push_back(bp);
    }
    
    return basepairs;
    
}

BasepairOPs
MotifFactory::_setup_basepair_ends(
    StructureOP const & structure,
    BasepairOPs const & basepairs) {
    
    ResidueOPs chain_ends;
    for(auto const & c : structure->chains()) {
        chain_ends.push_back(c->first());
        if(c->residues().size() > 1) { chain_ends.push_back(c->last()); }
    }
    
    BasepairOPs ends;
    for( auto const & bp : basepairs ) {
        for (auto const & ce1 : chain_ends) {
            for(auto const & ce2 : chain_ends) {
                if(bp->bp_type().compare("cW-W") == 0 && bp->res1() == ce1 && bp->res2() == ce2) {
                    ends.push_back(bp);
                }
            }
        }
    }
    
    return ends;
}

void
MotifFactory::_setup_secondary_structure(
    MotifOP & m) {
    
    auto ss = parser_.to_secondary_structure(m);
    Strings end_ids(m->ends().size());
    int i = 0;
    for(auto const & end : m->ends()) {
        auto res1 = ss->get_residue(end->res1()->uuid());
        auto res2 = ss->get_residue(end->res2()->uuid());
        auto ss_end = ss->get_bp(res1, res2);
        end_ids[i] = sstruct::assign_end_id(ss, ss_end);
    }
    
    ss->end_ids(end_ids);
    m->end_ids(end_ids);
    m->secondary_structure(ss);
    
    
}

void
MotifFactory::_align_chains(
    MotifOP & m) {
    
    auto chains = m->chains();
    ChainOP closest;
    float best = 1000, dist;
    auto c2 = center(ref_motif_->ends()[0]->atoms());
    Point c1;
    for(auto const & c : chains) {
        c1 = center(c->first()->atoms());
        dist = c1.distance(c2);
        if(dist < best) {
            best = dist;
            closest = c;
        }
    }
    
    auto updated_chains = ChainOPs{closest};
    for (auto const & c : chains) {
        if(c != closest) { updated_chains.push_back(c); }
    }
    
    m->structure(std::make_shared<Structure>(Structure(updated_chains)));
}

void
MotifFactory::_align_ends(
    MotifOP & m) {
    
    BasepairOP closest;
    float best = 1000, dist;
    auto c2 = center(ref_motif_->ends()[0]->atoms());
    Point c1;
    
    for(auto const & end : m->ends()) {
        c1 = center(end->atoms());
        dist = c1.distance(c2);
        if(dist < best) {
            best = dist;
            closest = end;
        }
    }
    
    auto updated_ends = BasepairOPs{closest};
    for(auto const & end : m->ends()) {
        if(end != closest) { updated_ends.push_back(end); }
    }
    
    int flip_res = 0;;
    for(auto & end : m->ends()) {
        flip_res = 0;
        for (auto const & c : m->chains()) {
            if(c->first() == end->res2()) {
                flip_res = 1;
                break;
            }
        }
        
        if(!flip_res) { continue; }
            
        auto temp = end->res1();
        end->res1(end->res2());
        end->res2(temp);
    }
    
    m->ends(updated_ends);
    
}

int
MotifFactory::_steric_clash(
    MotifOP const & m1,
    MotifOP const & m2) {
    
    float dist;
    for(auto const & c1 : m1->beads()) {
        for(auto const & c2 : m2->beads()) {
            if(c1.btype() == BeadType::PHOS || c2.btype() == BeadType::PHOS) {continue; }
            dist = c1.center().distance(c2.center());
            if( dist < clash_radius_) { return 1; }
        }
    }
    
    return 0;
    
}


















