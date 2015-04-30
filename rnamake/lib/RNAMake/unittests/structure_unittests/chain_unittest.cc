//
//  chain_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "chain_unittest.h"
#include "math/numerical.h"
#include "util/file_io.h"

ChainUnittest::ChainUnittest() {
    
    String path = unittest_resource_dir() + "/chain/test_str_to_chain.dat";
    Strings lines = get_lines_from_file(path);
    
    ResidueTypeSet rts;
    c_ = ChainOP(new Chain(str_to_chain(lines[0], rts)));
}

int
ChainUnittest::test_to_str() {

    String s = c_->to_str();
    ResidueTypeSet rts;
    Chain c2 = str_to_chain(s, rts);
    
    ResidueOPs org_res = c_->residues();
    ResidueOPs new_res = c2.residues();
    
    for(int i = 0; i < org_res.size(); i++) {
        if(org_res[i]->name().compare(new_res[i]->name()) != 0) { return 0; }
        if(org_res[i]->chain_id().compare(new_res[i]->chain_id())) { return  0; }
        if(org_res[i]->num() != new_res[i]->num()) { return 0; }
        
        AtomOPs atoms_1 = org_res[i]->atoms();
        AtomOPs atoms_2 = new_res[i]->atoms();
        
        for(int j = 0; j < atoms_1.size(); j++) {
            if(atoms_1[j].get() == NULL && atoms_2[j].get() != NULL) { return 0; }
            if(atoms_1[j].get() != NULL && atoms_2[j].get() == NULL) { return 0; }
            if(atoms_1[j].get() == NULL && atoms_2[j].get() == NULL) { continue; }
            
            if(!are_xyzVector_equal(atoms_1[j]->coords(), atoms_2[j]->coords())) { return 0; }
        }
    }
    
    return 1;
}

int
ChainUnittest::test_subchain() {
    
    try {
        ChainOP sc1 = c_->subchain(-1, 5);
        std::cout << "did not catch exception\n" << std::endl;
        exit(0);
    } catch (char const * e) {}
    
    try {
        ChainOP sc1 = c_->subchain(1, 1);
        std::cout << "did not catch exception\n" << std::endl;
        exit(0);
    } catch (char const * e) {}
    
    try {
        ChainOP sc1 = c_->subchain(1, 10000);
        std::cout << "did not catch exception\n" << std::endl;
        exit(0);
    } catch (char const * e) {}
    
    ChainOP sc = c_->subchain(1, 5);
    if(sc->residues().size() != 4) { return 0; }
    
    ResidueOP r1 = c_->residues()[1];
    ResidueOP r2 = c_->residues()[5];
    ChainOP sc2 = c_->subchain(r1, r2);
    if(sc2->residues().size() != 4) { return 0; }
    
    return 1;
}

int
ChainUnittest::test_to_pdb() {
    int i = 1;
    String s = c_->to_pdb_str(i);
    
    return 1;
}



int
ChainUnittest::run() {
    
    if (test_to_str() == 0)             { std::cout << "test_to_str failed" << std::endl; }
    if (test_subchain() == 0)           { std::cout << "test_subchain failed" << std::endl; }
    if (test_to_pdb() == 0)             { std::cout << "test_to_pdb failed" << std::endl; }
    
    return 1;
}