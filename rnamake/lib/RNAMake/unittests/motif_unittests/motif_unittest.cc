//
//  motif_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "unittest.h"
#include "motif_unittest.h"
#include "motif/motif_factory.h"
#include "motif/motif_to_secondary_structure.h"
#include "util/file_io.h"
#include "util/settings.h"

MotifUnittest::MotifUnittest() {
    
    MotifFactory mf;
    //String path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6";
    //m_ = mf.motif_from_file(path);
    
    String path2 = base_dir() + "/rnamake/resources/motifs/helices/HELIX.IDEAL";
    MotifOP m = mf.motif_from_file(path2);
    std::cout << m->dot_bracket() << std::endl;
    
    /*String path = unittest_resource_dir() + "/motif/test_str_to_motif.dat";
    Strings lines = get_lines_from_file(path);
    
    ResidueTypeSet rts;
    m_ = Motif(lines[0], rts);
    */
}
/*
int
MotifUnittest::test_copy() {
    Motif mcopy = m_.copy();
    Point p(50,0,0);
    m_.move(p);
    
    Point org_center = center(m_.atoms());
    Point new_center = center(mcopy.atoms());
    float dist = (50) - org_center.distance(new_center);
    if(dist > 0.001) { return 0; }
    
    return 1;
}

int
MotifUnittest::test_to_str() {
    ResidueTypeSet rts;
    Motif m2 = Motif(m_.to_str(), rts);
    
    AtomOPs new_atoms = m2.atoms();
    AtomOPs org_atoms = m_.atoms();
    
    float dist;
    for(int i = 0; i < org_atoms.size(); i++) {
        if(org_atoms[i].get() == NULL) { continue; }
        dist = org_atoms[i]->coords().distance(new_atoms[i]->coords());
    }
    return 1;
}

int
MotifUnittest::test_secondary_structure() {
    String org = "....((((((...(((.((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....)).)))..))))))(....((((..(((((((((.....)))))))).)))))......)";
    String new_ss = m_.secondary_structure();
    for(int i = 0; i < org.size(); i++) {
        if(org[i] != new_ss[i]) { return 0; }
    }
    return 1;
}

int
MotifUnittest::test_creation_from_dir() {
    //String m_path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/";
    //Motif m(m_path);
    
    return 1;
}

int
MotifUnittest::test_get_basepair_by_name() {
    MotifLibrary mlib(HELIX);
    MotifOP m = mlib.get_motif("HELIX.IDEAL");
    BasepairOP bp = m->get_basepair_by_name("A5-B7");

    return 1;
}
*/

int
MotifUnittest::run() {
    /*if (test_copy() == 0)                 { std::cout << "test_copy failed" << std::endl; }
    if (test_to_str() == 0)               { std::cout << "test_to_str failed" << std::endl; }
    if (test_secondary_structure() == 0)  { std::cout << "test_secondary_structure failed" << std::endl; }
    if (test_creation_from_dir() == 0)    { std::cout << "test_copy failed" << std::endl; }
    if (test_get_basepair_by_name() == 0) { std::cout << "test_get_basepair_by_name failed" << std::endl; }
     */
    return 0;
}

void
MotifUnittest::run_all() {
    /*String name = "MotifUnittest";
    typedef int (MotifUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_copy"   ] = &MotifUnittest::test_copy;
    func_map["test_to_str" ] = &MotifUnittest::test_to_str;
    func_map["test_secondary_structure"] = &MotifUnittest::test_secondary_structure;
    func_map["test_creation_from_dir"       ] = &MotifUnittest::test_creation_from_dir;
    func_map["test_get_basepair_by_name"       ] = &MotifUnittest::test_get_basepair_by_name;

    for(auto const & kv : func_map) {
        try {
            int result = (this->*kv.second)();
            if(result == 0) {
                std::cout << name << "::" << kv.first << " FAILED!" << std::endl;
            }
        }
        catch(...) {
            std::cout << name << "::" << kv.first << " returned ERROR!" << std::endl;
        }
        
    }*/
}
























