//
//  residue_type.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__residue_type__
#define __RNAMake__residue_type__

#include <stdio.h>
#include <vector>
#include "types.h"
#include "atom.h"

class ResidueType {
public:
    ResidueType() {}
    ResidueType(
        String const &,
        StringIntMap const &);
    
    ~ResidueType() {}

public:
    String
    get_correct_atom_name(
        Atom const &) const;
    
    int
    match_name(
        String const &) const;
    
    inline
    String
    const &
    name() const { return name_; }
    
    inline
    String
    const &
    short_name() const { return alt_names_[0]; }
    
    inline
    int
    atom_pos_by_name(
        String const & aname) const {
        
        StringIntMap::const_iterator iter( atom_map_.find(aname));
        if(iter != atom_map_.end()) {
            return iter->second;
        }
        else {
            return -1;
        }
    }
    
private:
    
    void
    extend_res_specific_altnames();
    
private:
    String name_;
    StringIntMap atom_map_;
    Strings alt_names_;
    StringStringMap atom_alt_names_;

};

typedef std::vector<ResidueType> ResidueTypes;


#endif /* defined(__RNAMake__residue_type__) */
