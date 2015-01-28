//
//  chain.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__chain__
#define __RNAMake__chain__

#include <stdio.h>
#include "residue.h"
#include "residue_type_set.h"
#include "types.h"

class Chain {
public:
    Chain() {}
    Chain(
        Residues const & residues):
        residues_ ( residues )
    {}
    
    Chain
    copy() const;
    
    ~Chain() {}

public:
    
    inline
    Residue const &
    first() { return residues_[0]; }
    
    inline
    Residue const &
    last() { return residues_.back(); }
    
    inline
    Chain
    subchain(int start, int end) {
        if(start < 0) { throw "start cannot be less then 0"; }
        return Chain(Residues(residues_.begin() + start, residues_.begin() + end));
    }
    
    inline
    Chain
    subchain(
        Residue const & r1,
        Residue const & r2) {
        int start = (int)(std::find(residues_.begin(), residues_.end(), r1) - residues_.begin());
        int end   = (int)(std::find(residues_.begin(), residues_.end(), r2) - residues_.begin());
        if( start > end) {
            int temp = start;
            start = end;
            end = temp;
        }
        return subchain(start, end);
    }
    
    String
    to_str() const;
    
    String
    to_pdb_str(int &) const;
    
    void
    to_pdb(String const) const;
        
public: //setters
    
    inline
    void
    residues(Residues const & nresidues) { residues_ = nresidues; }
    
public: //getters
    
    inline
    Residues
    residues() const { return residues_; }
    
private:
    Residues residues_;
};

Chain
str_to_chain(
    String const &,
    ResidueTypeSet const & );

typedef std::vector<Chain> Chains;

#endif /* defined(__RNAMake__chain__) */
