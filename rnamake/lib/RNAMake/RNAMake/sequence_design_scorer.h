//
//  sequence_design_scorer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 3/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sequence_design_scorer__
#define __RNAMake__sequence_design_scorer__

#include <stdio.h>
#include "vienna.h"

class SequenceDesignScorer {
public:
    SequenceDesignScorer():
    v_ ( Vienna() ),
    score_ ( 0 ),
    cofold_( 0 )  {}
    
    ~SequenceDesignScorer() {}
    
public:
    
    void
    setup(String const & ss) {
        ss_ = ss;
        for(auto const & s : ss_) {
            if(s == '&') { cofold_ = 1; }
        }
    }
    
public:
    float
    score(
        String const &);
    
private:
    Vienna v_;
    FoldResult fr_;
    String ss_;
    float score_;
    int cofold_;
    
};

typedef std::shared_ptr<SequenceDesignScorer> SequenceDesignScorerOP;

#endif /* defined(__RNAMake__sequence_design_scorer__) */