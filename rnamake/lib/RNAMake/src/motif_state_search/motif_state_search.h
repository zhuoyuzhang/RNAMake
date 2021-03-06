//
//  motif_state_search.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_search__
#define __RNAMake__motif_state_search__

#include <stdio.h>

//RNAMake Headers
#include "base/option.h"
#include "util/steric_lookup.hpp"
#include "motif/motif_state_aligner.h"
#include "motif_state_search/motif_state_selector.h"
#include "motif_state_search/motif_state_search_scorer.h"
#include "motif_state_search/motif_state_search_node.fwd.h"
#include "motif_state_search/motif_state_search_node.h"
#include "motif_state_search/motif_state_search_solution.h"


class MotifStateSearch : public OptionClass {
public:
    MotifStateSearch():
    queue_(MotifStateSearchNodeQueue()),
    selector_(default_selector()),
    scorer_(std::make_shared<MotifStateSearchScorer>(MSS_Astar())),
    solutions_(MotifStateSearchSolutionOPs()),
    beads_(Points()),
    aligner_(MotifStateAligner()),
    lookup_(StericLookup()),
    using_lookup_(0),
    options_(Options()) {
        setup_options();
    }
    
    ~MotifStateSearch() {}
    
public:
    
    void
    setup(
        BasepairStateOP const &,
        BasepairStateOP const &);
    
    MotifStateSearchSolutionOP
    next();
    
    int
    finished();
    
    void
    reset();
    
public: //option wrappers
    
    inline
    Options &
    options() { return options_; }
    
    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }
    
    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }
    
    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }
    
    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }
    
    inline
    bool
    has_option(String const & name) { return options_.has_option(name); }
    
    template<typename T>
    void
    set_option_value(
        String const & name,
        T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }
    
public:
    inline
    void
    beads(Points const & beads) { beads_ = beads; }
    
    inline
    void
    lookup(
        StericLookup const & sl) {
        using_lookup_ = 1;
        lookup_ = sl;
    }
    
    inline
    void
    scorer(MotifStateSearchScorerOP const & scorer) { scorer_ = scorer; }
    
    inline
    void
    selector(MotifStateSelectorOP const & selector) { selector_ = selector; }
    
    void
    update_var_options();
    
protected:
    
    void
    setup_options();
    
private:
    
    MotifStateSearchNodeOP
    _start_node(
        BasepairStateOP const &);
    
    MotifStateSearchSolutionOP
    _search();
    
    
private:
    MotifStateSearchNodeQueue queue_;
    MotifStateSelectorOP selector_;
    MotifStateSearchScorerOP scorer_;
    MotifStateSearchNodeOP test_node_;
    MotifStateSearchSolutionOPs solutions_;
    MotifStateAligner aligner_;
    MotifStateandTypes possible_children_;
    StericLookup lookup_;
    Points beads_;
    int no_more_solutions_;
    Options options_;
    //options
    bool sterics_, verbose_, helix_end_;
    int max_node_level_, min_size_, max_size_, max_solutions_;
    int sol_count_, min_node_level_;
    float accept_score_, min_ss_score_, max_steps_;
    int using_lookup_;
    
};

#endif /* defined(__RNAMake__motif_state_search__) */
