//
//  sequence_optimizer_3d.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/30/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "util/monte_carlo.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"


SequenceOptimizer3D::SequenceOptimizer3D():
scorer_(eternabot::Scorer()),
rng_(RandomNumberGenerator()) { setup_options(); }


void
SequenceOptimizer3D::_update_designable_bp(
    DesignableBPOP const & d_bp,
    MotifStateTreeOP & mst,
    Strings const & new_bp_state) {
    
    if(d_bp->m_id_bot != nullptr ) {
        auto n = mst->get_node(*d_bp->m_id_bot);
        auto spl = split_str_by_delimiter(n->data()->name(), "=");
        auto new_name = new_bp_state[0]+new_bp_state[1]+"="+spl[1];
        mst->replace_state(n->index(), RM::instance().motif_state(new_name));
    }
    if(d_bp->m_id_top != nullptr) {
        auto n = mst->get_node(*d_bp->m_id_top);
        auto spl = split_str_by_delimiter(n->data()->name(), "=");
        auto new_name = spl[0]+"="+new_bp_state[0]+new_bp_state[1];
        mst->replace_state(n->index(), RM::instance().motif_state(new_name));
    }
    
}

SequenceOptimizer3D::OptimizedSequenceOPs
SequenceOptimizer3D::get_optimized_sequences(
    MotifTreeOP const & mt,
    BasepairOP const & target_bp,
    int ni,
    int ei) {
    
    auto sols = OptimizedSequenceOPs();
    auto ss = mt->designable_secondary_structure();
    
    auto designable_bps = DesignableBPOPs();
    for(auto const & bp : ss->basepairs()) {
        auto bp_name = bp->res1()->name()+bp->res2()->name();
        if(bp_name == "NN") {
            designable_bps.push_back(
                std::make_shared<DesignableBP>(bp));
        }
    }
    
    auto possible_bps = std::vector<Strings>({{"A", "U"}, {"U", "A"}, {"G", "C"}, {"C", "G"}});
    
    for(auto const & d_bp : designable_bps) {
        for(auto const & m : ss->motifs()) {
            if(m->name() != "HELIX.IDEAL") { continue; }
            if(m->ends()[0] == d_bp->bp) { d_bp->m_id_bot = std::make_shared<Uuid>(m->id()); }
            if(m->ends()[1] == d_bp->bp) { d_bp->m_id_top = std::make_shared<Uuid>(m->id()); }
        }
        
        auto state = possible_bps[rng_.randrange(possible_bps.size())];
        d_bp->bp->res1()->name(state[0]);
        d_bp->bp->res2()->name(state[1]);
    }
    
    auto mst = std::make_shared<MotifStateTree>(mt);
    for(auto const & m : ss->motifs()) {
        if(m->name() != "HELIX.IDEAL") { continue; }
        auto n = mst->get_node(m->id());
        auto name = m->ends()[0]->res1()->name()+m->ends()[0]->res2()->name()+"=";
        name +=     m->ends()[1]->res1()->name()+m->ends()[1]->res2()->name();
        mst->replace_state(n->index(), RM::instance().motif_state(name));
        auto n2 = mst->get_node(m->id());
    }
    
    auto target_state = target_bp->state();
    auto last_score = target_state->diff(mst->get_node(ni)->data()->get_end_state(ei));
    auto new_score = 0.0f;
    auto eterna_score = 0.0f;
    
    scorer_.setup(ss);
    auto mc = MonteCarlo(1.0f);
    
    auto best = 1000.0f;
    
    int i = -1;
    auto d_bp = DesignableBPOP(nullptr);
    auto new_bp_state = Strings();
    while(i < 1000) {
        i++;

        d_bp = designable_bps[rng_.randrange(designable_bps.size())];
        new_bp_state = possible_bps[rng_.randrange(possible_bps.size())];
        d_bp->update_state(new_bp_state);
        
        _update_designable_bp(d_bp, mst, new_bp_state);
        
        new_score = target_state->diff(mst->get_node(ni)->data()->get_end_state(ei));

        if(mc.accept(last_score, new_score)) {
            last_score = new_score;
        }
        else {
            _update_designable_bp(d_bp, mst, d_bp->last_state);
            d_bp->revert_state();
            continue;
        }
                
        if(best > new_score) {
            best = new_score;
            if(verbose_) {
                std::cout << "SEQUENCE OPTIMIZER: best_score=" << best << std::endl;
            }
        }
        
        if(cutoff_ < new_score) { continue; }
        
        eterna_score = scorer_.score_secondary_structure(ss);
        if(eterna_score > eterna_cutoff_) {
            
            sols.push_back(std::make_shared<OptimizedSequence>(
                            OptimizedSequence{ss->sequence(), new_score, eterna_score}));
            
            if(verbose_) {
                std::cout << "SEQUENCE OPTIMIZER: found solution! score=" << new_score;
                std::cout << " eterna_score=" << eterna_score << std::endl;
            }
            
            if(sols.size() >= solutions_) { return sols; }
        }
        
    }
    
    
    
    
    
    return sols;
    
}


void
SequenceOptimizer3D::setup_options() {
    options_.add_option("cutoff", 5.0f, OptionType::FLOAT);
    options_.add_option("solutions", 10, OptionType::INT);
    options_.add_option("eterna_cutoff", -1.0f, OptionType::FLOAT);
    options_.add_option("verbose", true, OptionType::BOOL);
    options_.lock_option_adding();
    update_var_options();
}

void
SequenceOptimizer3D::update_var_options() {
    cutoff_         = options_.get_float("cutoff");
    solutions_      = options_.get_int("solutions");
    eterna_cutoff_  = options_.get_float("eterna_cutoff");
    verbose_        = options_.get_bool("verbose");
}
















