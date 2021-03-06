//
//  segmenter.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__segmenter__
#define __RNAMake__segmenter__

#include <stdio.h>
#include <queue>
#include <map>
#include "structure/residue.h"
#include "motif/motif.h"
#include "motif/motif_factory.h"

/*
 * Exception for Segmenter
 */
class SegmenterException : public std::runtime_error {
public:
    /**
     * Standard constructor for SegmenterException
     * @param   message   Error message for segmenter
     */
    SegmenterException(String const & message):
    std::runtime_error(message)
    {}
};

struct Pair {
    Pair(
        ResidueOP const & nres1,
        ResidueOP const & nres2,
        int ndist):
        res1(nres1),
        res2(nres2),
        dist(ndist)
    {}
        
    int
    inline
    contains(
        ResidueOP const & res) const {
        
        if(res == res1 || res == res2) { return 1; }
        else                           { return 0; }
    }
    
    ResidueOP res1, res2;
    int dist;
};

struct Segments {
    inline
    Segments(
        MotifOP const & nremoved,
        MotifOP const & nremaining):
        removed(nremoved),
        remaining(nremaining) {}
    
    MotifOP removed, remaining;
};

typedef std::shared_ptr<Pair> PairOP;
typedef std::vector<PairOP>   PairOPs;
typedef std::shared_ptr<Segments> SegmentsOP;

struct PairSearchNode {
    PairSearchNode(
        PairOPs const & npairs):
    pairs(npairs),
    score(0) {
        for(auto const & p : pairs) {
            score += p->dist;
        }
    
    }
    
    ~PairSearchNode() {}
    
    inline
    int
    contains(
        ResidueOP const & res) const {
        
        for(auto const & p : pairs) {
            if(p->contains(res)) { return 1; }
        }
        return 0;
    }
    
    inline
    int
    contains_pair(
        PairOP const & pair) {
        
        for(auto const & p : pairs) {
            if(pair == p) { return 1; }
        }
        return 0;

    }
    
    PairOPs pairs;
    int score;
};

typedef std::vector<PairSearchNode> PairSearchNodes;

struct PairSearchNodeCompare {
    bool
    operator () (
        PairSearchNode const & node1,
        PairSearchNode const & node2) {
        
        if (node1.score > node2.score) { return true;  }
        else     					   { return false; }
    }
};

typedef std::priority_queue<PairSearchNode, PairSearchNodes, PairSearchNodeCompare> PairSearchNodePriortyQueue;

class PairSearch {
public:
    PairSearch():
    queue_(PairSearchNodePriortyQueue()),
    solutions_(PairSearchNodes()),
    values_(std::map<ResidueOP, float>())
    {}
    
    ~PairSearch() {}
    
    
public:
    
    PairSearchNodes const &
    search(
        ResidueOPs const & res,
        PairOPs const & pairs,
        PairOPs const & end_pairs) {
        res_ = res;
        _get_default_values(res, pairs, end_pairs);
    
        auto all_pairs = PairOPs();
        for(auto const & p : pairs) { all_pairs.push_back(p); }
        for(auto const & p : end_pairs) { all_pairs.push_back(p); }
        
        for(auto const & p : pairs) {
            auto n = PairSearchNode(PairOPs{p});
            n.score += _get_estimated_score(n);
            queue_.push(n);
        }
        
        auto current = PairSearchNode(PairOPs{});
        while(!queue_.empty()) {
            current = queue_.top();
            queue_.pop();
            
            for(auto p : all_pairs) {
                if(current.contains(p->res1) || current.contains(p->res2)) { continue; }
                auto pairs = current.pairs;
                pairs.push_back(p);
                auto node = PairSearchNode(pairs);
                auto estimated_score = _get_estimated_score(node);
                if(estimated_score == 0) {
                    solutions_.push_back(node);
                    if(solutions_.size() > 10) {
                        return solutions_;
                    }
                }
                else {
                    node.score += _get_estimated_score(node);
                    queue_.push(node);
                }
            }
            
        }
        
        return solutions_;
    
    }
    
private:
    void
    _get_default_values(
        ResidueOPs const & res,
        PairOPs const & pairs,
        PairOPs const & end_pairs) {
        
        auto all_pairs = PairOPs();
        for(auto const & p : pairs) { all_pairs.push_back(p); }
        for(auto const & p : end_pairs) { all_pairs.push_back(p); }
        
        for(auto const & r : res) {
            auto total = 0.0, count = 0.0;
            for(auto const & p : all_pairs) {
                if(p->contains(r)) {
                    total += p->dist;
                    count += 1;
                }
            }
            values_[r] = (total / count) / 2;
        }
    }
    
    int
    _get_estimated_score(
        PairSearchNode const & n) {
        int score = 0;
        for(auto const & r : res_) {
            if(n.contains(r)) { continue; }
            score += values_[r];
        }
        return score;
        
    }
    
  
    
private:
    PairSearchNodePriortyQueue queue_;
    PairSearchNodes solutions_;
    std::map<ResidueOP, float> values_;
    ResidueOPs res_;
};

class Segmenter {
public:
    Segmenter():
    mf_(MotifFactory()){}
    
    ~Segmenter() {}
    
public:
    SegmentsOP
    apply(
        RNAStructureOP const &,
        BasepairOPs const &);
    
private:
    
    ChainOP
    _get_subchain(
        RNAStructureOP const &,
        PairOP const &);
    
    void
    _get_pairs(
        RNAStructureOP const &,
        ResidueOPs const &);
    
    SegmentsOP
    _get_segments(
        RNAStructureOP const &,
        ResidueOPs &,
        BasepairOPs const &,
        ResidueOPs const &,
        BasepairOPs const &);
    
private:
    
    MotifFactory mf_;
    PairOPs pairs_, end_pairs_;
    
};


#endif /* defined(__RNAMake__segmenter__) */
