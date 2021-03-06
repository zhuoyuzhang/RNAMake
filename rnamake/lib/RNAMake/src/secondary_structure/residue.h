//
//  residue.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_residue__
#define __RNAMake__ss_residue__

#include <stdio.h>
#include <sstream>
#include <memory>
#include <stdexcept>
#include <cassert>

#include "base/types.h"
#include "base/string.h"
#include "util/uuid.h"

namespace sstruct {
    
class SecondaryStructureException : public std::runtime_error {
public:
    SecondaryStructureException(
        String const & message) :
    std::runtime_error(message)
    {}
    
};
    
    
class Residue {
public:
    inline
    Residue(
        String const & name,
        String const & dot_bracket,
        int const & num,
        String const & chain_id,
        Uuid const & uuid,
        String const & i_code=""):
    name_(name),
    dot_bracket_(dot_bracket),
    num_(num),
    chain_id_(chain_id),
    uuid_(uuid),
    i_code_(i_code) {
        if     (name_ == "A") { res_type_ = 0; }
        else if(name_ == "C") { res_type_ = 1; }
        else if(name_ == "G") { res_type_ = 2; }
        else if(name_ == "U") { res_type_ = 3; }
        else if(name_ == "T") { res_type_ = 3; }
        else if(name_ == "N") { res_type_ = -1; }
        else {
            throw SecondaryStructureException(
                "in sstruct::Residue encountered a unknown name: "  + name_);
        }
        
    }
    
    inline
    Residue(
        Residue const & r):
    name_(r.name_),
    dot_bracket_(r.dot_bracket_),
    num_(r.num_),
    chain_id_(r.chain_id_),
    uuid_(r.uuid_),
    i_code_(r.i_code_),
    res_type_(r.res_type_)
    {}
    
    Residue(
        String const & s) {
        
        Strings spl = split_str_by_delimiter(s, ",");
        if(spl.size() < 4) {
            throw SecondaryStructureException("cannot build sstruct::Residue from str: " + s);
        }
        
        name_         = spl[0];
        dot_bracket_  = spl[1];
        num_          = std::stoi(spl[2]);
        chain_id_     = spl[3];
        uuid_         = Uuid();
        if(spl.size() == 5) {
            i_code_ = spl[4];
        }
        
        if     (name_ == "A") { res_type_ = 0; }
        else if(name_ == "C") { res_type_ = 1; }
        else if(name_ == "G") { res_type_ = 2; }
        else if(name_ == "U") { res_type_ = 3; }
        else if(name_ == "T") { res_type_ = 3; }
        else if(name_ == "N") { res_type_ = -1; }
        else {
            throw SecondaryStructureException(
                "in sstruct::Residue encountered a unknown name: " + name_);
        }

    }
    
    ~Residue() {}
    
public:

    inline
    String
    to_str() {
        std::stringstream ss;
        ss << name_ << "," << dot_bracket_ << "," << num_ << "," << chain_id_ << "," << i_code_;
        return ss.str();
    }
    
public: //getters
    
    inline
    String const &
    name() { return name_; }
    
    inline
    String const &
    dot_bracket() { return dot_bracket_; }
    
    inline
    int const &
    num() { return num_; }
    
    inline
    String const &
    chain_id() { return chain_id_; }
    
    inline
    String const &
    i_code() { return i_code_; }
    
    inline
    Uuid const &
    uuid() { return uuid_; }

    inline
    int
    res_type() { return res_type_; }
    
public: //setters
    
    inline
    void
    uuid(Uuid const & nuuid) { uuid_ = nuuid; }
    
    inline
    void
    name(String const & name) {
        name_ = name;
        if     (name_ == "A") { res_type_ = 0; }
        else if(name_ == "C") { res_type_ = 1; }
        else if(name_ == "G") { res_type_ = 2; }
        else if(name_ == "U") { res_type_ = 3; }
        else if(name_ == "T") { res_type_ = 3; }
        else if(name_ == "N") { res_type_ = -1; }
        else {
            throw SecondaryStructureException(
                "in sstruct::Residue encountered a unknown name: " + name_);
        }
    }
    

private:
    int num_;
    //A=0,C=1,G=2,U=3
    int res_type_;
    String name_, dot_bracket_, chain_id_, i_code_;
    Uuid uuid_;

};
    
typedef std::shared_ptr<Residue> ResidueOP;
typedef std::vector<ResidueOP> ResidueOPs;
    
struct res_less_than_key {
    inline
    bool
    operator() (ResidueOP const & r1, ResidueOP const & r2) {
        return (r1->num() < r2->num());
    }
};
    
} //sstruct

#endif /* defined(__RNAMake__ss_residue__) */
