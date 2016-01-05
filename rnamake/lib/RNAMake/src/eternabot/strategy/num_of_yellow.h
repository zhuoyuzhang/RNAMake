//
//  num_of_yellow.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__num_of_yellow__
#define __RNAMake__num_of_yellow__

#include <stdio.h>

#include "eternabot/strategy.h"

namespace eternabot {

class NumofYellowNucleotidesperLengthofString : public Strategy {
public:
    NumofYellowNucleotidesperLengthofString() {
        params_ = std::vector<float>(1);
        params_[0] = 791.641998291;
        upper_length_ = Ints(10);
        upper_length_[0] = 0; upper_length_[1] = 1; upper_length_[2] = 2;
        upper_length_[3] = 2; upper_length_[4] = 2; upper_length_[5] = 3;
        upper_length_[6] = 3; upper_length_[7] = 4; upper_length_[8] = 5;
        upper_length_[9] = 4;
        lower_length_ = Ints(10);
        lower_length_[0] = 0; lower_length_[1] = 0; lower_length_[2] = 0;
        lower_length_[3] = 1; lower_length_[4] = 1; lower_length_[5] = 1;
        lower_length_[6] = 2; lower_length_[7] = 3; lower_length_[8] = 2;
        lower_length_[9] = 1;
        mean_ = 91.2420911348;
        stdev_ = 12.5663926344;
    }
    
    ~NumofYellowNucleotidesperLengthofString() {}
    
    float
    score(FeaturesOP const & features) {
        
        float penalty = 0, count = 0;
        int stack_length;
        float yellow_count = 0;
        for(auto const & helix : features->helices) {
            if(helix->basepairs().size() < 2 || helix->basepairs().size() > 9) { continue; }
            stack_length = (int)helix->basepairs().size();
            count ++;
            for(auto const & bp : helix->basepairs()) {
                //is a bp of AU or UA
                if(sstruct::is_au_pair(bp)) { yellow_count ++; }
            }
            if     (upper_length_[stack_length] < yellow_count) {
                penalty += yellow_count - upper_length_[stack_length];
            }
            else if(lower_length_[stack_length] > yellow_count) {
                penalty += lower_length_[stack_length] - yellow_count;
            }
            
        }
        if(count == 0) { return 100; }
        
        return 100 - params_[0] * penalty/float(features->length);
        
    }
    
    
    
private:
    Ints upper_length_, lower_length_;
};
    
}



#endif /* defined(__RNAMake__num_of_yellow__) */
