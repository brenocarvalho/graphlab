/*
*  \author: Breno Carvalho
*/

#ifndef BC_BAYESIAN_NETWORK_HPP
#define BC_BAYESIAN_NETWORK_HPP

#include <cstdlib>
#include <cassert>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>

//#include <opencv2/opencv.hpp>

#include <graphlab.hpp>

#include <../../../toolkits/graphical_models/factors/factor_graph.hpp>
#include <../../../toolkits/graphical_models/factors/bp_vertex_program.hpp>

// Include the macro for each operation
#include <graphlab/macros_def.hpp>

//add functionalities to factor graph to allow evidence and change of the structure

class InferenceEngine{
    //Prepare the methods for a genral inference engine and consider moving that to another file
};

template<size_t MAX_DIM>
class BeliefNetwork {
    typedef graphlab::dense_table<MAX_DIM>       dense_table_t;
    typedef graphlab::discrete_domain<MAX_DIM>   domain_t;
    typedef graphlab::discrete_variable          variable_t;
private:
    belief_prop::factor_graph<MAX_DIM> fgraph;
    bool hasChanged;
    std::map<std::string, std::string> evidence;
public:
    BeliefNetwork(): hasChanged(false){};

    void add_variable(std::vector<std::string> values, std::vector<double> cpd, std::string variable_id){};
    void add_variable(std::vector<std::string> values, std::string variable_id){};

    void set_variable_parents(std::vector<std::string> parent_ids){};
    void set_variable_cpd(std::vector<double> cpd, variable_t id){};
    void set_variable_cpd(dense_table_t cpd, variable_t id){};
    dense_table_t get_variable_cpd(variable_t id){return NULL;};

    void set_evidence(std::map<std::string, std::string> evidence){this->evidence(evidence);};
    void add_evidence(std::string id, std::string value){};
    void clear_evidence(){};
    std::string get_evidence(std::string variable_id){
        return evidence;
    };
    std::map<std::string, std::string> get_evidence();

    dense_table_t infer_variable_cpd(InferenceEngine eng, std::string variable_id){return NULL;};
    dense_table_t infer_variables_cpd(InferenceEngine eng, std::vector<std::string> variable_ids){return NULL;};
    dense_table_t infer_all_variables_cpd(InferenceEngine eng){return NULL;};

    //Use SFrames to get info to set your parameters (complete data)
    //Use SFrames to get info to set your parameters (incomplete data) <- Use an engine

    //Use SFrames to get info to set your structure(complete data)
    //Use SFrames to get info to set your structures (incomplete data) <- Use an engine

    //Save belief network to file
    //Load belief network from file

    friend std::ostream& operator<<(std::ostream& os, const BeliefNetwork<MAX_DIM>& bn){
        return os<<"Belief Network:"<< std::endl <<
                    "has" << ((bn.hasChanged)? " ":" not ")<< "changed";
    };
};

int main(){
    BeliefNetwork<10> bnet;
    std::vector<std::string> values(2, "");
    values[0] = "true"; values[1] = "false";
    std::vector<double> cpd(4, 0.0);
    cpd[0] = 0.1; cpd[2] = 0.9;
    cpd[1] = 0.4; cpd[3] = 0.6;
    bnet.add_variable(values, cpd, "a");

    std::cout << bnet<< std::endl;
};

#endif //BC_BAYESIAN_NETWORK_HPP
