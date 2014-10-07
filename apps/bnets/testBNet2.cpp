/*
*  \author: Breno Carvalho
*/

//#ifndef BC_BAYESIAN_NETWORK_HPP
//#define BC_BAYESIAN_NETWORK_HPP

// #include <cstdlib>
// #include <cassert>
// #include <cstring>
// #include <fstream>
// #include <sstream>
// #include <iostream>
// #include <vector>
// #include <string>
// #include <algorithm>
// #include <limits>
// #include <cmath>

#include <graphlab.hpp>

#include <../../../toolkits/graphical_models/factors/factor_graph.hpp>
#include <../../../toolkits/graphical_models/factors/bp_vertex_program.hpp>

// Include the macro for each operation
#include <graphlab/macros_def.hpp>

//add functionalities to factor graph to allow evidence and change of the structure

class InferenceEngine{
    //Prepare the methods for a general inference engine and consider moving that to another file
};

template<size_t MAX_DIM> class Variable;
template<size_t MAX_DIM> class BeliefNetwork;


template<size_t MAX_DIM>
class BeliefNetwork {
public:
    typedef graphlab::dense_table<MAX_DIM>       dense_table_t;
    typedef graphlab::discrete_domain<MAX_DIM>   domain_t;
    typedef graphlab::discrete_variable          variable_t;
private:
    belief_prop::factor_graph<MAX_DIM> fgraph;
    std::map<std::string, Variable<MAX_DIM> > variables;
    bool hasChanged;
    std::map<std::string, std::string> evidence;
public:
    BeliefNetwork(): hasChanged(false){};

    void add_variable(std::vector<std::string> values, std::vector<double> cpd, std::string variable_id){
        add_variable(values, variable_id);
    //    variables[variable_id].set_cpd(dense_table_t(cpd));
    };
    void add_variable(std::vector<std::string> values, std::string variable_id){
        Variable<MAX_DIM> var(values, variable_id);
        variables[variable_id] = var;
    };

    void set_variable_parents(std::vector<std::string> parent_ids, std::string variable_id){
        std::vector< Variable<MAX_DIM> > parents(parent_ids.size());
        for(int i=0; parent_ids.size(); i++){
            parents[i] = variables[parent_ids[i]];
        }
        variables[variable_id].set_parents(parents);
    };

    void add_variable_parent(std::string parent_id, std::string variable_id){
        variables[variable_id].add_parent(variables[parent_id]);
    };
    void set_variable_cpd(std::vector<double> cpd, std::string id){};
    void set_variable_cpd(dense_table_t cpd, std::string id){};
    dense_table_t get_variable_cpd(variable_t id){return NULL;};

    void set_evidence(std::map<std::string, std::string> evidence){
        this->evidence(evidence);
        std::map<std::string, std::string >::const_iterator iterator = evidence.begin();
        for(; iterator != evidence.end(); iterator++){
            variables[(*iterator).first].set_evidence((*iterator).second);
        }
    };
    void add_evidence(std::string id, std::string value){
        evidence[id] = value;
        variables[id].set_evidence(value);
    };
    void clear_evidence(){};
    std::string get_evidence(std::string variable_id){
        return evidence;
    };
    std::map<std::string, std::string> get_evidence();

    dense_table_t infer_variable_cpd(InferenceEngine eng, std::string variable_id){dense_table_t d; return d;};
    dense_table_t infer_variables_cpd(InferenceEngine eng, std::vector<std::string> variable_ids){dense_table_t d; return d;};
    dense_table_t infer_all_variables_cpd(InferenceEngine eng){dense_table_t d; return d;};

    //Use SFrames to get info to set your parameters (complete data)
    //Use SFrames to get info to set your parameters (incomplete data) <- Use an engine

    //Use SFrames to get info to set your structure(complete data)
    //Use SFrames to get info to set your structures (incomplete data) <- Use an engine

    //Save belief network to file
    //Load belief network from file

    friend std::ostream& operator<<(std::ostream& os, const BeliefNetwork<MAX_DIM>& bn){
        os<<"Belief Network:"<< std::endl <<
            "has" << ((bn.hasChanged)? " ":" not ")<< "changed"<< std::endl;
        //Print each of its variables
        typedef typename std::map<std::string, Variable<MAX_DIM> >::const_iterator iter_map_t;
        iter_map_t iter = bn.variables.begin() ;
        for (; iter != bn.variables.end(); iter++) {
             os << (*iter).second << std::endl;
        }
        return os;
    };
};


template<size_t MAX_DIM>
class Variable{
    //TODO finish this class
private:
    typename BeliefNetwork<MAX_DIM>::dense_table_t distribution;
    typename BeliefNetwork<MAX_DIM>::variable_t f_g_var;
    std::vector<Variable> parents;
    std::vector<std::string> values;
    int evidence;
    std::string id;

public:
    Variable(std::vector<std::string> values, std::string variable_id): values(values), id(variable_id), evidence(-1){}
    Variable():id(""), evidence(-1){};

    std::string get_id(){ return id;}

    void set_parents(std::vector<Variable> parents){
        this->parents.clear();
        for(int i = 0; i < parents.size(); i++){
            this->parents.push_back(parents[i]);
        }
        distribution = typename BeliefNetwork<MAX_DIM>::dense_table_t();
    };

    void add_parent(Variable parent){
        parents.push_back(parent);
    };

    //void set_cpd(std::vector<double> cpd, std::string id){};
    void set_cpd(typename BeliefNetwork<MAX_DIM>::dense_table_t cpd){
    //    distribution = cpd;
    };
    //dense_table_t get_cpd(variable_t id){return NULL;};

    void set_evidence(std::string value){
        int i = 0;
        for(; i < values.size(); i++){
            if(!values[i].compare(value)){ // if we have a match with a value in the values list then:
                evidence = i;
                //TODO update the factor graph
            }
        }
    }

    std::string get_evidence(){
        if(evidence > 0){
            return values[evidence];
        }
        return "";
    }

    bool has_evidence(){ return evidence > 0; }

    friend std::ostream& operator<<(std::ostream& os, const Variable<MAX_DIM>& var){
        os << "Variable{id = '" << var.id << "'";

        os << ", values = {'" << var.values[0] << "'";
        for(int i = 1; i < var.values.size(); i++){
            os << ", '" << var.values[1] << "'";
        }
        os << "}";

        if(var.parents.size() > 0){
            os << ", parents = {'"<< var.parents[0].id<<"'";
            for(int i = 1; i < var.parents.size(); i++){
                os << ", '" << var.parents[i].id << "'";
            }
            os << "}";
        }
        if(var.evidence >= 0){
            os << ", evidence = '" << var.values[var.evidence] << "'";
        }

        //os << ", cpd = <" << var.distribution << ">";
        os << "}";
        return os;
    };
};

const size_t MAX_DIM1 = 10;

int main(){
    //Initiates a Bayesian Network
    BeliefNetwork<MAX_DIM1> bnet;

    //Set the values for the variables
    std::vector<std::string> values(2, "");
    values[0] = "true"; values[1] = "false";

    //Criate first variable, rain
    std::vector<double> cpd_r(4, 0.0);
    cpd_r[0] = 0.4; cpd_r[1] = 0.6;
    bnet.add_variable(values, cpd_r, "rain");

    //Criate second variable, wet
    std::vector<double> cpd_w(4, 0.0);
    bnet.add_variable(values, "wet");
    bnet.add_variable_parent("rain", "wet");
    cpd_w[0] = 0.9; cpd_w[2] = 0.1;
    cpd_w[1] = 0.4; cpd_w[3] = 0.6;
    bnet.set_variable_cpd(cpd_w, "wet");

    //Run inference
    InferenceEngine eng;
    bnet.infer_variable_cpd(eng, "wet");
    bnet.add_evidence("rain", "true");
    std::cout << 0.0 << std::endl;

    std::cout << bnet<< std::endl;
};

//#endif //BC_BAYESIAN_NETWORK_HPP
