/*
*  \author: Breno Carvalho
*/

#ifndef BC_BAYESIAN_NETWORK_HPP
#define BC_BAYESIAN_NETWORK_HPP

#include <graphlab.hpp>

#include <../../../toolkits/graphical_models/factors/factor_graph.hpp>
#include <../../../toolkits/graphical_models/factors/bp_vertex_program.hpp>
#include <../../../toolkits/graphical_models/factors/bp_graph_data.h>
#include <../../../apps/bnets/bnet_utils.h>
// Include the macro for each operation


const size_t MAX_DIM1 = 10;
template class InferenceEngine<MAX_DIM1>;

typedef graphlab::discrete_variable      variable_t;

void test1(int argc, char** argv, graphlab::distributed_control& dc){

    //TODO explain that this block is important for running the bp algorithm
    graphlab::command_line_options clopts("Run Loopy BP on a Network");
    clopts_vals clvals;
    if( setup_cli(clopts, clvals, argc, argv) != EXIT_SUCCESS ) return;


    //Initiates a Bayesian Network
    BeliefNetwork<MAX_DIM1> bnet;

    //Set the values for the variables
    std::vector<std::string> values(2, "");
    values[0] = "true"; values[1] = "false";

    //Criate first variable, rain
    std::vector<double> cpd_r(2, 0.0);
    cpd_r[0] = log(0.4); cpd_r[1] = log(0.6);
    std::string rain("rain");
    variable_t r =  bnet.add_variable(values, cpd_r, rain);

    //Criate second variable, wet
    std::cout << "Creating wet variable"<< std::endl;
    std::vector<double> cpd_w(4, 0.0);
    std::string wet("wet");
    variable_t w = bnet.add_variable(values, wet);
    bnet.add_variable_parent(rain, wet);
    cpd_w[0] = log(0.9); cpd_w[2] = log(0.1);
    cpd_w[1] = log(0.3); cpd_w[3] = log(0.7);
    bnet.set_variable_cpd(cpd_w, wet);

    //Run inference
    BeliefPropagationEngine<MAX_DIM1> eng;
    std::cout << std::endl << bnet << std::endl;

    //bnet.add_evidence("rain", "false");
    std::cout << "Inference of 'wet' " << bnet.infer_variable_cpd(eng, wet, dc, clvals, clopts) << std::endl;

    std::cout << std::endl << bnet << std::endl;

    std::cout << "Rain true : " << exp(bnet.fgraph.belief_for_variable(r).logP(0)) << std::endl;
    std::cout << "Rain false: " << exp(bnet.fgraph.belief_for_variable(r).logP(1)) << std::endl;

    std::cout << "Wet true : "  << exp(bnet.fgraph.belief_for_variable(w).logP(0)) << std::endl;
    std::cout << "Wet false: "  << exp(bnet.fgraph.belief_for_variable(w).logP(1)) << std::endl;
}

void test2(int argc, char** argv, graphlab::distributed_control& dc){
    //TODO explain that this block is important for running the bp algorithm
    graphlab::command_line_options clopts("Run Loopy BP on a Network");
    clopts_vals clvals;
    if( setup_cli(clopts, clvals, argc, argv) != EXIT_SUCCESS ) return;


    //Initiates a Bayesian Network
    BeliefNetwork<MAX_DIM1> bnet;

    //Set the values for the variables
    std::vector<std::string> values(2, "");
    values[0] = "true"; values[1] = "false";

    //Criate first variable, rain
    std::vector<double> cpd_r(2, 0.0);
    cpd_r[0] = 0; cpd_r[1] = 0;
    std::string rain("rain");
    variable_t r =  bnet.add_variable(values, cpd_r, rain);

    //Criate second variable, wet
    std::cout << "Creating wet variable"<< std::endl;
    std::vector<double> cpd_w(4, 0.0);
    std::string wet("wet");
    variable_t w = bnet.add_variable(values, wet);
    bnet.add_variable_parent(rain, wet);
    // cpd_w[0] = 0; cpd_w[4] = 0;
    // cpd_w[1] = 0; cpd_w[5] = 0;
    // cpd_w[2] = 0; cpd_w[6] = 0;
    // cpd_w[3] = 0; cpd_w[7] = 0;
    bnet.set_variable_cpd(cpd_w, wet);

    //Criate third variable, alarm
    std::cout << "Creating alarm variable"<< std::endl;
    std::vector<double> cpd_a(4, 0.0);
    std::string alarm("alarm");
    variable_t a = bnet.add_variable(values, alarm);
    bnet.add_variable_parent( wet, alarm);
    // cpd_a[0] = 0; cpd_a[2] = 0;
    // cpd_a[1] = 0; cpd_a[3] = 0;
    bnet.set_variable_cpd(cpd_a, alarm);

    DataSet data("test.csv");
    std::cout << data << std::endl;
    std::cout << bnet << std::endl;
    bnet.learn_parameters(dc, data);
    bnet.normalize(dc);

    std::cout << std::endl << bnet << std::endl;
}

int main(int argc, char** argv){
    graphlab::mpi_tools::init(argc, argv);
    // Create a distributed control object (must come after mpi_tools::init())
    graphlab::distributed_control dc;

    std::cout << std::endl << "TEST 1"<<std::endl << "=================" << std::endl;
    test1(argc, argv, dc);

    std::cout << std::endl << "TEST 2"<<std::endl << "=================" << std::endl;
    test2(argc, argv, dc);

    graphlab::mpi_tools::finalize();
};

#endif //BC_BAYESIAN_NETWORK_HP
