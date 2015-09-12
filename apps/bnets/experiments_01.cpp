/*
 *  \author: Breno Carvalho
 */

#include <../../../apps/bnets/csv_reader.h>
#include <../../../apps/bnets/fgraph_utils.h>

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>


// If a test runs ok it returns 0, it returns anything else otherwise
// Test functions start with test
// Auxiliary test functions starts with aux_test, they are used to run subtests

//%%%%%%%%%%%%%%%%%%%%%% Factor tests %%%%%%%%%%%%%%%%%%%%%%%%%%
std::vector<std::string> values_falseTrue(){
    std::vector< std::string > bin(2, "false");
    bin[1] = "true";
    return bin;
}

std::vector<std::string>values_lowMediumHigh(){
    std::vector< std::string > tri(3, "low");
    tri[1] = "medium"; tri[2] = "high";
    return tri;
}

// ( rain | false true )
RandomVariable new_rain_var(){
    return RandomVariable("rain", values_falseTrue());
};

// ( wet | false true maybe )
RandomVariable new_grass_wet_var(){
    std::vector<std::string> triple(3 ,"");
    triple[0] = "false"; triple[1] = "true"; triple[2] = "maybe";

    return RandomVariable("wet", triple);
};

// [< ( rain | false true ) > 0.5 0.5 ]
Factor new_rain_factor(){
    scope_type scope(1, new_rain_var());
    cpd_type cpd(2, 0.5);
    return Factor(scope, cpd);
};

// [< ( wet | false true maybe ) ( rain | false true ) > 0.3 0.5 0.1 0.9 0.9 0.1 ]
Factor new_grass_wet_factor(){
    scope_type scope(2, new_rain_var());
    scope[0] = new_grass_wet_var();
    cpd_type cpd(6, 0.0);
    //rain = false,  rain = true
    cpd[0] = 0.3; cpd[1] = 0.5; // wet = false
    cpd[2] = 0.1; cpd[3] = 0.9; // wet = true
    cpd[4] = 0.9; cpd[5] = 0.1; // wet = maybe
    return Factor(scope, cpd);
};

Factor new_a_factor(){
    RandomVariable A("A", values_falseTrue());
    scope_type a_scope  (1, A);
    cpd_type a_cpd(2, 0.0);
    a_cpd[0] = 0.2; a_cpd[1] = 0.8;

    return Factor(a_scope, a_cpd);
}

Factor new_b_factor(){
    RandomVariable B("B", values_falseTrue());
    scope_type b_scope  (1, B);
    cpd_type b_cpd(2, 0.0);
    b_cpd[0] = 0.4; b_cpd[1] = 0.6;

    return Factor(b_scope, b_cpd);
}

Factor new_cab_factor(){
    RandomVariable A("A", values_falseTrue());
    RandomVariable B("B", values_falseTrue());
    RandomVariable C("C", values_lowMediumHigh());

    scope_type cab_scope(3, C);

    cab_scope[1] = A;
    cab_scope[2] = B;
    cpd_type cab_cpd(2*2*3, 0.0);
    // (A = false B = false) (A = false B = true) (A = true B = false) (A = true B = true)
    cab_cpd[0] = 0.1;     cab_cpd[1] = 0.2;     cab_cpd[2]  = 0.3;     cab_cpd[3]  = 0.4; //C = low
    cab_cpd[4] = 0.2;     cab_cpd[5] = 0.28;    cab_cpd[6]  = 0.24;    cab_cpd[7]  = 0.28; //C = medium
    cab_cpd[8] = 0.2;     cab_cpd[9] = 0.2;     cab_cpd[10] = 0.2;     cab_cpd[11] = 0.4; //C = high

    return Factor(cab_scope, cab_cpd);
}

Factor new_db_factor(){
    RandomVariable B("B", values_falseTrue());
    RandomVariable D("D", values_lowMediumHigh());
    scope_type db_scope (2, D);

    db_scope[1] = B;
    cpd_type db_cpd(2*3, 0.0);

    //B = false         B = true
    db_cpd[0] = 0.9;    db_cpd[1] = 0.1; //D = low
    db_cpd[2] = 0.8;    db_cpd[3] = 0.2; //D = medium
    db_cpd[4] = 0.4;    db_cpd[5] = 0.6; //D = high

    return Factor(db_scope, db_cpd);
}

/*
 * Tests the Bnet creation and inference
 */
int test1(graphlab::distributed_control& dc){
    std::vector<std::string> binary(2 ,"");
    binary[0] = "false"; binary[1] = "true";

    std::vector<std::string> triple(3 ,"");
    triple[0] = "false"; triple[1] = "true"; triple[2] = "maybe";

    std::string rain_str = std::string("rain");
    RandomVariable rain(rain_str, binary);

    scope_type scope = scope_type(1, rain);
    cpd_type values(2, 0.0);
    values[0] = 0.4; values[1] = 0.6;
    Factor f1(scope, values);

    std::string grass_str = std::string("grass is wet");
    RandomVariable grass(grass_str, triple);

    scope_type scope_grass = scope_type(1, grass);
    cpd_type values_grass(3, 0.0);
    values_grass[0] = 0.3; values_grass[1] = 0.5; values_grass[2] = 0.2;

    Factor f2(scope_grass, values_grass);

    scope_type scope_rain_grass = scope_type(2, rain);
    scope_rain_grass[1] = grass;
    cpd_type values_rain_grass(6, 0.0);
    values_rain_grass[0] = 0.3; values_rain_grass[1] = 0.5;
    values_rain_grass[2] = 0.1; values_rain_grass[3] = 0.9;
    values_rain_grass[4] = 0.9; values_rain_grass[5] = 0.1;
    Factor f3(scope_rain_grass, values_rain_grass);

    Vertice v1(rain);
    Vertice v2(grass);
    Vertice a(f1);
    Vertice b(f2);
    Vertice c(f1*f2);
    std::cout << a     << std::endl \
    << v1    << std::endl \
    << f2    << std::endl \
    << f1*f2 << std::endl << (f1*f2).marginalizeAllBut(rain) << std::endl \
    << "V1: "<< v1.getInfoRandomVariable() << std::endl \
    << "a: " << a.getMessage() <<std::endl;

    Graph g(dc);
    g.addVertice(v1);
    g.addVertice(v2);
    g.addVertice(a);
    g.addVertice(b);
    g.addVertice(c);
    std::cout << "variables loaded" << std::endl;
    g.addEdge(a, v1);
    g.addEdge(b, v2);
    g.addEdge(c, v1);
    g.addEdge(c, v2);
    std::cout << "edges loaded" << std::endl;
    BeliefProp bp;
    std::cout << "bp constructed" << std::endl;
    g.startDGraph();
    g.runInference(bp);
    dc.barrier();

    return 0;
}

/*
 * Tests the Bnet creation and learning with complete data
 */
int test2(graphlab::distributed_control& dc){
    std::cout << "creating bn"<< std::endl;
    BeliefProp i_eng;
    BeliefNetwork bn(dc, &i_eng);
    std::cout << "bn created"<< std::endl;

    std::vector<std::string> binary(2 ,"");
    binary[0] = "false"; binary[1] = "true";

    std::vector<std::string> triple(3 ,"");
    triple[0] = "false"; triple[1] = "true"; triple[2] = "maybe";

    std::string rain_str = std::string("rain");
    RandomVariable rain(rain_str, binary);

    scope_type scope = scope_type(1, rain);
    cpd_type values(2, 0.0);
    values[0] = 0.1; values[1] = 0.1;
    Factor f1(scope, values);

    std::string grass_str = std::string("wet");
    RandomVariable grass(grass_str, triple);

    scope_type scope_rain_grass = scope_type(2, rain);
    scope_rain_grass[0] = grass;
    cpd_type values_rain_grass(6, 0.0);
    //values_rain_grass[0] = 0.3; values_rain_grass[1] = 0.5;
    //values_rain_grass[2] = 0.1; values_rain_grass[3] = 0.1;
    //values_rain_grass[4] = 0.9; values_rain_grass[5] = 0.1;
    Factor f2(scope_rain_grass, values_rain_grass);

    std::cout << "variables created"<< std::endl;

    BNet_vertice r = bn.addVariable(rain, f1);
    BNet_vertice g = bn.addVariable(grass, f2);

    std::cout << "variables inserted"<< std::endl;
    bn.addEdge(r,g);
    std::cout << "edge created"<< std::endl;

    DataSet data("/home/breno/Documents/git-repos/graphlab/apps/bnets/samples/test.csv");
    std::cout << "DATASET read" << std::endl;
    bn.learnParameters(data);
    std::cout << "learning done"<< std::endl;
    bn.save();
    std::cout << "saving done"<< std::endl;
    //std::cout << "P(wet) = " << bn.infer(grass) << std::endl;
    return 0;
}

/*
 * Tests if it loads a Bnet from file and learning with complete data
 */
int test3(graphlab::distributed_control& dc){
    dc.barrier();
    std::cout << "creating bn"<< std::endl;
    BeliefProp i_eng;
    BeliefNetwork bn(dc, &i_eng);
    std::cout << "bn created"<< std::endl;

    bn.load("/home/breno/Documents/git-repos/graphlab/apps/bnets/samples/test.txt");
    bn.save();
    std::vector<std::string> triple(3 ,"");
    triple[0] = "false"; triple[1] = "true"; triple[2] = "maybe";
    std::string grass_str = std::string("wet");
    RandomVariable grass(grass_str, triple);

    DataSet data("/home/breno/Documents/git-repos/graphlab/apps/bnets/samples/test.csv");
    std::cout << "DATASET read" << std::endl;
    bn.learnParameters(data);
    std::cout << "learning done"<< std::endl;
    //bn.save();
    //std::cout << "saving done"<< std::endl;
    //std::cout << "P(wet) = " << bn.infer(grass) << std::endl;
    return 0;
}

namespace exp_learning{
    /*
     * Testing for ONE network and ONE data set
     * Learns the parameters of the bayesian network specified in bnet_file_path using the data in data_path.
     * Returns the time elapsed **during learning** in seconds, not considerating network creation
     */
    double complete_data(graphlab::distributed_control& dc, BeliefNetwork& bn, DataSet& data){
        double begin, end;
        begin = clock();
        bn.learnParameters(data);
        end = clock();

        std::cout << "learning done"<< std::endl;
        dc.barrier();
        //bn.save(bnet_file_path+"_new");
        std::cout << "saving done"<< std::endl;
        return (end-begin) / CLOCKS_PER_SEC;
    }

    double complete_data(graphlab::distributed_control& dc, std::string bnet_file_path, DataSet& data){
        std::cout << "creating bn"<< std::endl;
        BeliefProp i_eng;
        BeliefNetwork bn(dc, &i_eng);
        bn.load(bnet_file_path);
        std::cout << "BNet read" << std::endl;
        return complete_data(dc, bn, data);
    }

    double complete_data(graphlab::distributed_control& dc, std::string bnet_file_path, std::string data_path){
        DataSet data(data_path);
        std::cout << "DATASET read" << std::endl;
        return complete_data(dc, bnet_file_path, data);
    }

    /*
     * Testing for many bayesian networks and many data sets, saving results in the file in out_str
     * bnet_str: bayesian network file's name
     * data_set_str: data set file's name
     * in_str: input file's name
     * out_str: output file's name
     * iterations: how many times the experiment is going to be run on the same data set using the same network
     * data_steps: the algorithm will repeat the experiment in data sets whose size grows as multiple of this variable
     * multiple_files: boolean variable that indicates the program to look for similar files to build the network
     */
    int complete_data(graphlab::distributed_control& dc, std::string bnet_str, std::string data_set_str, std::string in_str, std::string out_str, int iterations, int data_steps = 100, bool multiple_files = true){
        dc.barrier();

        //initialize: open the file
        double elapsed;
        std::ofstream output;
        output.open(out_str.c_str(), std::ios::out | std::ios::app);
        std::cout << "Data: "+in_str+data_set_str+".csv" << std::endl;
        DataSet data(in_str+data_set_str+".csv");
        int data_size = data_steps,
            max_data_size = data.getLinesQtd();
        std::cout << "FULL DATA READ " << max_data_size << "Lines";
        //Loop to run experiments
        for(; data_size <= max_data_size; data_size += data_steps){
            std::cout << "Start of iteration"<< std::endl;
            std::stringstream bnet_ss;
            //If the bayesian network is split across multiple files, take this information to GraphLab
            if(multiple_files){
                bnet_ss << bnet_str << dc.numprocs() << ".txt";
            }else{
                bnet_ss << bnet_str << ".txt";
            }

            //Repeat the complete data learning "iterations" times
            for(int i = 0;  i < iterations; i++){
                DataSet data_cut = data.head(data_size);
                elapsed = complete_data(dc, in_str+bnet_ss.str(), data_cut);
                std::cout << "complete_data learning done"<< std::endl;
                if(dc.procid() == 0){
                    output << bnet_str << "\t" << data_set_str << "\t" << data_size << "\t" << dc.numprocs() << "\t" << elapsed << std::endl;
                }
                output.flush();
            }
        }
        output.close();
        return 0;
    }
}

int main(int argc, char** argv){
    graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;
        std::string file_name = "exp_temp",
                    path = "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_net_fake_data/",
                    line, net_name;
        std::ifstream input_f("/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_list_complete_data.txt");
        if(input_f.is_open()){
            while(getline(input_f, line)){
                std::istringstream line_stream(line);
                if(line_stream){
                    line_stream >> net_name;
                    std::cout << net_name << "\n";
                    exp_learning::complete_data(dc, net_name, "fake_data_"+net_name, path,
                                                path+"results/"+net_name+".csv",
                        //100 iterations of the experiment for datasets that have at most 100000 lines in increases of 100
                                                100, 100, false);
                }
            }
            input_f.close();
        }

    std::cout << "closing mpi" << std::endl;
    dc.cout() << "num of procs: " << dc.numprocs() << std::endl;
    graphlab::mpi_tools::finalize();
    std::cout << "mpi closed" << std::endl;
}
