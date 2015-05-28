/*
 *  \author: Breno Carvalho
 */

#include <../../../apps/bnets/csv_reader.h>
#include <../../../apps/bnets/fgraph_utils.h>

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>


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
    std::cout << "saving done"<< std::endl;
    //std::cout << "P(wet) = " << bn.infer(grass) << std::endl;
    return 0;
}

namespace exp_learning{
    /*
     * Testing for ONE network and ONE data set
     * Learns the parameters of the bayesian network specified in bnet_file_path using the data in data_path.
     * Returns the time elapsed **during learning** in seconds, not considerating network creation
     */
    double complete_data(graphlab::distributed_control& dc, std::string bnet_file_path, std::string data_path){
        std::cout << "creating bn"<< std::endl;
        BeliefProp i_eng;
        BeliefNetwork bn(dc, &i_eng);
        bn.load(bnet_file_path);
        //bn.save();
        std::cout << "BNet read" << std::endl;
        DataSet data(data_path);
        std::cout << "DATASET read" << std::endl;

        double begin, end;
        begin = clock();
        bn.learnParameters(data);
        end = clock();

        std::cout << "learning done"<< std::endl;
        dc.barrier();
        //bn.save(bnet_file_path+"_new");
        std::cout << "saving done"<< std::endl;
        //std::cout << "P(wet) = " << bn.infer(grass) << std::endl;
        return (end-begin) / CLOCKS_PER_SEC;
    }

    /*
     * Testing for many bayesian networks and many data sets, saving results in the file in out_str
     */
    int complete_data(graphlab::distributed_control& dc, std::string bnet_str, std::string data_set_str, std::string in_str, std::string out_str, int iterations, int n_tests = 20, bool multiple_files = true){
        dc.barrier();

        //initialize: open the file
        double elapsed;
        std::ofstream output;
        output.open(out_str.c_str(), std::ios::out | std::ios::app);

        //Loop to run experiments
        for(int count = 0; count < n_tests; count++){
            std::cout << "start of iteration"<< std::endl;
            std::stringstream bnet_ss, data_set_ss;
            if(multiple_files){
                bnet_ss << bnet_str << dc.numprocs() << ".txt";
            }else{
                bnet_ss << bnet_str << ".txt";
            }
            data_set_ss << data_set_str << "-" << ((count+1)*100) << ".csv";

            for(int i = 0;  i < iterations; i++){
                elapsed = complete_data(dc, in_str+bnet_ss.str(), in_str+data_set_ss.str());
                std::cout << "Dataset: " << data_set_ss.str() << std::endl;
                std::cout << "complete_data learning done"<< std::endl;
                if(dc.procid() == 0){
                    output << bnet_str << "\t" << data_set_str << "\t" << ((count+1)*100) << "\t" << dc.numprocs() << "\t" << elapsed << std::endl;
                }
            }
        }
        output.close();
        return 0;
    }
}

/*
It mixes inference and learning
*/
int experiment1(graphlab::distributed_control& dc, std::string path, std::string net, std::string var_str, std::string data_str, std::string out_str){
    dc.barrier();
    std::cout << "creating bn"<< std::endl;
    BeliefProp i_eng;
    BeliefNetwork bn(dc, &i_eng);
    dc.barrier();
    double begin_inference, end_inference, end_learning;
    std::stringstream var_ss(var_str);
    RandomVariable var = RandomVariable::parseLine(var_ss);
    //std::cout << "bn created"<< std::endl;
    std::cout << path+net << " " << path+data_str;
    bn.load(path+net);
    bn.save();
    std::vector<std::string> triple(3 ,"");

    DataSet data(path+data_str);
    //std::cout << "DATASET read" << std::endl;
    std::cout << "---inference---" << std::endl;
    begin_inference = clock();
    //Factor grass_infer = bn.infer(var);
    end_inference = clock();
    std::cout << "---learning---" << std::endl;
    bn.learnParameters(data);
    end_learning = clock();
    std::cout << "learning done"<< std::endl;
    //bn.save();
    //std::cout << "saving done"<< std::endl;
    //std::cout << "P(wet) = " << bn.infer(grass) << std::endl;
    if(dc.procid() == 0){
        std::cout << "times" << (end_inference - begin_inference)*1000 / CLOCKS_PER_SEC << "ms\t" <<(end_learning - end_inference) *1000/ CLOCKS_PER_SEC << "ms";
        std::ofstream output;
        output.open(out_str.c_str(), std::ios::out | std::ios::app);
        output << net << "\t" << var.getName() << "\t" << dc.numprocs() << "\t" << (end_inference - begin_inference)*1000 / CLOCKS_PER_SEC << "\t" <<(end_learning - end_inference) *1000/ CLOCKS_PER_SEC << std::endl;
        output.close();
    }
    return 0;
}

int main(int argc, char** argv){
    graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;

    //double begin, end;

    /*
    std::cout << "test1" << std::endl;
    begin = clock();
    test1(dc); // Testing FactorGraph
    end = clock();
    std::cout << "Time Elapsed: " << (end-begin)/ CLOCKS_PER_SEC << "s." << std::endl;
    */
    //std::cout << "test2" << std::endl;
    //if you are the node 0
    //begin = clock();
    //test2(dc); // Testing BNet
    //end = clock();
    //std::cout << "Time Elapsed: " << (end-begin)/ CLOCKS_PER_SEC << "s." << std::endl;
    //}

    /*
    std::cout << "test3" << std::endl;
    //if you are the node 0
    begin = clock();
    test3(dc); // Testing BNet
    end = clock();
    std::cout << "Time Elapsed: " << (end-begin)/ CLOCKS_PER_SEC << "s." << std::endl;
    //}
    */

    /*for(int i = 0; i < 100; i++){
        experiment1(dc,
                        "/home/breno/Documents/git-repos/graphlab/apps/bnets/samples/",
                        "( rain |  false true )",
                        "test.gl",
                        "test.csv",
                        "/home/breno/pgm_test" );
    }*/

    /*
    for(int i = 0; i < 1; i++){
        experiment1(dc, "/home/breno/Documents/git-repos/graphlab/apps/bnets/samples/",
                        "BBNet-hailfinder.txt",
                        "P(CapInScen[LessThanAve,Average,MoreThanAve]|AMCINInScen,CapChange)",
                        "hailfinder_500_perc30_missing_data.txt",
                        "/home/breno/pgm_test" );
    }

    for(int i = 0; i < 1; i++){
        experiment1(dc,
                        "/home/breno/Documents/git-repos/graphlab/apps/bnets/samples/",
                        "BBNet-alarm.txt",
                        "P(CO[Low,Normal,High]|HR,StrokeVolume)",
                        "alarm_1000_perc30_missing_data.txt",
                        "/home/breno/pgm_test" );
    }
    */
    /*
    exp_learning::complete_data(dc, "bnet-sprinkler", "sprinkler-uniform-",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_net_fake_data/",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_net_fake_data/results/complete_data.txt", 10, 100);
    */
        exp_learning::complete_data(dc, "alarm", "fake_data_alarm/alarm",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_data/",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_data/results/complete_data.txt", 10, 100, false);

        exp_learning::complete_data(dc, "andes", "fake_data_andes/andes",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_data/",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_data/results/complete_data.txt", 10, 100, false);

        exp_learning::complete_data(dc, "barley", "fake_data_barley/barley",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_data/",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_data/results/complete_data.txt", 10, 100, false);

        exp_learning::complete_data(dc, "hepar2", "fake_data_hepar2/hepar2",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_data/",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_data/results/complete_data.txt", 10, 100, false);

        exp_learning::complete_data(dc, "munin", "fake_data_munin/munin",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_data/",
        "/home/breno/Documents/git-repos/graphlab/apps/bnets/exp_files/bayesian/fake_data/results/complete_data.txt", 10, 100, false);

    std::cout << "closing mpi" << std::endl;
    dc.cout() << "num of procs: " << dc.numprocs() << std::endl;
    graphlab::mpi_tools::finalize();
    std::cout << "mpi closed" << std::endl;
}
