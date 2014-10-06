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

#include <opencv2/opencv.hpp>

#include <graphlab.hpp>

#include <../../../toolkits/graphical_models/factors/factor_graph.hpp>
#include <../../../toolkits/graphical_models/factors/bp_vertex_program.hpp>

// Include the macro for each operation
#include <graphlab/macros_def.hpp>

typedef graphlab::dense_table<5>             dense_table_t;
typedef graphlab::discrete_domain<5>         domain_t;
typedef graphlab::discrete_variable          variable_t;


struct clopts_vals { 
  clopts_vals(double bound = 1E-4, double damping = 0.3, std::string exec_t="sync") : 
      BOUND(bound), DAMPING(damping), exec_type(exec_t) { }

  double BOUND;
  double DAMPING;
  std::string exec_type;
};


int setup_cli(graphlab::command_line_options& clopts, clopts_vals& clvals, 
    int argc, char** argv);
template<size_t MAX_DIM>
void run_engine(graphlab::distributed_control& dc, 
    typename belief_prop::graph_type<MAX_DIM>::type& graph, 
    const std::string& exec_type, 
    const graphlab::command_line_options& clopts);

// MAIN
// ============================================================================>
int main(int argc, char** argv) {
  std::cout << "This program solves the sum task."
            << std::endl;
  
  graphlab::mpi_tools::init(argc, argv);
  ///! Create a distributed control object (must come after mpi_tools::init())
  graphlab::distributed_control dc; 

  graphlab::command_line_options clopts("Run Loopy BP on a Belief Network");
  clopts_vals clvals;
  ///! Create a distributed graph object 
  belief_prop::graph_type<5>::type graph(dc, clopts);


  
  // Create the factor graph ------------------------------------------>
  std::cout << "Loading Factor Graph" << std::endl;
  belief_prop::factor_graph<5> fgraph; //5 = MAX_DIM

  std::vector<unsigned char> values(5);

  // create variables and prior factor nodes
  std::vector< std::vector< variable_t > > var_ids;
  var_ids.resize( 4, std::vector<variable_t>( 4, variable_t() ) );

  var_ids[0][0] = fgraph.add_variable(2, "c");
  var_ids[0][1] = fgraph.add_variable(2, "r");
  var_ids[1][0] = fgraph.add_variable(2, "s");
  var_ids[1][1] = fgraph.add_variable(2, "w");


  // add prior factor nodes
  double same_prob = 0.5f;
  double diff_prob = (1.0f - same_prob)/(5 - 1);
  for (unsigned i=0; i<2; ++i) {
    for (unsigned j=0; j<2; ++j) {
      variable_t var_id = var_ids[ j ][ i ];
      // NOTE the prior is always dense
      dense_table_t& prior = fgraph.prior_for_variable(var_id);

      domain_t::const_iterator end = prior.domain().end();
      for(domain_t::const_iterator asg = prior.domain().begin(); 
          asg != end; ++asg) {
        assert(asg->linear_index() < 5);
        double p;
        unsigned idx_og = 5;
        if(asg->linear_index() == idx_og)
          p = same_prob;
        else
          p = diff_prob;
        // Values are stored in log form
        prior.set_logP( *asg, log(p) );
      }
    }
  }


  // create factors and connect the variables
  float neighbor_same_prob = 0.5f;
  float neighbor_diff_prob = (1.0f - neighbor_same_prob)/(5 - 1);
  const float cost_scale = 10.0f;
  for (unsigned i=0; i<4; ++i) {
    for (unsigned j=0; j<4; ++j) {
      if (j != 0) {
        // Create the 2-way factor
        std::vector<variable_t> args;
        // connect vertical neighbors
        args.push_back(var_ids[ j-1 ][ i ]);
        args.push_back(var_ids[ j ][ i ]);

        // Construct the arguments (which will remap the domain)
        // std::cout << "domain: " << domain << std::endl;
        // Build the factor
        domain_t domain(args);
        dense_table_t factor(domain);

        // Set the weights
        domain_t::const_iterator end = factor.domain().end();
        for(domain_t::const_iterator asg = factor.domain().begin(); 
            asg != end; ++asg) {
          assert(asg->linear_index() < 5*5);
          double err;
          if(asg->asg(var_ids[j-1][i]) == asg->asg(var_ids[j][i]))
            err = neighbor_same_prob;
          else
            err = neighbor_diff_prob;

          // Values are stored in log form
          factor.set_logP( *asg, log(err*cost_scale) );
        }
        // Save the factor to the factor graph
        fgraph.add_factor(factor);
      }
      if (i != 0) {
        // Create the 2-way factor
        std::vector<variable_t> args;
        // connect horizontal neighbors
        args.push_back(var_ids[ j ][ i-1 ]);
        args.push_back(var_ids[ j ][ i ]);

        // Construct the arguments (which will remap the domain)
        // std::cout << "domain: " << domain << std::endl;
        // Build the factor
        domain_t domain(args);
        dense_table_t factor(domain);

        // Set the weights
        domain_t::const_iterator end = factor.domain().end();
        for(domain_t::const_iterator asg = factor.domain().begin(); 
            asg != end; ++asg) {
          assert(asg->linear_index() < 5*5);
          double err;
          if(asg->asg(var_ids[j][i-1]) == asg->asg(var_ids[j][i]))
            err = neighbor_same_prob;
          else
            err = neighbor_diff_prob;

          // Values are stored in log form
          factor.set_logP( *asg, log(err*cost_scale) );
        }
        // Save the factor to the factor graph
        fgraph.add_factor(factor);
      }
    }
  }
  assert(fgraph.num_factors() == 5*5 + (5-1)*5 + 5*(5-1));

  const size_t num_variables = fgraph.num_variables();
  const size_t num_factors = fgraph.num_factors();
  std::cout << "num_variables: " << num_variables << " "
            << "num_factors: " << num_factors << std::endl;
  std::cout << "Finished!" << std::endl;


  // Build the BP graph from the factor graph---------------------------------->
  std::cout << "Building BP graph from the factor graph" << std::endl;
  fgraph.make_bp_graph( graph, clvals.BOUND, clvals.DAMPING ); 
  run_engine<5>(dc, graph, clvals.exec_type, clopts);
  fgraph.pull_beliefs_for_variables( graph );


  // Saving the output -------------------------------------------------------->
  // NOTE: this can be done better. see loopybp_denoise.cpp
  std::cout << "Saving the predicted image" << std::endl;
  std::cout << "Collect the noisy image. " << std::endl;
//  merge_reduce_type pred_image = 
//    graph.map_reduce_vertices<merge_reduce_type>(pred_map_function);
  std::cout << "saving the pred image." << std::endl;
  if(dc.procid() == 0) {
    // Fill in output image------------------------------------------------------>
    cv::Mat_< uchar > output( 5, 5 );
    for (unsigned i=0; i<5; ++i) {
      for (unsigned j=0; j<5; ++j) {
        size_t ind = fgraph.belief_for_variable(var_ids[j][i]).max_index();
        output(j,i) = values[ind]; 
      }
    }
    cv::imwrite("denoised.png",output);

    cv::Mat_< uchar > gm = cv::imread("denoised_gm.png", 0); // force to grayscale with '0'
    cv::Scalar err = cv::sum(cv::abs(gm - output));
    ASSERT_LT(err(0), 5*5*1e-3);
  }

  std::cout << "All tests passed" << std::endl;
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main


// UTILS
// ============================================================================>
int setup_cli(graphlab::command_line_options& clopts, clopts_vals& clvals,
    int argc, char** argv) {

  clopts.attach_option("bound", clvals.BOUND,
                       "Residual termination bound");
  clopts.attach_option("damping", clvals.DAMPING,
                       "The amount of message damping (higher = more damping)");
//  clopts.attach_option("beliefs", &beliefs_filename,
//                       "The file to save the belief predictions"); 
  clopts.attach_option("engine", clvals.exec_type,
                       "The type of engine to use {async, sync}.");
  clopts.set_scheduler_type("fifo");


  bool success = clopts.parse(argc, argv);
  if(!success) {    
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    graphlab::mpi_tools::finalize();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template<size_t MAX_DIM>
void run_engine(graphlab::distributed_control& dc, 
    typename belief_prop::graph_type<MAX_DIM>::type& graph, 
    const std::string& exec_type, 
    const graphlab::command_line_options& clopts) 
{
  size_t num_vertices = graph.num_vertices();
  size_t num_edges = graph.num_edges();
  std::cout << "Loaded: " << num_vertices << " vertices "
            << "and " << num_edges << " edges." << std::endl;
  std::cout << "Finished!" << std::endl;

  // Create the engine -------------------------------------------------------->
  std::cout << "Creating the engine. " << std::endl;
  typedef graphlab::omni_engine<belief_prop::bp_vertex_program<5> > engine_type;
  engine_type engine(dc, graph, exec_type, clopts);

  std::cout << "Scheduling all vertices" << std::endl;
  engine.signal_all();
  std::cout << "Starting the engine" << std::endl;
  engine.start();
  const float runtime = engine.elapsed_seconds();
  size_t update_count = engine.num_updates();
  std::cout << "Finished Running engine in " << runtime 
            << " seconds with " << clopts.get_ncpus() << " cpus." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;
}


