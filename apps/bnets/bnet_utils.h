/*
*  \author: Breno Carvalho
*/

#ifndef BC_BAYESIAN_NETWORK_UTILS_HPP
#define BC_BAYESIAN_NETWORK_UTILS_HPP

#include <graphlab.hpp>

#include <../../../toolkits/graphical_models/factors/factor_graph.hpp>
#include <../../../toolkits/graphical_models/factors/bp_vertex_program.hpp>
#include <../../../toolkits/graphical_models/factors/bp_graph_data.h>

// Include the macro for each operation
#include <graphlab/macros_def.hpp>

//add functionalities to factor graph to allow evidence and change of the structure

struct clopts_vals {
  clopts_vals(double bound = 1E-4, double damping = 0.0,
      std::string exec_type="sync", int verbose = LOG_EMPH) :
      BOUND(bound), DAMPING(damping), EXEC_TYPE(exec_type), VERBOSE(verbose) { }

  double BOUND;
  double DAMPING;
  std::string EXEC_TYPE;
  int VERBOSE;
};

int setup_cli(graphlab::command_line_options& clopts, clopts_vals& opts,
    int argc, char** argv) {

  clopts.attach_option("bound", opts.BOUND,
                       "Residual termination bound");
  clopts.attach_option("damping", opts.DAMPING,
                       "The amount of message damping (higher = more damping)");
  clopts.attach_option("verbose", opts.VERBOSE,
                       "Verbosity of Printing: 0 (lots), 2 (default), 6 (no printing).");
//  clopts.attach_option("beliefs", &beliefs_filename,
//                       "The file to save the belief predictions");
  clopts.attach_option("engine", opts.EXEC_TYPE,
                       "The type of engine to use {async, sync}.");
  clopts.set_scheduler_type("fifo");

  bool success = clopts.parse(argc, argv);
  if(!success) {
//    std::cout << "Error parsing command line arguments!"
//              << std::endl;
    graphlab::mpi_tools::finalize();
    return EXIT_FAILURE;
  }

  //std::cout << "logging level: " << std::max(opts.VERBOSE, OUTPUTLEVEL) << std::endl;
  //global_logger().set_log_level(opts.VERBOSE);

  return EXIT_SUCCESS;
}

template<size_t MAX_DIM>
class InferenceEngine{
    //Prepare the methods for a general inference engine and consider moving that to another file
public:
    virtual void run(graphlab::distributed_control&dc ,
        graphlab::distributed_graph<belief_prop::vertex_data<MAX_DIM>, belief_prop::edge_data<MAX_DIM> >& graph,
        std::basic_string<char, std::char_traits<char>, std::allocator<char> > const& exec_type,
        graphlab::command_line_options const& clopts) =0; //virtual method with no implementation
};

template<size_t MAX_DIM>
class BeliefPropagationEngine : public InferenceEngine<MAX_DIM> {
    //Prepare the methods for a general inference engine and consider moving that to another file
public:
    void run(graphlab::distributed_control&dc ,
        graphlab::distributed_graph<belief_prop::vertex_data<MAX_DIM>, belief_prop::edge_data<MAX_DIM> >& graph,
        std::basic_string<char, std::char_traits<char>, std::allocator<char> > const& exec_type,
        graphlab::command_line_options const& clopts){
            // Create the engine -------------------------------------------------------->
            typedef graphlab::omni_engine<belief_prop::bp_vertex_program<MAX_DIM> > engine_type;
            engine_type engine(dc, graph, exec_type, clopts);

            engine.signal_all();
            engine.start();
            const float runtime = engine.elapsed_seconds();
            size_t update_count = engine.num_updates();

    }

};

template<size_t MAX_DIM> class Variable;
template<size_t MAX_DIM> class BeliefNetwork;


template<size_t MAX_DIM>
class BeliefNetwork {
public:
    typedef graphlab::dense_table<MAX_DIM>       dense_table_t;
    typedef graphlab::discrete_domain<MAX_DIM>   domain_t;
    typedef graphlab::discrete_variable          variable_t;
    typedef belief_prop::factor_graph<MAX_DIM>   factor_graph_t;
public:
    factor_graph_t fgraph;
    std::map<std::string, Variable<MAX_DIM> > variables;
    bool hasChanged;
    std::map<std::string, std::string> evidence;
public:
    BeliefNetwork(): hasChanged(false){};

    variable_t add_variable(std::vector<std::string>& values, std::vector<double>& cpd, std::string& variable_id){
        variable_t f_var = add_variable(values, variable_id);
        variables[variable_id].set_cpd(cpd);
        //hasChanged = true;
        return f_var;
    };

    variable_t add_variable(std::vector<std::string>& values, std::string& variable_id){
        Variable<MAX_DIM> var = Variable<MAX_DIM>(values, variable_id, fgraph);
        //std::cout << values.size() << " " << variable_id << std::endl;
        variable_t f_var = fgraph.add_variable(values.size(), variable_id);
        var.set_f_graph_node_num(f_var.id());
        //std::cout<<"varnum "<<f_var.id()<<std::endl;
        variables[variable_id] = var;
        hasChanged = true;
        return f_var;
    };

    void set_variable_parents(std::vector<std::string>& parent_ids, std::string& variable_id){
        std::vector< Variable<MAX_DIM> > parents(parent_ids.size());
        for(int i=0; parent_ids.size(); i++){
            parents[i] = variables[parent_ids[i]];
        }
        variables[variable_id].set_parents(parents);
        hasChanged = true;
    };

    void add_variable_parent(std::string& parent_id, std::string& variable_id){
        variables[variable_id].add_parent(variables[parent_id]);
        hasChanged = true;
    };
    void set_variable_cpd(std::vector<double>& cpd, std::string& id){
        variables[id].set_cpd(cpd);
    };
    //void set_variable_cpd(dense_table_t& cpd, std::string& id){};
    dense_table_t get_variable_cpd(variable_t& id){
        return *(variables[id]->get_cpd());
    };

    void set_evidence(std::map<std::string, std::string>& evidence){
        this->evidence(evidence);
        std::map<std::string, std::string >::const_iterator iterator = evidence.begin();
        for(; iterator != evidence.end(); iterator++){
            variables[(*iterator).first].set_evidence((*iterator).second);
        }
        hasChanged = true;
    };
    void add_evidence(std::string id, std::string value){
        evidence[id] = value;
        variables[id].set_evidence(value);
        hasChanged = true;
    };

    void clear_evidence(std::string& id){
        variables[id].clear_evidence();
        evidence.erase(id);
        hasChanged = true;
    }
    void clear_evidence(){
        std::map<std::string, std::string >::const_iterator iter = evidence.begin();
        for(; iter != evidence.end(); iter++){
            variables[(*iter).first].clear_evidence();
        }
        evidence.clear();
        hasChanged = true;
    };
    std::string get_evidence(std::string& variable_id) const{
        return evidence[variable_id];
    };
    std::map<std::string, std::string> get_evidence(){ return evidence;};

    dense_table_t infer_variable_cpd(InferenceEngine<MAX_DIM>& eng,
        std::string& variable_id,
        graphlab::distributed_control& dc,
        clopts_vals& clvals,
        const graphlab::command_line_options& clopts){
            typename belief_prop::graph_type<MAX_DIM>::type graph(dc, clopts);
            fgraph.make_bp_graph( graph, clvals.BOUND, clvals.DAMPING );
            eng.run(dc, graph, clvals.EXEC_TYPE, clopts);
            fgraph.pull_beliefs_for_variables( graph );

            return *(variables[variable_id].get_cpd());
        };
    dense_table_t infer_variables_cpd(InferenceEngine<MAX_DIM>& eng, std::vector<std::string>& variable_ids){
        dense_table_t d;
        hasChanged = false;
        return d;};

    //Use SFrames to get info to set your parameters (complete data)
    //Use SFrames to get info to set your parameters (incomplete data) <- Use an engine

    //Use SFrames to get info to set your structure(complete data)
    //Use SFrames to get info to set your structures (incomplete data) <- Use an engine

    //Save belief network to file
    //Load belief network from file

    friend std::ostream& operator<<(std::ostream& os, const BeliefNetwork<MAX_DIM>& bn){
        os<< "Belief Network:" << std::endl <<
            "has" << ((bn.hasChanged)? " ":" not ") << "changed" << std::endl;
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
    typedef  belief_prop::vertex_data<MAX_DIM> vertex_data_t;
    typename BeliefNetwork<MAX_DIM>::dense_table_t distribution;
    typename BeliefNetwork<MAX_DIM>::factor_graph_t *factor_graph;
    size_t f_graph_num;
    size_t f_graph_dist;
    std::vector<Variable*> parents;
    std::vector<typename BeliefNetwork<MAX_DIM>::variable_t> f_parents;
    std::vector<std::string> values;
    int evidence;
    std::string id;

public:
    Variable(std::vector<std::string>& values, std::string& variable_id,
            typename BeliefNetwork<MAX_DIM>::factor_graph_t &f_graph):
            values(values), id(variable_id), evidence(-1), factor_graph(&f_graph), f_graph_dist(0) {};

    Variable():
            id(""), evidence(-1), factor_graph(NULL), f_graph_dist(0) {};

    std::string get_id() const{ return id;}

    void set_parents(std::vector<Variable>& parents){
        this->parents.clear();
        f_parents.clear();
        for(int i = 0; i < parents.size(); i++){
            this->parents.push_back(&parents[i]);
            f_parents.push_back(parents[i]->get_f_graph_node());
        }
    };

    void add_parent(Variable& parent){
        parents.push_back(&parent);
        f_parents.push_back(parent.get_f_graph_node());
    };

    void set_f_graph_node_num(size_t num){
        f_graph_num = num;
    }

    typename BeliefNetwork<MAX_DIM>::variable_t get_f_graph_node() const{
        //std::cout << "id "<< id <<" f_graph_num" << f_graph_num<< std::endl;
        return factor_graph->get_variable(f_graph_num);}

    void set_cpd(std::vector<double>& cpd){
        f_parents.push_back(factor_graph->get_variable(f_graph_num));
        typename BeliefNetwork<MAX_DIM>::dense_table_t tmp(f_parents, cpd);
        distribution = tmp;
        f_graph_dist = factor_graph->add_factor(tmp, id+"_factor");
        f_parents.pop_back();
        //std::cout << "Setting cpd for "<< id << std::endl;
    };

    typename BeliefNetwork<MAX_DIM>::dense_table_t* get_cpd() const{
        std::vector<vertex_data_t>& facts = factor_graph->factors();
        vertex_data_t& vertex = facts[f_graph_dist];
        graphlab::table_factor<MAX_DIM> &belief = vertex.potential;
        typename BeliefNetwork<MAX_DIM>::dense_table_t *table = dynamic_cast<typename BeliefNetwork<MAX_DIM>::dense_table_t*>(belief.table());
        assert(table != NULL);

        return table;
    };

    void clear_evidence(){ evidence = -1;}

    void set_evidence(std::string value){
        int evid;
        typename BeliefNetwork<MAX_DIM>::dense_table_t* cpd = get_cpd();
        typename BeliefNetwork<MAX_DIM>::domain_t::const_iterator end = cpd->domain().end();
        typename BeliefNetwork<MAX_DIM>::domain_t::const_iterator asg = cpd->domain().begin();

        for(evid=0; evid < values.size(); evid++){
            if(values[evid].compare(value) == 0){ // if we have a match of value in the values list then:
                evidence = evid;
                break;
            }
        }
        for(; asg != end; ++asg) {
            if((*asg).asg(get_f_graph_node()) == evidence){
                //std::cout<< "Evidence met:"<< (*asg).asg(get_f_graph_node()) <<" "<< evidence;
                cpd->set_logP( *asg, 0);//distribution.logP(*asg) );
            }else{
                cpd->set_logP( *asg, -100000);//-std::numeric_limits<size_t>::infinity() );
            }
        }
    }

    std::string get_evidence() const{
        if(evidence >= 0){
            return values[evidence];
        }
        return "";
    }

    bool has_evidence() const{ return evidence > 0; }

    friend std::ostream& operator<<(std::ostream& os, const Variable<MAX_DIM>& var){
        os << "Variable{id = '" << var.id << "'";

        os << ", values = {'" << var.values[0] << "'";
        for(uint i = 1; i < var.values.size(); i++){
            os << ", '" << var.values[1] << "'";
        }
        os << "}";

        if(var.parents.size() > 0){
            os << ", parents = {'"<< var.parents[0]->id<<"'";
            for(int i = 1; i < var.parents.size(); i++){
                os << ", '" << var.parents[i]->id << "'";
            }
            os << "}";
        }
        if(var.evidence >= 0){
            os << ", evidence = '" << var.get_evidence() << "'";
        }

        os << ", cpd = <" << *(var.get_cpd()) << ">";
        os << "}";
        return os;
    };
};

#endif //BC_BAYESIAN_NETWORK_UTILS_HPP
