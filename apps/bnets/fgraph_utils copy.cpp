/*
*  \author: Breno Carvalho
*/

#ifndef FACTOR_GRAPH_NETWORK_UTILS_HPP
#define FACTOR_GRAPH_NETWORK_UTILS_HPP

#include <graphlab.hpp>

/*#include <../../../toolkits/graphical_models/factors/factor_graph.hpp>
#include <../../../toolkits/graphical_models/factors/bp_vertex_program.hpp>
#include <../../../toolkits/graphical_models/factors/bp_graph_data.h>*/
#include <../../../apps/bnets/csv_reader.h>

// Include the macro for each operation
#include <graphlab/macros_def.hpp>

class VerticeObject;
class RandomVariable;

typedef std::vector<RandomVariable> scope_type;
typedef std::vector<int> attribution_type;
typedef std::vector<double> cpd_type;

class VerticeObject{
public:
    virtual scope_type const getScope() =0;
};

class RandomVariable : public VerticeObject{
private:
    std::string name;
    std::vector<std::string> values;

public:
    RandomVariable(std::string const& var_name, std::vector<std::string> const& var_values){
        name   = var_name;
        values = var_values;
    }

    RandomVariable(){};

    std::string const& getName() const {return name;}

    std::vector<std::string> const& getValues() const {return values;}

    friend inline bool operator==(const RandomVariable& v1, const RandomVariable& v2){
        return v1.name == v2.name;
    }

    void save(graphlab::oarchive& oarc) const {
        oarc << name << values;
    }
    void load(graphlab::iarchive& iarc) {
        iarc >> name >> values;
    }

    friend std::ostream& operator<<(std::ostream& os, const RandomVariable& var){
        os << "(\"" << var.name << "\" |";
        std::vector<std::string>::const_iterator iter = var.values.begin();
        if( iter != var.values.end()){
            os << *iter; iter ++;
            for(; iter != var.values.end(); iter++){
                os << ", " << *iter;
            }
        }
        os << "|)";
        return os;
    };

    //TODO check if it is safe
    scope_type const getScope(){return scope_type(1, *this);};
};

class Factor : public VerticeObject{
private:
    scope_type scope;
    attribution_type scope_sizes;
    cpd_type cpd;

    static scope_type scopeUnion(scope_type const& s1,scope_type const& s2){
        scope_type out;
        scope_type const* bigger; scope_type const* smaller;
        if(s1.size() >= s2.size()){
            bigger = &s1;
            smaller = &s2;
        }else{
            bigger = &s2;
            smaller = &s1;
        }

        for(scope_type::const_iterator iter = bigger->begin();
        iter != bigger->end(); iter++){
            out.push_back(*iter);
        }
        for(scope_type::const_iterator iter = smaller->begin();
        iter != smaller->end(); iter++){
            if(!hasVariable(*bigger, *iter)){
                out.push_back(*iter);//out->push_back(tmp);
            }
        }
        return out;
    }
    static scope_type scopeMinus(scope_type const& s1, scope_type const& s2){
        scope_type out;
        for(scope_type::const_iterator iter = s1.begin(); iter != s1.end(); iter++){
            if(!hasVariable(s2, *iter)){
                out.push_back(*iter);//out->push_back(tmp);
            }
        }
        return out;
    }

    static attribution_type getAttribMap(scope_type const& from_scope, scope_type const& to_scope){
        attribution_type out(to_scope.size(), 0);
        int count, i = 0;
        for(scope_type::const_iterator to_iter = to_scope.begin(); to_iter != to_scope.end(); to_iter++){
            count = 0;
            for(scope_type::const_iterator from_iter = from_scope.begin(); from_iter != from_scope.end(); from_iter++){
                if(*to_iter == *from_iter){
                    out[i] = count;
                    break;
                }
                count++;
            }
            i++;
        }
        return out;
    }
    static attribution_type applyAttribMap(attribution_type const& attrib, attribution_type const& attrib_map){
        attribution_type out(attrib_map.size(),0);
        int count=0;
        //std::cout << "||";
        for(attribution_type::const_iterator i = attrib_map.begin(); i != attrib_map.end(); i++, count++){
            out[count] = attrib[*i];
            //std::cout << ":"<<out[count];
        }
        //std::cout << "||";
        return out;
    }

    static double& getValue(attribution_type const& scope_sizes, cpd_type& cpd, attribution_type const& assignment){
        int count = assignment.front();
        //TODO assert assignment.size() == scope.size();
        attribution_type::const_iterator i_assig = assignment.begin();
        attribution_type::const_iterator i_scope = scope_sizes.begin();
        i_assig++;
        for(;i_assig != assignment.end() && i_scope != scope_sizes.end();
        i_assig++, i_scope++){
            count += (*i_assig)*(*i_scope);
        }
        return cpd[count];
    }

public:
    scope_type const getScope(){return scope;};

    static bool hasVariable(scope_type const& scope, RandomVariable const& var){
        for(scope_type::const_iterator iter = scope.begin();
            iter != scope.end(); iter++){
            if(var == *iter){
                return true;
            }
        }
        return false;
    }
    bool hasVariable(RandomVariable const& var) const{
        return Factor::hasVariable(scope, var);
    }

    Factor(){}

    Factor(scope_type const& scope_list, cpd_type const& values){
        cpd   = cpd_type(values);
        scope = scope_type(scope_list);
        scope_sizes = attribution_type(scope_list.size(), 1);
        //int cpd_size = 1;
        for(scope_type::const_iterator i = scope.begin();
        i != scope.end(); i++){
            scope_sizes.push_back(i->getValues().size());
            //cpd_size *= i->getValues().size()
        }
    };

    void normalize(int variable_index){
        if(scope.size() == 1){
            double z = 0;
            for(cpd_type::const_iterator iter = cpd.begin(); iter != cpd.end(); iter++){
                z += *iter;
            }
            if(z != 0){
                for(cpd_type::iterator iter = cpd.begin(); iter != cpd.end(); iter++){
                    *iter /= z;
                }
            }
        }else{
            //TODO
        }
    }

    void save(graphlab::oarchive& oarc) const {
        oarc << scope << scope_sizes << cpd;
    }
    void load(graphlab::iarchive& iarc) {
        iarc >> scope >> scope_sizes >> cpd;
    }

    double& operator[](attribution_type const& idx) {
        return getValue(scope_sizes, cpd, idx);
    };
    double const& operator[](attribution_type const& idx) const {
        return const_cast<Factor&>(*this)[idx];
    };

    Factor& operator*=(Factor const& rhs){
        scope_type old_scope = this->scope;
        attribution_type old_scope_sizes = this->scope_sizes;
        cpd_type old_cpd = this->cpd;

        this->scope = scopeUnion(this->scope, rhs.scope);
        this->scope_sizes = attribution_type(this->scope.size(), 1);
        //~scope; ~scope_sizes;
        int cpd_size = 1; unsigned int count = 0;
        for(scope_type::const_iterator i = this->scope.begin();
        i != this->scope.end(); i++, count++){
            this->scope_sizes[count] = i->getValues().size();
            cpd_size *= i->getValues().size();
        }
        this->cpd = cpd_type(cpd_size, 0.0);
        attribution_type attrib  = attribution_type(this->scope.size(), 0);
        attribution_type attrib_map_left  = getAttribMap(this->scope, old_scope);
        attribution_type attrib_map_right = getAttribMap(this->scope, rhs.scope);

        count = 0;
        while(count < this->cpd.size()){
            for(int val = 0; val < this->scope_sizes[0]; val++){
                attrib[0] = val;
                //multiplication
                double mult = getValue(old_scope_sizes, old_cpd, applyAttribMap(attrib, attrib_map_left));
                mult *= rhs[applyAttribMap(attrib, attrib_map_right)];
                this->cpd[count] = mult;
                count++;
            }

            for(unsigned int var_pos = 0; var_pos < this->scope_sizes.size(); var_pos++){
                attrib[var_pos] += 1;
                if(attrib[var_pos] == this->scope_sizes[var_pos]){
                    attrib[var_pos] = 0;
                }else{
                    break;
                }
            }
        }
        return *this; // return the result by reference
    }

    inline friend Factor& operator*(Factor lhs, Factor const& rhs){
        return lhs*=rhs;
    }

    //created only for compatibility with graphlab::aggregator_type
    Factor& operator+=(Factor const& rhs){
        return this->operator*=(rhs);
    }

    inline Factor marginalize(RandomVariable const& var) const{
        return marginalize(scope_type(1,var));
    }

    inline Factor marginalizeAllBut(RandomVariable const& var) const{
        return marginalizeAllBut(scope_type(1, var));
    }

    Factor marginalize(scope_type const& vars) const{
        if (this->scope.size() == 0){
            return *this;
        }
        scope_type new_scope = scopeMinus(this->scope, vars);
        int cpd_size = 1;
        for(scope_type::const_iterator i = new_scope.begin();
        i != new_scope.end(); i++){
            cpd_size *= i->getValues().size();
        }
        cpd_type new_cpd = cpd_type(cpd_size, 0.0);
        Factor new_factor(new_scope, new_cpd);

        attribution_type attrib = attribution_type(this->scope.size(), 0);
        attribution_type new_f_attrib_map = getAttribMap(this->scope, new_scope);

        unsigned int count = 0;
        while(count < this->cpd.size()){
            for(int val = 0; val < this->scope_sizes[0]; val++){
                attrib[0] = val;
                new_factor[applyAttribMap(attrib, new_f_attrib_map)] += (*this)[attrib];
                count++;
            }

            for(unsigned int var_pos = 0; var_pos < this->scope_sizes.size(); var_pos++){
                attrib[var_pos] += 1;
                if(attrib[var_pos] == this->scope_sizes[var_pos]){
                    attrib[var_pos] = 0;
                }else{
                    break;
                }
            }
        }
        return new_factor;

    }

    inline Factor marginalizeAllBut(scope_type const& vars) const{
        return marginalize(scopeMinus(this->scope, vars));
    }

    inline Factor getUnitary() const{
        cpd_type tmp_cpd(cpd.size(), 1.);
        return Factor(this->scope, tmp_cpd);
    }

    friend std::ostream& operator<<(std::ostream& os, const Factor& fact){
        scope_type::const_iterator var = fact.scope.begin();
        os << "[<";
        if(var != fact.scope.end()){
            os << "'" << var->getName() << "'"; var++;
            for(;var != fact.scope.end(); var++){
                os << ", '" << var->getName()<<"'";
            }
        }
        os << "> Values: ";
        if(fact.cpd.size() > 0){
            os << fact.cpd[0];
            for(unsigned int i=1; i < fact.cpd.size(); i++){
                os<< ", " << fact.cpd[i];
            };
        }
        os << "]";
        return os;
    };
};

class Vertice {
private:
    Factor message;
    graphlab::vertex_id_type var_id;
    static size_t id;

public:
    enum {FACTOR, VARIABLE} type;

    explicit Vertice(int id): var_id(id){}
    explicit Vertice(): var_id(id++){}

    VerticeObject& getInfo();

    graphlab::vertex_id_type const& vid() const{ return var_id;}

    inline Factor const& getMessage() const{ return message;}

    //TODO check how I can change it to setMessage(Factor& msg)
    inline void setMessage(Factor msg){ message = msg;}

    inline int getType() const{ return type;}

    void save(graphlab::oarchive& oarc) const {
        oarc << message << var_id << id << type;
    }

    void load(graphlab::iarchive& iarc) {
        iarc >> message >> var_id >> id >> type;
    }
    friend inline bool operator==(const Vertice& v1, const Vertice& v2){
        return v1.vid() == v2.vid();
    }
};

size_t Vertice::id = 0;

class VariableVertice : public Vertice{
private:
    RandomVariable rVariable;

public:
    int n[1024];
    VariableVertice(RandomVariable& var):Vertice(){
        type      = VARIABLE;
        rVariable = var;
        cpd_type tmp_cpd(var.getValues().size(), 1.0);
        Factor tmp_fact(scope_type(1, var), tmp_cpd);
        this->setMessage(tmp_fact);

        n[0] = 1; n[1] = 0;
        std::cout << "vertice("<< var<<")" << n;
    }

    VariableVertice():Vertice(){
        type = VARIABLE;
        std::string tmp = "dog";
        n[0] = 2; n[1] = 0;
        std::cout << "vertice()" << n;
    }

    inline VerticeObject& getInfo(){
        return rVariable;
    }

    void save(graphlab::oarchive& oarc) const {
        ((Vertice &) *this).save(oarc);
        oarc << rVariable;
        std::cout << "saving " << rVariable << n;
    }

    void load(graphlab::iarchive& iarc) {
        ((Vertice &) *this).load(iarc);
        iarc >> rVariable;
        n[0] = 3; n[1] = 0;
        std::cout << "loading " << rVariable << n;
    }

    friend std::ostream& operator<<(std::ostream& os, const VariableVertice& vertice){
        os << "VARIABLE" << vertice.rVariable;
        return os;
    };
};

class FactorVertice : public Vertice{
private:
    Factor factor;
public:
    FactorVertice(Factor& fact):Vertice(){
        type = FACTOR;
        factor = fact;
        this->setMessage(fact);
    }

    FactorVertice():Vertice(){
        type = FACTOR;
    }

    inline VerticeObject& getInfo(){
        return factor;
    }

    void save(graphlab::oarchive& oarc) const {
        ((Vertice &) *this).save(oarc);
        oarc << factor;
        std::cout << "saving " << factor;
    }

    void load(graphlab::iarchive& iarc) {
        ((Vertice &) *this).load(iarc);
        iarc >> factor;
        std::cout << "loading " << factor;
    }

    friend std::ostream& operator<<(std::ostream& os, const FactorVertice& vertice){
        os << "FACTOR";
        os << vertice.factor;
        return os;
    };
};

typedef graphlab::distributed_graph<Vertice, graphlab::empty> GraphLab_type;
typedef std::pair<FactorVertice, VariableVertice> Edge_type;

class InferenceEngine:  public graphlab::ivertex_program<GraphLab_type, Factor>,
                        public graphlab::IS_POD_TYPE{
    public:
    virtual void runInference(graphlab::distributed_control& dist_ctrl, GraphLab_type & dgraph) =0;
};


class BeliefProp : public InferenceEngine{
    public:
        BeliefProp(){};
        void runInference(graphlab::distributed_control& dist_ctrl, GraphLab_type& dgraph){
            graphlab::omni_engine<BeliefProp> engine(dist_ctrl, dgraph, "sync");
            engine.signal_all();
            engine.start();
        }

        edge_dir_type gather_edges(icontext_type& context,
                                   const vertex_type& vertex) const {
            return graphlab::ALL_EDGES;
        }

        Factor gather(icontext_type& context, const vertex_type& vertex,
        edge_type& edge) const {
            if(vertex.data().getType() == Vertice::FACTOR){
                std::cout << "gather from factor " << edge.target().data().getMessage() <<std::endl;

                return edge.target().data().getMessage();
            }else{
                //if(vertex.data().getType() == Vertice::VARIABLE){
                std::cout << "gather from variable" << edge.source().data().getMessage() << std::endl;
                return edge.source().data().getMessage();
                //}
            }
        }

        // Use the total rank of adjacent pages to update this page
        void apply(icontext_type& context, vertex_type& vertex,
        const gather_type& total) {
            std::cout << "apply" << std::endl;
            if(vertex.data().getType() == Vertice::FACTOR){
                std::cout << "factor "<< std::endl;;
                FactorVertice& fact_vert = (FactorVertice&) (vertex.data());
                std::cout << "data "<< std::endl;
                Factor& fact = (Factor&) fact_vert.getInfo();
                std::cout << fact << total << std::endl; //<< total.marginalizeAllBut(fact.getScope()) << std::endl;
                vertex.data().setMessage(total);//.marginalizeAllBut(fact.getScope()));
                std::cout << "had message updated" << std::endl;
            }else{
                if(vertex.data().getType() == Vertice::VARIABLE){
                    std::cout << "variable " << std::endl;
                    VariableVertice& var_vert = (VariableVertice&) (vertex.data());
                    std::cout << "data "<< std::endl;
                    std::cout << var_vert.getMessage() << std::endl;
                    RandomVariable& var = (RandomVariable&) var_vert.getInfo();
                    std::cout << "-"<< var_vert.n[0] << std::string(((var_vert.n[0] == 0)?"!":"?")) <<"-";// << total << std::endl;
                    //vertex.data().setMessage(total);//.marginalizeAllBut(var));
                    std::cout << "had message updated" << std::endl;
                }
            }

            std::cout << " apply done" << std::endl;
        }

        // The scatter edges depend on whether the pagerank has converged
        edge_dir_type scatter_edges(icontext_type& context,
        const vertex_type& vertex) const {
        //    if (perform_scatter) return graphlab::OUT_EDGES;
            //else
            return graphlab::NO_EDGES;
        }
    };

//TODO
class Graph{
private:
    GraphLab_type dgraph;
    std::list< Edge_type > edges;
    std::list< FactorVertice > factors;
    std::list< VariableVertice > variables;
    graphlab::distributed_control & dist_ctrl;
    bool dgraph_done;

    void addEdge(Edge_type& edge){
        edges.push_back(edge);
        dgraph_done = false;
    }

public:
    explicit Graph(graphlab::distributed_control& dc): dist_ctrl(dc), dgraph_done(false), dgraph(GraphLab_type(dc)){};

    void addVertice(Vertice& vert){
        dgraph_done = false;
        if(vert.getType() == Vertice::FACTOR){
            factors.push_back((FactorVertice &) vert);
            return;
        }
        if(vert.getType() == Vertice::VARIABLE){
            variables.push_back((VariableVertice &) vert);
            return;
        }
    }

    inline void addEdge(VariableVertice & var, FactorVertice & fact){
        Edge_type edge(fact, var);
        addEdge(edge);
    }

    inline void addEdge(FactorVertice & fact, VariableVertice & var){
        Edge_type edge(fact, var);
        addEdge(edge);
    }

    void removeEdge(VariableVertice const& var, FactorVertice const& fact){
        for(std::list<Edge_type>::iterator edge = edges.begin(); edge != edges.end(); edge++){
            if((edge->first == fact) && (edge->second == var)){
                edges.erase(edge);
                dgraph_done = false;
            }
        }
    }

    void removeEdge(graphlab::vertex_id_type const& var_id, graphlab::vertex_id_type const& fact_id){
        for(std::list<Edge_type>::iterator edge = edges.begin(); edge != edges.end(); edge++){
            if((edge->first.vid() == fact_id) && (edge->second.vid() == var_id)){
                edges.erase(edge);
                dgraph_done = false;
            }
        }
    }

    //TODO make it parallel
    void startDGraph(){
        dgraph.clear();
        for(std::list< FactorVertice >::const_iterator vertice = factors.begin(); vertice != factors.end(); vertice++){
            dgraph.add_vertex(vertice->vid(), *vertice);
        }
        for(std::list< VariableVertice >::const_iterator vertice = variables.begin(); vertice != variables.end(); vertice++){
            dgraph.add_vertex(vertice->vid(), *vertice);
        }

        for(std::list< Edge_type >::const_iterator edge = edges.begin(); edge != edges.end(); edge++){
            dgraph.add_edge(edge->first.vid(), edge->second.vid());
        }

        dgraph.finalize();
        dgraph_done = true;
    }

    bool isDGraphUp(){ return dgraph.is_finalized() && dgraph_done;}

    //TODO
    cpd_type getPrior(VariableVertice const& var);

    void runInference(InferenceEngine& eng){
        if(!isDGraphUp()){
            startDGraph();
            eng.runInference(dist_ctrl, dgraph);
        }
    }

    //void setEvidence(VariableVertice& var, int value);
    //void setEvidence(RandomVariable& var, int value);

    friend std::ostream& operator<<(std::ostream& os, const Graph& vertice){
        return os << "Graph";
    }
};

int main(int argc, char** argv){
    graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;

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

    VariableVertice v1(rain);
    FactorVertice a(f1*f2);
    std::cout << a     << std::endl \
              << v1    << std::endl \
              << f2    << std::endl \
              << f1*f2 << std::endl << (f1*f2).marginalizeAllBut(rain) << std::endl \
              << "V1: "<< (RandomVariable&) v1.getInfo() << std::endl \
              << "a: " << a.getMessage() <<std::endl;

    Graph g(dc);
    g.addVertice(v1);
    g.addVertice(a);

    g.addEdge(a, v1);

    BeliefProp bp;
    g.runInference(bp);

    //dc::cout << g;

    graphlab::mpi_tools::finalize();
}

#endif //FACTOR_GRAPH_NETWORK_UTILS_HPP
