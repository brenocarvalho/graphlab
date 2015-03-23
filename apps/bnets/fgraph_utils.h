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
#include <math.h>

//used in experiments
#include <ctime>
#include <iostream>
#include <fstream>

// Include the macro for each operation
//#include <graphlab/macros_def.hpp>

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

        size_t positionOfValue(std::string const& value) const{
            for(int out = 0; out < values.size(); out++){
                if(values[out] == value){
                    return out;
                }
            }
            return -1;
        }

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
            os << "( " << var.name << " | ";
            std::vector<std::string>::const_iterator iter = var.values.begin();
            if( iter != var.values.end()){
                //os << *iter; iter ++;
                for(; iter != var.values.end(); iter++){
                    os << " " << *iter;
                }
            }
            os << ")";
            return os;
        };

        /*friend std::istream& operator>>(std::istream& input, RandomVariable& var){
            std::string tmp;
            input >> tmp;
            if(tmp != "("){
                //std::cerr << "RandomVariable parser: Invalid sintax" << std::endl;
                input.setstate(std::ios::failbit);
            }

            input >> var.name;
            input >> tmp;
            if(tmp != "|"){
                std::cerr << "RandomVariable parser: Invalid sintax"<< std::endl;
            }
            std::list<std::string> v_list;

            while(1){
                input >> tmp;
                if(tmp == ")" || input.fail()){
                    break;
                }
                v_list.push_back(tmp);
            }
            var.values = std::vector<std::string>(v_list.size());
            std::copy(v_list.begin(), v_list.end(), var.values.begin());
            return input;
        }*/

        static RandomVariable parseLine(std::stringstream& ss){
            std::string name, tmp;
            std::list<std::string> v_list;
            if( ss >> tmp && tmp == "("){
                ss >> name;
                if(ss >> tmp){
                    if (tmp == "|"){
                        while(ss >> tmp && tmp != ")"){
                            v_list.push_back(tmp);
                        }
                            v_list.push_back("None");
                    }
                }
            }
            std::vector<std::string> values(v_list.size());
            std::copy(v_list.begin(), v_list.end(), values.begin());
            return RandomVariable(name, values);
        }

        //TODO check if it is safe
        scope_type const getScope(){return scope_type(1, *this);};
};


/**
* Implements a discrete factor.
* It is used in the other graphical models to represent CPD's (Conditional Probability Distributions).
* Factor has a vector of Random Variables, called scope, and a vector of values, the cpd.
**/
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
        static scope_type scopeIntersection(scope_type const& s1, scope_type const& s2){
            scope_type out;
            for(scope_type::const_iterator iter = s1.begin(); iter != s1.end(); iter++){
                if(hasVariable(s2, *iter)){
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
            for(attribution_type::const_iterator i = attrib_map.begin(); i != attrib_map.end(); i++, count++){
                out[count] = attrib[*i];
            }
            return out;
        }

        static double& getValue(attribution_type const& scope_sizes, cpd_type& cpd, attribution_type const& assignment){
            int count;
            //TODO assert assignment.size() == scope.size();
            attribution_type::const_iterator i_assig = assignment.begin();
            attribution_type::const_iterator i_scope = scope_sizes.begin();

            count = *i_assig;
            i_assig++; i_scope++;
            for(;i_assig != assignment.end() && i_scope != scope_sizes.end();
            i_assig++, i_scope++){
                count = count*(*i_scope) + *i_assig;
            }
            return cpd[count];
        }

    public:
        scope_type const getScope(){return scope;};

        cpd_type const getCPD(){return cpd;};

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
            cpd   = values;
            scope = scope_list;
            scope_sizes = attribution_type(scope_list.size(), 1);
            //int cpd_size = 1;
            int count = 0;
            for(scope_type::const_iterator i = scope.begin();
            i != scope.end(); i++, count++){
                scope_sizes[count] = i->getValues().size();
                //cpd_size *= i->getValues().size()
            }
        };

        /**
        *Normalizes the factor towards it's first variable.
        **/
        void normalize(){
            double z = 0.0;
            if(scope.size() == 1){
                for(cpd_type::const_iterator iter = cpd.begin(); iter != cpd.end(); iter++){
                    z += *iter;
                }
                if(z != 0){
                    for(cpd_type::iterator iter = cpd.begin(); iter != cpd.end(); iter++){
                        *iter /= z;
                    }
                }
            }else{
                size_t cpd_size, var_size, sub_cpd_size;
                cpd_size = cpd.size();
                var_size = scope[0].getValues().size();
                sub_cpd_size = cpd_size / var_size;

                //For each value of the first variable in the cpd:
                for(size_t iter = 0; iter < sub_cpd_size; iter++){
                    z = 0.0;
                    //For each assignment in the cpd which the first variable is equals to 'first_var_iter':
                    for(size_t first_var_iter = 0; first_var_iter < var_size; first_var_iter++){
                        z += cpd[iter + sub_cpd_size*first_var_iter];
                        }
                    //For each assignment in the cpd which the first variable is equals to 'first_var_iter':
                    //If z is not zero or way too small do:
                    if(z > 1e-10 || -z > 1e-10){
                        for(size_t first_var_iter = 0; first_var_iter < var_size; first_var_iter++){
                            //Divide by the sum of the values, normalize
                            cpd[iter + sub_cpd_size*first_var_iter] /= z;
                        }
                    }
                }
            }
        }

        /**
        * Serialization method. Loads the factor into a graphlab file.
        **/
        void save(graphlab::oarchive& oarc) const {
            oarc << scope << scope_sizes << cpd;
        }

        /**
        * Serialization method. Reads the factor from a graphlab file.
        **/
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
        Factor operator+=(Factor const& rhs){
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

        Factor setEvidence(scope_type const& vars, std::vector<std::string> const& values) const{
            std::vector<size_t> int_values(values.size(), 0);

            for(int count = 0; count < vars.size() && count < values.size(); count++){
                int_values[count] = vars[count].positionOfValue(values[count]);
            }

            return setEvidence(vars, int_values);
        }

        inline Factor setEvidence(RandomVariable const& var, std::string const& value) const{
            return setEvidence(var, var.positionOfValue(value));
        }

        Factor setEvidence(scope_type const& vars, std::vector<size_t> const& values) const{
            Factor new_factor(this->scope, this->cpd);

            //* It seems a lot of fors chained.
            //* Don't worry, the complexity of it is still O(vars.size()*new_factor.scope.size()).
            //For each variable in the vars vector
            for(int var_count = 0; var_count < vars.size() && var_count < values.size(); var_count++){
                int pre_size = 1;
                //Find out if it is in the factor, suppose no matching ordenation
                for(size_t var_pos = 0; var_pos < new_factor.scope.size(); var_pos++){
                    if(new_factor.scope[var_pos] == vars[var_count]){
                        //* Here I break the loop thoughout the cpd in several smaller ones in order to evidence correctly
                        //* pre_partition breaks the cpd in chunks of the size of a cpd made of the variables that preced the variable to evidence
                        //* curr_part breaks the remaining in values[var_pos].size() chunks
                        //For each partion created using pre_partition
                        for(int pre_partition = 0; pre_partition < this->cpd.size()/pre_size; pre_partition++){
                            size_t var_values_size = new_factor.scope[var_pos].getValues().size();
                            //For each value of the evidenced variable
                            for(int value = 0; value < var_values_size; value++){
                                double& to_update = new_factor.cpd[pre_size*pre_partition + var_values_size*value];
                                if(value == values[var_pos]){
                                    //If it is the main variable set it to one, useful in bayesian nets
                                    if(var_pos == 0){
                                        to_update = (to_update == 0.)? 0. : 1.;
                                    }
                                }else{
                                    to_update = 0.;
                                }
                            }
                        }
                        break;
                    }
                    pre_size *= new_factor.scope[var_pos].getValues().size();
                }
            }
            return new_factor;
        }

        Factor setEvidence(RandomVariable const& var, size_t const& value) const{
            return setEvidence(scope_type(1, var), std::vector<size_t>(1, value));
        }

        //Variables that belong to both factors should appear in the same order
        double distanceTo(const Factor& other) const{
            double out = 0.0;
            //Reduce both to the same scope
            scope_type intersec_scopes = scopeIntersection(this->scope, other.scope);
            Factor left = this->marginalizeAllBut(intersec_scopes);
            Factor right = other.marginalizeAllBut(intersec_scopes);
            for(size_t iter = 0; iter < left.cpd.size() && iter < right.cpd.size(); iter++){
                out += pow(left.cpd[iter]*right.cpd[iter], 2);
            }
            return out;
        }

        friend std::ostream& operator<<(std::ostream& os, const Factor& fact){
            scope_type::const_iterator var = fact.scope.begin();
            os << "[< ";
            if(var != fact.scope.end()){
                //os << var->getName(); var++;
                for(;var != fact.scope.end(); var++){
                    os << " " << var->getName();
                }
            }
            os << "> ";
            if(fact.cpd.size() > 0){
                //os << fact.cpd[0];
                for(unsigned int i=0; i < fact.cpd.size(); i++){
                    os<< " " << fact.cpd[i];
                };
            }
            os << "]";
            return os;
        };

        static Factor parseLine(std::stringstream& ss){
            std::string name, tmp;
            std::list<RandomVariable> v_list;
            std::list<double> d_list;

            if(ss >> tmp && tmp == "[<"){
                int pos = ss.tellg();
                while(ss >> tmp && tmp != ">"){
                    ss.seekg(pos);
                    v_list.push_back(RandomVariable::parseLine(ss));
                    pos = ss.tellg();
                }

                while(ss >> tmp && tmp != "]"){
                    d_list.push_back(std::atof(tmp.c_str()));
                }

            }
            scope_type scope(v_list.size());
            std::copy(v_list.begin(), v_list.end(), scope.begin());

            cpd_type cpd(d_list.size());
            std::copy(d_list.begin(), d_list.end(), cpd.begin());
            return Factor(scope, cpd);
        }

        /*friend std::istream& operator>>(std::istream& input, Factor& fact){
            std::string tmp;
            input >> tmp;
            if(tmp != "[<"){
                //std::cerr << "Factor parser: Invalid sintax" << std::endl;
                input.setstate(std::ios::failbit);
            }

            RandomVariable tmp_var;
            std::list<RandomVariable> l_vars;
            while(1){
                input >> tmp_var;
                if(input.fail()){
                    break;
                }
                l_vars.push_back(tmp_var);
            }
            //TODO read variables and build scope
            fact.scope = scope_type(l_vars.size());
            std::copy(l_vars.begin(), l_vars.end(), fact.scope.begin());
            fact.scope_sizes = attribution_type(fact.scope.size(), 1);
            size_t count = 0;
            for(scope_type::const_iterator i = fact.scope.begin();
            i != fact.scope.end(); i++, count++){
                fact.scope_sizes[count] = i->getValues().size();
            }
            std::list<double> v_list;
            double value;
            while(1){
                input >> value;
                if(input.fail()){
                    break;
                }
                v_list.push_back(value);
            }
            fact.cpd = cpd_type(v_list.size());
            std::copy(v_list.begin(), v_list.end(),fact.cpd.begin());
            return input;
        }*/
};

/*Factor vertex and a variable vertex are in the same class
because graphlab does not work pretty well with polymorphism. It doesn't save
all the fields if the nodes are not homogeneously sized.*/
class Vertice {
    private:
        Factor message;
        Factor factor;
        RandomVariable rVariable;
        graphlab::vertex_id_type var_id;
        static size_t id;

    public:
        enum {FACTOR, VARIABLE} type;

        int iter_count;

        explicit Vertice(int id): var_id(id), iter_count(0){}
        explicit Vertice(): iter_count(0){}

        Vertice(Factor& fact, int vid = id++): var_id(vid), iter_count(0){
            type = FACTOR;
            factor = fact;
            this->setMessage(fact);
        }

        Vertice(RandomVariable& var, int vid = id++): var_id(vid), iter_count(0){
            type      = VARIABLE;
            rVariable = var;
            cpd_type tmp_cpd(var.getValues().size(), 1.0);
            Factor tmp_fact(scope_type(1, var), tmp_cpd);
            this->setMessage(tmp_fact);
        }

        static size_t getTotalId(){return id;}
        VerticeObject& getInfo(){
            if(getType() == FACTOR){
                return factor;
            }
            return rVariable;
        }

        void setInfoFactor(Factor const& info){
            factor = info;
        }

        void setInfoRandomVariable(RandomVariable const& info){
            rVariable = info;
        }

        VerticeObject const& getInfo() const{
            if(getType() == FACTOR){
                return factor;
            }
            return rVariable;
        }

        graphlab::vertex_id_type const& vid() const{ return var_id;}

        inline Factor const& getMessage() const{ return message;}

        //TODO check how I can change it to setMessage(Factor& msg)
        inline void setMessage(Factor const& msg){ message = msg;}

        inline int getType() const{ return type;}

        void save(graphlab::oarchive& oarc) const {
            oarc << message << var_id << id << type;
            if(type == FACTOR){
                oarc << factor;
            }else{
                oarc << rVariable;
            }
        }

        void load(graphlab::iarchive& iarc) {
            iarc >> message >> var_id >> id >> type;
            if(type == FACTOR){
                iarc >> factor;
            }else{
                iarc >> rVariable;
            }
        }
        friend inline bool operator==(const Vertice& v1, const Vertice& v2){
            return v1.vid() == v2.vid();
        }

        friend std::ostream& operator<<(std::ostream& os, const Vertice& vertice){
            if(vertice.getType() == FACTOR){
                os << "FACTOR" << vertice.var_id << vertice.factor;
            }else{
                os << "VARIABLE" << vertice.var_id << vertice.rVariable;
            }
            return os;
        };

        static Vertice parseLine(std::stringstream& ss){
            int id;
            std::string tmp;
            Vertice v;
            if(ss >> id){
                if(ss >> tmp){
                    if (tmp == "Factor"){

                        //parse factor and create factor
                        Factor fact = Factor::parseLine(ss);
                        v = Vertice(fact, id);
                        //add factor vertex
                    }else{
                        if(tmp == "Variable"){
                            //parse variable and create it
                            RandomVariable var = RandomVariable::parseLine(ss);
                            v = Vertice(var, id);
                            // add variable vertex
                        }
                    }
                }
            }else{
                //todo error
            }
            return v;
        }

        /*friend std::istream& operator>>(std::istream& input, Vertice& vertice){
            std::string tmp;
            input >> tmp;
            if(tmp == "FACTOR"){
                vertice.type = FACTOR;
                input >> vertice.var_id;
                input >> vertice.factor;
                vertice.setMessage(vertice.factor);
            }else{
                if(tmp == "VARIABLE"){
                    vertice.type = VARIABLE;
                    input >> vertice.var_id;
                    input >> vertice.rVariable;
                    cpd_type tmp_cpd(vertice.rVariable.getValues().size(), 1.0);
                    Factor tmp_fact(scope_type(1, vertice.rVariable), tmp_cpd);
                    vertice.setMessage(tmp_fact);
                }
            }
            return input;
        };*/
};

size_t Vertice::id = 0;

typedef graphlab::distributed_graph<Vertice, graphlab::empty> GraphLab_type;
typedef std::pair<Vertice*, Vertice*> Edge_type;

class convergence_message: public graphlab::IS_POD_TYPE{
    public:
        convergence_message():value(0.0){};
        convergence_message(double val):value(val){};
        double value;
        convergence_message& operator+=(convergence_message const& other){
            value += other.value;
            return *this;
        }
};

class InferenceEngine:  public graphlab::ivertex_program<GraphLab_type, Factor, convergence_message>,
                        public graphlab::IS_POD_TYPE{
    public:
        InferenceEngine(){};
        virtual void runInference(graphlab::distributed_control& dist_ctrl, GraphLab_type& dgraph) =0;
};

class BeliefProp : public InferenceEngine{
    public:
        static const double BOUND = 1e-4;
        static const int ITER_BOUND = 10000;//TODO USE THIS TO LIMIT ITERATIONS

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
                //std::cout << "gather from factor " << edge.target().data().getMessage() <<std::endl;
                return edge.target().data().getMessage();
            }else{
                //if(vertex.data().getType() == Vertice::VARIABLE){
                //std::cout << "gather from variable" << edge.source().data().getMessage() << std::endl;
                return edge.source().data().getMessage();
                //}
            }
        }

        // Use the total belief of adjacent factors to update this factor
        void apply(icontext_type& context, vertex_type& vertex,
        const gather_type& total) {
            vertex.data().iter_count++;
            if(vertex.data().getType() == Vertice::FACTOR){
                //std::cout << "Factor Apply: Updated to:" << std::endl;
                Vertice& fact_vert = vertex.data();
                Factor& fact = (Factor&) fact_vert.getInfo();
                //std::cout << "--*--" << std::endl<< std::endl << fact; std::cout << total;
                Factor new_fact = (fact*total).marginalizeAllBut(fact.getScope());
                //std::cout << "--*--" << std::endl;
                vertex.data().setMessage(new_fact);
                //vertex.data() = fact_vert;
                //std::cout << vertex.data().getMessage() << std::endl;
            }else{
                if(vertex.data().getType() == Vertice::VARIABLE){
                    //std::cout << "Variable Apply: Updated to:" << std::endl;
                    Vertice& var_vert = vertex.data();
                    RandomVariable& var = (RandomVariable&) var_vert.getInfo();
                    vertex.data().setMessage(total.marginalizeAllBut(var));
                    //vertex.data() = var_vert;
                    //std::cout << vertex.data().getMessage() <<std::endl;
                }
            }
        }
        //TODO finish bp
        // The scatter edges depend on whether the pagerank has converged
        edge_dir_type scatter_edges(icontext_type& context,
                                    const vertex_type& vertex) const {
        //    if (perform_scatter) return graphlab::OUT_EDGES;
            //else
            //return graphlab::NO_EDGES;
            return graphlab::ALL_EDGES;
        }

        void scatter(icontext_type& context, const vertex_type& vertex,
                     edge_type& edge) const {
            bool isAFactor = vertex.data().getType() == Vertice::FACTOR;
            vertex_type other_vertex = (isAFactor)? edge.target() : edge.source();
            cpd_type fact_cpd, var_cpd;

            Factor *fact, *var;

            if(isAFactor){
                fact = &(Factor&) vertex.data().getInfo();
                var = &(Factor&) edge.target().data().getMessage();
            }else{
                //if(vertex.data().getType() == Vertice::VARIABLE){
                fact = &(Factor&) edge.source().data().getInfo();
                var = &(Factor&) vertex.data().getMessage();
                //}
            }

            //std::cout << "BP: convergence: FACTOR: " << *fact << std::endl << "VAR: "<< *var << std::endl;
            double residual = fact->distanceTo(*var);
            if(residual > BeliefProp::BOUND &&
               other_vertex.data().iter_count < BeliefProp::ITER_BOUND){
                edge.source();
                context.signal(other_vertex, residual);
            }
        }
    };

class EvidenceEngine{
    private:
        static scope_type variables;
        static std::vector<std::string> values;
    public:

        EvidenceEngine(scope_type const& vars, std::vector<std::string> const& vals){
            variables = vars;
            values = vals;
        }

        // Use the total to update this vertex
        static void setEvidenceForVertex(GraphLab_type::vertex_type& vertex) {
            bool isAFactor = vertex.data().getType() == Vertice::FACTOR;
            Factor *fact;

            if(isAFactor){
                fact = &(Factor&) vertex.data().getInfo();
                vertex.data().setMessage(fact->setEvidence(variables, values));
            }else{
                fact = &(Factor&) vertex.data().getMessage();
                vertex.data().setMessage(Factor(fact->getScope(), cpd_type(fact->getCPD().size(), 1.)));
            }
        }

        void setEvidence(GraphLab_type& graph){
            graph.finalize();
            graph.transform_vertices(setEvidenceForVertex);
        }
};

class Graph{
    private:
        GraphLab_type dgraph;
        graphlab::distributed_control& dist_ctrl;
        bool dgraph_done;

        class graph_writer{
        public:
            std::string save_vertex(GraphLab_type::vertex_type v){
                std::stringstream sstream;
                if(v.data().getType() == Vertice::FACTOR){
                    Vertice& fact_vert = v.data();
                    Factor& fact = (Factor&) fact_vert.getInfo();
                    sstream << fact << std::endl;
                }else{
                    if(v.data().getType() == Vertice::VARIABLE){
                        Vertice& var_vert = v.data();
                        RandomVariable& var = (RandomVariable&) var_vert.getInfo();
                        sstream << var << std::endl;
                    }else{
                        sstream << "-*-" << std::endl;
                    }
                }
                return sstream.str();
            }

            std::string save_edge(GraphLab_type::edge_type e){
                std::stringstream sstream;
                sstream << e.source().data().vid() << " -> " << e.target().data().vid() << std::endl;
                return sstream.str();
            }
        };

        static std::pair<size_t, size_t> parseEdge(std::stringstream& ss){
            int idSource, idTarget;
            std::string tmp;

            if(ss >> idSource){
                if(ss >> tmp){
                    if (tmp == "->"){
                        ss >> idTarget;
                        return std::pair<size_t, size_t>(idSource, idTarget);
                    }
                }
            }else{
                return std::pair<size_t, size_t>(-1, -1);
            }
        }

        static bool line_parser(GraphLab_type& graph, std::string const& file_name, std::string const& line){
            std::stringstream sstream(line);
            int pos = sstream.tellg(), id;
            std::string tmp;
            sstream >> id;
            if(sstream >> tmp){
                sstream.seekg(pos);
                if(tmp == "->"){
                    std::pair<size_t, size_t> edge = parseEdge(sstream);
                    if(edge.first == -1){
                        return false;
                    }
                    graph.add_edge(edge.first, edge.second);
                }else{
                    Vertice v = Vertice::parseLine(sstream);
                    graph.add_vertex(v.vid(), v);
                }
            }
            return true;
        }

    public:
        explicit Graph(graphlab::distributed_control& dc): dist_ctrl(dc), dgraph_done(false), dgraph(GraphLab_type(dc)){};

        void addVertice(Vertice& vert){
            dgraph_done = false;
            dgraph.add_vertex(vert.vid(), vert);
        }

        inline void addEdge(Vertice& v1, Vertice& v2){
            Vertice* fact;
            Vertice* var;
            dgraph_done = false;
            if(v1.getType() == Vertice::FACTOR && v2.getType() == Vertice::VARIABLE){
                fact = &v1;
                var =  &v2;
                dgraph.add_edge(fact->vid(), var->vid());
            }else{
                if(v1.getType() == Vertice::VARIABLE && v2.getType() == Vertice::FACTOR){
                    fact = &v2;
                    var  =  &v1;
                    dgraph.add_edge(fact->vid(), var->vid());
                }//TODO raise exception if it is a wrong kind of edge
            }

        }

        void startDGraph(){
            if(!isDGraphUp()){
                dgraph.finalize(); // important
                dgraph_done = true;
            }
        }

        inline bool isDGraphUp(){return dgraph.is_finalized() && dgraph_done;}

        //TODO
        cpd_type getPrior(Vertice const& var);

        void runInference(InferenceEngine& eng){//, EvidenceEngine& evidEng){
            //TODO run evidence
            startDGraph();
            eng.runInference(dist_ctrl, dgraph);
            //TODO undo evidence
        }

        GraphLab_type& getDGraph(){ return dgraph;}

        void save(){
            // params: destination, writer, save as .gzip, save vertex, save edges, files per machine
            dgraph.save("/home/breno/BNETout", graph_writer(), false, true, true, 1);
        }

        void load(std::string prefix){
            dgraph.load(prefix, line_parser);
        }

        /*friend std::ostream& operator<<(std::ostream& os, const Graph& graph){
            os << "Graph<";//<< graph.getDGraph();
            //TODO
            os << ">";
            return os;
        }*/
};

//      std::pair< Factor, Variable >                Factor -> Variable
typedef std::pair<Vertice, Vertice> BNet_vertice;
class Learning_Complete_Data_Engine{
    private:
        static DataSet data_set;

        static Factor& absorb_data(Factor& fact, DataSet& data){
            data_table_type::iterator element;
            data_object_type const& header = data.getHeader();
            scope_type const& scope = fact.getScope();
            std::vector<int> scope_map(fact.getScope().size(),0);
            attribution_type attrib(scope.size(),0);

            bool missing_var;
            //Map the data header to the factor scope
            for(size_t attrib_pos = 0; attrib_pos < scope.size(); attrib_pos++){
                missing_var = true;
                //std::cout << std::endl<< "Header size: " << header.size();
                for(int map_pos = 0; map_pos < header.size(); map_pos++){
                    //bool var = (scope[attrib_pos].getName() == header[map_pos]);
                    //std::cout << std::endl<< "VAR: " << ((var)?"T":"F") << header.size();
                    if(scope[attrib_pos].getName() == header[map_pos]){
                        scope_map[attrib_pos] = map_pos;
                        missing_var = false;
                        break;
                    }
                }
                if(missing_var){
                    std::cerr << std::endl << "Invalid header!!"<< fact << std::endl;
                    //TODO raise exception and deal with it
                }
            }
            //Read each tuple and update the factor
            for(element = data.begin(); element != data.end(); element++){
                for(size_t attrib_pos = 0; attrib_pos < scope.size(); attrib_pos++){
                    attrib[attrib_pos] = scope[attrib_pos].positionOfValue((*element)[scope_map[attrib_pos]]);
                    //TODO raise an error if the value is invalid
                }
                fact[attrib] += 1;
                /*std::cout << "atrib: [";
                for(size_t attrib_pos = 0; attrib_pos < scope.size(); attrib_pos++){
                    std::cout << attrib[attrib_pos] << ", ";
                }
                std::cout << "\b]" << std::endl;*/
            }
            fact.normalize();
            return fact;
        }

    static void vertex_learn(GraphLab_type::vertex_type& vertex) {
        //std::cout << "Performing learning on a " << std::endl;
        if(vertex.data().getType() == Vertice::FACTOR){
            //std::cout << "factor";// << (Factor&) vertex.data();
            Vertice& fact_vert = vertex.data();
            Factor& fact = (Factor&) fact_vert.getInfo();
            //std::cout << std::endl << data_set;
            fact = absorb_data(fact, data_set);
            fact_vert.setInfoFactor(fact);
            //std::cout << "Learned!! = " << ((Factor&) (fact_vert.getInfo()));
        }else{
            //std::cout << "variable";
        }

        //std::cout << std::endl;
    }
    public:
        Learning_Complete_Data_Engine(DataSet& data){data_set = data;};

        void learn(GraphLab_type& graph){
            graph.finalize();
            graph.transform_vertices(vertex_learn);
        }
};

DataSet Learning_Complete_Data_Engine::data_set = DataSet();

class BeliefNetwork{
    private:
        Graph fgraph;
        graphlab::distributed_control& distributed_ctrl;
        InferenceEngine* i_engine;
        EvidenceEngine* e_engine;

        void learnParametersCompleteData(DataSet& data){
            Learning_Complete_Data_Engine l_eng(data);
            l_eng.learn(fgraph.getDGraph());
        };
        void learnParametersIncompleteData(DataSet& data){};

        class graph_writer{
            public:
                std::string save_vertex(GraphLab_type::vertex_type v){
                    std::stringstream sstream;
                    if(v.data().getType() == Vertice::FACTOR){
                        Vertice& fact_vert = v.data();
                        Factor& fact = (Factor&) fact_vert.getInfo();
                        sstream << fact.getScope()[0] << " " << fact << std::endl;
                    }
                    return sstream.str();
                }
                std::string save_edge(GraphLab_type::edge_type e){ return "a -> b"; }
        };

        class BeliefsSet{
            private:
                scope_type target_variables;
                graphlab::dht<std::string, Factor> belief_map;
                GraphLab_type& graph;
                static BeliefsSet* singleton;

                bool isTarget(GraphLab_type::vertex_type const& vertex){
                    Vertice vert = vertex.data();
                    if (vert.getType() == Vertice::VARIABLE){
                        for(scope_type::iterator var = target_variables.begin(); var != target_variables.end(); var++){
                            if(((RandomVariable&) vert.getInfo()) == *var){
                                return true;
                            }
                        }
                    }
                    return false;
                }

                /*void setBelief(GraphLab_type::vertex_type& vertex){
                    if(isTarget(vertex)){
                        Vertice& v = (Vertice&) vertex.data();
                        std::string name = ((RandomVariable&) v.getInfo()).getName();
                        belief_map.set(name, v.getMessage());
                    }
                }*/

                static void setBelief(GraphLab_type::vertex_type& vertex){
                    //The reason why this function uses the singleton pointer is described here
                    // http://stackoverflow.com/questions/15841338/c-unresolved-overloaded-function-type
                    //BeliefsSet.setBelief is a non-static-member-function, hence it has a special signature, that cannot be used in updateBeliefs()
                    if(singleton->isTarget(vertex)){
                        Vertice& v = (Vertice&) vertex.data();
                        std::string name = ((RandomVariable&) v.getInfo()).getName();
                        singleton->belief_map.set(name, v.getMessage());
                    }
                }

            public:
                BeliefsSet(GraphLab_type& g,
                                graphlab::distributed_control& d_ctrl,
                                scope_type const& target): graph(g), target_variables(target), belief_map(d_ctrl){
                    singleton = this;
                    d_ctrl.barrier();
                    updateBeliefs();
                    d_ctrl.barrier();
                }

                void updateBeliefs(){
                    graph.transform_vertices(setBelief);
                }

                Factor getBelief(RandomVariable const& var){
                    return belief_map.get(var.getName()).second;
                }
    };

    public:
        BeliefNetwork(graphlab::distributed_control& dc, InferenceEngine* eng): distributed_ctrl(dc), fgraph(dc), i_engine(eng){};

        BNet_vertice addVariable(RandomVariable& var, Factor& cpd){
            Vertice fact(cpd), rvar(var);
            fgraph.addVertice(fact);
            fgraph.addVertice(rvar);
            fgraph.addEdge(fact, rvar);
            return BNet_vertice(fact, rvar);
        };

        /**
        * The cpd is not being checked, so be careful when you add an edge
        **/
        void addEdge(BNet_vertice& source, BNet_vertice& target){
            //Linking the variable of the parent to the Factor of the child
            fgraph.addEdge(target.first, source.second);
        }
        //void addEdge(RandomVariable& source, RandomVariable& target);

        Factor infer(RandomVariable& variable){
            BeliefProp i_eng;
            fgraph.startDGraph();
            fgraph.runInference(i_eng);
            distributed_ctrl.barrier();
            //std::cout << "Inference done" << std::endl << "Runing inference engine" << std::endl;
            BeliefsSet bs(fgraph.getDGraph(), distributed_ctrl, scope_type(1, variable));
            return bs.getBelief(variable);
            //return Factor();
        }


        //TODO void evidence
        void setEvidence(scope_type const& vars, std::vector<std::string> const& vals){
            EvidenceEngine eng(vars, vals);
            eng.setEvidence(fgraph.getDGraph());
            distributed_ctrl.barrier();
        }

        void learnParameters(DataSet& data){
            //if(data.isComplete()){
                learnParametersCompleteData(data);
                distributed_ctrl.barrier();
            /*}else{
                learnParametersIncompleteData(data);
            }*/
        }

        //void learnStructure(DataSet data);

        void save(){
            //fgraph.getDGraph().save("/home/breno/BNETout", graph_writer(), false, true, true, 1);
            fgraph.save();
        }

        void load(std::string prefix){
            fgraph.load(prefix);
            distributed_ctrl.barrier();
        }

        /*friend std::ostream& operator<<(std::ostream& os, const BeliefNetwork& bnet){
            os << "Bnet<" << bnet.fgraph << ">";
            return os;
        }*/
};

BeliefNetwork::BeliefsSet* BeliefNetwork::BeliefsSet::singleton = NULL;

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
    values_rain_grass[2] = 0.1; values_rain_grass[3] = 0.1;
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
    << "V1: "<< (RandomVariable&) v1.getInfo() << std::endl \
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

    DataSet data("test.csv");
    std::cout << "DATASET read" << std::endl;
    bn.learnParameters(data);
    std::cout << "learning done"<< std::endl;
    bn.save();
    std::cout << "saving done"<< std::endl;
    std::cout << "P(wet) = " << bn.infer(grass) << std::endl;
    return 0;
}

int test3(graphlab::distributed_control& dc){
    dc.barrier();
    std::cout << "creating bn"<< std::endl;
    BeliefProp i_eng;
    BeliefNetwork bn(dc, &i_eng);
    std::cout << "bn created"<< std::endl;

    bn.load("/home/breno/Documents/git-repos/graphlab/apps/bnets/BNETout.1_of_1");
    bn.save();
    std::vector<std::string> triple(3 ,"");
    triple[0] = "false"; triple[1] = "true"; triple[2] = "maybe";
    std::string grass_str = std::string("wet");
    RandomVariable grass(grass_str, triple);

    DataSet data("test.csv");
    std::cout << "DATASET read" << std::endl;
    bn.learnParameters(data);
    std::cout << "learning done"<< std::endl;
    //bn.save();
    std::cout << "saving done"<< std::endl;
    std::cout << "P(wet) = " << bn.infer(grass) << std::endl;
    return 0;
}

int experiment1(graphlab::distributed_control& dc, std::string path, std::string net, std::string var_str, std::string data_str, std::string out_str){
    dc.barrier();
    std::cout << "creating bn"<< std::endl;
    BeliefProp i_eng;
    BeliefNetwork bn(dc, &i_eng);
    double begin_inference, end_inference, end_learning;
    std::stringstream var_ss(var_str);
    RandomVariable var = RandomVariable::parseLine(var_ss);
    //std::cout << "bn created"<< std::endl;

    bn.load(path+net);
    bn.save();
    std::vector<std::string> triple(3 ,"");

    DataSet data(path+data_str);
    //std::cout << "DATASET read" << std::endl;
    begin_inference = clock();
    Factor grass_infer = bn.infer(var);
    end_inference = clock();
    bn.learnParameters(data);
    end_learning = clock();
    //std::cout << "learning done"<< std::endl;
    //bn.save();
    //std::cout << "saving done"<< std::endl;
    //std::cout << "P(wet) = " << bn.infer(grass) << std::endl;
    if(dc.procid() == 0){
        std::cout << "times" << (end_inference - begin_inference)*1000 / CLOCKS_PER_SEC << "ms\t" <<(end_learning - end_inference) *1000/ CLOCKS_PER_SEC << "ms";
        std::ofstream output;
        output.open(out_str.c_str(), std::ios::out | std::ios::app);
        output << net << "\t" << var.getName() << "\t" << dc.numprocs() << "\t" << (end_inference - begin_inference)*1000 / CLOCKS_PER_SEC << "\t" <<(end_learning - end_learning) *1000/ CLOCKS_PER_SEC << std::endl;
        output.close();
    }
    return 0;
}


int main(int argc, char** argv){
    graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;

    double begin, end;

/*    std::cout << "test1" << std::endl;
    begin = clock();
    test1(dc); // Testing FactorGraph
    end = clock();
    std::cout << "Time Elapsed: " << (end-begin)/ CLOCKS_PER_SEC << "s." << std::endl;

    std::cout << "test2" << std::endl;
    //if you are the node 0
    begin = clock();
    test2(dc); // Testing BNet
    end = clock();
    std::cout << "Time Elapsed: " << (end-begin)/ CLOCKS_PER_SEC << "s." << std::endl;
    //}

    std::cout << "test3" << std::endl;
    //if you are the node 0
    begin = clock();
    test3(dc); // Testing BNet
    end = clock();
    std::cout << "Time Elapsed: " << (end-begin)/ CLOCKS_PER_SEC << "s." << std::endl;
    //} */

    /*for(int i = 0; i < 100; i++){
        experiment1(dc,
                        "/home/breno/Documents/git-repos/graphlab/apps/bnets/samples/",
                        "( rain |  false true )",
                        "test.gl",
                        "test.csv",
                        "/home/breno/pgm_test" );
    }*/

    for(int i = 0; i < 100; i++){
        experiment1(dc, "/home/breno/Documents/git-repos/graphlab/apps/bnets/samples/",
                        "BBNet-hailfinder.gl",
                        "P(CapInScen[LessThanAve,Average,MoreThanAve]|AMCINInScen,CapChange)",
                        "hailfinder_500_perc30_missing_data.txt",
                        "/home/breno/pgm_test" );
    }

    /*for(int i = 0; i < 100; i++){
        experiment1(dc,
                        "/home/breno/Documents/git-repos/graphlab/apps/bnets/samples/",
                        "BBNet-alarm.gl",
                        "P(CO[Low,Normal,High]|HR,StrokeVolume)",
                        "alarm_1000_perc30_missing_data.txt",
                        "/home/breno/pgm_test" );
    }*/

    std::cout << "closing mpi" << std::endl;
    graphlab::mpi_tools::finalize();
    std::cout << "mpi closed" << std::endl;
}

#endif //FACTOR_GRAPH_NETWORK_UTILS_HPP
