/*
 *  \author: Breno Carvalho
 */

#ifndef FACTOR_GRAPH_NETWORK_UTILS_HPP
#define FACTOR_GRAPH_NETWORK_UTILS_HPP

#include <graphlab.hpp>

/*
#include <../../../toolkits/graphical_models/factors/factor_graph.hpp>
#include <../../../toolkits/graphical_models/factors/bp_vertex_program.hpp>
#include <../../../toolkits/graphical_models/factors/bp_graph_data.h>
*/
#include <../../../apps/bnets/csv_reader.h>
#include <math.h>

//used in experiments
#include <ctime>
#include <iostream>
#include <fstream>

class VerticeObject;
class RandomVariable;

typedef std::vector<RandomVariable> scope_type;
typedef std::vector<int> attribution_type;
typedef std::vector<double> cpd_type;

class VerticeObject{
    public:
        virtual scope_type const& getScope() const =0;
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

        int positionOfValue(std::string const& value) const{
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
            os << "( " << var.name << " |";
            std::vector<std::string>::const_iterator iter = var.values.begin();
            if( iter != var.values.end()){
                //os << *iter; iter ++;
                for(; iter != var.values.end(); iter++){
                    os << " " << *iter;
                }
            }
            os << " )";
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
                            //v_list.push_back("None");
                    }
                }
            }
            std::vector<std::string> values(v_list.size());
            std::copy(v_list.begin(), v_list.end(), values.begin());
            return RandomVariable(name, values);
        }

        //TODO check if it is safe
        scope_type const& getScope() const{return scope_type(1, *this);}
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

    public:
        static scope_type scopeUnion(scope_type const& s1,scope_type const& s2){
            scope_type out;
            for(scope_type::const_iterator iter = s1.begin();
            iter != s1.end(); iter++){
                out.push_back(*iter);
            }
            for(scope_type::const_iterator iter = s2.begin();
            iter != s2.end(); iter++){
                if(!hasVariable(s1, *iter)){
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

        //If the second scope variables are in the first scope, it returns an
        //attribution_type object of the same size of the second scope such that
        //the numbers in each position correspond to the position of the variable in the first scope
        //E.g.: attribution_type([a,b,c,d], [b, c]) = [1,2]
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

        //Returns an attribution of values in attrib filtered by attrib_map
        static attribution_type applyAttribMap(attribution_type const& attrib, attribution_type const& attrib_map){
            attribution_type out(attrib_map.size(),0);
            int count=0;
            for(attribution_type::const_iterator i = attrib_map.begin(); i != attrib_map.end(); i++, count++){
                out[count] = attrib[*i];
            }
            return out;
        }

        /*
        * \param scope_sizes is a list of the sizes of the variables in the cpd
        * \param cpd is the cpd it-self
        * \param assignment is a position assignment
        *This method returns the value of the position determined by assignment
        *E.g.:
        *assignment = [ 1, 0 ]
        *Fact AB =
        * A | B | cpd()
        * 0 | 0 | v0
        * 0 | 1 | v1
        * 1 | 0 | v2
        * 1 | 1 | v3
        * Thus scope_sizes = [ 2, 2 ], cpd = [ v0, v1, v2, v3 ] and the return value should be v2
        *
        */
        static double& getValue(attribution_type const& scope_sizes, cpd_type& cpd, attribution_type const& assignment){
            int count;
            //TODO assert assignment.size() == scope.size();
            attribution_type::const_reverse_iterator i_assig = assignment.rbegin();
            attribution_type::const_reverse_iterator i_scope = scope_sizes.rbegin();

            count = *i_assig;
            i_assig++; i_scope++;
            for(;i_assig != assignment.rend() && i_scope != scope_sizes.rend();
            i_assig++, i_scope++){
                count = count*(*i_scope) + *i_assig;
            }
            return cpd[count];
        }

        scope_type const& getScope() const{return scope;}

        attribution_type const& getScopeSizes() const {return scope_sizes;}

        cpd_type const& getCPD() const{return cpd;}

        static bool hasVariable(scope_type const& scope, RandomVariable const& var){
            for(scope_type::const_iterator iter = scope.begin();
                iter != scope.end(); iter++){
                if(var == *iter){
                    return true;
                }
            }
            return false;
        }

        inline bool hasVariable(RandomVariable const& var) const{
            return Factor::hasVariable(scope, var);
        }

        Factor(){}

        Factor(scope_type const& scope_list, cpd_type const& values){
            cpd   = values;
            scope = scope_list;
            scope_sizes = attribution_type(scope.size(), 1);
            //int cpd_size = 1;
            int count = 0;
            for(scope_type::const_iterator i = scope.begin();
            i != scope.end(); i++, count++){
                scope_sizes[count] = i->getValues().size();
                //cpd_size *= i->getValues().size()
            }
        }

        /**
        *Returns the supremum distance between two factors
        * returns max_i(|cpd_1[i]- cpd_2[i]|)
        * This distance is the same of the n-distance when n tends to infinity
        **/
        double sup_distance(Factor const& other_fact) const{
            double max = 0.0, tmp;
            if(this->cpd.size() != other_fact.cpd.size()){
                throw "Incompatible factor comparation";
            }
            for(int pos = 0; pos < cpd.size(); pos++){
                tmp = std::fabs(this->cpd[pos] - other_fact.cpd[pos]);
                max = (max > tmp)? max : tmp;
            }
            return max;
        }

        /**
        *Returns the n-distance between two factors
        * returns pow(sum(pow(|cpd_1[i]- cpd_2[i]|, n)), 1./n)
        *n must be bigger than 0
        **/
        double n_distance(Factor const& other_fact, int n) const{
            double sum = 0.0;
            if(this->cpd.size() != other_fact.cpd.size()){
                throw "Incompatible factor comparation";
            }
            if( n > 0){
                if( n > 1){
                    for(int pos = 0; pos < cpd.size(); pos++){
                        sum += std::pow(std::fabs(this->cpd[pos]-other_fact.cpd[pos]), n);
                    }
                    return std::pow(sum, 1./n);
                }else{
                    //Manhattan distance
                    for(int pos = 0; pos < cpd.size(); pos++){
                        sum += std::fabs(this->cpd[pos]-other_fact.cpd[pos]);
                    }
                    return sum;
                }
            }//if n < 1 then n is a invalid value
            throw "Hey, bad choice of 'n' for your n-distance, n must be bigger than 0";
        }

        /**
        *Returns a distance between two factors, used as convergence critereon
        **/
        double distance(Factor const& other_fact) const{
            return n_distance(other_fact, 2);
            //return n_distance(other_fact, 3);
            //return n_distance(other_fact, 4);
            //return sup_distance(other_fact);
        }

        /**
        *Normalizes the factor towards it's first variable.
        *This method changes the original factor
        *E.g. this takes this factor:
        * A | B | C | cpd()
        * 0 | 0 | 0 | v0
        * 0 | 0 | 1 | v1
        * 0 | 1 | 0 | v2
        * 0 | 1 | 1 | v3
        * 1 | 0 | 0 | v4
        * 1 | 0 | 1 | v5
        * 1 | 1 | 0 | v6
        * 1 | 1 | 1 | v7
        * And returns this other one:
        * A | B | C | cpd()
        * 0 | 0 | 0 | v0/(v0+v4)
        * 0 | 0 | 1 | v1/(v1+v5)
        * 0 | 1 | 0 | v2/(v2+v6)
        * 0 | 1 | 1 | v3/(v3+v7)
        * 1 | 0 | 0 | v4/(v0+v4)
        * 1 | 0 | 1 | v5/(v1+v5)
        * 1 | 1 | 0 | v6/(v2+v6)
        * 1 | 1 | 1 | v7/(v3+v7)
        **/
        void normalize(){
            if(scope.size() < 1){return;}
            double z;
            if(scope.size() == 1){
                z = 0.0;
                for(cpd_type::const_iterator iter = cpd.begin(); iter != cpd.end(); iter++){
                    z += *iter;
                }
                if(z != 0){
                    for(cpd_type::iterator iter = cpd.begin(); iter != cpd.end(); iter++){
                        *iter /= z;
                    }
                }
            }else{
                int cpd_size, var_size, sub_cpd_size;
                cpd_size = cpd.size();
                var_size = scope[0].getValues().size();
                sub_cpd_size = cpd_size / var_size;
                for(int sec_vals_count = 0; sec_vals_count < sub_cpd_size; sec_vals_count++){
                    z = 0.0;
                    for(int first_val_count = 0; first_val_count < var_size; first_val_count++){
                        z += cpd[first_val_count*sub_cpd_size + sec_vals_count];
                    }
                    if(std::fabs(z) > 1e-10){ // if (z==0) then do z=1
                        for(int first_val_count = 0; first_val_count < var_size; first_val_count++){
                            cpd[first_val_count*sub_cpd_size + sec_vals_count] /= z;
                        }
                    }
                }
            }
            /*double z = 0.0;
            if(scope.size() < 1){return;}
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
                int cpd_size, var_size, sub_cpd_size;
                cpd_size = cpd.size();
                var_size = scope[0].getValues().size();
                sub_cpd_size = cpd_size / var_size;

                //For each value of the first variable in the cpd:
                for(int first_var_iter = 0; first_var_iter < var_size; first_var_iter++){
                    z = 0.0;
                    //For each assignment in the cpd which the first variable is equals to 'first_var_iter':
                    for(int iter = 0; iter < sub_cpd_size; iter++){
                        z += cpd[iter + sub_cpd_size*first_var_iter];
                        }
                    //For each assignment in the cpd which the first variable is equals to 'first_var_iter':
                    //If z is not zero or way too small do:
                    if(z > 1e-10 || -z > 1e-10){
                        for(int iter = 0; iter < sub_cpd_size; iter++){
                            //Divide by the sum of the values, normalize
                            cpd[iter + sub_cpd_size*first_var_iter] /= z;
                        }
                    }
                }
            }*/
        }

        /**
        *Normalizes the factorin a way its elements sum to one.
        *This method behaves like normalize() when scope.size() == 1.
        *It always consider scope.size() as being equal to 1, no matter the actual scope.size() (as long it is bigger than 0)
        *This method changes the original factor
        *E.g. this takes this factor:
        * A | B | C | cpd()
        * 0 | 0 | 0 | v0
        * 0 | 0 | 1 | v1
        * 0 | 1 | 0 | v2
        * 0 | 1 | 1 | v3
        * 1 | 0 | 0 | v4
        * 1 | 0 | 1 | v5
        * 1 | 1 | 0 | v6
        * 1 | 1 | 1 | v7
        * And returns this other one:
        * A | B | C | cpd()
        * 0 | 0 | 0 | v0/(v0+v1+v2+v3+v4+v5+v6+v7)
        * 0 | 0 | 1 | v1/(v0+v1+v2+v3+v4+v5+v6+v7)
        * 0 | 1 | 0 | v2/(v0+v1+v2+v3+v4+v5+v6+v7)
        * 0 | 1 | 1 | v3/(v0+v1+v2+v3+v4+v5+v6+v7)
        * 1 | 0 | 0 | v4/(v0+v1+v2+v3+v4+v5+v6+v7)
        * 1 | 0 | 1 | v5/(v0+v1+v2+v3+v4+v5+v6+v7)
        * 1 | 1 | 0 | v6/(v0+v1+v2+v3+v4+v5+v6+v7)
        * 1 | 1 | 1 | v7/(v0+v1+v2+v3+v4+v5+v6+v7)
        **/
        void simple_normalize(){
            if(scope.size() < 1){return;}
            double z;
            for(cpd_type::const_iterator iter = cpd.begin(); iter != cpd.end(); iter++){
                z += *iter;
            }
            if(z != 0){
                for(cpd_type::iterator iter = cpd.begin(); iter != cpd.end(); iter++){
                    *iter /= z;
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
        }

        double const& operator[](attribution_type const& idx) const {
            return const_cast<Factor&>(*this)[idx];
        }

        /**
        *Return the product of two factors.
        *If they have the same scope each cell of one factor is multiplied by the correspondent cell of the second factor.
        * If they have different scopes the resulting factor will have the union of the scopes and would function like a *join*.
        *This method changes the original factor
        *E.g. this takes this factor:
        * A | B | C | cpd()
        * 0 | 0 | 0 | v0
        * 0 | 0 | 1 | v1
        * 0 | 1 | 0 | v2
        * 0 | 1 | 1 | v3
        * 1 | 0 | 0 | v4
        * 1 | 0 | 1 | v5
        * 1 | 1 | 0 | v6
        * 1 | 1 | 1 | v7
        *
        *and
        *
        * C | D | cpd()
        * 0 | 0 | w0
        * 0 | 1 | w1
        * 1 | 0 | w2
        * 1 | 1 | w3
        *
        * Their product returns:
        * A | B | C | D | cpd()
        * 0 | 0 | 0 | 0 | v0 * w0
        * 0 | 0 | 0 | 1 | v0 * w1
        * 0 | 0 | 1 | 0 | v1 * w2
        * 0 | 0 | 1 | 1 | v1 * w3
        * 0 | 1 | 0 | 0 | v2 * w0
        * 0 | 1 | 0 | 1 | v2 * w1
        * 0 | 1 | 1 | 0 | v3 * w2
        * 0 | 1 | 1 | 1 | v3 * w3
        * 1 | 0 | 0 | 0 | v4 * w0
        * 1 | 0 | 0 | 1 | v4 * w1
        * 1 | 0 | 1 | 0 | v5 * w2
        * 1 | 0 | 1 | 1 | v5 * w3
        * 1 | 1 | 0 | 0 | v6 * w0
        * 1 | 1 | 0 | 1 | v6 * w1
        * 1 | 1 | 1 | 0 | v7 * w2
        * 1 | 1 | 1 | 1 | v7 * w3
        *
        *no normalization was done here
        *
        **/

        Factor& operator*=(Factor const& rhs){
            scope_type old_scope = this->scope;
            attribution_type old_scope_sizes = this->scope_sizes;
            cpd_type old_cpd = this->cpd;
            this->scope = scopeUnion(this->scope, rhs.scope);
            //std::cout << "SCOPE (*=) 0:" << this->scope[0];
            this->scope_sizes = attribution_type(this->scope.size(), 1);
            int cpd_size = 1; int count = 0;
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
                //This for loop updates the last value in the atributtion.
                for(int val = 0; val < this->scope_sizes[this->scope_sizes.size()-1]; val++){
                    attrib[this->scope_sizes.size()-1] = val;
                    //multiplication
                    double mult = getValue(old_scope_sizes, old_cpd, applyAttribMap(attrib, attrib_map_left));
                    mult *= rhs[applyAttribMap(attrib, attrib_map_right)];
                    this->cpd[count] = mult;
                    count++;
                }
                for(int var_pos = this->scope_sizes.size()-1; var_pos >= 0; var_pos--){
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

        friend Factor& operator*(Factor lhs, Factor const& rhs){
            return lhs*=rhs;
        }

        //created only for compatibility with graphlab::aggregator_type
        Factor operator+=(Factor const& rhs){
            std::cout << std::endl << "Printing product of " << *this << "by "<< rhs << " that is " << (*this)*rhs << std::endl << std::endl;
            return (*this) *= rhs;
        }

        inline Factor marginalize(RandomVariable const& var) const{
            return marginalize(scope_type(1,var));
        }

        inline Factor marginalizeAllBut(RandomVariable const& var) const{
            return marginalizeAllBut(scope_type(1, var));
        }

        Factor marginalize(scope_type const& vars) const{
            if(this->scope.size() == 0){
                throw "Cannot marginalize an empty Factor";
            }
            if(vars.size() == 0){
                return *this;
            }
            scope_type new_scope = scopeMinus(this->scope, vars);
            if(new_scope.size() == 0){
                return Factor(new_scope, cpd_type(1,0.0));
            }
            int cpd_size = 1;
            for(scope_type::const_iterator scope_iterator = new_scope.begin();
            scope_iterator != new_scope.end(); scope_iterator++){
                cpd_size *= (int) scope_iterator->getValues().size();
                //std::cout << *scope_iterator << " " << scope_iterator->getValues().size() << " ";
            }
            //std::cout << " CPD_SIZE: " << cpd_size << std::endl;
            cpd_type new_cpd = cpd_type(cpd_size, 0.0);
            Factor new_factor(new_scope, new_cpd);

            attribution_type attrib = attribution_type(this->scope.size(), 0);
            attribution_type new_f_attrib_map = getAttribMap(this->scope, new_scope);

            int count = 0;
            while(count < this->cpd.size()){
                for(int val = 0; val < this->scope_sizes[0]; val++){
                    attrib[0] = val;
                    new_factor[applyAttribMap(attrib, new_f_attrib_map)] += (*this)[attrib];
                    count++;
                }

                for(int var_pos = 0; var_pos < this->scope_sizes.size(); var_pos++){
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
            std::vector<int> int_values(values.size(), 0);

            for(int count = 0; count < vars.size() && count < values.size(); count++){
                int_values[count] = vars[count].positionOfValue(values[count]);
            }

            return setEvidence(vars, int_values); // Call proper function to evidenciate
        }

        inline Factor setEvidence(RandomVariable const& var, std::string const& value) const{
            return setEvidence(var, var.positionOfValue(value)); // Call proper function to evidenciate
        }

        Factor setEvidence(scope_type const& vars, std::vector<int> const& values) const{
            Factor new_factor(this->scope, this->cpd);

            //* It seems a lot of fors chained.
            //* Don't worry, the complexity of it is still O(vars.size()*new_factor.cpd.size()).
            //For each variable in the vars vector
            for(int var_count = 0; var_count < vars.size() && var_count < values.size(); var_count++){
                int prev_var_size = 1;
                for(int var_pos = 0; var_pos < new_factor.scope.size(); var_pos++){
                    if(new_factor.scope[var_pos] == vars[var_count]){
                        int value = values[var_count];
                        int curr_var_size = new_factor.scope[var_pos].getValues().size();
                        int next_var_size = new_factor.cpd.size()/(curr_var_size*prev_var_size);
                        for(int prev_var_count = 0; prev_var_count < prev_var_size; prev_var_count++){
                            for(int next_var_count = 0; next_var_count < next_var_size; next_var_count++){
                                for(int value_count = 0; value_count < curr_var_size; value_count++){
                                    double& to_update = new_factor.cpd[(prev_var_count*curr_var_size + value_count)*next_var_size + next_var_count];
                                    if(!(value_count == value)){
                                        to_update = 0.0;
                                    }
                                    //std::cout << "curr_var_size = " << curr_var_size;
                                }
                            }
                        }
                        break;
                    }
                    prev_var_size *= new_factor.scope[var_pos].getValues().size();
                }
            }
            new_factor.normalize();
            return new_factor;
        }

        Factor setEvidence(RandomVariable const& var, int const& value) const{
            return setEvidence(scope_type(1, var), std::vector<int>(1, value));
        }

        friend std::ostream& operator<<(std::ostream& os, const Factor& fact){
            scope_type::const_iterator var = fact.scope.begin();
            os << "[<";
            if(var != fact.scope.end()){
                //os << var->getName(); var++;
                for(;var != fact.scope.end(); var++){
                    os << " " << *var;
                    //os << " " << var->getName();
                }
            }
            os << " >";
            if(fact.cpd.size() > 0){
                //os << fact.cpd[0];
                for(int i=0; i < fact.cpd.size(); i++){
                    os<< " " << fact.cpd[i];
                };
            }
            os << " ]";
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

            if(scope.size() == 0){
                throw std::invalid_argument( "Factor: invalid line to parse" + name );
            }

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
            int count = 0;
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
        static int id; // global id count

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
            cpd_type tmp_cpd(var.getValues().size(), 1);
            Factor tmp_fact(scope_type(1, var), tmp_cpd);
            tmp_fact.normalize();
            this->setMessage(tmp_fact);
        }

        static int getTotalId(){return id;}
        VerticeObject& getInfo(){
            if(getType() == FACTOR){
                return factor;
            }
            return rVariable;
        }

        RandomVariable& getInfoRandomVariable(){ return rVariable; }

        RandomVariable const& getInfoRandomVariable() const{ return rVariable; }

        void setInfoRandomVariable(RandomVariable const& info){ rVariable = info; }

        Factor& getInfoFactor(){ return factor; }

        Factor const& getInfoFactor() const{ return factor; }

        void setInfoFactor(Factor const& info){ factor = info; }

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
            oarc << message << var_id << id << type << factor << rVariable;
        }

        void load(graphlab::iarchive& iarc) {
            iarc >> message >> var_id >> id >> type >> factor>> rVariable;
        }

        friend inline bool operator==(const Vertice& v1, const Vertice& v2){
            return v1.vid() == v2.vid();
        }

        friend std::ostream& operator<<(std::ostream& os, const Vertice& vertice){
            if(vertice.getType() == FACTOR){
                os << vertice.var_id << " Factor " << vertice.factor;
            }else{
                os << vertice.var_id << " Variable " << vertice.rVariable;
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
                        return v;
                        //add factor vertex
                    }else{
                        if(tmp == "Variable"){
                            //parse variable and create it
                            RandomVariable var = RandomVariable::parseLine(ss);
                            v = Vertice(var, id);
                            return v;
                            // add variable vertex
                        }
                    }
                }
            }
            throw "[parseline] could not read vertex!!";

        }
};

int Vertice::id = 0;

typedef graphlab::distributed_graph<Vertice, graphlab::empty> GraphLab_type;
typedef std::pair<Vertice*, Vertice*> Edge_type;

//Wrapper to help make the convergence message passing process happen
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
    private:
        double residual;

    public:
        static const double BOUND = 1e-4;
        static const int ITER_BOUND = 10;//TODO USE THIS TO LIMIT ITERATIONS

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
                //std::cout << "gather from factor " << edge.target().data().getMessage() << edge.target().data().getInfoFactor() << std::endl;
                return edge.target().data().getMessage();
            }else{
                if(vertex.data().getType() == Vertice::VARIABLE){
                //std::cout << "gather from variable" << edge.source().data().getMessage() << std::endl;
                    return edge.source().data().getMessage();
                }else{
                    return Factor();
                }
            }
            //return vertex.data().getMessage();
        }

        // Use the total belief of adjacent factors to update its own factor
        void apply(icontext_type& context, vertex_type& vertex,
        const gather_type& total) {
            std::cout << "APPLY" << std::endl;
            vertex.data().iter_count++;
            Factor received = total;
            received.normalize();

            if(vertex.data().getType() == Vertice::FACTOR){
                Vertice fact_vert = vertex.data();
                Factor fact = fact_vert.getInfoFactor();
                if(fact.getScope().size() == 0){
                    throw "READING BLANK FACTOR!";
                }
                Factor new_msg = (fact_vert.getMessage()*received).marginalizeAllBut(fact.getScope());
                //new_msg.normalize();

                std::cout
                << "Factor " << fact << std::endl
                << "    total : " << total << std::endl
                << "    received message " << received.marginalizeAllBut(fact.getScope()) << std::endl
                << "    sent message " << vertex.data().getMessage() << std::endl;

                residual = new_msg.distance(received);
                vertex.data().setMessage(new_msg);

                std::cout
                << "BP: convergence: " << std::endl
                << "Distance(fact, var) = " << residual << std::endl
                << "Bound = " << BeliefProp::BOUND << std::endl
                << "--*--" << std::endl << std::endl;

                if((residual > BeliefProp::BOUND) &&
                    (vertex.data().iter_count < BeliefProp::ITER_BOUND)){
                    context.signal(vertex, residual);
                }
            }else{
                if(vertex.data().getType() == Vertice::VARIABLE){

                    Vertice& var_vert = vertex.data();
                    RandomVariable var = (RandomVariable) var_vert.getInfoRandomVariable();
                    Factor old_msg = vertex.data().getMessage();
                    Factor new_msg = (old_msg*received).marginalizeAllBut(var);
                    //new_msg.normalize();

                    std::cout
                    << "Variable " << var << std::endl
                    << "    total : " << total << std::endl
                    << "    received message " << received.marginalizeAllBut(var) << std::endl
                    << "    sent message " << old_msg << std::endl;

                    residual = new_msg.distance(old_msg);
                    vertex.data().setMessage(new_msg);

                    std::cout
                    << "BP: convergence: " << std::endl
                    << "Distance(fact, var)= " << residual << std::endl
                    << "Bound = " << BeliefProp::BOUND << std::endl
                    << "--*--" << std::endl << std::endl;
                    /*if((residual > BeliefProp::BOUND) &&
                        (vertex.data().iter_count < BeliefProp::ITER_BOUND)){
                        context.signal(vertex, residual);
                    }*/
                }
            }
        }

        //TODO finish bp
        // The scatter edges depend on whether it has converged
        edge_dir_type scatter_edges(icontext_type& context,
                                    const vertex_type& vertex) const {
            if (residual > BeliefProp::BOUND && vertex.data().iter_count < BeliefProp::ITER_BOUND){
                return graphlab::ALL_EDGES;
            }
            else{
                return graphlab::NO_EDGES;
            }
        }

        void scatter(icontext_type& context, const vertex_type& vertex,
                     edge_type& edge) const {
            //bool isAFactor = vertex.data().getType() == Vertice::FACTOR;
            //vertex_type other_vertex = (isAFactor)? edge.target() : edge.source();
            context.signal(edge.source());
            context.signal(edge.target());
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
            if(isAFactor){
                Factor fact = vertex.data().getInfoFactor();
                vertex.data().setMessage(fact.setEvidence(variables, values));
            }else{
                Factor fact = vertex.data().getMessage();
                vertex.data().setMessage(Factor(fact.getScope(), cpd_type(fact.getCPD().size(), 1.)));
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
                        Factor& fact = fact_vert.getInfoFactor();
                        sstream << fact_vert.vid() << " Factor " << fact << std::endl;
                    }else{
                        if(v.data().getType() == Vertice::VARIABLE){
                            Vertice& var_vert = v.data();
                            RandomVariable& var = var_vert.getInfoRandomVariable();
                            sstream << var_vert.vid() << " Variable " << var << std::endl;
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

        static std::pair<int, int> parseEdge(std::stringstream& ss){
            int idSource, idTarget;
            std::string tmp;

            if(ss >> idSource){
                if(ss >> tmp){
                    if (tmp == "->"){
                        ss >> idTarget;
                        return std::pair<int, int>(idSource, idTarget);
                    }
                }
            }return std::pair<int, int>(-1, -1);
        }

        static bool line_parser(GraphLab_type& graph, std::string const& file_name, std::string const& line){
            std::stringstream sstream(line);
            int pos = sstream.tellg(), id;
            std::string tmp;
            sstream >> id;
            if(sstream >> tmp){
                sstream.seekg(pos);
                if(tmp == "->"){
                    std::pair<int, int> edge = parseEdge(sstream);
                    if(edge.first == -1){
                        return false;
                    }
                    graph.add_edge(edge.first, edge.second);
                }else{
                    Vertice v = Vertice::parseLine(sstream);
                    graph.add_vertex(v.vid(), v);
                    //std::cout << "Adding var" << v;
                }
            }
            //std::cout << "Read line";
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
                }//TODO raise exception if it is the wrong kind of edge
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

        void save(std::string dest_file){
            // params: destination, writer, save as .gzip, save vertex, save edges, files per machine
            dgraph.save(dest_file, graph_writer(), false, true, true, 1);
        }

        void save(){
            save("/home/breno/BNETout");
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

        static void absorb_data(Factor& fact, DataSet& data){
            //std::cout << "absorbing data - init" << std::endl;
            data_table_type::iterator element;
            data_object_type const& header = data.getHeader();
            scope_type const& scope = fact.getScope();
            std::vector<int> scope_map(scope.size(),0);
            attribution_type attrib(scope.size(),0);
            //std::cout << "absorbing data - create header map" << std::endl;
            bool missing_var;
            //Map the data header to the factor scope
            for(int attrib_pos = 0; attrib_pos < scope.size(); attrib_pos++){
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
                    std::cerr << std::endl << "Invalid header!! Missing variable: "<< scope[attrib_pos].getName() << std::endl;
                    //TODO raise exception and deal with it
                }
            }
            //std::cout << "absorbing data - iterate through data" << std::endl;
            //Read each tuple and update the factor
            for(element = data.begin(); element != data.end(); element++){
                for(int attrib_pos = 0; attrib_pos < attrib.size(); attrib_pos++){

                    attrib[attrib_pos] = scope[attrib_pos].positionOfValue((*element)[scope_map[attrib_pos]]);
                //    std::cout << "absorbing data - iterate through data - datum ";

                    //TODO raise an error if the value is invalid
                }
                 /*
                std::cout << "atrib: [";
                for(int att_pos = 0; att_pos < scope.size(); att_pos++){
                    std::cout << attrib[att_pos] << ", ";
                }
                std::cout << "\b]" << std::endl;
                 */
                fact[attrib] += 1;
            }
            //std::cout << "absorbing data - normalize" << std::endl;
            fact.normalize();
            //std::cout << "absorbing data - absorbed" << std::endl;
        }

    static void vertex_learn(GraphLab_type::vertex_type& vertex) {
        //std::cout << "Performing learning on a " << std::endl;
        if(vertex.data().getType() == Vertice::FACTOR){
            //std::cout << "factor" <<std::endl;
            std::cout << (Factor&) vertex.data() << std::endl;
            Factor fact = vertex.data().getInfoFactor();
            //std::cout << std::endl << data_set;
            absorb_data(fact, data_set);
            //std::cout  << "Data absorbed" << std::endl;
            vertex.data().setInfoFactor(fact);
            //std::cout << "Learned!! = " << vertex.data().getInfoFactor() << std::endl;
        }else{
            //std::cout << "variable" <<std::endl;
        }
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

        /*
        class graph_writer{
            public:
                std::string save_vertex(GraphLab_type::vertex_type v){
                    std::stringstream sstream;
                    if(v.data().getType() == Vertice::FACTOR){
                        Vertice& fact_vert = v.data();
                        Factor& fact = fact_vert.getInfoFactor();
                        sstream << fact.getScope()[0] << " " << fact << std::endl;
                    }
                    return sstream.str();
                }
                std::string save_edge(GraphLab_type::edge_type e){ return "a -> b"; }
        };
        */

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
                            if(vert.getInfoRandomVariable() == *var){
                                return true;
                            }
                        }
                    }
                    return false;
                }

                static void setBelief(GraphLab_type::vertex_type& vertex){
                    //The reason why this function uses the singleton pointer is described here
                    // http://stackoverflow.com/questions/15841338/c-unresolved-overloaded-function-type
                    //BeliefsSet.setBelief is a non-static-member-function, hence it has a special signature, that cannot be used in updateBeliefs()

                    if(singleton->isTarget(vertex)){
                        Vertice& v = (Vertice&) vertex.data();
                        std::string name = v.getInfoRandomVariable().getName();
                        singleton->belief_map.set(name, v.getMessage());
                    }
                }

            public:
                BeliefsSet(GraphLab_type& g, graphlab::distributed_control& d_ctrl, scope_type const& target):
                           graph(g), target_variables(target), belief_map(d_ctrl){
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
            //std::cout << "Adding var" <<  var << std::endl;
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
            //fgraph.startDGraph();
            fgraph.runInference(i_eng);
            distributed_ctrl.barrier();
            std::cout << "Inference done" << std::endl << "Getting value" << std::endl;
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

        void save(std::string dest_file_name){
            fgraph.save(dest_file_name);
        }

        void save(){
            //fgraph.getDGraph().save("/home/breno/BNETout", graph_writer(), false, true, true, 1);
            fgraph.save("/home/breno/BNETout");
        }

        void load(std::string prefix){
            fgraph.load(prefix);
            distributed_ctrl.barrier();
        }
};

BeliefNetwork::BeliefsSet* BeliefNetwork::BeliefsSet::singleton = NULL;

#endif //FACTOR_GRAPH_NETWORK_UTILS_HPP
