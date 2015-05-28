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

int aux_test_compare_factor_2_string(std::string header, Factor const& fact, std::string string){
    std::stringstream sstream;
    sstream << fact;
    if(string.compare(sstream.str())){
        std::cerr << header <<  std::endl <<
                        "\tExp:" << string << std::endl <<
                        "\tOut:" << sstream.str() << std::endl;
        return 1;
    }
    return 0;
}

int aux_test_compare_var_2_string(std::string header, RandomVariable const& var, std::string string){
    std::stringstream sstream;
    sstream << var;
    if(string.compare(sstream.str())){
        std::cerr << header <<  std::endl <<
                        "\tExp:" << string << std::endl <<
                        "\tOut:" << sstream.str() << std::endl;
        return 1;
    }
    return 0;
}

int aux_test_compare_vertice_2_string(std::string header, Vertice const& vertice, std::string string){
    std::stringstream sstream;
    sstream << vertice;
    if(string.compare(sstream.str())){
        std::cerr << header <<  std::endl <<
                        "\tExp:" << string << std::endl <<
                        "\tOut:" << sstream.str() << std::endl;
        return 1;
    }
    return 0;
}

int test_factor_initialization(){
    Factor rain = new_rain_factor();
    Factor wet  = new_grass_wet_factor();

    //Expected output
    std::string rain_str_exp = "[< ( rain | false true ) > 0.5 0.5 ]";
    std::string wet_str_exp  = "[< ( wet | false true maybe ) ( rain | false true ) > 0.3 0.5 0.1 0.9 0.9 0.1 ]";

    //Comparison
    //If any error print message
    std::string header = "TEST [initialization]:";
    int error_found = 0;

    error_found += aux_test_compare_factor_2_string(header, rain, rain_str_exp);
    error_found += aux_test_compare_factor_2_string(header, wet, wet_str_exp);

    return error_found;
}

int aux_test_attrib(Factor const& fact, cpd_type values){
    int cpd_size = fact.getCPD().size(), flag = 0;
    attribution_type attrib(fact.getScope().size(), 0);
    attribution_type attrib_sizes = fact.getScopeSizes();

    if(attrib_sizes.size() != attrib.size()){
        std::cerr << "TEST[attrib]: scope size does not match attribution size!" << std::endl;
        return -1;
    }

    if(fact.getScope().size() != attrib_sizes.size()){
        std::cerr << "TEST[attrib]: Scope size does not match attribution_sizes size! " <<
                     fact.getScope().size() << " != "<< attrib_sizes.size() << std::endl <<
                     "attrib_sizes = [";
        attribution_type::iterator iter = attrib_sizes.begin();
        if(iter != attrib_sizes.end()){
            std::cerr << *iter; iter++;
        }
        for(; iter != attrib_sizes.end(); iter++){
            std::cerr << ", " << *iter;
        }
        std::cerr << "]" << std::endl;
        return -2;
    }

    bool iterate = true;
    int count = 0;
    //Iterate through attrib checking if the cpd values match the expected
    while(iterate){
        attrib[attrib.size()-1] = 0;
        for(int& val = attrib[attrib.size()-1]; val < attrib_sizes[attrib_sizes.size()-1]; val++){
            if(fact[attrib] != values[count]){
                std::cerr << "TEST[attrib]: fact[";
                attribution_type::iterator i = attrib.begin();
                if (i != attrib.end()){
                    std::cerr << *i; i++;
                }
                for(; i != attrib.end(); i++){
                    std::cerr << ", " << *i;
                }
                std::cerr << "] = " << fact[attrib] << " != cpd[" << count << "] = " << values[count] << std::endl;
                std::cerr << "attrib_sizes = [ ";
                for(attribution_type::iterator att = attrib_sizes.begin(); att != attrib_sizes.end(); att++){
                    std::cerr << *att<<" ";
                }
                std::cerr << "]" << std::endl;
                flag++;
            }
            count++;
        }
        if(attrib.size() < 2){
            iterate = false;
        }
        for(int pos = attrib.size()-2; pos >= 0; pos--){
            attrib[pos]++;
            if(attrib[pos] < attrib_sizes[pos]){
                break;
            }else{
                attrib[pos] = 0;
                iterate = (pos != 0);
            }
        }
    }

    return flag;
}

int test_attribution(){
    Factor rain = new_rain_factor();
    Factor wet  = new_grass_wet_factor();
    Factor cab  = new_cab_factor();

    cpd_type cpd_rain(2, 0.5);
    cpd_type cpd_wet(2*3, 0.0);
    //rain = false,  rain = true
    cpd_wet[0] = 0.3; cpd_wet[1] = 0.5; // wet = false
    cpd_wet[2] = 0.1; cpd_wet[3] = 0.9; // wet = true
    cpd_wet[4] = 0.9; cpd_wet[5] = 0.1; // wet = maybe

    cpd_type cab_cpd(2*2*3, 0.0);
    // (A = false B = false) (A = false B = true) (A = true B = false) (A = true B = true)
    cab_cpd[0] = 0.1;     cab_cpd[1] = 0.2;     cab_cpd[2]  = 0.3;     cab_cpd[3]  = 0.4; //C = low
    cab_cpd[4] = 0.2;     cab_cpd[5] = 0.28;    cab_cpd[6]  = 0.24;    cab_cpd[7]  = 0.28; //C = medium
    cab_cpd[8] = 0.2;     cab_cpd[9] = 0.2;     cab_cpd[10] = 0.2;     cab_cpd[11] = 0.4; //C = high

    int output  = aux_test_attrib(rain, cpd_rain);
        output |= aux_test_attrib(wet, cpd_wet);
        output |= aux_test_attrib(cab, cab_cpd);
    return output;
}

int aux_test_scope_union(Factor const& fact1, Factor const& fact2){
    Factor factor2 = new_grass_wet_factor();
    bool error_found = false;

    scope_type scope_union = Factor::scopeUnion(fact1.getScope(), factor2.getScope());
    for(scope_type::const_iterator fact1_iter = fact1.getScope().begin(); fact1_iter != fact1.getScope().end(); fact1_iter++){
        bool flag_in = false;
        for(scope_type::const_iterator union_iter = scope_union.begin(); union_iter != scope_union.end(); union_iter++){
            flag_in |= (*fact1_iter == *union_iter);
        }
        if(!flag_in){
            std::cerr << "TEST[scope][union]: missing variable from scope union (" << *fact1_iter << ")" << std::endl;
            error_found = true;
        }
    }

    for(scope_type::const_iterator fact2_iter = fact2.getScope().begin(); fact2_iter != fact2.getScope().end(); fact2_iter++){
        bool flag_in = false;
        for(scope_type::const_iterator union_iter = scope_union.begin(); union_iter != scope_union.end(); union_iter++){
            flag_in |= (*fact2_iter == *union_iter);
        }
        if(!flag_in){
            std::cerr << "TEST[scope][union]: missing variable from scope union (" << *fact2_iter << ")" << std::endl;
            error_found = true;
        }
    }

    for(scope_type::const_iterator union_iter = scope_union.begin(); union_iter != scope_union.end(); union_iter++){
        for(scope_type::const_iterator sec_iter = union_iter+1; sec_iter != scope_union.end(); sec_iter++){
            if(*union_iter == *sec_iter){
                std::cerr << "TEST[scope][union]: variable is repeated (" << *union_iter << ")" << std::endl;
                error_found = true;
            }
        }
    }

    return error_found;
}

int test_scope_union(){
    Factor rain = new_rain_factor();
    Factor wet = new_grass_wet_factor();
    return aux_test_scope_union(rain, wet) || aux_test_scope_union(wet, rain);
}

int aux_test_scope_intersection(Factor const& fact1, Factor const& fact2){
    bool error_found = false;
    scope_type intersection = Factor::scopeIntersection(fact1.getScope(), fact2.getScope());
    for(scope_type::const_iterator intersect_iter = intersection.begin();
    intersect_iter != intersection.end(); intersect_iter++){
        if(!Factor::hasVariable(fact1.getScope(), *intersect_iter)){
            std::cerr << "TEST[scope][intersection]: missing variable (" << *intersect_iter << ") from first factor" << std::endl;
            error_found = true;
        }
        if(!Factor::hasVariable(fact2.getScope(), *intersect_iter)){
            std::cerr << "TEST[scope][intersection]: missing variable (" << *intersect_iter << ") from second factor" << std::endl;
            error_found = true;
        }
    }
    return error_found;
}

int test_scope_intersection(){
    Factor rain = new_rain_factor();
    Factor wet = new_grass_wet_factor();
    return aux_test_scope_intersection(rain, wet) || aux_test_scope_intersection(wet, rain);
}

int aux_test_scope_minus(Factor const& fact1, Factor const& fact2){
    bool error_found = false;
    scope_type f1_minus_f2 = Factor::scopeMinus(fact1.getScope(), fact2.getScope());

    int count_vars = 0; //this counter is used to check if there is extra variables not present in fact1
    //The following block checks:
    //  if all variables in fact1 that are in f1_minus_f2 are not in fact2
    //  if all variables in fact1 that are not in f1_minus_f2 are in fact2
    for(scope_type::const_iterator fact1_iter = fact1.getScope().begin();
    fact1_iter != fact1.getScope().end(); fact1_iter++){
        if(Factor::hasVariable(f1_minus_f2, *fact1_iter)){
            if(Factor::hasVariable(fact2.getScope(), *fact1_iter)){
                std::cerr << "TEST[scope][minus]: variable not removed (" << *fact1_iter << ")" << std::endl;
                error_found = true;
            }
            count_vars++;
        }else{
            if(!Factor::hasVariable(fact2.getScope(), *fact1_iter)){
                std::cerr << "TEST[scope][minus]: variable missing (" << *fact1_iter << ")" << std::endl;
                error_found = true;
            }
        }
    }

    //this "if" checks if all variables in f1_minus_f2 are in fact1 and not in fact2
    if(count_vars < f1_minus_f2.size()){
        std::cerr << "TEST[scope][minus]: " << f1_minus_f2.size() - count_vars << " extra variables!"<< std::endl;
    }
    return error_found;
}

int test_scope_minus(){
    Factor rain = new_rain_factor();
    Factor wet = new_grass_wet_factor();
    return aux_test_scope_minus(rain, wet) || aux_test_scope_minus(rain, wet);
}

int aux_test_attribution_map(scope_type const& scope1, scope_type const& scope2,
                             attribution_type const& attrib_map){
    if(scope1.size() < scope2.size()){
        std::cerr << "\tTEST[attribution][map]: Size of second scope is bigger than size of first scope." << std::endl;
        return -1;
    }
    attribution_type new_attrib = Factor::getAttribMap(scope1, scope2);

    int error_found = 0;

    if(attrib_map.size() != new_attrib.size()){
        std::cerr << "TEST[attribution][map]: attribution map of incorrect size, it should be "
        << attrib_map.size() << std::endl;
        return -2;
    }

    for(int pos = 0; pos < attrib_map.size(); pos++){
        if(attrib_map[pos] != new_attrib[pos]){
            std::cerr << "TEST[attribution][map]: attribution map mismatch exp:[";

            for(attribution_type::const_iterator i = attrib_map.begin(); i != attrib_map.end(); i++){
                std::cerr << " " << *i;
            }
            std::cerr << " ] != out:[";

            for(attribution_type::const_iterator i = new_attrib.begin(); i != new_attrib.end(); i++){
                std::cerr << " " << *i;
            }
            std::cerr << " ]" << std::endl;
            //std::cerr << scope1[attrib_map[pos]] << " == " << scope1[new_attrib[pos]];
            error_found = 1;
        }

    }
    return error_found;
}

int test_attribution_map(){
    //Creating Variables
    std::vector<std::string> values_humidity(2, "false");
    values_humidity[1] = "true";
    RandomVariable humidity("humidity", values_humidity);
    RandomVariable rain = new_rain_var();
    RandomVariable wet = new_grass_wet_var();

    //Creating scopes
    scope_type scope0(3, rain);
    scope0[0] = wet;            scope0[1] = humidity;
    scope_type scope1(3, rain);
    scope1[0] = humidity;       scope1[1] = wet;
    scope_type scope2(3, rain);
    scope2[0] = wet;            scope2[2] = humidity;
    scope_type scope3(3, rain);
    scope3[0] = humidity;       scope3[2] = wet;
    scope_type scope4(3, rain);
    scope4[1] = humidity;       scope4[2] = wet;
    scope_type scope5(3, rain);
    scope5[1] = wet;            scope5[2] = humidity;

    scope_type scope6(2, rain);
    scope6[1] = wet;
    scope_type scope7(2, rain);
    scope7[1] = humidity;
    scope_type scope8(2, rain);
    scope8[0] = humidity;
    scope_type scope9(2, rain);
    scope9[0] = humidity; scope9[1] = wet;

    scope_type scope10(1, rain);
    scope_type scope11(1, wet);
    scope_type scope12(1, humidity);

    //Creating attibution maps
    attribution_type map_0_1(3, 0),
                     map_0_2(3, 0),
                     map_0_3(3, 0),
                     map_0_4(3, 0),
                     map_0_5(3, 0),

                     map_0_6(2, 0),
                     map_0_7(2, 0),
                     map_0_8(2, 0),
                     map_0_9(2, 0),

                     map_0_10(1, 2),
                     map_0_11(1, 0),
                     map_0_12(1, 1);

    map_0_1[0] = 1; map_0_1[1] = 0; map_0_1[2] = 2;
    map_0_2[0] = 0; map_0_2[1] = 2; map_0_2[2] = 1;
    map_0_3[0] = 1; map_0_3[1] = 2; map_0_3[2] = 0;
    map_0_4[0] = 2; map_0_4[1] = 1; map_0_4[2] = 0;
    map_0_5[0] = 2; map_0_5[1] = 0; map_0_5[2] = 1;

    map_0_6[0] = 2; map_0_6[1] = 0;
    map_0_7[0] = 2; map_0_7[1] = 1;
    map_0_8[0] = 1; map_0_8[1] = 2;
    map_0_9[0] = 1; map_0_9[1] = 0;

    int error_found = aux_test_attribution_map(scope0, scope1, map_0_1) ||
                      aux_test_attribution_map(scope0, scope2, map_0_2) ||
                      aux_test_attribution_map(scope0, scope3, map_0_3) ||
                      aux_test_attribution_map(scope0, scope4, map_0_4) ||
                      aux_test_attribution_map(scope0, scope5, map_0_5) ||

                      aux_test_attribution_map(scope0, scope6, map_0_6) ||
                      aux_test_attribution_map(scope0, scope7, map_0_7) ||
                      aux_test_attribution_map(scope0, scope8, map_0_8) ||
                      aux_test_attribution_map(scope0, scope9, map_0_9) ||

                      aux_test_attribution_map(scope0, scope10, map_0_10) ||
                      aux_test_attribution_map(scope0, scope11, map_0_11) ||
                      aux_test_attribution_map(scope0, scope12, map_0_12);
    return error_found;
}

int test_normalization(){
    RandomVariable A("A", values_falseTrue()), B("B", values_falseTrue()),
                   C("C", values_lowMediumHigh()), D("D", values_lowMediumHigh());
    scope_type a_scope  (1, A),
               b_scope  (1, B),
               cab_scope(3, C),
               db_scope (2, D);
    //Setting scopes
    cab_scope[1] = A;
    cab_scope[2] = B;
    db_scope[1] = B;

    //Setting cpd's
    cpd_type a_cpd(2, 0.0), b_cpd(2, 0.0), cab_cpd(2*2*3, 0.0), db_cpd(2*3, 0.0);
    //A = false      A = true
    a_cpd[0] = 0.25; a_cpd[1] = 1;

    //B = false   B = true
    b_cpd[0] = 2; b_cpd[1] = 3;

    // (A = false B = false) (A = false B = true) (A = true B = false) (A = true B = true)
    cab_cpd[0] = 1;     cab_cpd[1] = 2;     cab_cpd[2]  = 3;     cab_cpd[3]  = 4; //C = low
    cab_cpd[4] = 5;     cab_cpd[5] = 7;     cab_cpd[6]  = 6;     cab_cpd[7]  = 7; //C = medium
    cab_cpd[8] = 1;     cab_cpd[9] = 1;     cab_cpd[10] = 1;     cab_cpd[11] = 2; //C = high

    //B = false        B = true
    db_cpd[0] = 9;    db_cpd[1] = 1; //D = low
    db_cpd[2] = 8;    db_cpd[3] = 2; //D = medium
    db_cpd[4] = 4;    db_cpd[5] = 6; //D = high

    //Creating factors
    Factor a_fact(a_scope,   a_cpd),
           b_fact(b_scope,   b_cpd),
           cab_fact(cab_scope, cab_cpd),
           db_fact(db_scope,  db_cpd);

    //Normalization
    a_fact.normalize(),
    b_fact.normalize(),
    cab_fact.normalize(),
    db_fact.normalize();

    int error_found = 0;
    std::string header = "TEST[normalization]: ";
    error_found += aux_test_compare_factor_2_string(header, a_fact,  "[< ( A | false true ) > 0.2 0.8 ]");
    error_found += aux_test_compare_factor_2_string(header, b_fact,  "[< ( B | false true ) > 0.4 0.6 ]");
    error_found += aux_test_compare_factor_2_string(header, cab_fact,"[< ( C | low medium high ) ( A | false true ) ( B | false true ) > 0.1 0.2 0.3 0.4 0.2 0.28 0.24 0.28 0.2 0.2 0.2 0.4 ]");
    error_found += aux_test_compare_factor_2_string(header, db_fact, "[< ( D | low medium high ) ( B | false true ) > 0.9 0.1 0.8 0.2 0.4 0.6 ]");

    return error_found;
}

int test_evidenciation(){
    //TODO remove it
    //return 1;

    RandomVariable A("A", values_falseTrue());
    RandomVariable B("B", values_falseTrue());
    RandomVariable C("C", values_lowMediumHigh());
    RandomVariable D("D", values_lowMediumHigh());

    Factor a_fact   = new_a_factor(),
           b_fact   = new_b_factor(),
           cab_fact = new_cab_factor(),
           db_fact  = new_db_factor();

    int error_found = 0;
    std::string header = "TEST[evidence]: ";
    error_found += aux_test_compare_factor_2_string(header, a_fact.setEvidence(A,  "true"),  "[< ( A | false true ) > 0 1 ]");
    error_found += aux_test_compare_factor_2_string(header, b_fact.setEvidence(B,  "false"), "[< ( B | false true ) > 1 0 ]");
    error_found += aux_test_compare_factor_2_string(header, cab_fact.setEvidence(B,"false"), "[< ( C | low medium high ) ( A | false true ) ( B | false true ) > 0.25 0 0.75 0 0.454545 0 0.545455 0 0.5 0 0.5 0 ]");
    error_found += aux_test_compare_factor_2_string(header, db_fact.setEvidence(D, "medium"),"[< ( D | low medium high ) ( B | false true ) > 0 0 0.8 0.2 0 0 ]");

    return error_found;
}

int test_marginalization(){
    RandomVariable A("A", values_falseTrue());
    RandomVariable B("B", values_falseTrue());
    RandomVariable C("C", values_lowMediumHigh());
    RandomVariable D("D", values_lowMediumHigh());

    Factor a_fact   = new_a_factor(),
           cab_fact = new_cab_factor(),
           db_fact  = new_db_factor();

    int error_found = 0;
    std::string header = "TEST[marginalization]: ";
    error_found += aux_test_compare_factor_2_string(header, a_fact.marginalize(A),   "[< > 0 ]");
    error_found += aux_test_compare_factor_2_string(header, a_fact.marginalize(D),   "[< ( A | false true ) > 0.2 0.8 ]");
    error_found += aux_test_compare_factor_2_string(header, cab_fact.marginalize(D), "[< ( C | low medium high ) ( A | false true ) ( B | false true ) > 0.1 0.2 0.3 0.4 0.2 0.28 0.24 0.28 0.2 0.2 0.2 0.4 ]");
    error_found += aux_test_compare_factor_2_string(header, cab_fact.marginalize(B), "[< ( C | low medium high ) ( A | false true ) > 0.3 0.7 0.48 0.52 0.4 0.6 ]");
    error_found += aux_test_compare_factor_2_string(header, db_fact.marginalize(D),  "[< ( B | false true ) > 2.1 0.9 ]");

    return error_found;
}

int test_getUnitary(){
    Factor cab_fact = new_cab_factor();
    std::string header = "TEST[getUnitary]: ";
    return aux_test_compare_factor_2_string(header, cab_fact.getUnitary(), "[< ( C | low medium high ) ( A | false true ) ( B | false true ) > 1 1 1 1 1 1 1 1 1 1 1 1 ]");
}

int aux_test_distance(Factor const& fact1, Factor const& fact2, double ref, double tolerance = 1e-5){
    double fact1_2_fact2 = fact1.distanceTo(fact2),
           fact2_2_fact1 = fact2.distanceTo(fact1);
    int error_found = 0;
    if(std::abs(fact1_2_fact2 - fact2_2_fact1) > tolerance){
        std::cerr << "TEST[distance]: this distance does not respect simetry." << std::endl;
        error_found++;
    }
    if(std::abs(fact1_2_fact2 - ref) > tolerance){
        std::cerr << "TEST[distance]: this distance is not the expected" <<  std::endl <<
                     "\tExp: " << ref << " Out: " << fact1_2_fact2 <<  std::endl;
        error_found++;
    }

    return error_found;
}

int test_distance(){
    Factor a_fact   = new_a_factor(),
           b_fact   = new_b_factor(),
           cab_fact = new_cab_factor(),
           db_fact  = new_db_factor();
    int error_found = 0;
    error_found += aux_test_distance(a_fact,  cab_fact, 1.4144964);
    error_found += aux_test_distance(b_fact,  cab_fact, 1.4322011);
    error_found += aux_test_distance(db_fact, cab_fact, 1.2162237);
    error_found += aux_test_distance(b_fact,  db_fact,  1.7262677);
    return error_found;
}

int test_multiplication(){
    Factor a_fact   = new_a_factor(),
           b_fact   = new_b_factor(),
           cab_fact = new_cab_factor(),
           db_fact  = new_db_factor();
    int error_found = 0;
    std::string header = "TEST[multiplication]: ";
    error_found += aux_test_compare_factor_2_string(header, a_fact*b_fact,    "[< ( A | false true ) ( B | false true ) > 0.08 0.12 0.32 0.48 ]");
    error_found += aux_test_compare_factor_2_string(header, b_fact*a_fact,    "[< ( B | false true ) ( A | false true ) > 0.08 0.32 0.12 0.48 ]");
    error_found += aux_test_compare_factor_2_string(header, a_fact*db_fact,   "[< ( A | false true ) ( D | low medium high ) ( B | false true ) > 0.18 0.02 0.16 0.04 0.08 0.12 0.72 0.08 0.64 0.16 0.32 0.48 ]");
    error_found += aux_test_compare_factor_2_string(header, b_fact*db_fact,   "[< ( B | false true ) ( D | low medium high ) > 0.36 0.32 0.16 0.06 0.12 0.36 ]");
    error_found += aux_test_compare_factor_2_string(header, a_fact*cab_fact,  "[< ( A | false true ) ( C | low medium high ) ( B | false true ) > 0.02 0.04 0.04 0.056 0.04 0.04 0.24 0.32 0.192 0.224 0.16 0.32 ]");

    return error_found;
}

int aux_test_parseLine_factor(std::string string){
    std::stringstream ss(string);
    Factor fact = Factor::parseLine(ss);
    std::string header = "TEST[parseLine]: ";
    return aux_test_compare_factor_2_string(header, fact, string);
}

int test_parseLine_factor(){
    std::string f1 = "[< ( A | false true ) ( B | false true ) > 1 1 1 1 ]",
                f2 = "[< ( C | low medium high ) ( A | false true ) ( B | false true ) > 0.1 0.2 0.3 0.4 0.2 0.28 0.24 0.28 0.2 0.2 0.2 0.4 ]",
                f3 = "[< ( D | low medium high ) ( B | false true ) > 0.9 0.1 0.8 0.2 0.4 0.6 ]";

    int error_found = 0;
    error_found += aux_test_parseLine_factor(f1);
    error_found += aux_test_parseLine_factor(f2);
    error_found += aux_test_parseLine_factor(f3);
    return error_found;
}

//%%%%%%%%%%%%%%%%%% RandomVariable tests %%%%%%%%%%%%%%%%%%%%%%

int aux_test_parseLine_var(std::string string){
    std::stringstream ss(string);
    RandomVariable var = RandomVariable::parseLine(ss);
    std::string header = "TEST[parseLine]: ";
    return aux_test_compare_var_2_string(header, var, string);
}

int test_parseLine_var(){
    std::string v1 = "( A | false true )",
                v2 = "( B | false true )",
                v3 = "( C | low medium high )",
                v4 = "( D | low medium high )";

    int error_found = 0;
    error_found += aux_test_parseLine_var(v1);
    error_found += aux_test_parseLine_var(v2);
    error_found += aux_test_parseLine_var(v3);
    error_found += aux_test_parseLine_var(v4);
    return error_found;
}

//%%%%%%%%%%%%%%%%%%%%%% Vertice tests %%%%%%%%%%%%%%%%%%%%%%%%%

int aux_test_parseLine_vertice(std::string string){
    std::stringstream ss(string);
    Vertice vertice = Vertice::parseLine(ss);
    std::string header = "TEST[parseLine]: ";
    return aux_test_compare_vertice_2_string(header, vertice, string);
}

int test_parseLine_vertice(){
    std::string v1 = "1 Factor [< ( rain | false true ) > 0.350746 0.649254 ]",
                v2 = "2 Variable ( rain | false true )",
                v3 = "3 Factor [< ( wet | false true maybe ) ( rain | false true ) > 0.642857 0.346154 0.357143 0.653846 0 0 ]",
                v4 = "4 Variable ( wet | false true maybe )";

    int error_found = 0;
    error_found += aux_test_parseLine_vertice(v1);
    error_found += aux_test_parseLine_vertice(v2);
    error_found += aux_test_parseLine_vertice(v3);
    error_found += aux_test_parseLine_vertice(v4);
    return error_found;
}

int test_vertice_initialization(){
    RandomVariable rain = new_rain_var();
    Factor fact = new_rain_factor();
    Vertice a(rain), b(fact);

    int error_found = 0;
    std::string header = "TEST[initialization]: ";
    error_found += aux_test_compare_vertice_2_string(header, a, "0 Variable ( rain | false true )");
    error_found += aux_test_compare_vertice_2_string(header, b, "1 Factor [< ( rain | false true ) > 0.5 0.5 ]");
    return error_found;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%Graph%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int graph_initialization(graphlab::distributed_control& dc){
    //build some basic graphs and find a way to print it
    return 1;
}

int graph_load(graphlab::distributed_control& dc){
    //load a graph from file
    // save it to another one
    //compare the two disconsidering order
    return 1;
}

int graph_save(graphlab::distributed_control& dc){
    return 1;
}

int graph_inference(graphlab::distributed_control& dc){
    //build a network
    //run inference
    // compare results
    return 1;
}

//%%%%%%%%%%%%%%%%%%%%%%%%Belief Network%%%%%%%%%%%%%%%%%%%%%%%

int bnet_initialization(graphlab::distributed_control& dc){
    //build some basic networks and find a way to print it
    return 1;
}

int bnet_load(graphlab::distributed_control& dc){
    //load a network from file
    // save it to another file
    //compare the two disconsidering order
    std::cout << "creating bn"<< std::endl;
    BeliefProp i_eng;
    BeliefNetwork bn(dc, &i_eng);
    //TODO Remove address constants
    bn.load("/home/breno/Documents/git-repos/graphlab/apps/bnets/test_files/bnet-test.txt");
    bn.save("/home/breno/Documents/git-repos/graphlab/apps/bnets/test_files/bnet-test-out.txt");
    //SEEMS OK
    //return compare_files_ special("/home/breno/Documents/git-repos/graphlab/apps/bnets/test_files/test.csv",\
    //                              "/home/breno/Documents/git-repos/graphlab/apps/bnets/test_files/bnet-test-out.txt");
    //std::cout << "saving done"<< std::endl;
    return 0;
}

int bnet_save(graphlab::distributed_control& dc){
    return 1;
}

/*
 * Tests the Bnet creation and learning with complete data
 */
int bnet_learning1(graphlab::distributed_control& dc){
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
    cpd_type values(2, 0.01);
    Factor f1(scope, values);

    std::string grass_str = std::string("wet");
    RandomVariable grass(grass_str, triple);

    scope_type scope_rain_grass = scope_type(2, rain);
    scope_rain_grass[0] = grass;
    cpd_type values_rain_grass(6, 0.01);
    Factor f2(scope_rain_grass, values_rain_grass);

    std::cout << "variables created"<< std::endl;

    BNet_vertice r = bn.addVariable(rain, f1);
    BNet_vertice g = bn.addVariable(grass, f2);

    std::cout << "variables inserted"<< std::endl;
    bn.addEdge(r,g);
    std::cout << "edge created"<< std::endl;

    DataSet data("/home/breno/Documents/git-repos/graphlab/apps/bnets/test_files/test.csv");
    std::cout << "DATASET read" << std::endl;
    bn.learnParameters(data);
    std::cout << "learning done"<< std::endl;
    bn.save();
    std::cout << "saving done"<< std::endl;
    //std::cout << "P(wet) = " << bn.infer(grass) << std::endl;
    return 0;
}

int main(int argc, char** argv){
    graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;

    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%TEST UNIT%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%Factor%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    dc.cout() << "TEST[initialization]      " << (test_factor_initialization()? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[attrib]              " << (test_attribution()          ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[scope][union]        " << (test_scope_union()          ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[scope][intersection] " << (test_scope_intersection()   ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[scope][minus]        " << (test_scope_minus()          ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[attribution][map]    " << (test_attribution_map()      ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[normalization]       " << (test_normalization()        ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[evidence]            " << (test_evidenciation()        ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[marginalization]     " << (test_marginalization()      ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[getUnitary]          " << (test_getUnitary()           ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[distance]            " << (test_distance()             ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[multiplication]      " << (test_multiplication()       ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[parseLine]           " << (test_parseLine_factor()     ? "Error" : "OK") << std::endl;
    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl << std::endl;

    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%TEST UNIT%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%Random Variable%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    dc.cout() << "TEST[parseLine]           " << (test_parseLine_var()        ? "Error" : "OK") << std::endl;
    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl << std::endl;

    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%TEST UNIT%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%%Vertice%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    dc.cout() << "TEST[parseLine]           " << (test_parseLine_vertice()    ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[initialization]      " <<(test_vertice_initialization()? "Error" : "OK") << std::endl;
    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl << std::endl;
/*
    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%TEST UNIT%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%Graph%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    dc.cout() << "TEST[initialization]      " << (graph_initialization(dc)    ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[load]                " << (graph_load(dc)              ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[save]                " << (graph_save(dc)              ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[inference]           " << (graph_inference(dc)         ? "Error" : "OK") << std::endl;
    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl << std::endl;

    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%TEST UNIT%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%Belief Network%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    dc.cout() << "TEST[initialization]      " << (bnet_initialization(dc)     ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[load]                " << (bnet_load(dc)               ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[save]                " << (bnet_save(dc)               ? "Error" : "OK") << std::endl;
    dc.cout() << "TEST[learning comp data]  " << (bnet_learning1(dc)          ? "Error" : "OK") << std::endl;
    dc.cout() << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl << std::endl;

*/
    dc.cout() << "closing mpi" << std::endl;
    graphlab::mpi_tools::finalize();
    dc.cout() << "mpi closed" << std::endl;
}
