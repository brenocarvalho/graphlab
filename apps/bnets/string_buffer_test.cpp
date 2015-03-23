#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <list>
#include <vector>
#include <iterator>

void parseVariable(std::stringstream& ss);

void parseFactor(std::stringstream& ss){
    std::string name, tmp;
    std::list<std::string> v_list;
    if( ss >> tmp && tmp == "[<"){
        int pos = ss.tellg();
        while(ss >> tmp && tmp != ">"){
            ss.seekg(pos);
            parseVariable(ss);
            pos = ss.tellg();
        }

        while(ss >> tmp && tmp != "]"){
            std::cout << " " << std::atof(tmp.c_str());
        }

    }
    //std::vector<std::string> values(v_list.size());
    //std::copy(v_list.begin(), v_list.end(), values.begin());

}

void parseVariable(std::stringstream& ss){
    std::string name, tmp;
    std::list<std::string> v_list;
    if( ss >> tmp && tmp == "("){
        ss >> name;
        if(ss >> tmp){
            if (tmp == "|"){
                while(ss >> tmp && tmp != ")"){
                    v_list.push_back(tmp);
                }
            }
        }
    }
    std::vector<std::string> values(v_list.size());
    std::copy(v_list.begin(), v_list.end(), values.begin());

    std::cout << "( " << name;
    std::vector<std::string>::iterator i = values.begin();
    if( i != values.end()){
        std::cout << " | " << *i;
        i++;
    }
    for(; i != values.end(); i++){
        std::cout << " " << *i;
    }
    std::cout << " )";
}

void parseVertice(std::stringstream& ss){
    int id;
    std::string tmp;
    if(ss >> id){
        if(ss >> tmp){
            if (tmp == "Factor"){
                std::cout << id << " Factor ";
                //parse factor and create factor
                parseFactor(ss);
                //add factor vertex
                std::cout << std::endl;
            }else{
                if(tmp == "Variable"){
                    std::cout << id << " Variable ";
                    //parse variable and create it
                    parseVariable(ss);
                    // add variable vertex
                    std::cout << std::endl;
                }
            }
        }
    }else{
        //todo error
    }

}

void parseEdge(std::stringstream& ss){
    int idSource, idTarget;
    std::string tmp;

    if(ss >> idSource){
        if(ss >> tmp){
            if (tmp == "->"){
                ss >> idTarget;
                std::cout << idSource << " -> " << idTarget << std::endl;
            }
        }
    }else{
        //todo error
    }
}

void parseLine(std::string const& st){
    std::stringstream sstream(st);
    int pos = sstream.tellg(), id;
    std::string tmp;
    sstream >> id;
    if(sstream >> tmp){
        sstream.seekg(pos);
        if(tmp == "->"){
            parseEdge(sstream);
        }else{
            parseVertice(sstream);
        }
    }
}


int main(){
    std::string s1 = "11 Factor [< ( HREKG | Low Normal High ) ( HR | Low Normal High ) ( ErrCauter | True False ) > 0.333333 0.98 0.333333 0.01 0.333333 0.01 0.333333 0.01 0.333333 0.98 0.333333 0.01 0.333333 0.01 0.333333 0.01 0.333333 0.98]",
                s2 = "00 Variable ( PAP | Low Normal High )",
                s3 = "01 -> 02";

    parseLine(s1);
    parseLine(s2);
    parseLine(s3);
}
