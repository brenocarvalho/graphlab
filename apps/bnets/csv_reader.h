/*
*  \author: Breno Carvalho
*/

#ifndef CSV_READER_HPP
#define CSV_READER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>
#include <cstring>
#include <iterator>
#include <algorithm>
#include <iomanip>

#include <graphlab.hpp>

typedef std::vector<std::string> data_object_type;
typedef std::vector< data_object_type > data_table_type;

class DataSet{
private:
    data_object_type header;
    //TODO print the data nicelly... std::list<int>* max_token_sizes; //os << std::setw(n) << values;
    data_table_type data;
    std::string sep;
    int row_size;

    static inline std::string &trim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    //Method taken from http://stackoverflow.com/questions/599989/is-there-a-built-in-way-to-split-strings-in-c
    void split(const std::string& str, const std::string& delimiters, std::vector<std::string>& tokens){
        // Skip delimiters at beginning.
        std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
        // Find first "non-delimiter".
        std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
        std::string token;
        while (std::string::npos != pos || std::string::npos != lastPos){
            token = str.substr(lastPos, pos - lastPos);
            tokens.push_back(trim(token));
            lastPos = str.find_first_not_of(delimiters, pos);
            pos = str.find_first_of(delimiters, lastPos);
        }
    }

    bool readLine(std::ifstream& file, std::string const& separator, data_object_type& elems){
        std::string line;//, item;
        bool out = std::getline(file, line);
        if(line.size() == 0){
            return false;
        }
        //std::string const line_const(line);
        std::vector<std::string> tokens;
        split(line.c_str(), separator, tokens);
        for(int i = 0; i < tokens.size(); i++){
            elems.push_back(tokens[i]);
            //std::cout << "adding: " << tokens[i] << elems.size() << std::endl;
        }
        return out;
    }

public:
    DataSet(){
        header = data_object_type();
        data = data_table_type();
        sep = "";
        row_size = -1;
    }

    data_object_type const& getHeader(){return header;}
    data_table_type  const& getDataTable(){return data;}

    int getRowSize(){return row_size;}
    int getLinesQtd(){return data.size();}
    std::string getSeparator(){return sep;}

    DataSet(std::string source_name, std::string separator = std::string(","), bool hasTitle = true){
        std::ifstream source;
        source.open(source_name.c_str(), std::ifstream::in);
        header = data_object_type();
        data = data_table_type();
        sep = separator;
        row_size = -1;
        if(hasTitle){
            readLine(source, sep, header);
            row_size = header.size();
            //std::cout << "Header size: " << header.size() << std::endl;
        }
        data_object_type line = data_object_type();
        while(readLine(source, sep, line)){
            if(line.size() != 0){
                if(row_size < 0){
                    row_size = line.size();
                }else{
                    if(line.size() != row_size){
                        //TODO Exception
                        std::cerr << "error: problem in the csv file: header.size(" << header.size() << ") != line.size(" << line.size() << ")" << std::endl;
                        std::cerr << "Header:"<< std::endl;
                        for(data_object_type::const_iterator i = header.begin(); i != header.end(); i++){
                            std::cerr << *i << " ";
                        }
                        std::cerr << "Line:"<< std::endl;
                        for(data_object_type::const_iterator i = line.begin(); i != line.end(); i++){
                            std::cerr << *i << ",";
                        }
                        std::cerr << std::endl;
                        throw "problem in the csv file";
                    }
                }
                data.push_back(line);
                line = data_object_type();
            }
        }
        source.close();
    }

    DataSet head(int max_size = 10){
        DataSet out;
        out.row_size = row_size;
        out.sep = sep;
        data_table_type::const_iterator i_rows= data.begin();
        data_object_type::const_iterator i_cells;
        data_object_type tmp;

        tmp = data_object_type();
        for (i_cells = header.begin(); i_cells != header.end(); ++i_cells) {
            tmp.push_back(*i_cells);
        }
        out.header = tmp;
        int to_read = (max_size > 0)? max_size : 0;
        for (i_rows = data.begin(); (to_read--) && i_rows != data.end(); ++i_rows) {
            tmp = data_object_type();
            for (i_cells = (*i_rows).begin(); i_cells != (*i_rows).end(); ++i_cells) {
                tmp.push_back(*i_cells);
            }
            out.data.push_back(tmp);
        }

        return out;
    }

    data_table_type::iterator begin(){
        return data.begin();
    }

    data_table_type::iterator end(){
        return data.end();
    }

    friend std::ostream& operator<<(std::ostream& os, const DataSet& d);

    void save(graphlab::oarchive& oarc) const {
        oarc << header << data << sep << row_size;
    }
    void load(graphlab::iarchive& iarc) {
        iarc >> header >> data >> sep >> row_size;
    }

};

std::ostream& operator<<(std::ostream& os, const data_object_type& l){
    data_object_type::const_iterator iterator = l.begin();
    data_object_type::const_iterator end = l.end();
    if( iterator != end){
        os << "'"<< *iterator <<"'"; iterator++;
        for (; iterator != end; iterator++) {
            os << ", '" << *iterator << "'";
        }
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const DataSet& d){
    os << "DataSet(" << d.data.size() << ", " << d.row_size << ")"<< std::endl;
    os << d.header << std::endl << std::endl;
    data_table_type::const_iterator iterator= d.data.begin();
    for (iterator; iterator != d.data.end(); ++iterator) {
        os << *iterator << std::endl;
    }
    return os;
}

//int main(int argc, char** argv){
//    DataSet data("test.csv");
//    std::cout << data.head();
//}

#endif //CSV_READER_HPP
