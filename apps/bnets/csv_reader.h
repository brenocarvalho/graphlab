/*
*  \author: Breno Carvalho
*/

#ifndef CSV_READER_HPP
#define CSV_READER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <cstring>
#include <iterator>
#include <algorithm>
#include <iomanip>

class DataSet{
public:
    std::list<std::string>* header;
    //TODO print the data nicelly... std::list<int>* max_token_sizes; //os << std::setw(n) << values;
    std::list< std::list<std::string> > *data;
    std::string sep;
    int row_size ;

    DataSet(){
        header = new std::list<std::string>();
        data = new std::list< std::list<std::string> >();
        sep = "";
        row_size = -1;
    }

    DataSet(std::string source_name, std::string separator = std::string(","), bool hasTitle = true){
        std::ifstream source;
        source.open(source_name.c_str(), std::ifstream::in);
        header = new std::list<std::string>();
        data = new std::list< std::list<std::string> >();
        sep = separator;
        row_size = -1;
        if(hasTitle){
            readLine(source, sep, header);
            row_size = header->size();
        }
        std::list<std::string>* line = new std::list<std::string>();
        while(readLine(source, sep, line)){
            if(row_size < 0){
                row_size = line->size();
            }else{
                if(line->size() != row_size){
                    //TODO Exception
                    std::cout << "error: problem in the csv file";
                }
            }
            data->push_back(*line);
            line = new std::list<std::string>();
        }
        delete line;
        source.close();
    }

    DataSet head(int max_size = 10){
        DataSet out;
        out.row_size = row_size;
        out.sep = sep; //TODO copy it
        std::list< std::list<std::string> >::const_iterator i_rows= data->begin();
        std::list<std::string>::const_iterator i_cells;
        std::list<std::string>* tmp;

        tmp = new std::list<std::string>();
        for (i_cells = header->begin(); i_cells != header->end(); ++i_cells) {
            tmp->push_back(*i_cells);
        }
        delete out.header;
        out.header = tmp;

        for (i_rows = data->begin(); i_rows != data->end(); ++i_rows) {
            tmp = new std::list<std::string>();
            for (i_cells = (*i_rows).begin(); i_cells != (*i_rows).end(); ++i_cells) {
                tmp->push_back(*i_cells);
            }
            out.data->push_back(*tmp);
        }

        return out;
    }

    std::list< std::list<std::string> >::const_iterator begin(){
        return data->begin();
    }

    std::list< std::list<std::string> >::const_iterator end(){
        return data->end();
    }

private:
    static inline std::string &trim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    //Method taken from http://stackoverflow.com/questions/599989/is-there-a-built-in-way-to-split-strings-in-c
    void split(const std::string& str, const std::string& delimiters , std::vector<std::string>& tokens)
    {
        // Skip delimiters at beginning.
        std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
        // Find first "non-delimiter".
        std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
        std::string token;
        while (std::string::npos != pos || std::string::npos != lastPos)
            {
                token = str.substr(lastPos, pos - lastPos);
                tokens.push_back(trim(token));
                lastPos = str.find_first_not_of(delimiters, pos);
                pos = str.find_first_of(delimiters, lastPos);
            }
        }

    bool readLine(std::ifstream& file, std::string& separator, std::list<std::string>* elems){
        std::string line;//, item;
        bool out = std::getline(file, line);
        if(line.size() <= 0){
            return false;
        }
        std::string const line_const(line);
        std::vector<std::string> tokens;
        split(line_const, separator, tokens);
        for(int i = 0; i < tokens.size(); i++){
            elems->push_back(tokens[i]);
            //std::cout << tokens[i] << std::endl;
        }
        return out;
    }

public:
    friend std::ostream& operator<<(std::ostream& os, const DataSet& d);

};

std::ostream& operator<<(std::ostream& os, const std::list<std::string>& l){
    std::list<std::string>::const_iterator iterator = l.begin();
    std::list<std::string>::const_iterator end = l.end();
    if( iterator != end){
        os << "'"<< *iterator <<"'"; iterator++;
        for (; iterator != end; iterator++) {
            os << ", '" << *iterator << "'";
        }
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const DataSet& d){
    os << "DataSet(" << d.data->size() << ", " << d.row_size << ")"<< std::endl;
    os << *(d.header) << std::endl << std::endl;
    std::list< std::list<std::string> >::const_iterator iterator= d.data->begin();
    for (iterator; iterator != d.data->end(); ++iterator) {
        os << *iterator << std::endl;
    }
    return os;
}

//int main(int argc, char** argv){
//    DataSet data("test.csv");
//    std::cout << data.head();
//}

#endif //CSV_READER_HPP
