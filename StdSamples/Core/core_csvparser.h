#ifndef __CSVPARSER_H
#define __CSVPARSER_H 

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "parser.hpp"

struct CSVParser
    {
        std::vector<std::vector<std::string>> m;

        CSVParser() = default;

        CSVParser(std::vector<std::vector<std::string>> & r) {
            m = r;
        }
        CSVParser(const std::string & file)
        {            
            std::ifstream f;
            size_t r = 0;
            f.open(file);
            if(f.is_open())
            {
                size_t row = 0;
                size_t col = 0;
                csv::CsvParser parser(f);
                while(!parser.empty())
                {                
                    std::vector<std::string> row;                    
                    while(1)
                    {
                        csv::Field field = parser.next_field();
                        if(field.type == csv::FieldType::ROW_END)
                        {
                            m.push_back(row);                            
                        }
                        else
                        {
                            row.push_back(*field.data);
                        }
                        if(field.type == csv::FieldType::ROW_END || parser.empty()) break;
                    }
                    
                }
            }            
        }

        size_t num_rows() { return m.size(); }
        size_t num_cols(size_t r) { return m[r].size(); }
        
        std::string operator()(size_t i, size_t j) { return m[i][j];  }        
        std::string get_col(size_t r, size_t c) { return m[r][c]; }
        
        void write_file(std::string& file)
        {
            FILE * f = fopen(file.c_str(),"w");
            for(size_t i = 0; i < m.size(); i++)
            {
                for(size_t j = 0; j < m[i].size()-1; j++)
                {
                    fprintf(f,"%s,",m[i][j].c_str());
                }
                fprintf(f,"\n");
            }
            fclose(f);
        }

        void print_row(size_t r)
        {
            for(size_t i = 0; i < m[r].size(); i++)
            {
                std::cout << m[r][i];
                if(i != m[r].size()-1) std::cout << ",";
            }
            std::cout << std::endl;
        }
    };

#endif