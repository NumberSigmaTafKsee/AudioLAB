%module Glob
%{
#include <dirent.h>
#include <sys/dir.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <string>
#include <iostream>
#include <glob.h>
%}
%include "std_string.i"
%include "std_vector.i"

%template(string_vector) std::vector<std::string>;

%inline %{
    struct Glob
    {
        std::vector<std::string> globs;

        Glob(const char * path) {
            glob_t globbuf;
            globbuf.gl_offs = 2;
            glob(path,GLOB_TILDE,NULL,&globbuf);            
            for(size_t i = 0; i < globbuf.gl_pathc; i++)
                if(globbuf.gl_pathv[i])                                
                    globs.push_back(std::string(globbuf.gl_pathv[i]));
            globfree(&globbuf);
        }

        void rescan(const char * path) {
            glob_t globbuf;
            globs.clear();
            globbuf.gl_offs = 2;
            glob(path,GLOB_TILDE,NULL,&globbuf);            
            for(size_t i = 0; i < globbuf.gl_pathc; i++)
                if(globbuf.gl_pathv[i])                                
                    globs.push_back(std::string(globbuf.gl_pathv[i]));
            globfree(&globbuf);
        }

        std::string operator[](size_t i) {
            return globs[i];
        }
        size_t size() const { return globs.size(); }

        std::string __getitem__(size_t i) { return globs[i]; }

    };
%}