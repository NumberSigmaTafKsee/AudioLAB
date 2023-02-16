%module dir
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

%inline %{
    struct Dir
    {
        DIR * dir;
        struct dirent * dent;
        size_t cnt;

        Dir(const char * path) {
            dir = opendir(path);        
            get_count();
        }
        ~Dir() {
            if(dir) closedir(dir);
        }

        long tell() const {
            return telldir(dir);
        }
        void get_count() {
            cnt = 0;
            while(dent = readdir(dir)) {            
                if( !strcmp(dent->d_name,".") || !strcmp(dent->d_name,"..")) continue;
                cnt++;
            }
            rewinddir(dir);
        }
        int  count() const { return cnt; }
        const char* entry() {
            dent = readdir(dir);
            return dent->d_name;
        }

        const char* getcwd() const {
            return getcwd();
        }
        void chdir(const char * path) {
            chdir(path);
        }
        void mkdir(const char * path) {
            mkdir(path);
        }
        void rename(const char * path, const char * new_path) {
            rename(path,new_path);
        }
        void system(const char * command) {
            system(command);
        }
    };
%}