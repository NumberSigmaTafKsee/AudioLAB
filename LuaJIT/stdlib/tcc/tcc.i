%module tcc
%{
extern "C" {
#include "libtcc.h"    
}
#include <cassert>
#include <iostream>
#include <vector>

#define ASSERT_TCC() assert(state != NULL)
%}

typedef struct TCCState TCCState;
typedef void (*TCCErrorFunc)(void *opaque, const char *msg);


TCCState *tcc_new(void);
void tcc_delete(TCCState *s);
void tcc_set_lib_path(TCCState *s, const char *path);
void tcc_set_error_func(TCCState *s, void *error_opaque, TCCErrorFunc error_func);
TCCErrorFunc tcc_get_error_func(TCCState *s);
void *tcc_get_error_opaque(TCCState *s);
void tcc_set_options(TCCState *s, const char *str);
int tcc_add_include_path(TCCState *s, const char *pathname);
int tcc_add_sysinclude_path(TCCState *s, const char *pathname);
void tcc_define_symbol(TCCState *s, const char *sym, const char *value);
void tcc_undefine_symbol(TCCState *s, const char *sym);
int tcc_add_file(TCCState *s, const char *filename);
int tcc_compile_string(TCCState *s, const char *buf);
int tcc_set_output_type(TCCState *s, int output_type);
int tcc_add_library_path(TCCState *s, const char *pathname);
int tcc_add_library(TCCState *s, const char *libraryname);
int tcc_add_symbol(TCCState *s, const char *name, const void *val);
int tcc_output_file(TCCState *s, const char *filename);
int tcc_run(TCCState *s, int argc, char **argv);
int tcc_relocate(TCCState *s1, void *ptr);
void *tcc_get_symbol(TCCState *s, const char *name);
//void tcc_list_symbols(TCCState *s, void *ctx, void (*symbol_cb)(void *ctx, const char *name, const void *val));


%inline %{

    static void _errorcb(void * opaque, const char * msg) {
        std::cout << "Error: " << msg << std::endl;        
    }

    struct Symbol
    {
        union {
            double d;
            char*  str;
            bool   b;
            void   *ptr;
        };
        bool is_string;
        

        Symbol(double value) {
            d = value;
            is_string = false;
        }
        Symbol(const char * s) {
            str = strdup(s);
            is_string = true;
        }
        Symbol(bool val) {
            b = val;
            is_string = false;
        }        
    };
    struct TinyCC
    {
        TCCState * state;

        TinyCC() {
            state = tcc_new();
            assert(state != NULL);
            tcc_set_error_func(state,NULL,_errorcb);
        }
        ~TinyCC() {
            if(state) tcc_delete(state);
        }

        void SetLibPath(const char * path) { ASSERT_TCC(); tcc_set_lib_path(state,path); }
        void SetOptions(const char * str)  { ASSERT_TCC(); tcc_set_options(state,str);  }
        int AddIncludePath(const char * path) { ASSERT_TCC(); return tcc_add_include_path(state,path); }
        int AddSysIncludePath(const char * path) { ASSERT_TCC(); return tcc_add_sysinclude_path(state,path); }
        void DefineSymbol(const char * symbol, const char * value) { ASSERT_TCC(); tcc_define_symbol(state, symbol, value); }
        
        void UndefineSymbol(const char *symbol) { ASSERT_TCC(); tcc_undefine_symbol(state, symbol); }

        int AddFile(const char * filename) { ASSERT_TCC(); return tcc_add_file(state, filename); }
        int CompileString(const char * buf) { ASSERT_TCC(); return tcc_compile_string(state,buf); }
        int SetOutputType(int type) { ASSERT_TCC(); return tcc_set_output_type(state,type); }
        int AddLibraryPath(const char * path) { ASSERT_TCC(); return tcc_add_library_path(state,path); }
        int AddLibrary(const char * name) { ASSERT_TCC(); return tcc_add_library(state, name); }
        int AddSymbol(const char * symbol, Symbol & sym) { ASSERT_TCC(); return tcc_add_symbol(state,symbol,sym.ptr); }
        int OutputFile(const char * file) { ASSERT_TCC(); return tcc_output_file(state,file); }
        int Run(int argc, std::vector<const char*> & args) {
            ASSERT_TCC();
            return tcc_run(state, argc, (char**)args.data() );
        }
     };
%}