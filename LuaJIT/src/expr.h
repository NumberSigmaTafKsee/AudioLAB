#ifndef __EXPR_H
#define __EXPR_H 

#include <cstdint>
#include <ctype.h>
using namespace std;
#include "expr/expr.h"

using namespace std;
typedef struct expr;
typedef struct ExprList;

struct Expr
{
    struct expr *e;
    ExprList * vars;
    
    Expr(const char *eq);
    ~Expr(); 

    double eval();

    void set_var(const char * symbol, double val);
    double get_var(const char * symbol);
    
};

#endif
