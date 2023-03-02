#ifndef __EXPR_H
#define __EXPR_H 


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