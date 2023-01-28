#include "expr.h"
#include "expr/expr.h"
#include <cassert>

struct ExprList
{
    struct expr_var_list v;    

    ExprList() { v = {0}; }
};

Expr::Expr(const char * eq) {         

        vars = new ExprList();
        e = expr_create(eq,strlen(eq), &vars->v,NULL);
        assert(e != NULL);                
    }

Expr::~Expr() { 
    if(e) expr_destroy(e,&vars->v); 
    delete vars;
}
double Expr::eval() { return expr_eval(e); }

void Expr::set_var(const char * symbol, double val)
{
    struct expr_var * var = expr_var(&vars->v, symbol, strlen(symbol));
    var->value = val;
}
double Expr::get_var(const char * symbol)
{
    struct expr_var * var = expr_var(&vars->v, symbol, strlen(symbol));
    return var->value;
}
