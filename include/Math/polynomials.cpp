
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <ostream>
#define eps 0.000001

using namespace std;

struct monom {  // monomial
    double coef; // coefficient
    int pow; // power
    monom * next;
    monom(double c, int p, monom * n = NULL): pow(p), coef(c), next(n) {};  // constructor
};
class polynomial {
public:
    polynomial(): head(NULL) {};
    polynomial(double c, int p) { head = new monom (c, p);};
    polynomial(const polynomial & p) { head = copy(p.head);};
    ~polynomial() { if (head != NULL) del_list(head);};
    
    polynomial operator +(const polynomial & p);  //
    polynomial & operator +=(const polynomial & p);
    polynomial & operator *(double a);
    polynomial & operator *=(double a);
    polynomial & operator *=(const polynomial p);
    polynomial operator *(const polynomial & p);
    polynomial & operator =(const polynomial & p);
    double operator ()(double a);  // polynomial value
    friend ostream & operator <<(ostream & out, const polynomial & p);  // write polynomial in stream
    double value(double a);  // value in point a
    void change(double a);// changing x0 in (x-x0)
private:
    monom * head;   // head of list
    void del_next(monom * m); // delete monom after h
    void add_after(double c, int p, monom * m);
    monom * del_head(monom * m); // delete head and return new head
    monom * copy(monom * m);
    monom * del_list(monom * m);  // deleting list
   
    monom * add(monom * m1, monom * m2);   // addition
    monom * mult_num(monom * m, double a);   // multiplying with number
    monom * mult_pol(monom * m1, monom * m2);   // multiplying with polynomial
};
int null(double a) { return (abs(a) < eps)?(1):(0);}
//
struct net {  // network
    net(double l, double r, int temp);  // l - left, r - right, temp - number of points
    double * x = NULL;
    double * y = NULL;
    int n; // number of points
    double func_value(double a);
};

net::net(double l, double r, int temp) {
    n = temp;
    x = new double [n];
    y = new double [n];
    double delta = (l + r)/(n - 1);  // delta of net
    for (int i = 0; i < n; i++) {
        x[i] = l + i * delta;
        y[i] = func_value(x[i]);
    }
}
double net::func_value(double a) {   // value of function in point a
    return a * a * sin(a);
    //return sin(a * a * a);
}
//
class lagrange   // class lagranje
{
public:
    lagrange(net f);
    double lag_value(double a) { return p.value(a);}
    friend ostream & operator <<(ostream & out, const lagrange & l) { return out << l.p;};
private:
    polynomial p;
};
//
class newton
{
public:
    newton(net f);
    double new_value(double a) { return p.value(a);}
    friend ostream & operator <<(ostream & out, const newton & n) { return out << n.p;};
private:
    polynomial p;
};
lagrange::lagrange(net f) {
    int i, j, n = f.n; // n - number of points
    double val, coef;
    polynomial c(1, 1);   // init of monom
    for (i = 0; i < n; i++) {
        polynomial temp(1, 0);
        val = f.x[i];
        for (j = 0; j < n; j++) {
            if (i != j) {
                c.change(-f.x[j]);
                temp *= c;
                coef = (1/(val - f.x[j]));
                temp *= coef;
            }
        }
        if (f.y[i] != 0) {
            temp *= f.y[i];
            p += temp;
            cout << temp << endl;
        }
    }
}
newton::newton(net f) {
    int i, j, n = f.n;
    double * val1 = new double [n];
    double * val2 = new double [n];
    polynomial one(1, 0);
    polynomial temp(1, 0);   // polynomial (1*x^0)
    polynomial c(1, 1);     // polynomial (1*x^1)
    for (i = 0; i < n; i++) {
        val1[i] = f.y[i];
        val2[i] = 0;
    }
    if (!null(val1[0])) {
        p = one;
        p *= val1[0];
    }
    for (i = 1; i < n; i++) {
        for (j = 0; j < n - i; j++)  // counting the difference
            val2[j] = (val1[j] - val1[j + 1])/(f.x[j] - f.x[j + i]);
        for (j = 0; j < i; j++) {    // forming of polynomial in power i
            c.change(-f.x[j]);
            temp *= c;
        }
        if (!null(val2[0])) {  // if number != 0
            temp *= val2[0];
            p += temp;
        }
        double * t = val1;
        val1 = val2;
        val2 = t;
        temp = one;
    }
}
double polynomial::value(double a) {  // gorner's sheme for polynomial
    int n = head->pow;
    double s = 0;
    monom * m = head;
    while(m != NULL) {   // while there elements in list
        s *= a;
        if (n == m->pow) {
            s += m->coef;
            m = m->next;
        }
        n--;
    }
    while(n >= 0) {  // if the last polynomial had the power more than null
        s *= a;
        n--;
    }
    return s;  // value of polynomial in point a
}
void polynomial::change(double a) {
    monom * m = head->next;
    if (!null(a)) {     // if (a != 0)
        if (m != NULL)
            m->coef = a;
        else
            head->next = new monom(a, 0);
    }
    else
        del_next(head);
}
void polynomial::add_after(double c, int p, monom * m) {   // addition after current monom
    monom * m1 = new monom(c, p, m->next);
    m->next = m1;
}
void polynomial::del_next(monom * m) {   //  deleting next monom
    if (m->next != NULL) {
        monom * m1 = m->next;
        m->next = m1->next;
        delete m1;
    }
}
monom * polynomial::del_head(monom * m) {  // deleting head
    monom * m1 = m;
    m1 = m->next;
    delete m;
    return m1;
}
monom * polynomial::copy(monom * m) { // copying list to another one
    monom * m1 = m;
    monom * m2 = m->next;
    monom * p = new monom(m->coef, m->pow);
    monom * temp = p;
    while (m2 != NULL) {
        m1 = m2;
        m2 = m2->next;
        add_after(m1->coef, m1->pow, temp);
           temp = temp->next;
    }
    return p;
}
monom * polynomial::del_list(monom * m) {
    monom * m1 = m;
    while (m1 != NULL) {
        m = m->next;
        delete m1;
        m1 = m;
    }
    return m;
}
monom * polynomial::add(monom * m1, monom * m2) { // polinoms addition
    monom * p1 = m1;  // keep result in p1 list
    monom * p2 = m2;  // pointer to the previous monom
    monom * p = p1;  // pointer to the head of list
    m1 = p;
    while(p1 != NULL && p2 != NULL) { // going through the list to the end
        if (p1->pow < p2->pow) {
            if (m1 == p1) {
                p = new monom(p2->coef, p2->pow, m1);
                m1 = p;
                p2 = p2->next;
                p1 = m1;   //
            }
            else {
                add_after(p2->coef, p2->pow, p);  // addition of new elements in proper places
                p = p->next;
                p1 = p->next;
                p2 = p2->next;
            }
        }
        else {
            if (p1->pow == p2->pow) {
                p1->coef += p2->coef;
                if (null(p1->coef)) {
                    if (p1 == m1) { // deleting head
                        p1 = del_head(m1);
                        m1 = p1;
                        p = m1;
                        p2 = p2->next;
                    }
                    else {
                        del_next(p);
                        p1 = p->next;
                        p2 = p2->next;
                    }
                }
                else {
                    p = p1;
                    p1 = p1->next;
                    p2 = p2->next;
                }
            }
            else {
                p = p1;
                p1 = p1->next;
            }
        }
    }
    if (p2 != NULL)  // copying the last part of p2 list
        p1 = copy(p2);
    return m1;
}
monom * polynomial::mult_num(monom * m, double a) // multiplying polynomial on the number
{
    monom * p1 = m;
    monom * p = p1; //the previous element of p1
    m = p;
    while(p1 != NULL) { // deleting of null-elements
        p1->coef *= a;
        if (null(p1->coef)) {  // if the element less than epsilon
            if (p1 == m) { // if the element in head
                m = del_head(m);
                p1 = m;
                p = p1;
            }
            else {
                del_next(p);
                p = p1;
                p1 = p1->next;
            }
        }
        else {
            p = p1;
            p1 = p1->next;
        }
    }
    return m;
}
monom * polynomial::mult_pol(monom * m1, monom * m2) {  // multiplying of monoms
    monom* p1 = m1;
    monom* p2 = m2;
    int power;
    monom * temp = new monom(p1->coef * p2->coef, p1->pow + p2->pow, NULL); // initit the head of list
    p2 = p2->next;
    monom* p_temp = temp;
    monom* p_temp_next = p_temp->next;
    monom* p_temp_prev = p_temp;
    while (p1 != NULL) { // going through the first list
        while (p2 != NULL) {
            power = p1->pow + p2->pow;
            if (power == p_temp->pow) {
                p_temp->coef += p1->coef * p2->coef;
                if (null(p_temp->coef)) {
                    if (temp == p_temp) {
                        temp = del_head(temp);
                        p_temp = temp;
                        p_temp_prev = p_temp;
                    }
                    else
                        del_next(p_temp_prev);
                }
                p2 = p2->next;
            }
            else {
                if ((p_temp_next == NULL) || (power > p_temp_next->pow)) {
                    add_after(p1->coef * p2->coef, power, p_temp);
                    p2 = p2->next;
                }
                p_temp_prev = p_temp;
                p_temp = p_temp->next;
                p_temp_next = p_temp->next;
            }
        }
        p2 = m2;
        p_temp = temp; // returning to the beginning of the list
        p_temp_prev = p_temp;
        p_temp_next = p_temp->next;
        p1 = p1->next;  // transition to the next monom
    }
    m1 = del_list(m1);
    m1 = temp;
    return m1;
}
ostream & operator << (ostream& out, const polynomial & p) {
    monom * h = p.head;
    while (h != NULL) {
        out << setw(15) << scientific << h->coef << "*x^" << h->pow << " ";
        h = h->next;
    }
    out << endl;
    return out;
}
double polynomial::operator ()(double a) {  //  value in point a
    if (head == NULL)
        return 0;
    return value(a);
}

polynomial & polynomial::operator +=(const polynomial & p) {
    if (head == NULL) // if the fist polynom is empty, copy the second to the first
        head = copy(p.head);
    else
        head = add(head, p.head);
    return * this;
}
polynomial polynomial::operator +(const polynomial & p) {
    polynomial h;
    h.head = copy(head);
    h += p;
    return h;
}
polynomial & polynomial::operator *(double a) {
    polynomial new_p;
    new_p.head = copy(head);
    new_p *= a;
    return * this;
}
polynomial & polynomial::operator *=(double a) {
    if (head != NULL)
        head = mult_num(head, a);
    return * this;
}
polynomial & polynomial::operator *=(const polynomial p) { // multiplication with polynomial
    head = mult_pol(head, p.head);
    return * this;
}
polynomial polynomial::operator *(const polynomial & p) {  // multiplication with polynomial
    polynomial h;
    h.head = copy(head);
    h *= p;
    return h;
}
polynomial & polynomial::operator =(const polynomial & p) {  //  assignment
    if (p.head != head) {
        head = del_list(head);
        head = copy(p.head);
    }
    return * this;
}

int main() {
    net s(0, 1.5, 4); // (a, b, n) [a, b] otrezok, n - kolichestvo
    lagrange l(s);
    cout << "lagrange: " << l << endl;
    newton n(s);
    cout << "newton: " << n << endl;
    int i = 0;
    long double x = 0;
    cout << setw(15) << "x" << setw(15) << "y" << setw(15) << "lagrange" << setw(15) << "newton" << endl;
    cout << "№0" << "  ";
    cout << scientific << s.x[i] << setw(15) << s.y[i];
    cout << scientific << setw(15) << l.lag_value(s.x[i]) << setw(15) << n.new_value(s.x[i]) << endl;
    for (i = 1; i < s.n; i++) {
        x = (s.x[i - 1] + s.x[i])/2;
        cout << endl << "    ";
        cout << scientific << x << setw(15) << s.func_value(x) << setw(15); //
        cout << scientific << setw(15) << l.lag_value(x) << setw(15) << n.new_value(x) << endl;
        cout << endl << "№" << i << "  ";
        cout << scientific << s.x[i] << setw(15) << s.y[i] << setw(15); // net pointsz`
        cout << scientific << setw(15) << l.lag_value(s.x[i]) << setw(15) << n.new_value(s.x[i]) << endl;
    }
    return 0;
}