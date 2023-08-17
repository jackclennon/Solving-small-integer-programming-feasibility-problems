//#include "Highs.h"
//#include "lp_data/HighsModelUtils.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream> 
#include <tuple>
#include <bitset>
#include <algorithm>
#include <time.h>

using namespace std;

struct IP {

  int nVars;
  int nIneqs;
  int nEqs;
  vector<vector<int>> ineqs;
  vector<vector<int>> eqs;
  int presolveRows;
  int presolveCols;
  size_t hash_code;

  bool searchHash(vector<size_t> hks) {
    for (int i = 0; i < hks.size(); i++)
      if (hks[i] == hash_code)
        return true;
    return false;
  }
};

vector<vector<long double>> transpose(vector<vector<long double>> b){
  if (b.size() == 0)
    return b;
  vector<vector<long double>> t(b[0].size(), vector<long double>(b.size(), 0));
  for (int i = 0; i < b.size(); i++){
    for (int j = 0; j < b[i].size(); j++){
      t[j][i] = b[i][j];
    }
  }
  return t;
}

template <typename S> 
ostream& operator<<(ostream& os, const vector<S>& vector) {
    os << " ";
    for (auto element : vector) {
        os << element << " ";
    }
    os << "\n";
    return os;
}

struct tableauStruct {
  vector<vector<long double>> A;
  vector<long double> b;
  vector<long double> c;
  vector<long double> z;
  vector<long double> reducedCosts;
  vector<long double> div;
  vector<int> basis;
  int nCols;
  int nRows;
  int nEqs; //this is used to make sure we don't put an augmentation variable back into the basis
  string status;
};

void divide(tableauStruct *tab, long double a, int rowInd) {
  if (a == 0) {
    cin.clear();
    //cout << "\na == 0";
    cin.clear();
    cin.get();
  }
  //tab->div[rowInd] /= a;
  tab->b[rowInd] /= a;
  for (int i= 0; i < tab->nCols; i++)
    tab->A[rowInd][i] /= a;
}

void multiply(tableauStruct *tab, long double a, int rowInd) {
  for (int i= 0; i < tab->nCols; i++)
    tab->A[rowInd][i] *= a;
  tab->b[rowInd] *= a;
}

void fraction(tableauStruct *tab, long double num, long double den, int rowInd) {
  // divides the row "rowInd" by num/den
  // best to use fraction(tab, 1, a, row) over multiply(tab, a, row) as it maintains normalisation
  multiply(tab, den, rowInd);
  divide(tab, num, rowInd);
}

enum presolveType {EmptyColumn=0, RowSingleton=1, ColumnSingleton=2};

struct presolveStep {
  presolveType type;
  vector<long double> row;
  long double b;
  bool eq; //true means equality, false means inequality
  long double x; //value of the x component that has been removed
  int columnIndex;
};

struct problem {
  bool infeasible = false;
  vector<long double> lb;
  vector<long double> ub;
  vector<presolveStep> presolveStack;
};


void pivotRow(tableauStruct *tab, int pr, int pivotCol, int row) {
  /*
  //get pivot element
  auto rpc = tab->A[row][pivotCol];
  //auto x = tab->div[row];
  //auto y = tab->div[pivotRow];
  //tab->div[row] = x*y;
  for (int i = 0; i < tab->nCols; i++) {
    tab->A[row][i] = y*tab->A[row][i] - rpc*tab->A[pivotRow][i];
  }
  tab->b[row] =  y*tab->b[row] - rpc*tab->b[pivotRow];
  */
  //reduce by the gcd
  //divide(tab, 1, row);

  /* THIS IS REDUNDANT!!!
  auto v = tab->A[pivotRow][pivotCol];
  long double q = 1.;
  for (int i=0; i<tab->nCols; i++) {
    tab->A[pivotRow][i]/=v;
    tab->b[pivotRow]/=v;
  }
  */
  //cout << "\nPivot Row: " << tab->A[pivotRow];
  //cout << "\nPivot Col: " << transpose(tab->A)[pivotCol];
  long double q = 1.;
  q = tab->A[row][pivotCol];
  //cout << "\nMultiplication Factor: " << q;
  for (int i=0; i<tab->nCols; i++)
    tab->A[row][i] -= q*tab->A[pr][i];
  tab->b[row]-= q*tab->b[pr];
}

void print(tableauStruct *tab, int row) {
  vector<double> toPrint;
  for (auto e: tab->A[row]) {
    toPrint.push_back((double)e/(double)tab->div[row]);
  }
  toPrint.push_back((double)tab->b[row]/(double)tab->div[row]);
  cout << toPrint;
}

tableauStruct *makeTableau(vector<vector<long double>> A, vector<long double> b, vector<long double> c, int nEqs) {
  tableauStruct *tab = new tableauStruct{};
  tab->A = A;
  tab->b = b;
  tab->div = vector<long double>(A.size(), 1);
  tab->nCols = A[0].size();
  tab->nRows = A.size();
  tab->basis;
  tab->c = c;
  tab->z = vector<long double> (tab->nCols, 0.);
  tab->reducedCosts = vector<long double> (tab->nCols, 0.);

  //for (int i = tab->nCols-tab->nRows; i < tab->nCols; i++)
  /*
  for (int i = 0; i < tab->nRows; i++) {
    if (i < tab->nRows-nEqs)
      tab->basis.push_back(tab->nCols-(tab->nRows-nEqs)+i);
    else
      tab->basis.push_back(i);
  }
  */
  for (int i = 0; i < tab->nRows; i++)
    tab->basis.push_back(tab->nCols-tab->nRows+i);
  tab->nEqs = nEqs;
  tab->status = "unknown";
  return tab;
}

inline void hash_combine(std::size_t& seed) { }

template <typename T, typename... Rest>
inline void hash_combine(std::size_t& seed, const T& v, Rest... rest) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    hash_combine(seed, rest...);
}

struct IPHash {
  size_t operator()(IP const& ip) const noexcept {
    string str = "";
    str += to_string(ip.nVars);
    str += ".";
    str += to_string(ip.nIneqs); 
    str += ".";
    str += to_string(ip.nEqs);
    str += ".";
    for (auto l: ip.ineqs) {
      for (auto e: l) {
        str += to_string(e);
        str += ",";
      }
      str += ".";
    }
    for (auto l: ip.eqs) {
      for (auto e: l) {
        str += to_string(e);
        str += ",";
      }
      str += ".";
    }
    return hash<string>{}(str);
  }  
};

vector<vector<long double>> findMinus(vector<vector<long double>> A) {
  for (int row=0; row<A.size(); row++)
    for (int col=0; col<A[0].size(); col++) 
      A[row][col] = -A[row][col];
  return A;
}

vector<vector<int>> findMinus(vector<vector<int>> A) {
  for (int row=0; row<A.size(); row++)
    for (int col=0; col<A[0].size(); col++) 
      A[row][col] = -A[row][col];
  return A;
}

bool IPcomp(IP a, IP b) {
  if (a.nVars == b.nVars) {
    if (a.nIneqs == b.nIneqs) {
      if (a.nEqs == b.nEqs) {
        bool flag1 = true;
        for (int i = 0; i < a.ineqs.size(); i++) {
          if (a.ineqs[i] != b.ineqs[i]) {
            flag1 = false;
          }
        }
        bool flag2 = true;
        for (int i = 0; i < a.eqs.size(); i++) {
          if (a.eqs[i] != b.eqs[i]) {
            flag2 = false;
          }
        }
        if (flag1 && flag2)
          return true;
      }
    }
  }
  return false;
}

IP cases[nCases];
int recent = 0;
int emptyCount = 0;
int emptyCountInt = 0;

int countPresolveRuleApplications[14];

vector<string> split(string str) {
  vector<string> strings;
  int startIndex = 0, endIndex = 0;
  for (int i = 0; i <= str.size(); i++) {
      if (str[i] == ' ') {
          endIndex = i;
          string temp;
          temp.append(str, startIndex, endIndex - startIndex);
          strings.push_back(temp);
          startIndex = endIndex + 1;
      }
  }
  return strings;
}

void createIP(struct IP *ip, int nV, vector<string> in, vector<string> eq) {
  ip->nVars = nV;
  ip->nIneqs = in.size();
  ip->nEqs = eq.size();

  for (int i = 0; i < in.size(); i++) {
    vector<int> constr;
    vector<string> strConstr = split(in[i]);
    for (int j = 0; j < strConstr.size(); j++)
      constr.push_back(stoll(strConstr[j]));
    ip->ineqs.push_back(constr);
  }

  for (int i = 0; i < eq.size(); i++) {
    vector<int> constr;
    vector<string> strConstr = split(eq[i]);
    for (int j = 0; j < strConstr.size(); j++)
      constr.push_back(stoll(strConstr[j]));
    ip->eqs.push_back(constr);
  }

  ip->ineqs = findMinus(ip->ineqs);
}

void saveTestCase(IP *ip) {
  cases[recent] = *ip;
  recent++;
}

vector<vector<long double>> transposeD(vector<vector<long double>> b){
  if (b.size() == 0)
    return b;
  vector<vector<long double>> t(b[0].size(), vector<long double>());
  for (int i = 0; i < b.size(); i++){
    for (int j = 0; j < b[i].size(); j++){
      t[j].push_back(b[i][j]);
    }
  }
  
  return t;
}


vector<vector<long double>> glue(vector<vector<long double>> A, vector<vector<long double>> B) {
  if (A.size() == 0) 
    return B;
  if (B.size() == 0)
    return A;
  if (A.size() == B.size()) {
    for (int row=0; row<A.size(); row++) {
      for (auto e: B[row])
        A[row].push_back(e);
    }
  } /*else {
    cout << "\n\n\tWRONG DIMENSIONS\n\n";
    cout << "\nA: \n";
    for (auto e: A)
      cout << e;
    cout << "\nB: \n";
    for (auto e: B)
      cout << "row: " << e;
  }*/
  return A;
}

vector<long double> glueVec(vector<long double> A, vector<long double> B) {
  for (auto e: B)
    A.push_back(e);
  return A;
}

vector<vector<long double>> I(int x, int y) {
  vector<vector<long double>> A;
  for (int row=0; row<y; row++) {
    vector<long double> r;
    for (int col=0; col<x; col++) {
        r.push_back((long double)(col==row));
    }
    A.push_back(r);
  }
  return A;
}

tuple<vector<vector<long double>>, vector<long double>> extractAb(vector<vector<int>> p) {
  vector<vector<long double>> A;
  vector<long double> b;

  for (int i=0; i<p.size(); i++){
    vector<long double> r;
    for (int j=1; j<p[i].size(); j++) {
      r.push_back(p[i][j]);
    }
    A.push_back(r);
    b.push_back(-p[i][0]);
  }
  return make_tuple(A, b);
}

vector<vector<long double>> upToRow(vector<vector<long double>> A, int row) {
  vector<vector<long double>> B;
  for (int i=0; i<row; i++)
    B.push_back(A[i]);
  return B;
}

long double maxVec(vector<long double> vec) {
  long double max = -100000000000.;
  for (auto e: vec)
    if (e > max)
      max = e;
  return max;
}

bool isIn(int i, vector<int> s) {
  for (auto e: s)
    if (e == i)
      return true;
  return false;
}

bool isInLd(int i, vector<long double> s) {
  for (auto e: s)
    if (e == i)
      return true;
  return false;
}

tuple<tuple<int, int>, bool> getPivot(vector<vector<long double>> tableau, vector<long double> r, vector<int> B) {
  int row = -1;
  int col = -1;

  vector<long double> b;
  for (auto row: tableau)
    b.push_back(row.back());

  vector<int> rows;

  //pick largest non-zero index of the last row, first part of bland's rule
  vector<long double> dr;
  for (int i = 0; i < r.size()-1; i++) 
    if (!isIn(i, B))
      dr.push_back(r[i]);
  auto l = maxVec(dr);
  for (int i = 0; i < r.size(); i++) {
    if (l == r[i]) {
      col = i;
      break;
    }
  }

  //cout << "\ndr: " << dr;

  //cout << "\ncol: " << col << "\n";
  //find the min row, taking the smallest index if we have a tie, second part if bland's rule
  long double min = 10000000;
  auto t = transposeD(tableau);
  auto A = transposeD(upToRow(t, t.size()-1));
  for (int j=0; j<tableau.size(); j++) {
    auto v = b[j]/A[j][col];
    if (v < min) {
      min = v;
      rows = vector<int> {j};
    } else if (min == v) {
      rows.push_back(j);
    }
  }

  /*
  bool unbounded = true;
  for (auto e: rows)
    if (e >= 0)
      unbounded = false;

  if (unbounded)
    return make_tuple(make_tuple(row, col), unbounded);
  */
  //cout << "\nRow length: " << rows.size();
  row = 100000000;
  for (int i = 0; i < rows.size(); i++)
    if (rows[i] < row) {
      row = rows[i];
    }
  //cout << "\nrow: " << row << "\n";

  //I suspect this means that it is infeasible????
  if (row == 100000000)
    return make_tuple(make_tuple(row, col), true);

  return make_tuple(make_tuple(row, col), false);
}

long double get_b(tableauStruct *tab, int i) {
  if (i < tab->nRows)
    return (long double)tab->b[i]/(long double)tab->div[i];
  //cout << "\nIndex too large";
  return 0.;
}

vector<vector<long double>> pivot(vector<vector<long double>> tableau, tuple<int, int> p) {
  int row, col;
  tie(row, col) = p;

  auto v = tableau[row][col];

  /*
  cout << "\nv: " << v;
  cin.clear();
  if (v==0.) {
    cout << "\nv is 0";
    cin.get();
  }
  */

  //this is where we're getting the error from
  auto gamma = transposeD(tableau)[col];

  for (int i=0; i<tableau.size(); i++){
    for (int j=0; j<tableau[i].size(); j++){
      if (i==row) {
        tableau[i][j] /= v;
      }
    }
  }

  for (int i=0; i<tableau.size(); i++){
    for (int j=0; j<tableau[i].size(); j++){
      if (i!=row) {
        tableau[i][j] -= gamma[i]*tableau[row][j];
      }
    }
  }

  return tableau;
}

void pivotDual(tableauStruct *tableau, int row, int col) {

  auto v = tableau->A[row][col];

  //cout << "\nrow: " << row << "\tcol: " << col;

  //auto gamma = transpose(tableau->A)[col];
  divide(tableau, v, row);

  for (int i=0; i<tableau->nRows; i++){
    if (i!=row) {
      pivotRow(tableau, row, col, i);
    }
  }
}


bool canImprove(vector<long double> r) {
  for (auto e: r) {
    //cout << "CAN IMPROVE\n";
    if (e > 0)
      return true;
  }
  return false;
}

bool canImproveDual(tableauStruct *tab) {
  //if we have any augmentation variables in the basis then we should remove them
  /*
  for (auto e: tab->basis) {
    if (e >= tab->nCols - tab->nEqs)
      return true;
  }
  */

  //check that the reduced costs are all greater than or equal to 0
  /*
  for (auto e: tab->reducedCosts) {
    if (e < 0) {
      tab->status = "infeasible";
      return false;
    }
  }
  */

  /*
  cout << "\nb: " << tab->b;
  cout << "\nnRows: " << tab->nRows;
  cout << "\nb.size(): " << tab->b.size();
  */
  //cout << "\nThis far";
  bool flag2 = false;
  for (int i = 0; i < tab->nRows; i++) {
    if (tab->b[i] < 0.) {
      //cout << "\nCAN IMPROVE: " << get_b(*tab, i);
      flag2 = true;
    }
  }
  
  //if (!flag2)
  //  tab->status = "can't improve"; // cant improve

  //cout << "\nFLAG 2: " << flag2;

  return flag2;
}

vector<long double> updateZ(vector<long double> c, vector<int> B, vector<vector<long double>> tab) {
  vector<long double> z(tab[0].size(), 0.);
  for (int row=0; row<tab.size(); row++) {
    for (int col=0; col<tab[row].size(); col++) {
      /*
      if (isnan(tab[row][col])) {
        cin.clear();
        cout << "\n v is nan" << tab[row][col] << " col " << col << " row " << row;
        cin.get();
      }
      if (isnan(c[B[row]])) {
        cin.clear();
        cout << "\n c is nan " << " row: " << row << " B[row] " << B[row] << " cost " << c;
        cout << c[B[row]];
        cin.get();
      }
      */
      z[col] += c[B[row]] * tab[row][col];
    }
  }
  return z;
}

vector<int> updateZInt(vector<int> c, vector<int> B, vector<vector<int>> A, vector<int> b) {
  vector<vector<int>> tab;
  for (int i =0; i < b.size(); i++) {
    auto row = A[i];
    row.push_back(b[i]);
    tab.push_back(row);
  }
  vector<int> z(tab[0].size(), 0.);
  for (int row=0; row<tab.size(); row++) {
    for (int col=0; col<tab[row].size(); col++) {
      /*
      if (isnan(tab[row][col])) {
        cin.clear();
        cout << "\n v is nan" << tab[row][col] << " col " << col << " row " << row;
        cin.get();
      }
      if (isnan(c[B[row]])) {
        cin.clear();
        cout << "\n c is nan " << " row: " << row << " B[row] " << B[row] << " cost " << c;
        cout << c[B[row]];
        cin.get();
      }
      */
      z[col] += c[B[row]] * tab[row][col];
    }
  }
  return z;
}

void updateZTab(tableauStruct *tab) {
  tab->z = vector<long double> (tab->nCols, 0.);
  for (int row=0; row<tab->nRows; row++) {
    for (int col=0; col<tab->nCols; col++) {
      //cout << "\nSize of Cost Vector: " << tab->c.size();
      //cout << "\nBasis Value: " << tab->basis[row] << "\nIndex: " << row;
      auto cost = tab->c[tab->basis[row]];
      auto val = (long double)tab->A[row][col]/(long double)tab->div[row];
      tab->z[col] += cost * val;
    }
  }
}

void updateReducedCostsTab(tableauStruct *tab) {
  for (int i=0; i < tab->nCols; i++)
    tab->reducedCosts[i] = tab->c[i]-tab->z[i];
}

vector<long double> updateR(vector<long double> z, vector<long double> c) {
  vector<long double> u;
  for (int i=0; i < z.size(); i++)
    u.push_back(c[i]-z[i]);
  return u;
}

vector<int> updateRInt(vector<int> z, vector<int> c) {
  vector<int> u;
  for (int i=0; i < z.size(); i++)
    u.push_back(c[i]-z[i]);
  return u;
}

long double simplex(vector<long double> cost, vector<vector<long double>> A, vector<long double> b, int nonSlack) {

  /*
  
  Last task is to figure out how to keep track of the divisions to keep everything representable
  */


  //make the tableau
  //this will only work well if we use the dual where we don't need an initial solution
  vector<vector<long double>> tableau;
  for (int i=0; i < A.size(); i++) {
    auto row = A[i];
    row.push_back((long double)b[i]);
    tableau.push_back(row);
  }

  //add a 1 to the end of each row, where we will keep track of the divisor for the row ?

  bool q = false;
  bool unbounded = false;
  tuple<int, int> p;
  string toPrint = "";
  vector<int> B;

  for (int i=0; i<A.size(); i++) {
    B.push_back(nonSlack + i);
  }

  /*
  cout << "\n B: " << B;
  cout << "\n # of vars " << A[0].size();
  cout << "\n # of constraints " << A.size();
  cout << "\nLength of cost " << cost.size();
  */
  auto z = updateZ(cost, B, tableau);
  auto r = updateR(z, cost);
  /*
  cout << "\nz: " << z;
  cout << "\ncost: " << cost;

  cout << "\n INITIAL R: " << r.size();
  cout << "\n INITIAL R: " << r;
  */
  int its = 0;
  
  while (canImprove(r) && !unbounded && its < 10) {
    its++;
    tie(p, unbounded) = getPivot(tableau, r, B);
    if (unbounded) {
      cout << "\nunbouded\n";
      break;
    }

    int row, col;
    tie(row, col) = p;
    if (row < 0 || col < 0)
      break;
    B[row] = col;

    tableau = pivot(tableau, p);
    toPrint += "\n";
    for (auto row: tableau) {
      for (auto e: row)
        toPrint += to_string((int)floor(e)) + " ";
      toPrint += "\n";
    }

    z = updateZ(cost, B, tableau);
    r = updateR(z, cost);

    for (auto e: r)
    if (isnan(e)) {
      cout << "\n\tContinuing\n";
      break;
    }

    cout << "\nAfter updating\n";
    for (auto row: tableau)
      cout << row;
    cout << "\n R: " << r << "\n" << z;

    q = true;
  }
  //extract the data we want and return
  return z.back();
}

void printCase(int ind) {
  cout << "\n";
    cout << cases[ind].nVars << "\n";
    cout << cases[ind].nIneqs << "\n";
    cout << cases[ind].ineqs << "\n";
    cout << cases[ind].nEqs << "\n";
    cout << cases[ind].eqs << "\n";
}

vector<vector<long double>> makeZeroes(int h, int w) {
  return vector<vector<long double>> (h, vector<long double> (w, 0));
}

vector<long double> getSplitSolution(tableauStruct *tableau) {
  vector<long double> splitSolution(tableau->nCols, 0.);
  //vector<long double> temp(tableau->nRows+1, 0.);
  //vector<vector<long double>> solMatrix (tableau->nRows, temp);

  //int bc = 0;
  //cout << "\nA\tdiv\tb\tindex\n";
  for (int i=0; i < tableau->nRows; i++) {
    if (tableau->A[i][tableau->basis[i]] == 0) {
      //cout << "\nSplit sol size: " << splitSolution.size() << "\tnCols: " << tableau->nCols;
      //cout << "\ni: " << i << "\tbasis[i]: " << tableau->basis[i] << "\n";
      splitSolution[tableau->basis[i]] = 0;
    }
    else
      splitSolution[tableau->basis[i]] = tableau->b[i]/tableau->A[i][tableau->basis[i]];
  }
  return splitSolution;
}

vector<long double> doPostsolve(problem *p, vector<long double> x) {
  //cout << "\nStart postsolve: " << p->presolveStack.size();
  //This should go after we recombine x
  //There are only 3 postsolves that I need to apply so far
  for (int i=p->presolveStack.size()-1; i>=0; i--) {
    //cout << "\n\nType: " << p->presolveStack[i].type << "\tx0: " << x;
    switch (p->presolveStack[i].type) {
      case (EmptyColumn): {
        //cout << "\nEmpty Row";
        //In this case we should put back in a variable in x set to 0
        //cout << "\nx1" << x;
        //cout << "\ninsert index: " << p->presolveStack[i].columnIndex;
        x.insert(x.begin()+p->presolveStack[i].columnIndex, 0.);
        //cout << "\nx2: " << x;
      }break;
      case (ColumnSingleton): {
        //cout << "\nColumn Singleton";
        //Here we need to use the row we stored to calculate a value for a new x variable, which we then need to add back to the problem
        //first add back in x, set it to 0 for now
        //cout << "\nx1: " << x;
        //cout << "\ninsert index: " << p->presolveStack[i].columnIndex;
        x.insert(x.begin()+p->presolveStack[i].columnIndex, 0.);
        //then find (row)x
        long double q = 0.;
        //segfault here. Somehow, even though its an EmptyColumn, we're entering into this case>??
        //cout << "\n\nRow: " << p->presolveStack[i].row;
        for (int j=0; j<x.size(); j++)
          q += p->presolveStack[i].row[j]*x[j];
        //now we find the new value for x
        x[p->presolveStack[i].columnIndex] = (p->presolveStack[i].b - q)/p->presolveStack[i].row[p->presolveStack[i].columnIndex];
        //cout << "\nx1: " << x;
        //test to see if x is an integer, if the constraint was an equality constraint
        if (p->presolveStack[i].eq) {
          if (x[p->presolveStack[i].columnIndex] - floor(x[p->presolveStack[i].columnIndex]) != 0.) {
             p->infeasible = true;
             return x;
          }
        }
      }break;
    }
    if (p->infeasible)
      return x;
  }

  //Now check the bounds
  for (int i=0; i<x.size(); i++) {
    if (x[i] > p->ub[i])
      x[i] = p->ub[i];
    if (x[i] < p->lb[i])
      x[i] = p->lb[i];
  }

  return x;
}

int isPrimalFeasible(problem *p, tableauStruct *tab, vector<long double> splitSolution, vector<vector<long double>> A1, vector<vector<long double>> A2, vector<long double> b, int nVars) {
  
  //cout << "\nSplit Solution:\n";
  //cout << splitSolution;

  // reconstruct x
  vector<long double> x;
  /*
  // I don' think this is correct
  for (int i=0; i < nVars-nEqs; i++) {
    x.push_back(splitSolution[i] - splitSolution[i+nVars]);
  }
  for (int i=nVars-nEqs; i<nVars; i++)
    x.push_back(splitSolution[i]);
  */

  /*
  cout << "\nSplitsol: " << splitSolution;
  cout << "\nSize of A: " << tab->A.size();
  if (tab->A.size() > 0)
    cout << "\nnVars: " << nVars << "\tnumber of ineqs: " << tab->nRows-tab->nEqs << "\tActual nVars of A: " << tab->A[0].size()  << "\nNEqs: " << tab->nEqs;
  */

  for (int i=0; i<nVars; i++)
    x.push_back(splitSolution[i] - splitSolution[i+nVars]);

  //cout <<  "\nX entering the postsolve: " << x;

  if (p->presolveStack.size() == 0) {
    if (tab->A.size() == 0)
        return 1;
  } else {
    if (tab->A.size() == 0)
      return 1;
    x = doPostsolve(p, x);
    if (p->infeasible) {
      //cout << "\nPOSTSOLVE INFEASIBLE";
      return 0;
    }
  }

  /*
  cout << "\nSplit Solution: " << splitSolution;
  cout << "\nnVars: " << nVars << "\tx:\n";
  cout << x;
  */

  //impose upper and lower bounds on the variables
  for (int i=0; i<x.size(); i++) {
    if (x[i] < p->lb[i])
      x[i] = p->lb[i];
    if (x[i] > p->ub[i])
      x[i] = p->ub[i];
  }

  //check if it is a primal feasible solution
  vector<long double> v1, v2;
  for (auto row: A1) {
    long double e = 0;
    for (int i=0; i<row.size(); i++) {
      e += (long double)row[i]*x[i];
    }
    v1.push_back(e);
  }
  for (auto row: A2) {
    long double e = 0;
    for (int i=0; i<row.size(); i++) {
      e += (long double)row[i]*x[i];
    }
    v2.push_back(e);
  }

  /*
  cout << "\nSolutions\n";
  cout << x;
  cout << "\nSolution Values\n";
  cout << v1 << v2;
  cout << "\nbounds\n";
  cout << b;
  */

  for (auto e: v1)
    if (isnan(e))
        cin.get();

  for (auto e: v2)
    if (isnan(e))
        cin.get();

  //check if inequalities are valid
  for (int i=0; i<v1.size(); i++)
    if (v1[i] > b[i]) {
      //cout << "\nWHERE IT FAILS 1: " << v1[i] << "  >  " << b[i];
      return 0;
    }
  //check if the equalities are valid
  long double eps = 0.01;
  for (int i=v1.size(); i<b.size(); i++)
    if (v2[i-v1.size()] < b[i] - eps ||  v2[i-v1.size()] > b[i] + eps){
      //cout << "\nWHERE IT FAILS 2: " << v2[i-v1.size()] << "  >  " << b[i];
      return 0;
    }
  
  //delete(tableau);
  return 1;
}

int optimalityCheck(problem *p, tableauStruct *tab, vector<long double> x) {
  /*for (int i=0; i<x.size(); i++) 
    if (p->lb[i] > x[i] || p->ub[i] < x[i])
      return 0;*/
  //check if it is a primal feasible solution
  vector<long double> v1, v2;
  vector<vector<long double>> A1, A2;
  for (int i=0; i<tab->nRows-tab->nEqs; i++)
    A1.push_back(tab->A[i]);
  for (int i=tab->nRows-tab->nEqs; i<tab->nRows; i++)
    A2.push_back(tab->A[i]);
  for (auto row: A1) {
    long double e = 0;
    for (int i=0; i<row.size(); i++) {
      e += (long double)row[i]*x[i];
    }
    v1.push_back(e);
  }
  for (auto row: A2) {
    long double e = 0;
    for (int i=0; i<row.size(); i++) {
      e += (long double)row[i]*x[i];
    }
    v2.push_back(e);
  }

  //check if inequalities are valid
  for (int i=0; i<v1.size(); i++)
    if (v1[i] > tab->b[i]) {
      //cout << "\nWHERE IT FAILS 1: " << v1[i] << "  >  " << b[i];
      return 0;
    }
  //check if the equalities are valid
  long double eps = 0.01;
  for (int i=v1.size(); i<tab->b.size(); i++)
    if (v2[i-v1.size()] < tab->b[i] - eps ||  v2[i-v1.size()] > tab->b[i] + eps){
      //cout << "\nWHERE IT FAILS 2: " << v2[i-v1.size()] << "  >  " << b[i];
      return 0;
    }
  

  return 1;
}

void doPresolve(problem *p, tableauStruct *tab) {
  //This should go before we split x to get lower bounds on it
  vector<vector<long double>> At;
  vector<vector<long double>> Tt;
  vector<long double> bt;

  bool changed = true;
  //loop the presolve until it makes no changes
  while (changed) {
    /*
    cout << "\nA:\n";
    for (int y = 0; y<tab->A.size(); y++)
      cout << tab->b[y] << "\t" << tab->A[y];
    cout << "\n" << p->lb;
    cout << p->ub;
    */
    changed = false;

    //remove any columns with closed bounds
    auto T = transpose(tab->A);
    int I = -1;
    if (tab->A.size() != 0) {
      for (int i=0; i<tab->nCols; i++) {
        if (p->lb[i] == p->ub[i]) {
          for (int j=0; j<tab->A.size(); j++) {
            tab->b[j] -= tab->A[j][i]*p->lb[i];
          }
          I  = i;
          changed = true;
          break;
        }
        if (changed)
          break;
      }
      if (changed) {
        T.erase(T.begin()+I);
        tab->A = transpose(T);
        p->ub.erase(p->ub.begin()+I);
        p->lb.erase(p->lb.begin()+I);
        if (tab->A.size() == 0)
          return;
        if (tab->nRows > 0)
          tab->nCols = tab->A[0].size();
        else
          tab->nCols = 0;
        if (tab->nRows == 0)
          return;
        //cout << "\nRemoved Columns with closed bounds";
        continue;
      }
    }

    //Rule 1
    //remove all empty rows
    //In this rule we simply remove all empty rows
    Tt.clear();
    At.clear();
    bt.clear();
    int newNRows = 0;
    int newNEqs = 0;
    for (int r=0; r<tab->A.size(); r++) {
      auto row = tab->A[r];
      bool flag = false;
      for (auto e: row) {
        if (e != 0.) {
          At.push_back(row);
          bt.push_back(tab->b[r]);
          newNRows++;
          if (r >= tab->nRows - tab->nEqs)
            newNEqs++;
          flag = true;
          break;
        }
      }
      if (!flag) {
        if (r < tab->nRows - tab->nEqs) {
          if (tab->b[r] < 0) {
            //then it's infeasible
            p->infeasible = true;
            return;
          } 
        } else {
          if (tab->b[r] != 0) {
            //then it's infeasible
            p->infeasible = true;
            return;
          }
        }
        //cout << "\nEmpty Row: " << row;
        changed = true;
      }
    }
    if (changed) {
      tab->A = At;
      tab->b = bt;
      tab->nRows = newNRows;
      tab->nEqs = newNEqs;
      //cout << "\nEmpty Row";
      continue;
    }

    //Check if we reach a state with no possible solution
    bool triviallyInfeasible = false;
    for (int i=0; i < tab->nRows - tab->nEqs; i++) {
      auto eq = tab->A[i];
      bool nonZero = false;
      if (eq[0] > 0) {
        for (int i=0; i<tab->A[0].size(); i++) {
          if (eq[i] != 0)
            nonZero = true;
        }
        if (!nonZero) {
          triviallyInfeasible = true;
          break;
        }
      }
    }
    for (int i=tab->nRows - tab->nEqs; i<tab->nRows; i++) {
      auto eq = tab->A[i];
      bool nonZero = false;
      if (eq[0] != 0) {
        for (int i=0; i<tab->A[0].size(); i++) {
          if (eq[i] != 0)
            nonZero = true;
        }
        if (!nonZero) {
          triviallyInfeasible = true;
          break;
        }
      }
    }
   if (triviallyInfeasible) {
      //cout << "\nTrivially Infeasible";
      p->infeasible = true;
      return;
   }


    //Check if we reach a state with a trivial solution
    bool triviallyFeasible = true;
    //first, make sure that the trivial solution is feasible
    for (int i=0; i<tab->A[0].size(); i++) {
      if (p->lb[i] <= 0 && p->ub[i] >= 0) {
      } else {
        triviallyFeasible = false;
      }
    }
    if (triviallyFeasible) {
      //check if inequalities are valid
      for (int i=0; i < tab->nRows - tab->nEqs; i++)
        if (0 > tab->b[i])
          triviallyFeasible = false;
    }
    if (triviallyFeasible) {
      //check if the equalities are valid
      long double eps = 0.000001;
      for (int i=tab->nRows - tab->nEqs; i<tab->nRows; i++)
        if (0 < tab->b[i] - eps ||  0 > tab->b[i] + eps)
          triviallyFeasible = false;
    }
    if (triviallyFeasible) {
      //This one has a trivial solution, so we treat it as if it has reduced to empty
      /*cout << "\nA\n";
      for (int i=0; i<tab->A.size(); i++)
        cout << tab->b[i] << "\t" << tab->A[i];*/
      tab->A = vector<vector<long double>> {};
      tab->b = vector<long double> {};
      //cout << "\nTRIVIALLY FEASIBLE";
      return;
    }

    //Rule 2
    //empty column
    //In this rule we simply remove variables in x
    int newNCols = 0;
    T = transpose(tab->A);
    Tt.clear();
    bt.clear();
    I = -1;
    for (int c=0; c<T.size(); c++) {
      int It;
      auto col = T[c];
      bool flag = false;
      for (int r = 0; r<col.size(); r++) {
        auto e = col[r];
        if (e != 0.) {
          flag = false;
          break;
        }
        flag = true;
        It = c;
      }
      if (flag) {
        //we've found an empty column
        I = It;
        //cout << "\nEmpty Column: " << I;
        changed = true;
        break;
      }
    }
    
    if (changed) {
       //Add to the presolve stack
      p->presolveStack.push_back(presolveStep{type: EmptyColumn, eq: false, x:0., columnIndex: I});
      //remove the Ith column
      //cout << "\nWill remove column: " << I << ": " << T[I];
      T.erase(T.begin()+I);
      tab->A = transpose(T);

      p->lb.erase(p->lb.begin()+I);
      p->ub.erase(p->ub.begin()+I);
    
      tab->nRows = tab->A.size();
      if (tab->nRows > 0)
        tab->nCols = tab->A[0].size();
      else
        tab->nCols = 0;
      if ( tab->nRows == 0)
        tab->nEqs = 0;
      
      //cout << "\nEmpty Column";
      continue;
    }

    //Forcing Constraint
    vector<long double> h, g;
    for (int row=0; row<tab->nRows; row++) {
      long double eh = 0.;
      long double eg = 0.;
      for (int col=0; col<tab->A[0].size(); col++) {
        if (tab->A[row][col] < 0.) {
          eg += tab->A[row][col]*p->ub[col];
          eh += tab->A[row][col]*p->lb[col];
        } else if (tab->A[row][col] > 0.) {
          eh += tab->A[row][col]*p->ub[col];
          eg += tab->A[row][col]*p->lb[col];
        }
      }
      if (tab->b[row] > eh) {
        p->infeasible = true;
        return;
      }
      if (row >= tab->nRows - tab->nEqs) {
        if (tab->b[row] < eg) {
          p->infeasible = true;
          return;
        }
      }
      if ((eg == tab->b[row] || eh == tab->b[row]) && row >= tab->nRows - tab->nEqs) {
        //we have a forcing constraint
        T = transpose(tab->A);
        for (int col=tab->A[0].size()-1; col>=0; col--) {
          if (tab->A[row][col] != 0) {
            //remove these variables
            T.erase(T.begin() + col);
          }
        }
        tab->A  = transpose(T);
        //remove this row
        if (tab->A.size() == 0)
          return;
        tab->A.erase(tab->A.begin() + row);
        //update the nRows and nEqs
        if (row >= tab->nRows-tab->nEqs)
          tab->nEqs--;
        tab->nRows--;

        //cout << "\nForcing Constraint";

        changed = true;
      } else {
        //we can do nothing. However, we can record these bounds for use later
        h.push_back(eh);
        g.push_back(eg);
      }
    }

    //Rule 4
    //fixed column
    //In this rule we remove a column and a row
    //get the column index of the first column with a singleton 
    I = -1;
    int R = -1;
    T = transpose(tab->A);
    Tt.clear();
    bt.clear();
    for (int c=0; c<T.size(); c++) {
      int It;
      int Rt;
      auto col = T[c];
      bool flag = false;
      bool isFree = false;
      int n = 0;
      for (int r = 0; r<col.size(); r++) {
        if (col[r] != 0.) {
          n++;
          I = c;
          R = r;
        }
      }
      if (n==1) {
        //check if it is free
        if (p->ub[c] == inf && p->lb[c] == -inf) {
          changed = true;
          break;
        }
        //check if it is implied free
        for (int i=0; i<tab->nCols; i++) {
          if (i != I) {
            long double lp, up;
            if (tab->A[R][i] != 0.) {
              if (tab->A[R][i] > 0) {
                lp = (tab->b[R]-h[R])/tab->A[R][i] + p->ub[i];
                up = (tab->b[R]-g[R])/tab->A[R][i] + p->lb[i];
              } else if (tab->A[R][i] < 0) {
                lp = (tab->b[R]-g[R])/tab->A[R][i] + p->ub[i];
                up = (tab->b[R]-h[R])/tab->A[R][i] + p->lb[i];
              }
              if (p->lb[i] <= lp && lp <= up && up <= p->ub[i]) {
                changed = true;
                break;
              }
            }
          }
        }
        if (changed)
          break;
      }
    }

    if (changed) {
      //Add to the presolve stack
      //cout << "\nR: " << R;
      p->presolveStack.push_back(presolveStep{type: ColumnSingleton, row: tab->A[R], b: tab->b[R], eq: false, x:0, columnIndex: I});
      //remove the Ith column
      //cout << "\nWill remove column: " << I << ": " << T[I];
      /*cout << "\n\nStart Removing column and row\nA\n";
      for (auto e: transpose(T))
        cout << e;
      cout << "\nb: " << tab->b;
      */
      T.erase(T.begin()+I);
      p->lb.erase(p->lb.begin()+I);
      p->ub.erase(p->ub.begin()+I);
      if (T.size() > 0) {
        At = transpose(T);
        At.erase(At.begin()+R);
        //remove the Rth row
        tab->b.erase(tab->b.begin()+R);
      }
      else
        At = vector<vector<long double>> {};
      
      //cout << "\nWill remove row: " << R << ": " << At[R];
      /*
      cout << "\n\nAfter Removing column and row\nA\n";
      for (auto e: At)
        cout << e;
      cout << "\nb: " << tab->b;
      */
      
      tab->A = At;

      if (tab->nRows > 0)
        tab->nCols = tab->A[0].size();
      else
        tab->nCols = 0;
      //if (R >= tab->nCols - tab->nRows && tab->nEqs > 0)
      if (R >= tab->nRows - tab->nEqs && tab->nEqs > 0)
        tab->nEqs -= 1;

      tab->nRows = tab->A.size();
      //cout << "\nFixed Column";
      continue;
    } else if (changed) {
      changed = false;
    }


    changed = true;
    while (changed) {
      changed = false;
      //Rule 3
      //singleton row
      //In this one we keep the variable for x and we will fix it to the value that this technique finds in the postsolve
      I = -1;
      R = -1;
      At.clear();
      bt.clear();
      for (int r=0; r<tab->A.size(); r++) {
        int It;
        int Rt;
        auto row = tab->A[r];
        bool flag = false;
        for (int c = 0; c<row.size(); c++) {
          auto e = row[c];
          if (e != 0.) {
            if (!flag) {
              flag = true;
              It = c;
              Rt = r;
            } else {
              flag = false;
              It = -1;
              Rt = -1;
              break;
            }
          }
        }
        if (flag) {
          //we've found a row singleton
          I = It;
          R = Rt;
          //cout << "\nRow Singleton: " << I << ", " << R;
          changed = true;
          break;
        }
      }

      if (changed && R >= 0) {
        long double x = tab->b[R]/tab->A[R][I];
        //is it an equality constraint or an inequality constraint
        if (R < tab->nRows-tab->nEqs) {
          //inequality constraint
          if (tab->A[R][I]/abs(tab->A[R][I]) == 1) {
            //new upper bound
            if (p->lb[I] <= x) {
              if (p->ub[I] >= x)
                p->ub[I] = x;
            } else {
              //cout << "\n1";
              p->infeasible = true;
              return;
            }
          } else {
            //new lower bound
            if (p->ub[I] >= x) {
              if (p->lb[I] <= x)
                p->lb[I] = x;
            } else {
              //cout << "\n2";
              p->infeasible = true;
              return;
            }
          } 
        } else {
          //equational constraint
          //new upper bound
          if (p->lb[I] <= x)
            //new lower bound
            if (p->ub[I] >= x) {
              p->ub[I] = x;
              p->lb[I] = x;
            } else {
              //cout << "\n3";
              p->infeasible = true;
              return;
          }
          else {
            //cout << "\n4";
            p->infeasible = true;
            return;
          }
        }

        //remove the now redundant row
        for (int r=0; r<tab->A.size(); r++) {
          if (r != R) {
            At.push_back(tab->A[r]);
            bt.push_back(tab->b[r]);
          }
        }
        tab->A = At;
        tab->b = bt;

        //update the nRows and nEqs
        if (R >= tab->nRows-tab->nEqs)
          tab->nEqs--;
        tab->nRows--;

        if (tab->nRows == 0) 
          return;

        //Here would be a good place to do weakly dominated columns
        //1. check that it is on an equality constraint
        //2. check if the variable I is a column singleton
        //3. check that there is a column singleton in the row
        //R >= tab->nRows-tab->nEqs && 
        /*
        if (tab->nCols > 1 && (p->ub[I] <= inf || p->lb[I] >= -inf)) {
          bool Y = true;
          bool F = false;
          for (int csc = 0; csc < tab->A[0].size(); csc++) {
            for (int row=0; row<tab->A.size(); row++) {
              if (I != csc && R != row) {
                if (tab->A[row][csc] == 0.) {
                  F = true;
                }
                if (F) {
                  if (tab->A[row][csc] == 0.) {
                    F = false;
                    //we do not have a column singleton
                  } else {
                    break;
                  }
                }
              }
            }
          }
          for (int row=0; row<tab->A.size(); row++) {
            if (row != R && tab->A[row][I] != 0) {
              //we don't have a column singleton
              Y = false;
            }
          }
          if (!Y && F) {
            //now we can remove it (I think)
            p->ub.erase(p->ub.begin()+I);
            p->lb.erase(p->lb.begin()+I);
            T = transpose(tab->A);
            cout << "\nWeakly Dominated";
            T.erase(T.begin() + I);
            if (T.size() == 0)
              return;
            tab->A = transpose(T);
            if (tab->A.size() == 0)
              return;
            tab->nCols = tab->A[0].size();
          }

          cout << "\nA:\n";
          for (auto e: tab->A)
            cout << e;
          cout << "\nb: " << tab->b;
          cout << "\nlb: " << p->lb;
          cout << "\nub: " << p->ub;
          cout << "\nSingleton Row";
          continue;
        }*/
      }
    }

    /*
    cout << "\nA:\n";
    for (auto y: tab->A)
      cout << y;
    cout << "\nChanged: " << changed;
    */

    //Rule 5
    //redundant row
    //In this rule we remove all constraints that are dominated by a different constraint
    At.clear();
    bt.clear();
    newNRows = 0;
    /*
    for (int r1=0; r1<tab->nRows; r1++) {
      auto row1 = tab->A[r1];
      bool add = true;
      for (int r2=0; r2<tab->nRows; r2++) {
        auto row2 = tab->A[r2];
        if (row1 == row2 && r1 != r2) {
          if (tab->b[r1] < tab->b[r2]) {
            At.push_back(row1);
            bt.push_back(tab->b[r1]);
            newNRows++;
          } else if (tab->b[r1] == tab->b[r2]) {
            if (r1 < r2) {
              At.push_back(row1);
              bt.push_back(tab->b[r1]);
              newNRows++;
            }
          }
        }
      }
    }
    tab->nRows = newNRows;
    tab->A = At;
    tab->b = bt;
    */
    vector<vector<int>> likeRows;
    bool flag = false;
    for (int j=0;j<tab->A.size();j++) {
        for (int i=0;i<likeRows.size();i++) {
            if (tab->A[j] == tab->A[likeRows[i][0]]) {
                likeRows[i].push_back(j);
                flag = true;
                break;
            }
        }
        if (!flag)
            likeRows.push_back(vector<int> {j});
    }
    At.clear();
    bt.clear();
    int decEqs = 0;
    for (auto I: likeRows) {
        if (I.size() == 1) {
            At.push_back(tab->A[I[0]]);
            bt.push_back(tab->b[I[0]]);
        } else {
            changed = true;
            //find which rows dominate
            //find the minimum b for each set in likeRows, and this is the dominant
            //in the process, identify if we get an equality being greater than an inequality
            int mini = -1;
            long double mine = 1000000000000;
            for (auto i: I) {
                if (tab->b[i] <= mine) {
                    if (i > -1) {
                        //check if the old dominant row is an equality
                        //cout << "\nThis far 1";
                        //cout <<  "\nnRows: " << tab->nRows << "\tnEqs:  " << tab->nEqs <<  "\teqs starting index: " << tab->nRows - tab->nEqs;
                        if (tab->b[i] > mine && mini >= tab->nRows - tab->nEqs) {
                            //cout << "\nThis far 2";
                            //this must be infeasible because it means an equality is dominated
                            //cout << "\n5";
                            p->infeasible = true;
                            return;
                        }
                    }
                    mine = tab->b[i];
                    mini = i;
                } else {
                    if (tab->b[i] > mine && i >= tab->nRows - tab->nEqs) {
                        //this must be infeasible because it means an equality is dominated
                        //cout << "\n6";
                        p->infeasible = true;
                        return;
                    }
                }
            }
            //now that we know the dominant row we can make the new A matrix
            At.push_back(tab->A[mini]);
            bt.push_back(tab->b[mini]);
        }
    }
    if (changed) {
        tab->A = At;
        tab->b = bt;
        tab->nRows = tab->A.size();
        if (tab->nRows > 0)
          tab->nCols = tab->A[0].size();
        else
          tab->nCols = 0;
        //if (R >= tab->nCols - tab->nRows && tab->nEqs > 0)
        if (R >= tab->nRows - tab->nEqs && tab->nEqs > 0)
          tab->nEqs -= 1;
        
        //cout << "\nRedundant Row";
        continue;
    }
  }

  return;
}


void zip(vector<int> &a, vector<long double> &b, vector<pair<int, int>> &zipped) {
  for(size_t i=0; i<a.size(); ++i){
      zipped.push_back(std::make_pair(a[i], b[i]));
  }
}

void unzip(vector<int> &a, vector<long double> &b, vector<pair<int, int>> &zipped) {
  for(size_t i=0; i<a.size(); i++) {
      a[i] = zipped[i].first;
      b[i] = zipped[i].second;
  }
}



int dualSimplex(vector<vector<long double>> A1, vector<vector<long double>> A2, vector<long double> b, vector<long double> c) {
  // A1 should be of the form <=
  // A2 shoud be of the form =

  // return 7 : Integer solution
  // return 6 : reduced to empty
  // return 5 : Artificial variables in the basis
  // return 4 : Method Fails
  // return 3 : found a solution, but it is infeasible
  // return 2 : Got stuck at a point somehow
  // return 1 : found a feasible solution
  // return 0 : infeasible

  //store the original A1 and A2
  auto A1o = A1;
  auto A2o = A2;
  auto bo = b;

  int nEqs = A2.size();

  //combine the 2 constraint matrices
  for (auto row: A2)
    A1.push_back(row);

  //do the presolve
  auto tableau = makeTableau(A1, b, c, nEqs);
  //cout << "\n\nSize of reduced costs: " << tableau->reducedCosts.size() << "\n";
  //cout << "\n\nnCols: " << tableau->nCols << "\n";
  problem *p = new problem{infeasible: false, lb: vector<long double> (tableau->nCols, -inf), ub: vector<long double> (tableau->nCols, inf)};
  
  doPresolve(p, tableau);
  if (p->infeasible){
    //cout << "\nInfeasible on Presolve";
    return 0;
  }

  /*
  cout << "\n\nlb: " << p->lb;
  cout << "\n\nub: " << p->ub;
  cout << "\nA:\n";
  for (auto e: tableau->A)
    cout << e;
  cout << "\nb: " << tableau->b;
  */
  //get the number of columns in A1 and A2 so I can reconstruct the result at the end
  int nVars = tableau->nCols;

  if (tableau->A.size() > 0) {

    /*
    //add A-
    tableau->A = glue(tableau->A, findMinus(tableau->A));
    //add the cost of x-
    int h = tableau->c.size();
    for (int i=0; i<h;i++)
      tableau->c.push_back(0.);

    tableau->nCols = tableau->A[0].size();

    */

    tableau->nCols = tableau->A[0].size();

    for (int i=0; i<tableau->nCols; i++) {
      if (p->lb[i] >= 0.) {
        for (int j=0; j<tableau->A.size(); j++) {
          tableau->b[j] -= tableau->A[j][i]*p->lb[i];
        }
        p->ub[i] -= p->lb[i];
        p->lb[i] = 0.;
      } else if (p->ub[i] <= 0.) {
        for (int j=0; j<tableau->A.size(); j++) {
          tableau->b[j] -= tableau->A[j][i]*p->ub[i];
          tableau->A[j][i] *= -1.;
        }
        p->ub[i] -= p->lb[i];
        p->lb[i] = 0.;
      } else {
        for (int j=0; j<tableau->A.size(); j++) {
          tableau->A[j].push_back(-tableau->A[j][i]);
        }
        p->lb.push_back(0.);
        p->ub.push_back(abs(p->lb[i]));
        p->lb[i] = 0.;
        c.push_back(0.);
      }
    }

    tableau->nCols = tableau->A[0].size();

    //remove any columns with closed bounds
    /*
    bool changed = true;
    while (changed) {
      changed = false;
      auto T = transpose(tableau->A);
      int I = -1;
      if (tableau->A.size() != 0) {
        for (int i=0; i<tableau->nCols; i++) {
          if (p->lb[i] == p->ub[i]) {
            for (int j=0; j<tableau->A.size(); j++) {
              tableau->b[j] -= tableau->A[j][i]*p->lb[i];
            }
            I  = i;
            changed = true;
            break;
          }
          if (changed)
            break;
        }
        if (changed) {
          T.erase(T.begin()+I);
          tableau->A = transpose(T);
          p->ub.erase(p->ub.begin()+I);
          p->lb.erase(p->lb.begin()+I);
          if (tableau->A.size() == 0)
            return 6;
          if (tableau->nRows > 0)
            tableau->nCols = tableau->A[0].size();
          else
            tableau->nCols = 0;
          if (tableau->nRows == 0)
            return 6;
          continue;
        }
      }
    }
    */
    
    //make sure we have enough columns for the number of rows, if not, add artificial variables
    //add the cost of the slack variables
    //c = glueVec(c, vector<long double>(tableau->nRows-tableau->nEqs, 0.))

    //This one includes the slacks on the equalities
    c = glueVec(c, vector<long double>(tableau->nRows+tableau->nEqs, 0.));
    //if (tableau->A.size() <= tableau->A[0].size()) {
        //find the number of rows of zeroes for the slack variables and add them
    //auto X = I(tableau->nRows-tableau->nEqs, tableau->nRows);

    /*
    //get X1
    auto X1 = glue(I(tableau->nRows-tableau->nEqs, tableau->nRows-tableau->nEqs), makeZeroes(tableau->nRows-tableau->nEqs, 2*tableau->nEqs));
    auto X2 = glue(makeZeroes(tableau->nEqs, tableau->nRows-tableau->nEqs), glue(I(tableau->nEqs, tableau->nEqs), findMinus(I(tableau->nEqs, tableau->nEqs))));
    auto X3 = findMinus(X2);

    vector<vector<long double>> X;
    for (auto x: X1) {
      if (x.size() != 0)
        X.push_back(x);
    }
    for (auto x: X2){
      if (x.size() != 0)
        X.push_back(x);
    }
    for (auto x: X3){
      if (x.size() != 0)
        X.push_back(x);
    }

    //Add the extra constraints
    auto Ae = makeZeroes(tableau->nEqs, tableau->nCols);
    for (auto e: Ae)
      tableau->A.push_back(e);
    //update b
    b = glueVec(b, vector<long double> (tableau->nEqs, 0.));
    */

    auto Am = findMinus(tableau->A);
    for (int i=tableau->nRows - tableau->nEqs; i<tableau->nRows; i++) {
      tableau->A.push_back(Am[i]);
      tableau->b.push_back(-tableau->b[i]);
    }

    tableau->nRows = tableau->A.size();

    auto X = I(tableau->A.size(), tableau->A.size());

    /*
    cout << "\n\n\tX\n";
    for (auto e: X)
      cout << e;
    */
    
    tableau->A = glue(tableau->A, X);

    //add the upper and lower bounds for the new slacks
    for (auto q: X) {
      p->lb.push_back(0.);
      p->ub.push_back(inf);
    }

    int basisBegin = tableau->nCols;

    tableau->nCols = tableau->A[0].size() - tableau->nEqs;
    /*} else {
        //add the augmentation vars
        tableau->A = glue(tableau->A, I(tableau->nRows, tableau->nRows));
        tableau->nCols = tableau->A[0].size();
        //add the cost of the artificial variables
        c = glueVec(c, vector<long double>(tableau->nEqs, -10000000.));
    }*/

    //find a new basis
    tableau->basis.clear();
    for (int i = 0; i < tableau->nRows; i++)
        tableau->basis.push_back(basisBegin + i);
        //tableau->basis.push_back((i+tableau->nRows-tableau->nEqs)%tableau->nCols);


    //split x so that it is bound
    //A1 = glue(A1, findMinus(A1));

    //A1 = glue(A1, I(A1.size(), 0));

    //A2 = glue(A2, findMinus(A2));

    //add some zeroes where the slack variables would go
    //A2 = glue(A2, makeZeroes(A2.size(), A1.size()));

    //combine the 2 constraint matrices
    //for (auto row: A2)
    //  A1.push_back(row);

    //add the cost of the augmentation variables
    //c = glueVec(c, vector<long double>(nEqs, -1000000000000));
    //add the coefficients for the slack and augmentation variables
    //A1 = glue(A1, I(A1.size(), 0));

    int leaving = 0;          // this means the row
    int entering = 1;         // this means the column
    //auto tableau = makeTableau(A1, b, c, nEqs);

    //Need to update the dimensions of reduced costs and z
    //tableau->z = vector<long double> (tableau->nCols, 0.);
    //tableau->reducedCosts = vector<long double> (tableau->nCols, 0.);

    //updateZTab(tableau);
    //updateReducedCostsTab(tableau);

    /*
    cout << "\n\nSize of reduced costs: " << tableau->reducedCosts.size() << "\n";
    cout << "\n\nnCols: " << tableau->nCols << "\n";
    */

    int its = 0;

    /*
    cout << "\n\nIteration: " << its << "\n";
    for (auto row: tableau->A)
      cout << row;
    cout << "\nb:";
    cout << tableau->b << "\n";
    cout << "\nDiv: " << tableau->div;
    cout << "\nB: " << tableau->basis << "\n";
  
    cout << "\nlb: " << p->lb;
    cout << "\nub: " << p->ub;
    */

    while (canImproveDual(tableau)) {
      /*
      cout << "\n\nIteration: " << its << "\n";
      for (auto row: tableau->A)
        cout << row;
      cout << "\nb:";
      cout << tableau->b << "\n";
      cout << "\nDiv: " << tableau->div;
      cout << "\nB: " << tableau->basis << "\n";
      */
      /*if (its > 40) {
        return 2;
      }*/
      //get pivot row, this gets the index with the most negative value of b, implementing the first part of blands rule
      vector<long double> v;
      vector<int> indSet;
      for (int i = 0; i < tableau->basis.size(); i++) {
        /*
        if (tableau->basis[i] >= tableau->nCols-tableau->nEqs) {
          indSet.push_back(i);
          v.push_back(-100000000.);
          continue;
        }
        */
        long double theta = tableau->b[i];
        if (theta < 0 && !isInLd(theta, v)) {
          indSet.push_back(i);
          v.push_back(theta);
        }
      }

      if (v.size() == 0) {
        //cout << "\n\n\t2 Exiting the dual simplex here leaving: " << leaving << "\n\n";
        break;
        //return 0;
      }
      //this gets the index of the variable corresponding to the smallest theta
      //int leaving = indSet[distance(begin(v), min_element(begin(v), end(v)))];
      /*
      auto e = *min_element(begin(v), end(v));
      for (int i=0; i<v.size(); i++)
        if (e==v[i])
          leaving = indSet[i];
      */

      //sort v
      //then go through v from smallest to largest to find the best pivot element
      vector<pair<int, int>> zipped;
      zip(indSet, v, zipped);
      std::sort(std::begin(zipped), std::end(zipped), [&](const auto& a, const auto& b) {
            return a.second > b.second;
      });
      unzip(indSet, v, zipped);

      auto u = indSet;

      bool foundPivot = false;
      for (auto l: u) {
        leaving = l;
        //get pivot col
        //they will all have the same reduced cost, so we're always going to choose the entry with the most negative value with the smallest index
        v.clear();
        indSet.clear();
        //cout << "\nBasis: " << tableau->basis;
        for (int i = 0; i < tableau->nCols; i++) {
          if (tableau->A[leaving][i] == 0)
            continue;
          long double e = tableau->A[leaving][i];
          if (!isIn(i, tableau->basis) && e < 0 && !isInLd(e, v) && tableau->b[leaving]/e <= p->ub[i]) {
            indSet.push_back(i);
            v.push_back(e);
            foundPivot = true;
          }
        }
        if (foundPivot)
          break;
      }

      //if (tableau->status=="failure")
      //  break;
      
      //if (tableau->status == "infeasible") {
      //  return 0;
      //}

      if (v.size() == 0) {
        //cout << "\n\n\t1 Exiting the dual simplex here leaving: " << leaving << "\n\n";
        //break;
        return 0;
      }

      //entering = indSet[distance(begin(v), max_element(begin(v), end(v)))];
      auto e = *min_element(begin(v), end(v));
      for (int i=0; i<v.size(); i++)
        if (e==v[i])
          entering = indSet[i];
      /*
      cout << "\nLEAVING: " << leaving << "\n";
      cout << "\nENTERING: " << entering << "\n";
      cout << "\nPivot Element: " << tableau->A[leaving][entering];
      */

      tableau->basis[leaving] = entering;

      //pivot
      pivotDual(tableau, leaving, entering);

      /*
      updateZTab(tableau);
      updateReducedCostsTab(tableau);
      */

      /*
      cout << "\nEndBasis: " << tableau->basis;
      for (auto row: tableau->A)
        cout << row;
      cout << "\nb:";
      cout << tableau->b << "\n";
      */
      its++;
    } 
  } else {
    return 6;
  }
  //if (!canImproveDual(tableau))
  //  cout << "\nCAN'T IMPROVE\n";

  for (int i=0; i<tableau->b.size(); i++) {
    //cout << "\nb" << i << " :\t" << tableau->b[i];
    if (tableau->b[i] < 0)
      return 0;
  }

  for (auto e: tableau->b)
    if (e - floor(e) != 0.)
      return 7;
  
  //make sure reduced costs are all viable
  /*cout << "\nReduced Costs: " << tableau->reducedCosts;
  for (auto e: tableau->reducedCosts)
    if (e < 0)
      return 4;*/
  

  //make sure we don't have any artificial variables in the basis
  /*
  for (auto e: tableau->basis)
    if (e > tableau->nRows - 2*tableau->nEqs)
      return 5;
  */

  //cout << "\nbasis:" << tableau->basis << "\nStart of arti: " << tableau->nCols-A2o.size()-1;

  /*
  
  View output
  for (auto row: tableau->A)
      cout << row;
  cout << "\nb";
  cout << tableau->b;
  cout << "\ndiv";
  cout << tableau->div;
  cout << "\nb after div\n";
  for (int i = 0; i < tableau->nRows; i++)
    cout << get_b(*tableau, i)<< " " ;
  getchar();
  */

  // extract solution
  //auto s = getSplitSolution(tableau);
  //auto r = isPrimalFeasible(p, tableau, s, A1o, A2o, bo, nVars);
  //auto r = optimalityCheck(p, tableau, s);

  return 1;
}

/*
HighsModel makeModel(IP *tc, bool isInt) {
  //make a highs model with A, b taken from cases[id]
  HighsModel model;
  model.lp_.num_col_ = tc->nVars;
  model.lp_.num_row_ = tc->nIneqs + tc->nEqs;
  if (isInt)
    model.lp_.integrality_ = vector<HighsVarType> (tc->nVars, HighsVarType::kInteger);

  vector<double> lb;
  vector<double> ub;
  for (auto row: tc->ineqs) {
    ub.push_back(-(double)row[0]);
    lb.push_back(-1.0e30);
  }
  for (auto row: tc->eqs) {
    lb.push_back(-(double)row[0]);
    ub.push_back(-(double)row[0]);
  }


  //cout << "\nlb: " << lb;
  //cout << "\nub: " << ub;

  /*
  cout << "\nIneqs\n";
  for (auto row: tc->ineqs)
    cout << row;
  */
 /*
  model.lp_.col_lower_ = vector<double> (tc->nVars, -1.0e30);
  model.lp_.col_upper_ = vector<double> (tc->nVars, 1.0e30);
  model.lp_.col_cost_ = vector<double> (tc->nVars, 0.0);
  model.lp_.row_lower_ = lb;
  model.lp_.row_upper_ = ub;

  //count the number of non-zeros in the matrix A
  int nonZeros = 0;
  for (auto row: tc->ineqs){
    int k = 0;
    for (auto e: row) {
      //skip the first column because it is the constant value already taken across to the lower bound
      if (k > 0) {
        if (e != 0)
          nonZeros++;
      }
      k++;
    }
  }
  for (auto row: tc->eqs){
    int k = 0;
    for (auto e: row) {
      //skip the first column because it is the constant value already taken across to the lower bound
      if (k > 0) {
        if (e != 0)
          nonZeros++;
      }
      k++;
    }
  }
//index each non zero element along and down columns
  vector<int> start;
  vector<int> index;
  vector<double> values;

  int c = 0;
  for (int col = 1; col <= tc->nVars; col++) {
    int cN = true;
    for (int row = 0; row < tc->nIneqs; row++)
      if (tc->ineqs[row][col] != 0) {
        index.push_back(row);
        values.push_back((long double)tc->ineqs[row][col]);
        if (cN) {
          start.push_back(c);
          cN = false;
        }
        c++;
      }
    for (int row = 0; row < tc->nEqs; row++)
      if (tc->eqs[row][col] != 0) {
        index.push_back(row+tc->nIneqs);
        values.push_back((double)tc->eqs[row][col]);
        if (cN) {
          start.push_back(c);
          cN = false;
        }
        c++;
      }
    if (cN)
      start.push_back(c);
  }

  start.push_back(nonZeros);

  model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
  model.lp_.a_matrix_.start_ = start;
  model.lp_.a_matrix_.index_ = index;
  model.lp_.a_matrix_.value_ = values;

  return model;
}
*/

/*
int getPresolve(IP *ip, bool isInt) {
  //presolve the model using highs
  HighsModel model = makeModel(ip, isInt);
  Highs highs;
  highs.passModel(model);
  highs.setOptionValue("presolve_rule_logging", true);
  //11110110111111 = 0x3FFF disallows all
  //0x2FFF this is just the aggregator and the basic rules
  //11001110111111
  //00110101000000 = 0xC40
  highs.setOptionValue("presolve_rule_off", 0x3DBF);
  highs.presolve();
  const HighsPresolveStatus& model_status = highs.getModelPresolveStatus();
  ip->presolveCols = model.lp_.num_col_;
  ip->presolveRows = model.lp_.num_row_;

  HighsPresolveLog log = highs.getPresolveLog();
  for (int i=0; i<log.rule.size(); i++)
    if (model_status == HighsPresolveStatus::kReducedToEmpty)
      countPresolveRuleApplications[i] += log.rule[i].call;

  return model_status == HighsPresolveStatus::kReducedToEmpty;
}
*/

int main() {
  clock_t tStart = clock();
  // vector<size_t> hashKeys;
  //vector<int> uniqueProblems;

  //vector<IP> found;

  //int dupCount = 0;
  
  fstream file;
  file.open("uniqueTestCases.txt",ios::in); //open a file to perform read operation using file object
  if (file.is_open()){   //checking whether the file is open
    string elem;

    int nVars = -1;
    int nIneqs = -1;
    vector<string> ineqs;
    int nEqs = -1;
    vector<string> eqs;
    bool readIneqs = false;

    while(getline(file, elem) && recent < nCases) {
      stringstream ss(elem);
      if (nVars < 0) {
        ss >> nVars;
      } else if (!readIneqs) {
        ss >> nIneqs;
        readIneqs = true;
      } else if (nIneqs > 0) {
        ineqs.push_back(elem);
        nIneqs--;
      } else if (nIneqs == 0) {
        nIneqs--;
        ss >> nEqs;
      } else if (nEqs > 0) {
        eqs.push_back(elem);
        nEqs--;
      } else if (nEqs == 0) {
        nEqs--;

        IP *ip = new IP{};
        createIP(ip, nVars, ineqs, eqs);

        
        //hash it
        /*
        ip->hash_code = IPHash{}(*ip);
        if (!ip->searchHash(hashKeys)) {
          hashKeys.push_back(ip->hash_code);
          uniqueProblems.push_back(recent);
        }
        */

        /*
        //comparing without hashing
        bool flag = false;
        for (auto c: found){
          if (IPcomp(c, *ip)) {
            dupCount++;
            flag = true;
            break;
          }
        }
        if (!flag)
          found.push_back(*ip);

        if (recent%1000 == 0)
          cout << recent << "\t" << found.size() << "\t" << dupCount << "\n";
        */
       
        saveTestCase(ip);
        nVars = -1;
        nIneqs = -1;
        nEqs = -1;
        ineqs.clear();
        eqs.clear();
        readIneqs = false;
        delete ip;
      }
    }
    file.close(); 
  }

  //write to a file all the unique cases
  //writeNewFile(uniqueProblems, "uniqueTestCases");
  //cout << hashKeys.size();
  //cout << "\n" << uniqueProblems.size() << "\n";

  //Finds how many are reduced to empty by presolve
  /*
  int both = 0;
  for (int i = 0; i<nCases; i++) {
    bool ic = getPresolve(&cases[i], true);
    bool c = getPresolve(&cases[i], false);
    emptyCountInt += ic;
    emptyCount += c;
    both += ic*c;
  }

  cout  << "\n\n" << "Not Ints reduced to empty:" << "\t\t" << emptyCount << "\n";
  cout  << "\n\n" << "Ints reduced to empty:" << "\t\t" << emptyCountInt << "\n\n";
  cout << "\n\n" << "Both reduced to empty:" << "\t\t" << both;

  //Ints reduced to empty: 2276
  //Not ints reduced to empty 2679
  //Both together at the same time reduced to empty 2260

  cout << "\nPresolve Rule Application Count:\n";
  for (auto e: countPresolveRuleApplications)
    cout << e << "\t";
  cout << "\n";


  cin.clear();
  cin.get();
  */

  /*

  // It looks like this dual is accurate (probably) as it has optimal value 0 for all the test cases

  cout << "\nWorking: " << working;
  cout << "\nFailed: " << notWorking;
  cout << "\nNo Constraints: " << noConstraints;
  cout << "\n";

  */

  /*

  //test that the simplex alg I wrote works
  int working = 0;
  int notWorking = 0;
  int failed = 0;
  for (int ind=0; ind<nCases; ind++){
    cout << "\n \t IND:: " << ind << "\n";
    cout << "\n";
    cout << cases[ind].nVars << "\n";
    cout << cases[ind].nIneqs << "\n";
    cout << cases[ind].ineqs << "\n";
    cout << cases[ind].nEqs << "\n";
    cout << cases[ind].eqs << "\n";

    //there are cases with no variables, only a constant term, this will skip over those cases
    if (cases[ind].nVars == 1)
      continue;

    if (cases[ind].nEqs==0 && cases[ind].nIneqs==0)
      continue;

    auto dual = findDual(ind);
    cout << "this far";

    IP tc = cases[ind];
    //find transposes of the 2 parts of the A matrix
    vector<vector<long double>> A1, A2;
    vector<long double> b, b2;
    tie(A1, b) = extractAb(tc.ineqs);
    tie(A2, b2) = extractAb(tc.eqs);
    auto A1t = transpose(A1);
    auto A2t = transpose(A2);
    //find the new A matrix
    
    auto row1 = glue(A1, glue(findMinus(A1), I(0, A1.size())));
    auto row2 = glue(A2, glue(findMinus(A2), findMinus(I(A2.size(), 0))));
    for (auto row: row2)
      row1.push_back(row);
    auto A = glue(transpose(row1), findMinus(I(row1[0].size(), 0)));

    //make b vector
    b.insert(b.end(), b2.begin(), b2.end());

    //convert b to long double, and change the problem to a maximisation problem
    vector<long double> bd;
    for (auto e: b)
      bd.push_back(-(long double)e);
    //convert A to long double
    vector<vector<long double>> Ad;
    for (auto row: A) {
      vector<long double> newRow;
      for (auto e: row) {
        newRow.push_back((long double)e);
      }
      Ad.push_back(newRow);
    }
    while (bd.size() <= Ad[0].size()) {
      bd.push_back(0);
    }
    if (Ad.size() > Ad[0].size())
      continue;
    auto s = simplex(bd, Ad, vector<long double> (Ad.size(), 0), 0);
    if (info.objective_function_value == s)
      working++;
    else
      notWorking++;
    if (s > 0)
      failed++;
    cout << "\n WORKING: " << working << "\t" << working + notWorking << "\t" << failed << "\n";
  }

  // It looks like this algorithm is accurate, at least when it works
  // sometimes we get segfaults, which I can't explain yet
  // There appears to be a pattern in how many iterations it takes to solve these problems
  // It would be interesting to see where these problems are coming from
  // If I could then I could maybe find why they seem to have similar behaviour in the simplex alg

  cout << "\nWorking: " << working;
  cout << "\nFailed: " << notWorking;

  /*
  vector<vector<long double>> A = {{2., 4., 1., 0.}, {3., 2., 0., 1.}};
  vector<long double> b = {16., 12.};
  vector<long double> c = {7., 6., 0., 0., 0.};
  cout << "\n\nIS IT WORKING? SHOULD BE 32!!!    " << simplex(c, A, b, 2) << "\n\n";
  */
  
  //Test the dual simplex
  int infeasible = 0;
  int working = 0;
  int incorrect = 0;
  int failed = 0;
  int methodFails = 0;
  int artInBasis = 0;
  int trivFea = 0;
  int trivInf = 0;
  int highsInf = 0;
  int highsFea = 0;
  int disagree = 0;
  int gap = 0;
  int empty = 0;
  int mineEmpty = 0;
  int integer = 0;
  for (int ind=0; ind<nCases; ind++){
    /*
    cout << "\n \t IND:: " << ind << "\n";
    cout << "\n";
    cout << cases[ind].nVars << "\n";
    cout << cases[ind].nIneqs << "\n";
    cout << cases[ind].ineqs << "\n";
    cout << cases[ind].nEqs << "\n";
    cout << cases[ind].eqs << "\n";
    */

    // See what the highs solver finds
    /*
    auto primal = makeModel(&cases[ind], false);
    Highs highs;
    HighsStatus return_status;
    return_status = highs.passModel(primal);
    assert(return_status==HighsStatus::kOk);
    const HighsLp& lp = highs.getLp();
    return_status = highs.run();
    assert(return_status==HighsStatus::kOk);
    const HighsModelStatus& model_status = highs.getModelStatus();
    //assert(model_status==HighsModelStatus::kUnbounded);
    const HighsInfo& info = highs.getInfo();
    if (model_status == HighsModelStatus::kInfeasible)
      highsInf++;
    if (model_status == HighsModelStatus::kOptimal || model_status == HighsModelStatus::kModelEmpty)
      highsFea++;
    if (model_status == HighsModelStatus::kModelEmpty)
      empty++;
    */

    /*  code to count the trivially feasible cases
    bool triviallyFeasible = true;
    //check if inequalities are valid
    for (int i=1; i<cases[ind].nIneqs; i++)
      if (0 > -cases[ind].ineqs[i][0])
        triviallyFeasible = false;
    //check if the equalities are valid
    long double eps = 0.000001;
    for (int i=0; i<cases[ind].nEqs; i++)
      if (0 < -cases[ind].eqs[i][0] - eps ||  0 > -cases[ind].eqs[i][0] + eps)
        triviallyFeasible = false;
    if (triviallyFeasible) {
      working++;
      trivFea++;
      continue;
    }
    */

   /*  code to count the trivially infeasible cases

    bool triviallyInfeasible = false;
    for (auto eq: cases[ind].ineqs) {
      bool nonZero = false;
      if (eq[0] > 0) {
        for (int i=1; i<cases[ind].nVars; i++) {
          if (eq[i] != 0)
            nonZero = true;
        }
        if (!nonZero) {
          triviallyInfeasible = true;
          break;
        }
      }
    }
    for (auto eq: cases[ind].eqs) {
      bool nonZero = false;
      if (eq[0] != 0) {
        for (int i=1; i<cases[ind].nVars; i++) {
          if (eq[i] != 0)
            nonZero = true;
        }
        if (!nonZero) {
          triviallyInfeasible = true;
          break;
        }
      }
    }


   if (triviallyInfeasible) {
      infeasible++;
      trivInf++;
      continue;
   }
   */

    if (cases[ind].nVars == 1) {
      infeasible++;
      trivInf++;
      continue;
    }

    if (cases[ind].nEqs==0 && cases[ind].nIneqs==0) {
      working++;
      continue;
    }

    vector<vector<long double>> A1;
    vector<vector<long double>> A2;
    vector<long double> b;
    vector<long double> b2;
    tie(A1, b) = extractAb(cases[ind].ineqs);
    tie(A2, b2) = extractAb(cases[ind].eqs);
    b.insert(b.end(), b2.begin(), b2.end()); 

    auto x = dualSimplex(A1, A2, b, vector<long double>(cases[ind].nVars, 0.));
    switch (x) {
    case (0):
      infeasible++;
    break;
    case (1):
      working++;
    break;
    case (2):
      failed++;
    break;
    case (3):
      incorrect++;
    break;
    case (4):
      methodFails++;
    break;
    case (5):
      artInBasis++;
    break;
    case (6):
      mineEmpty++;
    break;
    case (7):
      integer++;
    break;
    }
    

    /*
    cout << "\n\ninfeasible: " << infeasible << "\tHighs Infeasible: " << highsInf;
    cout << "\nWorking: " << working << "\tHighs Solved to Optimality: " << highsFea;
    cout << "\nIncorrect: " << incorrect << "\tHighs reduced to Empty: " << empty;
    cout << "\nMethod Fails: " << methodFails;
    cout << "\nArtificial Variable in the Basis: " << artInBasis;
    cout << "\nfailed: " << failed << "\n";
    cout << "\nEmpty: " << mineEmpty + trivFea << "\n";

    cout << "\n\nHighs Disagrees With my dual simplex: " << disagree << " times\n";

    cout << "\nTrivially Feasible: " << trivFea << "\n";
    */

    /*
    bool c = getPresolve(&cases[ind], false);
    if (x == 6 && !c && cases[ind].nEqs + cases[ind].nIneqs > 1) {
      cout << "\n\nGOES WRONG SOMEWHERE highs: nonEmpty, mine: empty\n\n";
      //cin.clear();
      //cin.get();
        if (!((model_status == HighsModelStatus::kOptimal || model_status == HighsModelStatus::kModelEmpty))) {
          cout << "\n\nGOES WRONG SOMEWHERE highs: optimal, mine: infeasible\n\n";
          cin.clear();
          cin.get();
          disagree++;
      }
    }
    */
    /*
    if ((model_status == HighsModelStatus::kOptimal || model_status == HighsModelStatus::kModelEmpty) && x == 0) {
        cout << "\n\nGOES WRONG SOMEWHERE highs: optimal, mine: infeasible\n\n";
        cin.clear();
        cin.get();
        cin.get();
        cin.clear();
        disagree++;
    }
    */
    /* 
    if ((model_status == HighsModelStatus::kInfeasible && x == 1) || (model_status == HighsModelStatus::kInfeasible && x == 6)) {
        cout << "\n\nGOES WRONG SOMEWHERE highs: infeasible, mine: feasible\n\n";
        cin.clear();
        cin.get();
        disagree++;
    }
    */
    
    //cin.clear();
    //cin.get();
    //cin.clear();

  }


  //cout  << "\n\n" << "Not Ints reduced to empty:" << "\t\t" << emptyCount << "\n";
  //cout  << "\n\n" << "Ints reduced to empty:" << "\t\t" << emptyCountInt << "\n\n";
  //cout << "\n\n" << "Both reduced to empty:" << "\t\t" << both;

  cout << "\n\ninfeasible: " << infeasible << "\tHighs Infeasible: " << highsInf;
  cout << "\nWorking: " << working << "\tHighs Solved to Optimality: " << highsFea;
  cout << "\nIncorrect: " << incorrect << "\tHighs reduced to Empty: " << empty;
  cout << "\nMethod Fails: " << methodFails;
  cout << "\nArtificial Variable in the Basis: " << artInBasis;
  cout << "\nfailed: " << failed << "\n";
  cout << "\nEmpty: " << mineEmpty + trivFea << "\n";
  cout << "\nInteger: " << integer << "\n";

  cout << "\n\nHighs Disagrees With my dual simplex: " << disagree << " times\n";

  cout << "\nTrivially Feasible: " << trivFea << "\n";

  cout << "\nTIME\n";
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);


  /*
  // Test case for the dual simplex, it works
  vector<vector<long double>> A = {{-3, -5}, {-5, -2}};
  vector<vector<long double>> A2;
  vector<long double> b = {-5, -3};
  int x = dualSimplex(A, A2, b);
  cout << "\nOUTPUT: " << x << "\n";

  */

  return 0;
}


