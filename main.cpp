
#include "Highs.h"
#include "lp_data/HighsModelUtils.h"
#include <fstream>
#include <iostream>
#include <sstream> 
#include <tuple>
#include <bitset>

using namespace std;


//non-ints reduced to empty : 55013
//                            54923
//                            54594
//                            52204   01111101000000 = 0x1F40 This is within 1% (actually 0.9073307771%), so I consider it acceptable
//ints reduced to empty :     52682
//it looks like it reduces more to empty when we don't consider the integers
//but not many

/*

Let's say anything more than 1.5% of the cases that reduce to empty is too much to lose
as we are already losing around this much from ignoring the integer property
so we will not implement anything more costly than 1%

Presolve Rule Application Count:
2446    375516  709     426345  6803    71674   66      16890   1011    44438   43      0       60281   55486


keep    Keep    keep?   keep    keep    keep    lose    keep    keep?   keep    lose    lose    keep    keep
1,      1,      1,      1,      1,      1,      0,      1,      1,      1,      0,      0,      1,      1
0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0

0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0
I will now test this rule configuration, and see how much difference the questionable rules make

making this change reduces to 55013, which is the exact same


Total ints reduced to empty: 58344
                             56937  2680
                             56937  2EE0

7828    414446  1107    471118  7941    60222   0       20404   48556   0       0       0       69314   0
7658    87306   158     429134  1281    5655    0       14212   14546   0       0       0       64918   0
7858    423231  1111    460229  8185    68641   0       1394    49184   0       0       0       69506   0
0       1       2       3       4       5       6       7       8       9       10      11      12      13
0       0       0       0       0       0       1       0       0       1       1       1       0       1 
          
this seems to be the best outcome with the mimimal rules, so I'll implement this one
We're only allowing 0-5, 8 and 12


With all the rules, non-ints reduced to empty: 62456
2, 3, 7


There are 37 times that my alg finds a solution and HiGHS does not
Investigate this

The 2 algs dissagree 300 times, roughly 5% of the time
Investigate this


*/

const int nCases = 5980; // <- unique set, full set -> 150218;
// #duplicates
//                 145564    hashing    4624
//                 144238    IPcomp     5980
//                 144238    better hashing 5980, also much more efficient

//Not very many is it

/*

Of these 5980 2301 reduce to empty if you relax the integer constraint
it would be interesting to see what the frequency of each problem is within the set

*/


struct IP {

  int nVars;
  int nIneqs;
  int nEqs;
  vector<vector<HighsInt>> ineqs;
  vector<vector<HighsInt>> eqs;
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
    cout << "\na == 0";
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

struct postsolveStep {
  string type;
  vector<long double> row;
  long double b;
  bool eq; //true means equality, false means inequality
  int columnIndex;
};


//DOUBLE CHECK THIS PART
void pivotRow(tableauStruct *tab, int pivotRow, int pivotCol, int row) {
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
    tab->A[row][i] -= q*tab->A[pivotRow][i];
  tab->b[row]-= q*tab->b[pivotRow];
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

vector<vector<HighsInt>> findMinus(vector<vector<HighsInt>> A) {
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
    vector<HighsInt> constr;
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

  cout << "\nlb: " << lb;
  cout << "\nub: " << ub;

  cout << "\nIneqs\n";
  for (auto row: tc->ineqs)
    cout << row;

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

int getPresolve(IP *ip, bool isInt) {
  //presolve the model using highs
  HighsModel model = makeModel(ip, isInt);
  Highs highs;
  highs.passModel(model);
  highs.setOptionValue("presolve_rule_logging", true);
  //11111111111111 = 0x3FFF disallows all
  //11001010101011
  //11110111010100
  //11001110111111
  //00110101000000 = 0xC40
  highs.setOptionValue("presolve_rule_off", 0x3DD4);
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

void writeNewFile(vector<int> dataPos, string name) {
  ofstream newFile;
  newFile.open (name);
  for (auto uId: dataPos) {
    auto u = cases[uId];
    newFile << "~~~~~\n";
    newFile << to_string(u.nVars) << "\n";
    newFile << to_string(u.nIneqs) << "\n";
    for (auto row: u.ineqs) {
      for (auto e: row)
        newFile << to_string(e) << " ";
      newFile << "\n";
    }
    newFile << to_string(u.nEqs) << "\n";
    for (auto row: u.eqs) {
      for (auto e: row)
        newFile << to_string(e) << " ";
      newFile << "\n";
    }
  }
  newFile.close();
}

/*
vector<vector<long double>> transpose(vector<vector<long double>> A) {
  vector<long double> r(A.size(), 0);
  vector<vector<long double>> At(A[0].size(), r);
  for (int row=0; row<A.size(); row++){
    for (int col=0; col<A[0].size(); col++) {
      At[col][row] = A[row][col];
    }
  }
  return At;
}
*/


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
  }
  return A;
}

vector<long double> glueVec(vector<long double> A, vector<long double> B) {
  for (auto e: B)
    A.push_back(e);
  return A;
}

vector<vector<long double>> I(int n, int zeroes) {
  vector<vector<long double>> A;
  for (int row=0; row<zeroes; row++)
    A.push_back(vector<long double> (n+zeroes, 0));
  for (int row=0; row<n; row++) {
    vector<long double> r;
    for (int col=0; col<n; col++)
        r.push_back((int)(row==col));
    A.push_back(r);
  }    
  return A;
}

tuple<vector<vector<long double>>, vector<long double>> extractAb(vector<vector<HighsInt>> p) {
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

HighsModel findDual(int ind) {
  IP tc = cases[ind];

  vector<vector<long double>> A1, A2;
  vector<long double> b, b2;
  tie(A1, b) = extractAb(tc.ineqs);
  tie(A2, b2) = extractAb(tc.eqs);
  //find the new A matrix
  auto row1 = glue(A1, glue(findMinus(A1), I(0, A1.size())));
  auto row2 = glue(A2, glue(findMinus(A2), findMinus(I(A2.size(), 0))));
  for (auto row: row2)
    row1.push_back(row);
  auto A = glue(transpose(row1), findMinus(I(row1[0].size(), 0)));

  //make b vector
  b.insert(b.end(), b2.begin(), b2.end());

  //convert b to long double
  vector<double> bd;
  for (auto e: b)
    bd.push_back((double)e);

  HighsModel model;
  model.lp_.num_col_ = b.size();
  model.lp_.num_row_ = A.size();

  vector<double> lb;
  vector<double> ub;
  for (auto row: A) {
    lb.push_back(0.);
    ub.push_back(1.0e30);
  }

  model.lp_.col_lower_ = vector<double> (b.size(), -1.0e30);
  model.lp_.col_upper_ = vector<double> (b.size(), 1.0e30);
  model.lp_.col_cost_ = bd;
  model.lp_.row_lower_ = lb;
  model.lp_.row_upper_ = ub;

  //count the number of non-zeros in the matrix A
  int nonZeros = 0;
  for (auto row: A){
    for (auto e: row) {
      if (e != 0)
        nonZeros++;
    }
  }

  //index each non zero element along and down columns
  vector<int> start;
  vector<int> index;
  vector<double> values;

  int c = 0;
  for (int col = 0; col <= b.size(); col++) {
    int cN = true;
    for (int row = 0; row < A.size(); row++)
      if (A[row][col] != 0) {
        index.push_back(row);
        values.push_back((double)A[row][col]);
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

  cout << "\ndr: " << dr;

  cout << "\ncol: " << col << "\n";
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
  cout << "\nRow length: " << rows.size();
  row = 100000000;
  for (int i = 0; i < rows.size(); i++)
    if (rows[i] < row) {
      row = rows[i];
    }
  cout << "\nrow: " << row << "\n";

  //I suspect this means that it is infeasible????
  if (row == 100000000)
    return make_tuple(make_tuple(row, col), true);

  return make_tuple(make_tuple(row, col), false);
}

long double get_b(tableauStruct tab, int i) {
  if (i < tab.nRows)
    return (long double)tab.b[i]/(long double)tab.div[i];
  cout << "\nIndex too large";
  return 0.;
}

vector<vector<long double>> pivot(vector<vector<long double>> tableau, tuple<int, int> p) {
  int row, col;
  tie(row, col) = p;

  auto v = tableau[row][col];

  cout << "\nv: " << v;
  cin.clear();
  if (v==0.) {
    cout << "\nv is 0";
    cin.get();
  }

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
  for (auto e: tab->basis) {
    if (e >= tab->nCols - tab->nEqs)
      return true;
  }

  //check that the reduced costs are all greater than or equal to 0
  for (auto e: tab->reducedCosts) {
    if (e < 0) {
      tab->status = "infeasible";
      return false;
    }
  }

  //cout << "\nThis far";
  bool flag2 = false;
  for (int i = 0; i < tab->nRows; i++) {
    if (get_b(*tab, i) < 0.) {
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

  cout << "\n B: " << B;
  cout << "\n # of vars " << A[0].size();
  cout << "\n # of constraints " << A.size();
  cout << "\nLength of cost " << cost.size();

  auto z = updateZ(cost, B, tableau);
  auto r = updateR(z, cost);

  cout << "\nz: " << z;
  cout << "\ncost: " << cost;

  cout << "\n INITIAL R: " << r.size();
  cout << "\n INITIAL R: " << r;

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
    if (tableau->A[i][tableau->basis[i]] == 0)
      splitSolution[tableau->basis[i]] = 0;
    else
      splitSolution[tableau->basis[i]] = get_b(*tableau, i)/tableau->A[i][tableau->basis[i]];
  }
  return splitSolution;
}

int isPrimalFeasible(vector<long double> splitSolution, vector<vector<long double>> A1, vector<vector<long double>> A2, vector<long double> b, int nVars) {
  
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

  for (int i=0; i<nVars; i++)
    x.push_back(splitSolution[i] - splitSolution[i+nVars]);

  cout << "\nSplit Solution: " << splitSolution;
  cout << "\nnVars: " << nVars << "\tx:\n";
  cout << x;

  //check if it is a primal feasible solution
  vector<long double> v1, v2;
  for (auto row: A1) {
    long double e = 0;
    for (int i=0; i<nVars; i++) {
      e += (long double)row[i]*x[i];
    }
    v1.push_back(e);
  }
  for (auto row: A2) {
    long double e = 0;
    for (int i=0; i<nVars; i++) {
      e += (long double)row[i]*x[i];
    }
    v2.push_back(e);
  }

  cout << "\nSolutions\n";
  cout << x;
  cout << "\nSolution Values\n";
  cout << v1 << v2;
  cout << "\nbounds\n";
  cout << b;

  //check if inequalities are valid
  for (int i=0; i<v1.size(); i++)
    if (v1[i] > b[i]) {
      cout << "\nWHERE IT FAILS: " << v1[i] << "  >  " << b[i];
      return 0;
    }
  //check if the equalities are valid
  long double eps = 0.01;
  for (int i=v1.size(); i<b.size(); i++)
    if (v2[i-v1.size()] < b[i] - eps ||  v2[i-v1.size()] > b[i] + eps){
      cout << "\nWHERE IT FAILS: " << v2[i-v1.size()] << "  >  " << b[i];
      return 0;
    }
  
  //delete(tableau);
  return 1;
}


vector<postsolveStep> presolveStack;

tuple<tuple<tuple<vector<vector<long double>>, vector<vector<long double>>>, vector<long double>>, bool> doPresolve(vector<vector<long double>> A1, vector<vector<long double>> A2, vector<long double> b) {
  // return (((A1, A2), b), could still be feasible)
  vector<vector<long double>> A1t;
  vector<vector<long double>> A2t;
  vector<long double> bt;

  //Rule 1
  //empty column
  auto T1 = transpose(A1);
  auto T2 = transpose(A2);
  vector<vector<long double>> T1t;
  for (int r=0; r<A1.size(); r++) {
    auto row = A1[r];
    for (auto e: row) {
      if (e != 0) {
        T1t.push_back(row);
        bt.push_back(b[r]);
        break;
      }
    }
  }
  vector<vector<long double>> T2t;
  for (int r=0; r<A2.size(); r++) {
    auto row = A2[r];
    for (auto e: row) {
      if (e != 0) {
        T2t.push_back(row);
        bt.push_back(b[r+A1.size()]);
        break;
      }
    }
  }
  A1 = transpose(T1t);
  A2 = transpose(T2t);
  b = bt;

  //singleton row
  /*
    Where I don't know what to do for this one, maybe remove the constraint here, then add it back in when it comes time 
    to solve it. This makes sense. It means I'd have to pass the original problem to the dual simplex, then after the 
    dual simplex receives the problem and remembers it, then I'll apply the presolve
  */
  for (int r=0; r<A1.size(); r++) {
    auto row = A1[r];
    bool flag = false;
    int I = 0;
    for (int i=0; i<row.size(); i++) {
      auto e = row[i];
      if (e != 0) {
        if (flag) {
          flag = false;
          I = i;
          break;
        }
        flag = true;
      }
    }
    if (flag) {
      //we have a singleton row. What to put here?
      //as this is in A1 it is unbounded
      //I should calculate the value of x
      long double x = b[r]/A1[r][I];
      //apply this as an upper bound on x[I], if I decide to implement dominated columns
    }
  }
  for (int r=0; r<A2.size(); r++) {
    auto row = A2[r];
    bool flag = false;
    int I = 0;
    for (int i=0; i<row.size(); i++) {
      auto e = row[i];
      if (e != 0) {
        if (flag) {
          flag = false;
          I = i;
          break;
        }
        flag = true;
      }
    }
    if (flag) {
      //we have a singleton row. What to put here?
      //as this is in A2 we have a fixed variable in x
      long double x = b[r]/A1[r][I];
      if (x - floor(x) != 0.) {
        return make_tuple(make_tuple(make_tuple(A1t, A2t), bt), false);
      }
    }
  }

  //Rule 2
  //fixed column
  T1 = transpose(A1);
  T2 = transpose(A2);
  T1t.clear();
  T2t.clear();
  bt.clear();
  int I = -1;
  for (int r=0; r<T1.size(); r++) {
    auto row = T1[r];
    bool flag = false;
    for (int i=0; i<row.size(); i++) {
      auto e = row[i];
      if (e != 0) {
        if (flag) {
          flag = false;
          I = i;
          break;
        }
        flag = true;
      }
    }
    if (!flag) {
      T1t.push_back(row);
    } else {
      // this is where we have a column singleton
      A1[I][r] = 0.;
      presolveStack.push_back(postsolveStep{type: "Column Singleton", row: A1[I], b: b[I], eq: false, columnIndex: r});
    }
  }
  if (I > 0) {
    //Now row I should be removed
    for (int i=0; i<A1.size(); i++) {
      if (i != I) {
        A1t.push_back(A1[i]);
        bt.push_back(b[i]);
      }
    }
    A1 = A1t;
    b = bt;
  }
  I = -1;
  for (int r=0; r<T2.size(); r++) {
    auto row = T2[r];
    bool flag = false;
    for (int i=0; i<row.size(); i++) {
      auto e = row[i];
      if (e != 0) {
        if (flag) {
          flag = false;
          I = i;
          break;
        }
        flag = true;
      }
    }
    if (!flag) {
      T2t.push_back(row);
    } else {
      // this is where we have a column singleton
      A2[I][r] = 0.;
      presolveStack.push_back(postsolveStep{type: "Column Singleton", row: A1[I], b: b[I+A1.size()], eq: true, columnIndex: r});
    }
  }
  A1 = transpose(T1t);
  A2 = transpose(T2t);
  if (I > 0) {
    A2t.clear();
    bt.clear();
    //Now row I should be removed
    for (int i=0; i<A2.size(); i++) {
      if (i != I) {
        A2t.push_back(A2[i]);
        bt.push_back(b[i+A1.size()]);
      }
    }
    A2 = A2t;
    b = bt;
  }

  //remove all empty rows
  A1t.clear();
  bt.clear();
  for (int r=0; r<A1.size(); r++) {
    auto row = A1[r];
    for (auto e: row) {
      if (e != 0) {
        A1t.push_back(row);
        bt.push_back(b[r]);
        break;
      }
    }
  }
  A2t.clear();
  for (int r=0; r<A2.size(); r++) {
    auto row = A2[r];
    for (auto e: row) {
      if (e != 0) {
        A2t.push_back(row);
        bt.push_back(b[r+A1.size()]);
        break;
      }
    }
  }
  A1 = A1t;
  A2 = A2t;
  b = bt;

  //redundant row
  A1t.clear();
  A2t.clear();
  bt.clear();
  for (int r1=0; r1<A1.size(); r1++) {
    auto row1 = A1[r1];
    bool add = true;
    for (int r2=0; r2<A1.size(); r2++) {
      auto row2 = A1[r2];
      if (row1 == row2 && r1 != r2) {
        if (b[r1] < b[r2]) {
          A1t.push_back(row1);
          bt.push_back(b[r1]);
        } else if (b[r1] == b[r2]) {
          if (r1 < r2) {
            A1t.push_back(row1);
            bt.push_back(r1);
          }
        }
      }
    }
  }
  for (int r1=0; r1<A2.size(); r1++) {
    auto row1 = A2[r1];
    bool add = true;
    for (int r2=0; r2<A2.size(); r2++) {
      auto row2 = A2[r2];
      if (row1 == row2 && r1 != r2) {
        if (b[r1] < b[r2]) {
          A1t.push_back(row1);
          bt.push_back(b[r1]);
        } else if (b[r1] == b[r2]) {
          if (r1 < r2) {
            A1t.push_back(row1);
            bt.push_back(r1);
          }
        }
      }
    }
  }
  A1 = A1t;
  A2 = A2t;
  b = bt;

  return make_tuple(make_tuple(make_tuple(A1t, A2t), bt), true);
}

int dualSimplex(vector<vector<long double>> A1, vector<vector<long double>> A2, vector<long double> b, vector<long double> c) {
  // A1 should be of the form <=
  // A2 shoud be of the form =

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

  
  //get the number of columns in A1 and A2 so I can reconstruct the result at the end
  int nVars;
  if (A1.size() > 0)
    nVars = A1[0].size();
  else if (A2.size() > 0)
    nVars = A2[0].size();
  else
    return 1;

  //split x so that it is bound
  A1 = glue(A1, findMinus(A1));

  //A1 = glue(A1, I(A1.size(), 0));


  A2 = glue(A2, findMinus(A2));

  //add some zeroes where the slack variables would go
  //A2 = glue(A2, makeZeroes(A2.size(), A1.size()));

  int nEqs = A2.size();

  //add the cost of the slack variables
  c = glueVec(c, vector<long double>(A1.size(), 0));

  //combine the 2 constraint matrices
  for (auto row: A2)
    A1.push_back(row);

  //add the cost of the augmentation variables
  c = glueVec(c, vector<long double>(nEqs, -1000000000000));
  //add the coefficients for the slack and augmentation variables
  A1 = glue(A1, I(A1.size(), 0));

  int leaving = 0;          // this means the row
  int entering = 1;         // this means the column
  auto tableau = makeTableau(A1, b, c, nEqs);

  updateZTab(tableau);
  updateReducedCostsTab(tableau);

  int its = 0;
  /*
  cout << "\n\nIteration: " << its << "\n";
  for (auto row: tableau->A)
    cout << row;
  cout << "\nb:";
  cout << tableau->b << "\n";
  cout << "\nDiv: " << tableau->div;
  cout << "\nB: " << tableau->basis << "\n";
  */

  while (canImproveDual(tableau)) {
    cout << "\n\nIteration: " << its << "\n";
    for (auto row: tableau->A)
      cout << row;
    cout << "\nb:";
    cout << tableau->b << "\n";
    cout << "\nDiv: " << tableau->div;
    cout << "\nB: " << tableau->basis << "\n";
    if (its > 40) {
      return 2;
    }
    //get pivot row, this gets the index with the most negative value of b, implementing the first part of blands rule
    vector<long double> v;
    vector<int> indSet;
    for (int i = 0; i < tableau->basis.size(); i++) {
      if (tableau->basis[i] >= tableau->nCols-tableau->nEqs) {
        indSet.push_back(i);
        v.push_back(-1000000.);
      }
      long double theta = get_b(*tableau, i);
      if (theta < 0 && !isInLd(theta, v)) {
        indSet.push_back(i);
        v.push_back(theta);
      }
    }

    if (v.size() == 0) {
      break;
    }
    //this gets the index of the variable corresponding to the smallest theta
    //int leaving = indSet[distance(begin(v), min_element(begin(v), end(v)))];
    auto e = *min_element(begin(v), end(v));
    for (int i=0; i<v.size(); i++)
      if (e==v[i])
        leaving = indSet[i];

    //get pivot col
    //they will all have the same reduced cost, so we're always going to choose the entry with the most negative value with the smallest index
    v.clear();
    indSet.clear();
    //cout << "\nBasis: " << tableau->basis;
    for (int i = 0; i < tableau->nCols-nEqs; i++) {
      if (tableau->A[leaving][i] == 0)
        continue;
      long double e = tableau->A[leaving][i];
      if (!isIn(i, tableau->basis) && e < 0 && !isInLd(e, v)) {
        indSet.push_back(i);
        v.push_back(e);
      }
    }

    if (tableau->status=="failure")
      break;
    
    if (tableau->status == "infeasible")
      return 0;

    if (v.size() == 0) {
      break;
    }

    //entering = indSet[distance(begin(v), max_element(begin(v), end(v)))];
    e = *max_element(begin(v), end(v));
    for (int i=0; i<v.size(); i++)
      if (e==v[i])
        entering = indSet[i];

    //cout << "\nLEAVING: " << leaving << "\n";
    //cout << "\nENTERING: " << entering << "\n";
    //cout << "\nPivot Element: " << tableau->A[leaving][entering];

    tableau->basis[leaving] = entering;

    //pivot
    pivotDual(tableau, leaving, entering);

    updateZTab(tableau);
    updateReducedCostsTab(tableau);

    its++;
  } 
  //if (!canImproveDual(tableau))
  //  cout << "\nCAN'T IMPROVE\n";
  cout << "\n\nIteration: " << its << "\n";
  for (auto row: tableau->A)
    cout << row;
  cout << "\nb:";
  cout << tableau->b << "\n";
  cout << "\nDiv: " << tableau->div;
  cout << "\nB: " << tableau->basis << "\n";
  cout << "\nStatus: " << tableau->status;
  /*
  for (int i=0; i<tableau->nRows; i++)
    if (get_b(*tableau, i) < 0)
      return 4;
  */

  /*
  //make sure reduced costs are all viable
  for (auto e: tableau->reducedCosts)
    if (e < 0)
      return 4;
  */

  //make sure we don't have any artificial variables in the basis
  /*
  for (auto e: tableau->basis)
    if (e > nEqs-1)
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
  return isPrimalFeasible(getSplitSolution(tableau), A1o, A2o, bo, nVars);
}

int main() {

  vector<size_t> hashKeys;
  vector<int> uniqueProblems;

  vector<IP> found;

  int dupCount = 0;
  
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
      }
    }
    file.close(); 
  }

  //write to a file all the unique cases
  //writeNewFile(uniqueProblems, "uniqueTestCases");
  //cout << hashKeys.size();
  //cout << "\n" << uniqueProblems.size() << "\n";

  //Finds how many are reduced to empty by presolve
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

  /*

  // Test the duality
  // 1. find the optimal value of the dual
  // 2. if its <= 0 for all of the cases then the dual is probably ok
  int working = 0;
  int notWorking = 0;
  int noConstraints = 0;
  for (int ind=0; ind<nCases; ind++){
    
    printCase(ind);

    if (cases[ind].nIneqs == 0 && cases[ind].nEqs == 0) {
      noConstraints++;
      continue;
    }

    auto dual = findDual(ind);
    Highs highs;
    HighsStatus return_status;
    return_status = highs.passModel(dual);
    assert(return_status==HighsStatus::kOk);
    const HighsLp& lp = highs.getLp();
    return_status = highs.run();
    assert(return_status==HighsStatus::kOk);
    const HighsModelStatus& model_status = highs.getModelStatus();
    //assert(model_status==HighsModelStatus::kUnbounded);
    const HighsInfo& info = highs.getInfo();
    cout << "\n\n" << highs.solutionStatusToString(info.primal_solution_status) << "\n";
    if (info.objective_function_value == 0)
      working++;
    else
      notWorking++;
    cout << "\n" << working << "\t" << notWorking << "\n";
  }

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

    Highs highs;
    HighsStatus return_status;
    return_status = highs.passModel(dual);
    assert(return_status==HighsStatus::kOk);
    const HighsLp& lp = highs.getLp();
    return_status = highs.run();
    assert(return_status==HighsStatus::kOk);
    const HighsModelStatus& model_status = highs.getModelStatus();
    //assert(model_status==HighsModelStatus::kUnbounded);
    const HighsInfo& info = highs.getInfo();
    
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
  
  /*
  //Test the dual simplex
  int infeasible = 0;
  int working = 0;
  int incorrect = 0;
  int failed = 0;
  int methodFails = 0;
  int artInBasis = 0;
  int trivFea = 0;
  int highsInf = 0;
  int highsFea = 0;
  int disagree = 0;
  int gap =0;
  for (int ind=0; ind<nCases; ind++){
    cout << "\n \t IND:: " << ind << "\n";

    cout << "\n \t IND:: " << ind << "\n";
    cout << "\n";
    cout << cases[ind].nVars << "\n";
    cout << cases[ind].nIneqs << "\n";
    cout << cases[ind].ineqs << "\n";
    cout << cases[ind].nEqs << "\n";
    cout << cases[ind].eqs << "\n";
    

    // See what the highs solver finds
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
    

    bool triviallyFeasible = true;
    //check if inequalities are valid
    for (int i=1; i<cases[ind].nIneqs; i++)
      if (0 > -cases[ind].ineqs[i][0])
        triviallyFeasible = false;
    //check if the equalities are valid
    long double eps = 0.01;
    for (int i=0; i<cases[ind].nEqs; i++)
      if (0 < -cases[ind].eqs[i][0] - eps ||  0 > -cases[ind].eqs[i][0] + eps)
        triviallyFeasible = false;
    if (triviallyFeasible) {
      trivFea++;
    }

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
          infeasible++;
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
          infeasible++;
          break;
        }
      }
    }

    if (triviallyInfeasible)
      continue;

    if (cases[ind].nVars == 1) {
      infeasible++;
      continue;
    }

    if (cases[ind].nEqs==0 && cases[ind].nIneqs==0) {
      working++;
      continue;
    }

    vector<vector<long double>> A1, A2;
    vector<long double> b, b2;
    tie(A1, b) = extractAb(cases[ind].ineqs);
    tie(A2, b2) = extractAb(cases[ind].eqs);
    b.insert(b.end(), b2.begin(), b2.end());

    /*
    cout << "\nA1\n";
    for (auto row: A1)
      cout << row;
    cout << "\nA2\n";
    for (auto row: A2)
      cout << row;
      */

    /*

    cout << "\nDo Presolve";
    tuple<vector<vector<long double>>, vector<vector<long double>>> tup;
    presolveStack.clear();
    tie(tup, b) = doPresolve(A1, A2, b);
    tie(A1, A2) = tup;

    /*
    cout << "\nA1\n";
    for (auto row: A1)
      cout << row;
    cout << "\nA2\n";
    for (auto row: A2)
      cout << row;
    */

    /*
    //testing a different approach
    auto A2m = findMinus(A2);
    vector<long double> b2m;
    for (auto e: b2)
      b2m.push_back(-e);

    for (auto row: A2)
      A1.push_back(row);
    for (auto row: A2m)
      A1.push_back(row);
    b.insert(b.end(), b2m.begin(), b2m.end());
    b.insert(b.end(), b2.begin(), b2.end());
    A2.clear();
    */

    /*

    auto x = dualSimplex(A1, A2, b, vector<long double>(cases[ind].nVars*2, 0.));
    if (x == 0)
      infeasible++;
    if (x == 1) 
      working++;
    if (x == 2)
      failed++;
    if (x == 3)
      incorrect++;
    if (x == 4)
      methodFails++;
    if (x == 5)
      artInBasis++;

    cout << "\n\ninfeasible: " << infeasible << "\tHighs Infeasible: " << highsInf;
    cout << "\nWorking: " << working << "\tHighs Solved to Optimality: " << highsFea;
    cout << "\nIncorrect: " << incorrect;
    cout << "\nMethod Fails: " << methodFails;
    cout << "\nArtificial Variable in the Basis: " << artInBasis;
    cout << "\nfailed: " << failed << "\n";

    cout << "\n\nHighs Disagrees With my dual simplex: " << disagree << " times\n";

    cout << "\nTrivially Feasible: " << trivFea << "\n";

    if (x==2) {
      cout << "\nHIGHS STATUS: " << utilModelStatusToString(model_status);
      cout << "\nDual Simplex status: " << x;
      cout << "\nA1 size: " << A1.size();
      cout << "\nA2 size: " << A2.size();
      cin.clear();
      //cin.get();
    }

  }

  cout << "\n\ninfeasible: " << infeasible << "\tHighs Infeasible: " << highsInf;
  cout << "\nWorking: " << working << "\tHighs Solved to Optimality: " << highsFea;
  cout << "\nIncorrect: " << incorrect;
  cout << "\nMethod Fails: " << methodFails;
  cout << "\nArtificial Variable in the Basis: " << artInBasis;
  cout << "\nfailed: " << failed << "\n";

  cout << "\n\nHighs Disagrees With my dual simplex: " << disagree << " times\n";

  cout << "\nTrivially Feasible: " << trivFea << "\n";

  */


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


