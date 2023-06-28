
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

struct vectorHash {
  std::size_t operator()(std::vector<HighsInt> const& vec) const {
      std::hash<uint32_t> h;
      std::size_t ret = vec.size();
      for(auto& i : vec) {
          ret ^= h(i) | i;
      }
      return ret;
  }
};

struct rowStruct {
  vector<int64_t> e;
  int64_t d;
};

int64_t gcd(rowStruct row) {
  int64_t x = row.d;
  for (auto y: row.e) {
    x = __gcd(x, y);
    if (x == 1)
      return x;
  }
  return x;
}

void divide(rowStruct *row, int64_t a) {
  row->d *= a;
  int64_t q = gcd(*row);
  row->d /= q;
  for (int i= 0; i < row->e.size(); i++)
    row->e[i] /= q;
}

rowStruct *makeRow(vector<int64_t> vec) {
  rowStruct* row = new rowStruct();
  row->d = 1;
  row->e = vec;
  return row;
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
    lb.push_back(-(double)row[0]);
    ub.push_back(1.0e30);
  }
  for (auto row: tc->eqs) {
    lb.push_back(-(double)row[0]);
    ub.push_back(-(double)row[0]);
  }

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

  //11001110111111
  //00110101000000 = 0xC40
  highs.setOptionValue("presolve_rule_off", 0x2EE0);
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
vector<vector<HighsInt>> transpose(vector<vector<HighsInt>> A) {
  vector<HighsInt> r(A.size(), 0);
  vector<vector<HighsInt>> At(A[0].size(), r);
  for (int row=0; row<A.size(); row++){
    for (int col=0; col<A[0].size(); col++) {
      At[col][row] = A[row][col];
    }
  }
  return At;
}
*/

template <typename S> 
ostream& operator<<(ostream& os, const vector<S>& vector) {
    os << " ";
    for (auto element : vector) {
        os << element << " ";
    }
    os << "\n";
    return os;
}

vector<vector<HighsInt>> transpose(vector<vector<HighsInt>> b){
  if (b.size() == 0)
    return b;
  vector<vector<HighsInt>> t(b[0].size(), vector<HighsInt>());
  for (int i = 0; i < b.size(); i++){
    for (int j = 0; j < b[i].size(); j++){
      t[j].push_back(b[i][j]);
    }
  }
  return t;
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


vector<vector<HighsInt>> glue(vector<vector<HighsInt>> A, vector<vector<HighsInt>> B) {
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

vector<vector<HighsInt>> findMinus(vector<vector<HighsInt>> A) {
  for (int row=0; row<A.size(); row++)
    for (int col=0; col<A[0].size(); col++) 
      A[row][col] = -A[row][col];
  return A;
}

vector<vector<HighsInt>> I(int n, int zeroes) {
  vector<vector<HighsInt>> A;
  for (int row=0; row<zeroes; row++)
    A.push_back(vector<HighsInt> (n+zeroes, 0));
  for (int row=0; row<n; row++) {
    vector<HighsInt> r;
    for (int col=0; col<n; col++)
        r.push_back(-(int)(row==col));
    A.push_back(r);
  }    
  return A;
}

tuple<vector<vector<HighsInt>>, vector<HighsInt>> extractAb(vector<vector<HighsInt>> p) {
  vector<vector<HighsInt>> A;
  vector<HighsInt> b;

  for (int i=0; i<p.size(); i++){
    vector<HighsInt> r;
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

  vector<vector<HighsInt>> A1, A2;
  vector<HighsInt> b, b2;
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

vector<vector<long double>> pivot(vector<vector<long double>> tableau, tuple<int, int> p) {
  int row, col;
  tie(row, col) = p;

  auto v = tableau[row][col];

  cout << "\nv: " << v;
  cin.clear();
  if (v==0.)
    cin.get();

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

bool canImprove(vector<long double> r) {
  for (auto e: r) {
    //cout << "CAN IMPROVE\n";
    if (e > 0)
      return true;
  }
  return false;
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

vector<long double> updateR(vector<long double> z, vector<long double> c) {
  vector<long double> u;
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

  /*
  //Finds how many are reduced to empty by presolve
  
  for (int i = 0; i<nCases; i++) {
    //emptyCountInt += getPresolve(&cases[i], true);
    emptyCount += getPresolve(&cases[i], false);
  }

  cout  << "\n\n" << "Not Ints reduced to empty:" << "\t\t" << emptyCount << "\n";
  //cout  << "\n\n" << "Ints reduced to empty:" << "\t\t" << emptyCountInt << "\n\n";

  cout << "\nPresolve Rule Application Count:\n";
  for (auto e: countPresolveRuleApplications)
    cout << e << "\t";
  cout << "\n";

  */

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
    vector<vector<HighsInt>> A1, A2;
    vector<HighsInt> b, b2;
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

  //  Fixed the segfaulting
  //  it was because I was getting a negative index
  //  Must implement the thing to normalise the rows
  //  this will avoid the root cause of the segfaults
  //  It doesn't look ike my dual is working perfecty
  //  learn about simd

  //  The alg works now
  //  The problem is with the dual
  //  Find the dual

  return 0;
}


