
#include "Highs.h"
#include "lp_data/HighsModelUtils.h"
#include <fstream>
#include <iostream>
#include <sstream> 
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

const int nCases = 150218;
// #duplicates
//                 145564    hashing    4624
//                 144238    IPcomp     5980

//Not very many is it


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

inline void hash_combine(std::size_t& seed) { }

template <typename T, typename... Rest>
inline void hash_combine(std::size_t& seed, const T& v, Rest... rest) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    hash_combine(seed, rest...);
}

struct IPHash {
  size_t operator()(IP const& ip) const noexcept {
    size_t h1 = hash<int>{}(ip.nVars);
    size_t h2 = hash<int>{}(ip.nIneqs);
    size_t h3 = hash<int>{}(ip.nEqs);
    vector<size_t> h4l;
    for (int i=1; i < ip.ineqs.size(); i++)
      h4l.push_back(vectorHash{}(ip.ineqs[i]));
    vector<size_t> h5l;
    for (int i=1; i < ip.eqs.size(); i++)
      h5l.push_back(vectorHash{}(ip.eqs[i]));
    size_t h4;
    if (ip.nIneqs > 0)
      h4 = vectorHash{}(ip.ineqs[0]);
    else
      h4 = vectorHash{}(vector<HighsInt>(1, 0));
    size_t h5;
    if (ip.nEqs > 0)
      h5 = vectorHash{}(ip.eqs[0]);
    else
      h5 = vectorHash{}(vector<HighsInt>(1, 0));
    for (auto h: h4l)
      hash_combine(h4, h);
    for (auto h: h5l)
      hash_combine(h5, h);
    hash_combine(h1, h2, h3, h4, h5);
    return h1;
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
        values.push_back((double)tc->ineqs[row][col]);
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

int main() {

  vector<size_t> hashKeys;
  vector<int> uniqueProblems;

  vector<IP> found;

  int dupCount = 0;
  
  fstream file;
  file.open("feasibility_testcases.txt",ios::in); //open a file to perform read operation using file object
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

        /*
        //hash it
        ip->hash_code = IPHash{}(*ip);
        if (!ip->searchHash(hashKeys)) {
          hashKeys.push_back(ip->hash_code);
          uniqueProblems.push_back(recent);
        }
        */

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

  /*

  Finds how many are reduced to empty by presolve

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

  //cout << hashKeys.size();
  cout << "\n" << dupCount << "\n";

  return 0;
}


