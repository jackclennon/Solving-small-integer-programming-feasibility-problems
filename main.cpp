
#include "Highs.h"
#include <fstream>
#include <iostream>
#include <sstream> 

using namespace std;

struct IP {
  int id;
  int nVars;
  int nIneqs;
  int nEqs;
  vector<vector<int>> ineqs;
  vector<vector<int>> eqs;
};

typedef struct {
  IP *cases[100]; //[150218];
  int recent = 0;
} testcases;

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

struct IP* createIP(struct IP *ip, int id, int nV, vector<string> in, vector<string> eq) {
  ip->id = id;
  ip->nVars = nV;
  ip->nIneqs = in.size();
  ip->nEqs = eq.size();

  for (int i = 0; i < in.size(); i++) {
    vector<int> constr;
    vector<string> strConstr = split(in[i]);
    for (int j = 0; j < strConstr.size(); j++) {
      constr.push_back(stoi(strConstr[j]));
    }
    ip->ineqs.push_back(constr);
  }

  for (int i = 0; i < eq.size(); i++) {
    vector<int> constr;
    vector<string> strConstr = split(eq[i]);
    for (int j = 0; j < strConstr.size(); j++) {
      constr.push_back(stoi(strConstr[j]));
    }
    ip->eqs.push_back(constr);
  }

  return ip;
}

void saveTestCase(testcases *test, IP *ip) {
  test->cases[test->recent] = ip;
  test->recent++;
}

int main() {
  testcases *Test = new testcases{};

  fstream file;
  file.open("feasibility_testcases.txt",ios::in); //open a file to perform read operation using file object
  if (file.is_open()){   //checking whether the file is open
    string elem;
    int id = 0;
    int nVars = -1;
    int nIneqs = -1;
    vector<string> ineqs;
    int nEqs = -1;
    vector<string> eqs;
    bool readIneqs = false;

    while(getline(file, elem) && id < 100) {
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
        ip = createIP(ip, id, nVars, ineqs, eqs);
        saveTestCase(Test, ip);

        id++;
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

  for (auto ip: Test->cases) {
    cout << ip->nIneqs;
    for (auto row: ip->ineqs){
      cout << "\n";
      for (auto e: row)
        cout << e << " ";
    }
    cout << "\n\n\n";
  }

  return 0;
}


