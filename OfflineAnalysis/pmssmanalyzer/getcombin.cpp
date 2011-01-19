#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include "TMath.h"

using namespace std;
void getcombin(int n, int r, vector<vector<int> >& vcomb)
{
  vector<int> a;
  for(int i=0; i < n; i++) a.push_back(i);

  vector<int> c(r);

  map<string, int> comb;
  int count = 0;

  copy(a.begin(),a.begin()+r, c.begin());

  sort(c.begin(), c.end());
  ostringstream out;
  for(int i=0; i < r; i++) out << c[i] << " ";
  comb[out.str()] = count;

  while (TMath::Permute(n, &a[0])) {
    count++;
    //cout << "Count: " << count << "\t";
    for (int i=0; i<n; i++) {
      //cout << a[i] << " ";
    }
    //cout << endl;

    vector<int> c(r);

    copy(a.begin(),a.begin()+r, c.begin());

    sort(c.begin(), c.end());
    ostringstream out;
    for(int i=0; i < r; i++) out << c[i] << " ";
    comb[out.str()] = count;
  }

  //cout << endl;
  vcomb.resize(comb.size());

  int index = 0;
  for(map<string, int>::iterator it=comb.begin(); it != comb.end(); it++)
    {
      //cout << it->second << "\t" << it->first << endl;
      istringstream inp(it->first);
      int j=0;
      vcomb[index] = vector<int>(r);
      while (inp >> vcomb[index][j]) j++;
      index++;
    }

/*
  //cout << endl << "just checking..." << endl;
  cout << "vcomb size: " << vcomb.size() << endl;
  for(unsigned int i=0; i < vcomb.size(); i++)
    {
      cout << i+1;
      for(unsigned int j=0; j < vcomb[i].size(); j++) 
	cout << "\t" << vcomb[i][j];
      cout << endl;
    }
  cout << endl;
*/
}
