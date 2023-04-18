//#define _STANDALONE_TEST_
#ifndef _STANDALONE_TEST_
#include <chemicalSymbols.h>
#endif

// C++ include files
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>

using namespace std;

// return the atomic number for a chemical element
int symbol_to_z(string sym)
{
  static vector<string> elements;
  if (elements.size()==0)
  {
    string elementsString;
    elementsString.append("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V ");
    elementsString.append("Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc ");
    elementsString.append("Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu ");
    elementsString.append("Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi ");
    elementsString.append("Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr ");
    elementsString.append("Rf Db Sg Bh Hs Mt Ds Rg");
    istringstream iss(elementsString);
    copy(istream_iterator<string>(iss),
         istream_iterator<string>(),
         back_inserter<vector<string> >(elements));
  }

  int z;
  z = 1;

  for (unsigned int i=0; i<elements.size(); ++i)
  {
    if (sym.compare(elements[i])==0)
      return z;
    else
      ++z;
  }

  return -1;
}

#ifdef _STANDALONE_TEST_
int main()
{
  string element = "Ar";
  cout << symbol_to_z(element) << endl;
  element = "foo";
  cout << symbol_to_z(element) << endl;
  element = "Pt";
  cout << symbol_to_z(element) << endl;
  return 0;
}
#endif
