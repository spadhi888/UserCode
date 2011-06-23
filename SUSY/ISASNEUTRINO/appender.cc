/*
appender.cc

Running: The executable takes two arguments. For example, if the executable is called massextractor, then run by writing
massextractor slhafile.slha output.slha
where output.slha just contains everything after (and including) the masses and before Decays
*/

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main(int argc, char* argv[])
{
	ifstream ISAJET(argv[1]); //this is the input slha file containing the ISAJET data
	ifstream BR(argv[2]); //this is the input slha file containing the BRs from SDECAY
	ofstream out(argv[3]); //
	
	string line; //this stores each line from the input slha file
	
	//Copy contents of ISAJET file into new file
	while(getline(ISAJET, line)) //while loop runs till end of file (not really efficient)
	{
		out << line << endl; //copy all contents into out stream
	}
	
	//Copy contents of BR file into new file
	while(getline(BR, line))
	{
		out << line << endl;
	}
	
	//Done
	BR.close();
	ISAJET.close();
	out.close();
}
