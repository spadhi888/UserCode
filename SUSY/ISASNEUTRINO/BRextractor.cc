/*
massextractor.cc

Running: The executable takes two arguments. For example, if the executable is called massextractor, then run by writing
massextractor slhafile.slha output.slha
where output.slha just contains everything after (and including) the masses and before Decays
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

int main(int argc, char* argv[])
{
	ifstream in(argv[1]); //this is the input slha file 
	ofstream out(argv[2], ios::app); //this is the output slha file

	stringstream ss;
	ss << argv[3];
	string decayid;
	ss >> decayid;
	
	string line; //this stores each line from the input slha file
	
	while(getline(in, line)) //while loop runs till end of file (not really efficient)
	{
		if(line.find(decayid)!=string::npos) //run if find the phrase "Block MASS" in line	
		{
			out << line << endl;
			break;
		}
	}
	
	while(getline(in, line))
	{
		//NOTE: The end identifier used below depends on whether the output is from ISAJET or from SDECAY.
		//It's easier to use just Width as the identifier.
		//if(line.find("#         PDG            Width")==string::npos) //if don't find string, write to file
		if(line.find("Width")==string::npos)
		{
			out << line << endl;
		}
		else //if string in if statement found, end program
		{
			return 0;
		}

	}
}
