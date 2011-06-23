/*
decaycutter.cc

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
	ifstream in1(argv[1]); //this is the input slha file 
	ofstream out1("cut1.slha"); //this is the output slha file

	string thing1 = "EL-   decays";
	string thing2 = "MUL-  decays";
	string thing3 = "TAU1- decays";
	
	string line; //this stores each line from the input slha file
	
	while(getline(in1, line)) //while loop runs till end of file (not really efficient)
	{
		if(line.find(thing1)!=string::npos) //run if find the phrase "Block MASS" in line	
		{
			break;
		}
		out1 << line << endl;
	}
	
	while(getline(in1, line))
	{
		if(line.find("Width")==string::npos)
		{
			//Keep going till find end of particular decay
		}
		else //find end of decay
		{
			break; //end while loop
		}
	
	}

	while(getline(in1, line))
	{
		out1 << line << endl;
	}

	in1.close();

	//-------------------------------------

	ifstream in2("cut1.slha"); //this is the input slha file 
	ofstream out2("cut2.slha"); //this is the output slha file

	while(getline(in2, line)) //while loop runs till end of file (not really efficient)
	{
		if(line.find(thing2)!=string::npos) //run if find the phrase "Block MASS" in line	
		{
			break;
		}
		out2 << line << endl;
	}
	
	while(getline(in2, line))
	{
		if(line.find("Width")==string::npos)
		{
			//Keep going till find end of particular decay
		}
		else //find end of decay
		{
			break; //end while loop
		}
	
	}
	
	while(getline(in2, line))
	{
		out2 << line << endl;
	}

	in2.close();

	//-------------------------------------

	ifstream in3("cut2.slha"); //this is the input slha file 
	ofstream out3(argv[2]); //this is the output slha file
	
	while(getline(in3, line)) //while loop runs till end of file (not really efficient)
	{
		if(line.find(thing3)!=string::npos) //run if find the phrase "Block MASS" in line	
		{
			break;
		}
		out3 << line << endl;
	}
	
	while(getline(in3, line))
	{
		if(line.find("Width")==string::npos)
		{
			//Keep going till find end of particular decay
		}
		else //find end of decay
		{
			break; //end while loop
		}
	
	}
	
	while(getline(in3, line))
	{
		out3 << line << endl;
	}

	in3.close();


}
