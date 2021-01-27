#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <iterator>
#include <time.h>
#include "bits/stdc++.h"

using namespace std;

bool is_file_exist(const string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

vector<string> splitstr2string(const string str, const char delim)
{
	stringstream ss(str);
	string s;
	vector<string> res;
	while (getline(ss, s, delim))
	{
		res.push_back(s);
	}
	return res;
}

double dist(double x1, double y1, double x2, double y2) {
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

string degenerate(vector<string> barcodes) {
//R = A,G; Y = C,T; M = A,C; K = G,T; S = C,G; W = A,T; H = A,C,T; B = C,G,T; V = A,C,G; D = A,G,T; N = A,C,G,T

	string res = "";
	unordered_set<char> s;
	for (int j = 0; j < barcodes[0].size(); j++) {
		s.clear();
		for (int i = 0; i < barcodes.size(); i++) {
			s.insert(barcodes[i][j]);
		}
		
		if (s.find('N') != s.end() || (s.find('A') != s.end() && s.find('C') != s.end() && s.find('G') != s.end() && s.find('T') != s.end()))
			res += "N";
		else if (s.find('A') != s.end() && s.find('C') != s.end() && s.find('T') != s.end())
			res += "H";
		else if (s.find('C') != s.end() && s.find('G') != s.end() && s.find('T') != s.end())
			res += "B";
		else if (s.find('A') != s.end() && s.find('C') != s.end() && s.find('G') != s.end())
			res += "V";
		else if (s.find('A') != s.end() && s.find('G') != s.end() && s.find('T') != s.end())
			res += "D";
		else if (s.find('A') != s.end() && s.find('G') != s.end())
			res += "R";
		else if (s.find('A') != s.end() && s.find('C') != s.end())
			res += "M";
		else if (s.find('A') != s.end() && s.find('T') != s.end())
			res += "W";
		else if (s.find('C') != s.end() && s.find('T') != s.end())
			res += "Y";
		else if (s.find('G') != s.end() && s.find('T') != s.end())
			res += "K";
		else if (s.find('C') != s.end() && s.find('G') != s.end())
			res += "S";
		else
			res += barcodes[0].substr(j, 1);
	}
	return res;
}

bool isRaw(string barcode) 
{
	for (auto c : barcode) {
		if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
			return false;
		}
	}
	return true;
}

// degenerate in situ bead barcodes
int main(int argc, char const *argv[]) 
{	
	time_t my_time = time(NULL); 
    printf("%s", ctime(&my_time));
	
	if (argc != 6)
	{
		cout << "Please provide 5 parameters!" << endl;
		return 0;
	}
	
	string beadBarcodeFile = argv[1];
	string beadLocationFile = argv[2];
    string matchingFile = argv[3];
	string matchingFile2 = argv[4];
	string outputFile = argv[5];
	
	if (!is_file_exist(beadBarcodeFile))
	{
		cout << beadBarcodeFile + " not found" << endl;
		return 0;
	}
	
	if (!is_file_exist(beadLocationFile))
	{
		cout << beadLocationFile + " not found" << endl;
		return 0;
	}
	
	if (!is_file_exist(matchingFile))
	{
		cout << matchingFile + " not found" << endl;
		return 0;
	}
	
	if (!is_file_exist(matchingFile2))
	{
		cout << matchingFile2 + " not found" << endl;
		return 0;
	}
	
	std::ifstream infile(beadBarcodeFile);
	std::string line;
	vector<string> beadBarcodes;
	if (infile.is_open())
	{
		while (std::getline(infile, line)) 
		{
			line.erase(std::remove(line.begin(), line.end(), ','), line.end());
			beadBarcodes.push_back(line);
		}
		infile.close();
	}
	
	unordered_map<string, double> mpx;
	unordered_map<string, double> mpy;
	std::ifstream infile1(beadLocationFile);
	if (infile1.is_open())
	{
		std::getline(infile1, line);
		string corrx = line;
		std::getline(infile1, line);
		string corry = line;

		infile1.close();
		
		vector<string> xs = splitstr2string(corrx, ',');
		vector<string> ys = splitstr2string(corry, ',');
		
		for (int i = 0; i < beadBarcodes.size(); i++)
		{
			mpx[beadBarcodes[i]] = stod(xs[i]);
			mpy[beadBarcodes[i]] = stod(ys[i]);
		}
	}

	unordered_map<string, vector<string>> neighbors;
	vector<string> myBarcodes;
	vector<int> myDistances;
	
	std::ifstream infile2(matchingFile);
	if (infile2.is_open())
	{
		while (std::getline(infile2, line)) 
		{
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			vector<string> items = splitstr2string(line, '\t');
			myBarcodes.push_back(items[0]);
			myDistances.push_back(stoi(items[2]));
			neighbors[items[0]].push_back(items[1]);
		}
		infile2.close();
	}
	
	unordered_set<string> used;
	queue<string> Q;
	
	ofstream outfile;
	outfile.open (outputFile);
	
	for(int i = 0; i < myBarcodes.size(); i++)
	{
		string bc = myBarcodes[i];
		if (used.find(bc) != used.end())
			continue;
		
		vector<string> res;
		Q.push(bc);
		used.insert(bc);
		
		while (!Q.empty()) {
			string cur = Q.front();
			Q.pop();
			res.push_back(cur);
			
			for (string nei : neighbors[cur]) {
				if (used.find(nei) != used.end())
					continue;
				
				if (dist(mpx[cur], mpy[cur], mpx[nei], mpy[nei]) <= 10.0) {
					Q.push(nei);
					used.insert(nei);
				}
			}
		}
		
		if (res.size() == 1) {
			if (isRaw(res[0]))
				outfile << res[0] << "\t" << myDistances[i] << "\tY" << "\t" << mpx[res[0]] << "\t" << mpy[res[0]] << endl;
			else
				outfile << res[0] << "\t" << myDistances[i] << "\tN" << "\t" << mpx[res[0]] << "\t" << mpy[res[0]] << endl;
		}
		else {
			double x = 0.0, y = 0.0;
			for (string s : res) {
				x += mpx[s];
				y += mpy[s];
			}
			x /= res.size();
			y /= res.size();
			outfile << degenerate(res) << "\t" << myDistances[i] << "\tN" << "\t" << x << "\t" << y << endl;
		}
	}
	
	std::ifstream infile3(matchingFile2);
	if (infile3.is_open())
	{
		while (std::getline(infile3, line)) 
		{
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			outfile << line << endl;
		}
		infile3.close();
	}
	
	outfile.close();
	
	my_time = time(NULL); 
    printf("%s", ctime(&my_time));
	
    return 0;
}
