#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_map>
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

int hammingDist(string seq1, string seq2)
{
    int cou = 0;
    for (int i = 0; i < min(seq1.size(), seq2.size()); i++)
	{
		if (seq1[i] != seq2[i] && seq1[i] != 'N' && seq2[i] != 'N')
            cou++;
	}
    return cou;
}

string seqrcomplement(string barcode)
{
	string res = "";
	for (int i = barcode.size() - 1; i >= 0; --i)
	{
		if (barcode[i] == 'G')
			res += 'C';
		if (barcode[i] == 'C')
			res += 'G';
		if (barcode[i] == 'T')
			res += 'A';
		if (barcode[i] == 'A')
			res += 'T';
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

int main(int argc, char const *argv[]) 
{	
	time_t my_time = time(NULL); 
    printf("%s", ctime(&my_time));
	
	if (argc != 7)
	{
		cout << "Please provide 6 parameters!" << endl;
		return 0;
	}
	
    string beadBarcodeFile1 = argv[1];
	string beadBarcodeFile2 = argv[2];
	string beadLocationFile = argv[3];
	string outputFile = argv[4];   // barcode_matching_01.txt
	string outputFile2 = argv[5];  // barcode_matching_2.txt
	string beadtype = argv[6];
	
	if (!is_file_exist(beadBarcodeFile1))
	{
		cout << beadBarcodeFile1 + " not found" << endl;
		return 0;
	}
	
	if (!is_file_exist(beadBarcodeFile2))
	{
		cout << beadBarcodeFile2 + " not found" << endl;
		return 0;
	}
	
	if (!is_file_exist(beadLocationFile))
	{
		cout << beadLocationFile + " not found" << endl;
		return 0;
	}
	
	std::ifstream infile(beadBarcodeFile1);
	std::string line;
	vector<string> beadBarcodes1;
	if (infile.is_open())
	{
		while (std::getline(infile, line)) 
		{
			line.erase(std::remove(line.begin(), line.end(), ','), line.end());
			string barcode = line;
			if (beadtype == "180402") // x(1:8),x(9:14)
				barcode = barcode.substr(0, 14);
			if (beadtype == "SLACBeads_1_RC")
				barcode = seqrcomplement(barcode);
			if (beadtype == "BobMistake")
				barcode = barcode.substr(0, 7) + barcode.substr(8, 6) + "N";
			beadBarcodes1.push_back(barcode);
		}
		infile.close();
	}
	
	std::ifstream infile2(beadBarcodeFile2);
	vector<string> beadBarcodes2;
	if (infile2.is_open())
	{
		while (std::getline(infile2, line)) 
		{
			line.erase(std::remove(line.begin(), line.end(), ','), line.end());
			string barcode = line;
			if (beadtype == "180402") // x(1:8),x(9:14)
				barcode = barcode.substr(0, 14);
			if (beadtype == "SLACBeads_1_RC")
				barcode = seqrcomplement(barcode);
			if (beadtype == "BobMistake")
				barcode = barcode.substr(0, 7) + barcode.substr(8, 6) + "N";
			beadBarcodes2.push_back(barcode);
		}
		infile2.close();
	}

	unordered_map<string, string> mpx;
	unordered_map<string, string> mpy;
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
		
		for (int i = 0; i < beadBarcodes2.size(); i++)
		{
			mpx[beadBarcodes2[i]] = xs[i];
			mpy[beadBarcodes2[i]] = ys[i];
		}
	}
	
    ofstream outfile;
	outfile.open (outputFile);
	
	ofstream outfile2;
	outfile2.open (outputFile2);
	
	for(int i=0; i<beadBarcodes1.size(); i++)
	{
		string barcode1 = beadBarcodes1[i];
		int minDist = 1000;
		int minCou = 0;
		string matchedb = "";
		
		for(int j=0; j<beadBarcodes2.size(); j++) 
		{
			string barcode2 = beadBarcodes2[j];
			if (barcode1 == barcode2) continue;
			
			int dist = hammingDist(barcode1, barcode2);
			if (dist <= 1)
			{
				outfile << barcode1 << "\t" << barcode2 << "\t" << dist << endl;
			}
			
			if (dist < minDist)
			{
				minDist = dist;
			}
		}
		
		if (minDist > 1)
		{
			if (isRaw(barcode1))
				outfile2 << barcode1 << "\t" << minDist << "\tY" << "\t" << mpx[barcode1] << "\t" << mpy[barcode1] << endl;
			else
				outfile2 << barcode1 << "\t" << minDist << "\tN" << "\t" << mpx[barcode1] << "\t" << mpy[barcode1] << endl;
		}
	}
	
	outfile.close();
	outfile2.close();
	
	ofstream outfile3;
	outfile3.open (outputFile.substr(0, outputFile.length() - 3) + "finished");
	outfile3 << "finished" << endl;
	outfile3.close();
	
	my_time = time(NULL); 
    printf("%s", ctime(&my_time));
	
    return 0; 
}

