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
	unordered_map<char, string> mp = {{'A', "A"}, {'C', "C"}, {'G', "G"}, {'T', "T"}, {'R', "AG"}, {'Y', "CT"}, {'M', "AC"}, {'K', "GT"}, {'S', "CG"}, {'W', "AT"}, {'H', "ACT"}, {'B', "CGT"}, {'V', "ACG"}, {'D', "AGT"}};
	
    int cou = 0;
    for (int i = 0; i < min(seq1.size(), seq2.size()); i++)
	{
		if (seq1[i] == seq2[i] || seq1[i] == 'N' || seq2[i] == 'N')
            continue;
		bool f = false;
		for (char c1 : mp[seq1[i]]) {
			for (char c2 : mp[seq2[i]]) {
				if (c1 == c2) {
					f = true;
					break;
				}
			}
			if (f) break;
		}
		if (!f) cou++;
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

// return unique matched illumina barcodes
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
	string illuminaBarcodeFile = argv[2];
	string outputFile = argv[3];   // barcode_matching_distance.txt
	string outputFile3 = argv[4];  // barcode_matched_details.txt
	string beadtype = argv[5];
	
	if (!is_file_exist(beadBarcodeFile))
	{
		cout << beadBarcodeFile + " not found" << endl;
		return 0;
	}
	
	if (!is_file_exist(illuminaBarcodeFile))
	{
		cout << illuminaBarcodeFile + " not found" << endl;
		return 0;
	}
	
	std::ifstream infile(beadBarcodeFile);
	std::string line;
	vector<string> beadBarcodes;
	vector<int> hd;
	vector<string> original;
	vector<string> coordx;
	vector<string> coordy;
	if (infile.is_open())
	{
		while (std::getline(infile, line)) 
		{
			//res[0] << "\t" << myDistances[i] << "\tY" << "\t" << mpx[res[0]] << "\t" << mpy[res[0]]
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			vector<string> items = splitstr2string(line, '\t');
			beadBarcodes.push_back(items[0]);
			hd.push_back(stoi(items[1]));
			original.push_back(items[2]);
			coordx.push_back(items[3]);
			coordy.push_back(items[4]);
		}
		infile.close();
	}
	
	// Run 3: Real barcodes, HD=0 matching threshold for the min HD1 barcodes, HD=1 matching threshold for all other barcodes, don’t treat degenerate barcodes differently.
	//{1, 0}, {2, 1}, {3, 1}, {4, 1}
	unordered_map<int, int> distmap = {{1, 0}};
	for (int i = 2; i <= 100; i++)
		distmap[i] = 1;
	
    std::ifstream infile2(illuminaBarcodeFile);
	vector<string> illuminaBarcodes;
	if (infile2.is_open())
	{
		while (std::getline(infile2, line)) 
		{
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			illuminaBarcodes.push_back(line);
		}
		infile2.close();
	}
	
	ofstream outfile;
	outfile.open (outputFile);
	outfile << "IlluminaBarcodes" << "\t" << "ProcessedIlluminaBarcodes" << "\t" << "BeadBarcodes" << "\t" << "Distance" << "\tX\tY" << endl;
	
	ofstream outfile3;
	outfile3.open (outputFile3);
	outfile3 << "IlluminaBarcodes" << "\t" << "ProcessedIlluminaBarcodes" << "\t" << "BeadBarcodes" << "\t" << "Distance" << "\tX\tY" << endl;
	
	string logfile = outputFile3 + ".log";
	ofstream outfile4;
	outfile4.open (logfile);
	
	outfile4 << "Num of Illumina barcodes: " << illuminaBarcodes.size() << endl;
	int unique_matched = 0;
	int multi_matched = 0;
	
	for(int i=0; i<illuminaBarcodes.size(); i++)
	{
		string barcode = illuminaBarcodes[i];
		if (beadtype == "180402") // x(1:8),x(9:14)
			barcode = barcode.substr(0, 14);
		if (beadtype == "SLACBeads_1_RC")
			barcode = seqrcomplement(barcode);
		if (beadtype == "BobMistake")
			barcode = barcode.substr(0, 7) + barcode.substr(8, 6) + "N";
		
		// calculate hamming distance		
		int minDist = 1000;
		int minCou = 0;
		string matchedb = "";
		string matchedx = "";
		string matchedy = "";
		int minDist2 = 1000;
		string matchedb2 = "";
		string matchedx2 = "";
		string matchedy2 = "";
		for(int j=0; j<beadBarcodes.size(); j++) 
		{
			int dist = hammingDist(barcode, beadBarcodes[j]);
			
			if (dist < minDist2) {
				minDist2 = dist;
				matchedb2 = beadBarcodes[j];
				matchedx2 = coordx[j];
				matchedy2 = coordy[j];
			}
			
			if (dist > distmap[hd[j]])
				continue;
			
			if (dist == minDist)
			{
				minCou++;
			}
			if (dist < minDist)
			{
				minDist = dist;
				minCou = 1;
				matchedb = beadBarcodes[j];
				matchedx = coordx[j];
				matchedy = coordy[j];
			}
		}
		
		outfile << illuminaBarcodes[i] << "\t" << barcode << "\t" << matchedb2 << "\t" << minDist2 << "\t" << matchedx2 << "\t" << matchedy2  << endl;
		
		if (minDist < 1000)
		{
			if (minCou == 1)
			{
				unique_matched++;
				if (matchedb != "")
					outfile3 << illuminaBarcodes[i] << "\t" << barcode << "\t" << matchedb << "\t" << minDist << "\t" << matchedx << "\t" << matchedy << endl;
			}
			else
				multi_matched++;
		}
	}
	
	outfile.close();
	outfile3.close();
	
	outfile4 << "Num of unique matched Illumina barcodes: " << unique_matched << endl;
	outfile4 << "Num of multiple matched Illumina barcodes: " << multi_matched << endl;
	outfile4.close();
	
	ofstream outfile2;
	outfile2.open (outputFile3.substr(0, outputFile3.length() - 3) + "finished");
	outfile2 << "finished" << endl;
	outfile2.close();
	
	my_time = time(NULL); 
    printf("%s", ctime(&my_time));
	
    return 0; 
}

