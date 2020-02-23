#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <time.h>
#include "bits/stdc++.h"

// Calculate hamming distance between each Illumina barcode and all of bead barcodes
// and output unique matched Illumina barcodes

using namespace std;

// Check if file exists
bool is_file_exist(const string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

// Split a string into an array of substrings
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

// Calculate hamming distance between two strings
int hammingDist(string seq1, string seq2)
{
    int cou = 0;
    for (int i = 0; i < min(seq1.size(), seq2.size()); i++)
	{
        if (seq1[i] != seq2[i])
            cou++;
	}
    return cou;
}

// Process Illumina barcode
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
	
	if (argc != 7)
	{
		cout << "Please provide six parameters" << endl;
		return 0;
	}
	
    string beadBarcodeFile = argv[1];
	string illuminaBarcodeFile = argv[2];
	string outputFile = argv[3];	// barcode_matching_distance.txt
	string outputFile1 = argv[4];	// barcode_matched_details.txt
	string beadtype = argv[5];
	int dist_threshold = stoi(argv[6]);
	
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
	unordered_map<string, double> mpx;
	unordered_map<string, double> mpy;
	if (infile.is_open())
	{
		while (std::getline(infile, line)) 
		{
			vector<string> items = splitstr2string(line, '\t');
			beadBarcodes.push_back(items[0]);
			mpx[items[0]] = stod(items[1]);
			mpy[items[0]] = stod(items[2]);
		}
		infile.close();
	}
	
    std::ifstream infile2(illuminaBarcodeFile);
	vector<string> illuminaBarcodes;
	if (infile2.is_open())
	{
		while (std::getline(infile2, line)) 
		{
			illuminaBarcodes.push_back(line);
		}
		infile2.close();
	}
	
	ofstream outfile;
	outfile.open (outputFile);
	outfile << "IlluminaBarcodes" << "\t" << "ProcessedIlluminaBarcodes" << "\t" << "BeadBarcodes" << "\t" << "Distance" << endl;
	
	ofstream outfile1;
	outfile1.open (outputFile1);
	outfile1 << "IlluminaBarcodes" << "\t" << "ProcessedIlluminaBarcodes" << "\t" << "BeadBarcodes" << "\t" << "Distance" << "\tX\tY" << endl;
	
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
		int minDist = 100;
		int minCou = 0;
		string matchedb = "";
		for(int j=0; j<beadBarcodes.size(); j++) 
		{
			int dist = hammingDist(barcode, beadBarcodes[j]);
			if (dist == minDist)
			{
				minCou++;
			}
			if (dist < minDist)
			{
				minDist = dist;
				minCou = 1;
				matchedb = beadBarcodes[j];
				
				if (minDist == 0)
					break;
			}
		}
		
		outfile << illuminaBarcodes[i] << "\t" << barcode << "\t" << matchedb << "\t" << minDist << endl;
		
		if (minDist <= dist_threshold && minCou == 1 && matchedb != "")
		{
			outfile1 << illuminaBarcodes[i] << "\t" << barcode << "\t" << matchedb << "\t" << minDist << "\t" << mpx[matchedb] << "\t" << mpy[matchedb] << endl;
		}
	}
	
	outfile.close();
	outfile1.close();
	
	ofstream outfile2;
	outfile2.open (outputFile1.substr(0, outputFile1.length() - 3) + "finished");
	outfile2 << "finished" << endl;
	outfile2.close();
	
	my_time = time(NULL); 
    printf("%s", ctime(&my_time));
	
    return 0; 
}

