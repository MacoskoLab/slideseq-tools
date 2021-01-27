#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
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

int get_matched_len(string cigar) {
	int res = 0, v = 0;
	for (char c : cigar) {
		if (c >= '0' && c <= '9')
			v = v * 10 + c - '0';
		else {
			if (c == 'M')
				res += v;
			v = 0;
		}
	}
	return res;
}

int main(int argc, char const *argv[]) 
{	
	if (argc != 2)
	{
		cout << "Please provide one parameter: input sam file!" << endl;
		return 0;
	}
	
	string inputSamFile = argv[1];
	
	if (!is_file_exist(inputSamFile))
	{
		cout << inputSamFile + " not found" << endl;
		return 0;
	}
	
	unordered_map<string, int> mp;
	std::ifstream infile(inputSamFile);
	std::string line;
	if (infile.is_open())
	{
		while (std::getline(infile, line))
		{
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			if (line.size() <= 1 || line[0] == '@')
				continue;
			
			vector<string> items = splitstr2string(line, '\t');
			string read = items[0];
			mp[read]++;
		}
		infile.close();
	}
	
	std::ifstream infile1(inputSamFile);
	map<int, int> unique_score;
	map<int, int> unique_mismatch;
	map<int, int> unique_ratio;
	if (infile1.is_open())
	{
		while (std::getline(infile1, line))
		{
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			if (line.size() <= 1 || line[0] == '@')
				continue;
			
			vector<string> items = splitstr2string(line, '\t');
			string read = items[0];
			
			if (mp[read] == 1) {
				//AS:i:87	nM:i:0
				string quality_score = "";
				string mismatch = "";
				for (string is : items) {
					if (is.substr(0, 5) == "AS:i:") {
						quality_score = is.substr(5);
						unique_score[stoi(quality_score)]++;
					}
					if (is.substr(0, 5) == "nM:i:") {
						mismatch = is.substr(5);
						unique_mismatch[stoi(mismatch)]++;
					}
				}
				
				int l1 = get_matched_len(items[5]);
				int l2 = items[9].size();
				int mapping_ratio = l1 * 100 / l2;
				unique_ratio[mapping_ratio]++;
			}
		}
		infile1.close();
	}
	
	unordered_set<string> used;
	std::ifstream infile2(inputSamFile);
	map<int, int> multi_score;
	map<int, int> multi_mismatch;
	map<int, int> multi_ratio;
	if (infile2.is_open())
	{
		while (std::getline(infile2, line))
		{
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			if (line.size() <= 1 || line[0] == '@')
				continue;
			
			vector<string> items = splitstr2string(line, '\t');
			string read = items[0];
			
			if (mp[read] > 1) {
				if (used.find(read) == used.end()) {
					used.insert(read);
					
					//AS:i:87	nM:i:0
					string quality_score = "";
					string mismatch = "";
					for (string is : items) {
						if (is.substr(0, 5) == "AS:i:") {
							quality_score = is.substr(5);
							multi_score[stoi(quality_score)]++;
						}
						if (is.substr(0, 5) == "nM:i:") {
							mismatch = is.substr(5);
							multi_mismatch[stoi(mismatch)]++;
						}
					}
					
					int l1 = get_matched_len(items[5]);
					int l2 = items[9].size();
					int mapping_ratio = l1 * 100 / l2;
					multi_ratio[mapping_ratio]++;
				}
			}
		}
		infile2.close();
	}
	
	string outputFile1 = inputSamFile + ".unique.score";
	string outputFile2 = inputSamFile + ".multi.score";
	string outputFile3 = inputSamFile + ".unique.mismatch";
	string outputFile4 = inputSamFile + ".multi.mismatch";
	string outputFile5 = inputSamFile + ".unique.ratio";
	string outputFile6 = inputSamFile + ".multi.ratio";
	
	ofstream outfile1;
	outfile1.open (outputFile1);
	for (auto it = unique_score.begin(); it != unique_score.end(); it++)
		outfile1 << it->first << "\t" << it->second << endl;
	outfile1.close();
	
	ofstream outfile2;
	outfile2.open (outputFile2);
	for (auto it = multi_score.begin(); it != multi_score.end(); it++)
		outfile2 << it->first << "\t" << it->second << endl;
	outfile2.close();
	
	ofstream outfile3;
	outfile3.open (outputFile3);
	for (auto it = unique_mismatch.begin(); it != unique_mismatch.end(); it++)
		outfile3 << it->first << "\t" << it->second << endl;
	outfile3.close();
	
	ofstream outfile4;
	outfile4.open (outputFile4);
	for (auto it = multi_mismatch.begin(); it != multi_mismatch.end(); it++)
		outfile4 << it->first << "\t" << it->second << endl;
	outfile4.close();
	
	ofstream outfile5;
	outfile5.open (outputFile5);
	for (auto it = unique_ratio.begin(); it != unique_ratio.end(); it++)
		outfile5 << it->first << "\t" << it->second << endl;
	outfile5.close();
	
	ofstream outfile6;
	outfile6.open (outputFile6);
	for (auto it = multi_ratio.begin(); it != multi_ratio.end(); it++)
		outfile6 << it->first << "\t" << it->second << endl;
	outfile6.close();
	
    return 0; 
}
