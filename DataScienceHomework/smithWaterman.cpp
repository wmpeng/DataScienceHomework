/**
* Implement Smith-Waterman Algorithm
* author: ÍõÃ÷Åô 1511212
* date: 2017/7/1 23:11
*/
#include<iostream>
#include<stdio.h>
#include<queue>
#include<algorithm>
#include<fstream>
#include<string>
#define singleGapPenalty
#define cerr(x) cout<<(#x)<<": "<<(x)<<endl;
using namespace std;

const int A = 0, T = 1, G = 2, C = 3;

/**
* Store the result of Smith-Waterman algorithm
*/
struct Result
{
	deque<int> resultA, resultB;
	int startPosA = 0, startPosB = 0, endPosA = 0, endPosB = 0;
	double maxScore;
};

/**
* Smith-Waterman algorithm
* @param seqA			first DNA number sequence
* @param lengthA		length of first sequence
* @param seqB			second DNA number sequence
* @param lengthB		length of second sequence
* @param scoreM			substitution-score matrix
* @param gapPenalty		gap penalty
* @return		the result of Smith-Waterman algorithm
*/
void SmithWaterman(
	const int lengthA, int *seqA, const int lengthB, int *seqB,
	double scoreM[4][4], double gapPenalty[4],
	Result *result
)
{
	deque<int> &resultA = result->resultA, &resultB = result->resultB;
	int &startPosA = result->startPosA, &startPosB = result->startPosB, &endPosA = result->endPosA, &endPosB = result->endPosB;
	double &maxScore = result->maxScore;

	/*
	* Construct a scoring matrix f and initialize
	*/
	double **f = new double*[lengthA + 1];
	for (int i = 0; i <= lengthA; i++)
		f[i] = new double[lengthB + 1]{ 0 };

	/*
	* Computation of scores
	*/
	for (int i = 1; i <= lengthA; i++)
		for (int j = 1; j <= lengthB; j++)
		{
			double seq[4] = { 0,f[i - 1][j - 1] + scoreM[seqA[i]][seqB[j]], f[i - 1][j] + gapPenalty[seqA[i]],f[i][j - 1] + gapPenalty[seqB[j]] };
			f[i][j] = *max_element(seq, seq + 4);
		}

	/*
	* Find the position of maximum value
	*/
	double max = -1;
	int maxi, maxj;
	for (int i = 1; i <= lengthA; i++)
		for (int j = 1; j <= lengthB; j++)
			if (f[i][j] > max)
				maxi = i, maxj = j, max = f[i][j];
	endPosA = maxi;
	endPosB = maxj;
	maxScore = max;

	/*
	* Trace back
	*/
	for (int i = maxi; i >= 1; )
		for (int j = maxj; j >= 1; )
		{
			startPosA = i;
			startPosB = j;
			if (f[i][j] == f[i - 1][j - 1] + scoreM[seqA[i]][seqB[j]])
			{
				resultA.push_front(seqA[i]);
				resultB.push_front(seqB[j]);
				i--;
				j--;
			}
			else if (f[i][j] == f[i - 1][j] + gapPenalty[seqA[i]])
			{
				resultA.push_front(seqA[i]);
				resultB.push_front(-1);
				i--;
			}
			else if (f[i][j] == f[i][j - 1] + gapPenalty[seqB[j]])
			{
				resultA.push_front(-1);
				resultB.push_front(seqB[j]);
				j--;
			}
			if (f[i][j] == 0)
				return;
		}
}

/*
* check if the arguments are correct or not
* @param argc		same as function main()
* @param argv		same as function main()
* @param ifSeqA		ifstream instance of the first DNA sequence file
* @param ifSeqB		ifstream instance of the second DNA sequence file
* @param ifScore	ifstream instance of the score matrix and gap penalty file
* @return		if the arguments are correct
*/
bool check(const int argc, char** const argv, ifstream &ifSeqA, ifstream &ifSeqB, ifstream &ifScore)
{
	/**
	* Give an example of arguments
	*/
	if (argc != 4)
	{
		printf(
			"usage: %s <1stFile> <2ndFile> <3rdFile>\n"
			"  <1stFile> : the first DNA sequence file\n"
			"  <2ndFile> : the second DNA sequence file\n"
			"  <3rdFile> : score matrix and gap penalty file\n"
			"\n"
			"example:\n"
			"  %s ./A.txt ./B.txt ./S.txt\n\n", argv[0], argv[0]);
		return false;
	}

	/**
	* Open the files and check
	*/
	ifSeqA.open(argv[1]);
	ifSeqB.open(argv[2]);
	ifScore.open(argv[3]);
	if (!ifSeqA.is_open())
	{
		cout << "\"" << argv[1] << "\" not exist." << endl << endl;
		return false;
	}
	if (!ifSeqB.is_open())
	{
		cout << "\"" << argv[2] << "\" not exist." << endl << endl;
		return false;
	}
	if (!ifScore.is_open())
	{
		cout << "\"" << argv[3] << "\" not exist." << endl << endl;
		return false;
	}
	return true;
}

int main(int argc, char** argv)
{
	string seqA, seqB;
	double scoreM[4][4], gapPenalty[4];
	ifstream ifSeqA, ifSeqB, ifScore;

	if (!check(argc, argv, ifSeqA, ifSeqB, ifScore))
		return 0;

	/**
	* input data
	*/
	string temp;
	while (ifSeqA >> temp)
		seqA += temp;
	while (ifSeqB >> temp)
		seqB += temp;

	for (int i = 1; i <= 5; i++)
		ifScore >> temp;
	for (int i = 0; i < 4; i++)
	{
		ifScore >> temp;
		for (int j = 0; j < 4; j++)
			ifScore >> scoreM[i][j];
	}
	ifScore >> temp;

#ifdef singleGapPenalty
	ifScore >> gapPenalty[0];
	gapPenalty[3] = gapPenalty[2] = gapPenalty[1] = gapPenalty[0];
#else
	for (int i = 0; i<4; i++)
		ifScore >> gapPenalty[i];
#endif

	/**
	* Translate DNA character sequences to number sequences by rules:A -> 0, T -> 1, G -> 2, C -> 3;
	*/
	int *numSeqA = new int[seqA.length()], *numSeqB = new int[seqB.length()];
	for (int i = 0; i<seqA.length(); i++)
		switch (seqA[i])
		{
		case 'A': numSeqA[i + 1] = A; break;
		case 'T': numSeqA[i + 1] = T; break;
		case 'G': numSeqA[i + 1] = G; break;
		case 'C': numSeqA[i + 1] = C; break;
		}
	for (int i = 0; i<seqB.length(); i++)
		switch (seqB[i])
		{
		case 'A': numSeqB[i + 1] = A; break;
		case 'T': numSeqB[i + 1] = T; break;
		case 'G': numSeqB[i + 1] = G; break;
		case 'C': numSeqB[i + 1] = C; break;
		}

	/**
	* Call the SmithWaterman function
	*/
	Result result;
	string resultA, resultB;
	SmithWaterman(seqA.length(), numSeqA, seqB.length(), numSeqB, scoreM, gapPenalty, &result);

	/**
	* Translate result sequences to character sequences: -1 -> '-', 0 -> 'A', 1 -> 'T', 2 -> 'G', 3 -> 'C';
	*/
	for (int i = 0; i<result.resultA.size(); i++)
		switch (result.resultA[i])
		{
		case -1: resultA += '-'; break;
		case A: resultA += 'A'; break;
		case T: resultA += 'T'; break;
		case G: resultA += 'G'; break;
		case C: resultA += 'C'; break;
		}
	for (int i = 0; i<result.resultB.size(); i++)
		switch (result.resultB[i])
		{
		case -1: resultB += '-'; break;
		case A: resultB += 'A'; break;
		case T: resultB += 'T'; break;
		case G: resultB += 'G'; break;
		case C: resultB += 'C'; break;
		}

	printf(
		"highest value is %f\n"
		"sequences 1 from %5d to %5d : %s\n"
		"sequences 2 from %5d to %5d : %s\n", result.maxScore, result.startPosA, result.endPosA, resultA.c_str(), result.startPosB, result.endPosB, resultB.c_str()
	);

	return 0;
}