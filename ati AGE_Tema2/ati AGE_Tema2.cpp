// ati AGE_Tema0.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <ctime> 
#include <cstdlib>
#include <cmath>
#include <math.h>
#define sizeOfSample 30
#define numberOfCromosomes 100 // 1000 best results1
#define MIN_VAL 9999999.0
#define MAX_VAL 0.0
#define PI 3.1415926535897
#define EPSILON 0.1
#define PROB_MUT 0.001
#define PROB_CROSS 0.7
#define STOP 500

std::ofstream fout;


using namespace std;

int ***M, ***NewM;
double bestChromosome = MIN_VAL, bestValue = MIN_VAL;
int numberOfBits;
double C = -sizeOfSample * 418.9829, A, B;
double chromosomeRes[numberOfCromosomes + 1];
double nOfBits;
enum Name { DeJong, Schwefels, Rastrigins } name;

int Calculate_noOfBits(double a, double b) {

	//10^d, d = 2;
	double N = (b - a) * 100;
	nOfBits = log2(N);
	int result = nOfBits;
	if (nOfBits == result)
		return result;
	return result + 1;
}

double DeJong1(double sample[sizeOfSample]) {
	double sum = 0;
	for (int index = 0; index < sizeOfSample; index++) {
		sum += sample[index] * sample[index];
	}
	return sum;
}

double Schwefel(double sample[sizeOfSample]) {
	double sum = 0;
	for (int index = 0; index < sizeOfSample; index++) {
		sum += (-1) * sample[index] * sin(sqrt(abs(sample[index])));
	}
	return sum;
}

double Rastrigin(double sample[sizeOfSample]) {
	double sum = 0;
	for (int index = 0; index < sizeOfSample; index++) {
		sum += (sample[index] * sample[index] - 10 * cos(2 * PI*sample[index]));
	}
	return sum;
}

double CalculateResult(Name functionName, double sample[sizeOfSample]) {

	double result;
	switch (functionName) {
	case DeJong:
		result = DeJong1(sample);
		break;
	case Schwefels:
		result = Schwefel(sample);
		break;
	case Rastrigins:
		result = 10 * sizeOfSample + Rastrigin(sample);
		break;
	default:
		break;
	}
	return result;
}

double RandomValue(double lowValue, double highValue) {
	return lowValue + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (highValue - lowValue)));
}

double ToBase10(int vector[], int length) {

	double number = 0;
	for (int index = length - 1; index >= 1; index--) {
		if (vector[index] == 1) {
			number += pow(2, length - index - 1);
		}
	}
	//return number / 100;
	return A + number*(B - A) / (pow(2, nOfBits) - 1);
}

double FitnessFunction(Name name, double result) {

	switch (name) {
		case Schwefels:
			return abs(result + C);
			break;
		default:
			return 1 / (result + EPSILON);
			break;
	}
}

void InitialPopulation() {

	int crom = numberOfCromosomes, rows = sizeOfSample, cols = numberOfBits;
	M = new int**[crom], NewM = new int**[crom];
	for (int i = 0; i < crom; i++) {
		M[i] = new int*[rows];
		NewM[i] = new int*[rows];
		for (int j = 0; j < rows; j++) {
			M[i][j] = new int[cols];
			NewM[i][j] = new int[cols];
		}
	}

	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		for (int index2 = 0; index2 < sizeOfSample; index2++) {
			for (int bit = 0; bit < numberOfBits; bit++) {
				M[index1][index2][bit] = rand() % 2;
			}
		}
	}//am generat o 3-matrice random de biti / populatie de cromozomi
}

void CopyData(int survivors[]) {

	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		for (int index2 = 0; index2 < sizeOfSample; index2++) {
			for (int bit = 0; bit < numberOfBits; bit++) {
				NewM[index1][index2][bit] = M[survivors[index1]][index2][bit];
			}
		}
	}
}

void CopyMatrix() {

	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		for (int index2 = 0; index2 < sizeOfSample; index2++) {
			for (int bit = 0; bit < numberOfBits; bit++) {
				M[index1][index2][bit] = NewM[index1][index2][bit];
			}
		}
	}
}

void RouletteWheel(Name name) {

	double result, cromosomeSum = 0, random;
	double sample[sizeOfSample];
	double prob[numberOfCromosomes + 1];
	double qprob[numberOfCromosomes + 1];

	int survivorsChrom[numberOfCromosomes];

	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		for (int index2 = 0; index2 < sizeOfSample; index2++) {
			sample[index2] = ToBase10(M[index1][index2], numberOfBits);
			if (M[index1][index2][0] == 1)
				sample[index2] *= -1;
		}
		result = CalculateResult(name, sample);
		fout << result << '\n';
		if (result < bestValue) {
			bestValue = result;
		}
		chromosomeRes[index1] = FitnessFunction(name, result);
	}

	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		cromosomeSum += chromosomeRes[index1];
		if (chromosomeRes[index1] > bestChromosome) {
			bestChromosome = chromosomeRes[index1];
		}
	}

	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		prob[index1] = chromosomeRes[index1] / cromosomeSum;
	}

	qprob[0] = 0;
	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		qprob[index1 + 1] = qprob[index1] + prob[index1];
	}

	for (int index = 0; index < numberOfCromosomes; index++) {
		random = RandomValue(0, 1);
		survivorsChrom[index] = 0;
		for (int index1 = 0; index1 < numberOfCromosomes; index1++)
			if (qprob[index1] <= random && random <= qprob[index1 + 1]) {
				survivorsChrom[index] = index1;
			}
	}
	CopyData(survivorsChrom);
}

void Mutation() {
	
	double random;
	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		for (int index2 = 0; index2 < sizeOfSample; index2++) {
			for (int bit = 0; bit < numberOfBits; bit++) {
				random = RandomValue(0, 1);
				if (random < PROB_MUT) {
					NewM[index1][index2][bit] = 1 - NewM[index1][index2][bit];
				}
			}
		}
	}
}

void Cross1() {

	double random;
	for (int index1 = 0; index1 < numberOfCromosomes; index1 += 2) {
		for (int index2 = 0; index2 < sizeOfSample; index2++) {
			for (int bit = 0; bit < numberOfBits; bit++) {
				random = RandomValue(0, 1);
				if (random < 0.5) {
					NewM[index1][index2][bit] = NewM[index1 + 1][index2][bit];
				}
			}
		}
	}

}


void Cross() {

	double random;
	int crossSurv[numberOfCromosomes];
	int k = 0;
	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		random = RandomValue(0, 1);
		if (random < PROB_CROSS) {
			crossSurv[k++] = index1;
		}
	}

	if (k % 2 == 1) {
		k--;
	}

	int randomLin = rand() % sizeOfSample;
	int randomCol = rand() % numberOfBits;
	int aux, j;

	for (int index = 0; index < k; index += 2) {
		for (int i = randomLin; i < sizeOfSample; i++) {
			if (i == randomLin) {
				j = randomCol;
			}
			else {
				j = 0;
			}
			for (; j < numberOfBits; j++) {

				aux = NewM[crossSurv[index]][i][j];
				NewM[crossSurv[index]][i][j] = NewM[crossSurv[index + 1]][i][j];
				NewM[crossSurv[index + 1]][i][j] = aux;
			}
		}
	}
}


double EvaluateOffSprings(Name name) {

	double best = MIN_VAL, result;
	double sample[sizeOfSample];
	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		for (int index2 = 0; index2 < sizeOfSample; index2++) {
			sample[index2] = ToBase10(NewM[index1][index2], numberOfBits);
			if (NewM[index1][index2][0] == 1)
				sample[index2] *= -1;
		}
		result = CalculateResult(name, sample);
		if (result < bestValue) {
			bestValue = result;
		}
		chromosomeRes[index1] = FitnessFunction(name, result);
		if (chromosomeRes[index1] > best) {
			best = chromosomeRes[index1];
		}
	}
	return best;
}

void GeneticAlgorithm(Name functionName) {

	int counter = 0;
	double result;

	InitialPopulation();

	while (counter < STOP) {

		RouletteWheel(functionName);
		Mutation();
		Cross();
		result = EvaluateOffSprings(functionName);
		//cout << bestValue << '\n';// " --> " << result << '\n';
		if (result > bestChromosome) {
			counter = 0;
			bestChromosome = result;
		}
		CopyMatrix();
		counter++;
	}
	//cout << "Best Chromosome: " << bestChromosome << '\n';
	cout << "Best Value: " << bestValue << '\n';
	

	
	fout << bestValue<<'\n';
}

void SelectFunction(int option) {

	switch (option) {
	case 1:
		name = DeJong;
		A = -5.12; B = 5.12;
		numberOfBits = Calculate_noOfBits(A, B);
		break;
	case 2:
		A = -500.0; B = 500.0;
		name = Schwefels;
		numberOfBits = Calculate_noOfBits(A, B);
		break;
	case 3:
		A = -5.12; B = 5.12;
		name = Rastrigins;
		numberOfBits = Calculate_noOfBits(A, B);
		break;
	default:
		break;
	}
}

void ReadOption() {
	int option;
	cout << "Introduceti optiunea (1 - DeJong, 2 - Schwefels, 3 - Rastrigin): ";
	cin >> option;
	SelectFunction(option);
}

int main() {
	fout.open("result.txt", std::ios_base::app);
	srand((unsigned int)time(NULL));
	ReadOption();
	GeneticAlgorithm(name);
	return 0;
}