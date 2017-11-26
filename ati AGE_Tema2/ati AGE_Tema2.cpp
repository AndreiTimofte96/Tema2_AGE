// ati AGE_Tema0.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <ctime> 
#include <cstdlib>
#include <cmath>
#include <math.h>
#define sizeOfSample 10
#define numberOfRuns 10000
#define numberOfCromosomes 10
#define maxNoOfParam 3
#define MIN_VAL 9999999.0
#define MAX_VAL 0.0
#define PI 3.1415926535897
#define DMAX 1000
#define EPSILON 1


using namespace std;

int numberOfBits;
int ***M;
double function[numberOfCromosomes + 1];
double prob[numberOfCromosomes + 1];
double qprob[numberOfCromosomes + 1];

enum Name {Rastrigins } name;

int Calculate_noOfBits(double a, double b) {

	//10^d, d = 2;
	double N = (b - a) * 100;
	double n = log2(N);
	int result = n;
	if (n == result)
		return result;
	return result + 1;
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
	case Rastrigins:
		result = 10 * sizeOfSample + Rastrigin(sample);
		break;
	default:
		break;
	}
	return result;
}

double ToBase10(int vector[], int length) {

	double number = 0;
	for (int index = length - 1; index >= 1; index--) {
		if (vector[index] == 1) {
			number += pow(2, length - index - 1);
		}
	}
	return number / 100;
}

void InitialPopulation() {
	
	int crom = numberOfCromosomes, rows = sizeOfSample, cols = numberOfBits;
	M = new int**[crom];
	for (int i = 0; i < crom; i++){
		M[i] = new int*[rows];
		for (int j = 0; j < rows; j++)
			M[i][j] = new int[cols];
	}

	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		for (int index2 = 0; index2 < sizeOfSample; index2++) {
			for (int bit = 0; bit < numberOfBits; bit++) {
				M[index1][index2][bit] = rand() % 2;
			}
		}
	}//am generat o 3-matrice random de biti / populatie de cromozomi
}

double FitnessFunction(double result) {

	return 1 / (result + EPSILON);
}
void RouletteWheel(Name name) {
	
	double result, cromosomeSum = 0;
	double sample[sizeOfSample], minSample[sizeOfSample];
	double function[numberOfCromosomes + 1];
	
		for (int index1 = 0; index1 < numberOfCromosomes; index1++){
		for (int index2 = 0; index2 < sizeOfSample; index2++) {
			sample[index2] = ToBase10(M[index1][index2], numberOfBits);
			if (M[index1][index2][0] == 1)
				sample[index2] *= -1;
		}
		result = CalculateResult(name, sample);
		function[index1] = FitnessFunction(result);
	}

	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		cromosomeSum += function[index1];
	}

	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		prob[index1] = function[index1] / cromosomeSum;
	}
	
	qprob[0] = 0;
	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		qprob[index1 + 1] = qprob[index1] + prob[index1];
	}
	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
		cout << qprob[index1] << ' ';
	}

	for (int index1 = 0; index1 < numberOfCromosomes; index1++) {
	
	}
}

void GeneticAlgorithm(Name functionName) {
	
	InitialPopulation();
	RouletteWheel(functionName);
	
}


void Rastrigins_Init() {
	name = Rastrigins;
	numberOfBits = Calculate_noOfBits(-5.12, 5.12);
}

int main() {

	srand((unsigned int)time(NULL));
	Rastrigins_Init();
	GeneticAlgorithm(Rastrigins);
	return 0;
}