#include <iostream>
#include <random>
#include <vector>
#include <cmath>
using namespace std;
#define Nsample 100000
#define M 100

int main(){
	default_random_engine generator;
	uniform_real_distribution<double> dist(0.0,1.0);
	vector<double> res;
	double coin, l;
	int count;
	res.clear();
	for (int isample=0; isample<Nsample; isample++){
		count = 0;
		for (int i=0; i<M; i++){
			coin = dist(generator);
			if (coin<=0.8) count++;
		}
		if (count>=80){
			l = pow(0.5/0.8,count*1.0)*pow(0.5/0.2,(M-count)*1.0);
			res.push_back(l);
		}
		else res.push_back(0);
	}
	double mean = accumulate(begin(res),end(res),0.0) / res.size();
	cout << "Probability of exceeding 80: ";
	cout << mean << endl;
	double var = 0.0;
	for (vector<double>::iterator i=begin(res); i!=end(res); i++)
		var += pow(*i-mean,2.0);
	cout << "Standard deviation: ";
	cout << sqrt(var/Nsample) << endl;
	cout << "Confidence interval: [" << mean-1.96*sqrt(var)/Nsample
		<< "," << mean+1.96*sqrt(var)/Nsample << "]" << endl;
	return 0;
}

