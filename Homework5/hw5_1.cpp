#include <iostream>
#include <random>
#include <numeric>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace std;
#define Nsample 100000

double solve_mu(double T, double m, double K, double S0){
	double l=0, dt=T/m, r, mid, func;
	r = m*K/S0;
	while (1){
		mid = (l+r)/2;
		func = 0;
		for (int i=1; i<=m; i++)
			func += exp(i*mid*dt);
		if (abs(func-m*K/S0)<=0.01)
			return mid;
		if (func-m*K/S0>0)
			r = mid;
		else l = mid;
	}
	return mid;
}

double price_Asian(double mu, double m, double K, double S0, double r,
	double sigma, double T){
	double dt = T/m, rand_X, St, mean, l;
	vector<double> series, X, res;
	default_random_engine generator;
	normal_distribution<double> dist(mu*dt, sigma*sqrt(dt));
	for (int isample=0; isample<Nsample; isample++){
		St = S0;
		series.clear();
		X.clear();
		for (int i=1; i<=m; i++){
			rand_X = dist(generator);
			X.push_back(rand_X);
			St *= exp(rand_X);
			series.push_back(St);
		}
		mean = accumulate(begin(series),end(series),0.0) / series.size();
		l = 1;
		for (int i=0; i<m; i++)
			l *= exp((2*X[i]-mu*dt-r*dt+sigma*sigma*dt/2)*(r*dt-sigma*sigma
						*dt/2-mu*dt)/(2.0*sigma*sigma*dt));
		res.push_back(max(mean-K,0.0)*l);
	}
	double price = accumulate(begin(res),end(res),0.0)/Nsample;
	cout << "The price of Asian option by importance sampling: ";
	cout << price*exp(-r*T) << endl;
	double var = 0.0;
	for (vector<double>::iterator i=begin(res); i!=end(res); i++)
		var += pow(*i-price,2.0);
	cout << "Standard deviation: ";
	cout << sqrt(var/Nsample) << endl;
	cout << "Confidence interval: [" << price-1.96*sqrt(var)/Nsample
		<< "," << price+1.96*sqrt(var)/Nsample << "]" << endl;
	return sqrt(var);
}

double price_Asian_naive(double m, double K, double S0, double r,
	double sigma, double T){
	double dt = T/m, rand_X, St, mean, l, X;
	vector<double> series, res;
	default_random_engine generator;
	normal_distribution<double> dist((r-sigma*sigma/2)*dt, sigma*sqrt(dt));
	for (int isample=0; isample<Nsample; isample++){
		St = S0;
		series.clear();
		for (int i=1; i<=m; i++){
			rand_X = dist(generator);
			St *= exp(rand_X);
			series.push_back(St);
		}
		mean = accumulate(begin(series),end(series),0.0) / series.size();
		res.push_back(max(mean-K,0.0));
	}
	double price = accumulate(begin(res),end(res),0.0)/Nsample;
	cout << "The price of Asian option by naive Monte-Carlo: ";
	cout << price*exp(-r*T) << endl;
	double var = 0.0;
	for (vector<double>::iterator i=begin(res); i!=end(res); i++)
		var += pow(*i-price,2.0);
	cout << "Standard deviation: ";
	cout << sqrt(var/Nsample) << endl;
	cout << "Confidence interval: [" << price-1.96*sqrt(var)/Nsample
		<< "," << price+1.96*sqrt(var)/Nsample << "]" << endl;
	return sqrt(var);
}

int main(){
	double S0=50, sigma=0.2, r=0.05, T=0.25;
	double m=5, K=60, mu;
	mu = solve_mu(T, m, K, S0);
	cout << "We choose mu as: " << mu << endl;
	double var1 = price_Asian(mu,m,K,S0,r,sigma,T);
	for (int i=0; i<80; i++) cout << "-";
	cout << endl;
	double var2 = price_Asian_naive(m,K,S0,r,sigma,T);
	cout << "Reduction ratio: " << var2/var1 << endl;
	return 0;
}