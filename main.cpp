#include<iostream>
#include<random>
#include<chrono>
#include<fstream>
#include"opt_alg.h"
#include"ode_solver.h"

using namespace std;

int main()
{
	try
	{
		cout << "LAB NUMBER " << LAB_NO << endl;
		cout << "LAB PART " << LAB_PART << endl << endl;
#if LAB_NO==1
		matrix Y0 = matrix(new double[3]{ 5,1,10 }, 3);
		matrix *Y = solve_ode(0, 1, 1000, Y0);
		ofstream S("sim_t.csv");
		S << Y[0];
		S.close();
		S.open("sim_t.csv");
		S << Y[1];
		S.close();
		
#elif LAB_NO==2

#if LAB_PART == 1

		double x0, d = 1, alfa = 1.5, epsilon = 1e-5, gamma = 1e-200;
		int Nmax = 1000;
		random_device R;

		x0 =200.0 * R() / R.max() - 100; // trzeba przejsc na double
		//x0 = 80;
		cout << x0 << endl << endl;

		double *ab = expansion(x0, d, alfa, Nmax);
		cout << ab[0] << "\t" << ab[1] << endl;
		cout << solution::f_calls << endl << endl;
		solution::clear_calls();

		solution opt_F = fib(ab[0], ab[1], epsilon);
		cout << opt_F << endl << endl;
		solution::clear_calls();

		solution opt_L = lag(ab[0], ab[1], epsilon, gamma, Nmax);
		cout << opt_L << endl << endl;
		solution::clear_calls();



#elif LAB_PART == 2

		double x0, d = 1e-4, alfa = 1.5, epsilon = 1e-10, gamma = 1e-200;
		int Nmax = 1000;
		random_device R;

		x0 = 200.0 * R() / R.max() - 100; // trzeba przejsc na double
										  //x0 = 80;
		x0 = (1e-2 - 1e-4)*R() / R.max() + 1e-4;
		cout << x0 << endl << endl;

		double *ab = expansion(x0, d, alfa, Nmax);
		cout << ab[0] << "\t" << ab[1] << endl;
		cout << solution::f_calls << endl << endl;
		solution::clear_calls();

		solution opt_F = fib(ab[0], ab[1], epsilon);
		cout << opt_F << endl << endl;
		solution::clear_calls();

		solution opt_L = lag(ab[0], ab[1], epsilon, gamma, Nmax);
		cout << opt_L << endl << endl;
		solution::clear_calls();

#endif

#elif LAB_NO==3

#if LAB_PART==1
		matrix x0(2, 1), s0(2, 1);
		double s = 0.1, epsilon = 1e-3, alfa, beta;
		int Nmax = 1000;
		random_device R;

		x0(0) = 2.0*R() / R.max() - 1;
		x0(1) = 2.0*R() / R.max() - 1;

		cout << x0 << endl << endl;
		alfa = 0.5;
		solution opt_HJ = HJ(x0, s, alfa, epsilon, Nmax);
		cout << opt_HJ << endl << endl;
		solution::clear_calls();
		alfa = 2;
		beta = 0.5;
		s0(0) = s;
		s0(1) = s;
		solution opt_R = Rosen(x0, s0, alfa, beta, epsilon, Nmax);
		cout << opt_R << endl << endl;
		solution::clear_calls();

#elif LAB_PART==2

		matrix x0(2, 1), s0(2, 1);
		double s = 0.5, epsilon = 1e-3, alfa, beta;
		int Nmax = 1000;
		random_device R;

		x0(0) = 10.0*R() / R.max() - 1;
		x0(1) = 10.0*R() / R.max() - 1;

		cout << x0 << endl << endl;
		alfa = 0.5;
		solution opt_HJ = HJ(x0, s, alfa, epsilon, Nmax);
		cout << opt_HJ << endl << endl;
		solution::clear_calls();
		alfa = 2;
		beta = 0.5;
		s0(0) = s;
		s0(1) = s;
		solution opt_R = Rosen(x0, s0, alfa, beta, epsilon, Nmax);
		cout << opt_R << endl << endl;
		solution::clear_calls();

#endif // LAB_PART==1 i LAB_PART == 2


#elif LAB_NO==4
	#if LAB_PART == 1

		//kara zewnetrzna
		matrix x0(2);

		double a, c = 1, dc = 2, epsilon = 1e-4;
		int Nmax = 5000;
		random_device R;

		//a przyjmuje 3 wartosci
		//wyznaczamy punkt poczatkowy
		
		//w przypadku kary zewnetrznej, caly simpleks moze lezec poza obszarem
		//w karze wewnetrznej nie moze byc poza, (przynajmniej jeden musi byc w srodku?)

		//wartosc funkcji zmodyfikowanej wynosi nieskonczonosc, dlatego nie moze (? xD)

		a = 4;
		do {
			x0(0) = 4.0 * R() / R.max() + 1;
			x0(1) = 4.0 * R() / R.max() + 1;
		} while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) > a);

			cout << x0 << endl;
			cout << sqrt(pow(x0(0), 2) + pow(x0(1), 2)) << endl;

			solution opt = pen(x0, c, dc, epsilon, Nmax, a);

			cout << opt << endl;
			cout << sqrt(pow(x0(0), 2) + pow(x0(1), 2)) << endl;

			solution::clear_calls();



	#elif LAB_PART == 2

		//wspolczynnikk c przy karze wewnterznej 

			matrix x0(2);

			double a, c = 10, dc = 0.5, epsilon = 1e-4;
			int Nmax = 10000;
			random_device R;

			//a przyjmuje 3 wartosci
			//wyznaczamy punkt poczatkowy

			//w przypadku kary zewnetrznej, caly simpleks moze lezec poza obszarem
			//w karze wewnetrznej nie moze byc poza, (przynajmniej jeden musi byc w srodku?)

			//wartosc funkcji zmodyfikowanej wynosi nieskonczonosc, dlatego nie moze (? xD)

			a = 4;
			do {
				x0(0) = 4.0 * R() / R.max() + 1;
				x0(1) = 4.0 * R() / R.max() + 1;
			} while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) > a);

			cout << x0 << endl;
			cout << sqrt(pow(x0(0), 2) + pow(x0(1), 2)) << endl;

			solution opt = pen(x0, c, dc, epsilon, Nmax, a);

			cout << opt << endl;
			cout << sqrt(pow(x0(0), 2) + pow(x0(1), 2)) << endl;

			solution::clear_calls();


	#elif LAB_PART == 3
		//problem rzeczywisty
		

#endif // LAB_PART == 1


#elif LAB_NO==5
matrix x0(2, 1), limits(2, 2);
double epsilon = 1e-3, h0;
int Nmax = 5000;
random_device R;
limits(0, 0) = limits(1, 0) = -10;
limits(0, 1) = limits(1, 1) = 10;
x0(0) = (limits(0, 1) - limits(0, 0)*R() / R.max() + limits(0, 0));
x0(1) = (limits(1, 1) - limits(1, 0)*R() / R.max() + limits(1, 0));
cout << x0 << endl << endl;

//dodatnie to stalokrokowe; ujemne - zmiennokrokowe
h0 = 0.05;
solution opt_SD = SD(x0, h0, epsilon, Nmax, limits);
cout << opt_SD << endl << endl;
solution::clear_calls();

solution opt_CG = CG(x0, h0, epsilon, Nmax, limits);
cout << opt_CG << endl << endl;
solution::clear_calls();

solution opt_Newt = Newton(x0, h0, epsilon, Nmax, limits);
cout << opt_Newt << endl << endl;
solution::clear_calls();

//wykres 4-6 jak jedna metoda z rozna dlugoscia kroku
		
#elif LAB_NO==6

#if LAB_PART == 1
	matrix x0(2, 1);
	matrix O(2, 3);
	int Nmax = 5000;
	double epsilon = 1e-3;
	random_device R;
	double w = 0;
	//w = 1;
	//ograniczenia
	O(0, 0) = O(1, 0) = -10;
	O(0, 1) = O(1, 1) = 10;

	double x0_array[101] = { 7.96527,
2.39591,
- 5.17126,
8.92462,
6.4856,
1.47486,
- 2.26168,
0.560467,
6.72734,
0.988852,
- 8.45072,
- 1.09074,
- 9.62573,
2.86431,
9.15457,
- 1.92729,
0.994045,
- 8.28892,
0.48236,
- 4.67534,
- 9.89904,
9.4151,
0.124189,
- 5.68619,
4.14937,
- 0.0740952,
- 2.29028,
5.27272,
2.8773,
- 6.50731,
8.55928,
5.81242,
- 0.436744,
0.0501724,
1.76049,
7.2017,
- 4.60202,
0.140079,
- 7.52943,
6.4443,
- 9.26969,
- 2.56712,
- 2.68219,
- 5.37014,
- 5.1124,
7.94999,
- 2.28765,
- 2.21282,
- 5.58407,
0.910586,
- 8.77296,
- 9.67457,
0.927007,
- 6.05435,
3.89766,
- 3.97773,
- 7.89917,
3.77392,
5.72313,
1.61512,
- 8.3785,
- 5.8107,
- 3.83416,
- 6.29355,
7.7401,
2.61836,
2.50513,
- 5.97435,
9.91325,
- 6.65652,
4.83507,
- 1.11553,
9.80829,
0.0439439,
- 6.2415,
- 9.8036,
5.76168,
3.67674,
3.24383,
2.6554,
- 7.73964,
5.11435,
8.47141,
7.00235,
4.00175,
- 3.76135,
- 7.69588,
- 9.66801,
0.59108,
- 2.13897,
- 4.12153,
8.19349,
9.89293,
3.40943,
0.261936,
0.48127,
- 2.22218,
- 7.51516,
- 5.02161,
- 6.43464,
- 8.96859
	};
	double x1_array[101] = { -0.874441,
9.12389,
6.05134,
- 2.41976,
4.70723,
4.94828,
- 6.76496,
- 2.71745,
6.80103,
1.05123,
5.51448,
- 7.6754,
4.12699,
4.32642,
2.14499,
- 5.30446,
- 0.730289,
- 2.53116,
7.86208,
8.47123,
9.9369,
4.55587,
4.13882,
2.11048,
3.79239,
9.63813,
7.94845,
- 3.55623,
- 9.09906,
- 6.99104,
2.06283,
8.16564,
- 4.414,
7.30409,
4.80947,
- 2.89541,
- 4.62761,
- 4.55823,
- 8.44859,
6.43147,
1.50753,
6.88343,
8.94914,
5.66581,
- 2.30908,
- 7.32836,
- 9.68588,
- 2.9695,
- 0.466833,
- 8.61815,
7.82512,
- 3.16861,
- 7.82401,
1.67261,
- 4.34793,
- 6.22094,
- 8.78577,
6.6132,
1.02491,
- 2.31942,
4.56971,
3.58901,
- 7.78739,
- 3.23362,
- 5.22658,
- 1.51339,
- 4.24152,
9.66144,
2.08578,
- 6.98224,
- 4.06927,
- 4.24527,
2.93689,
6.96404,
- 6.85191,
- 5.33478,
- 4.68382,
- 9.78456,
7.68567,
1.49867,
- 7.27811,
- 0.76402,
- 4.15111,
- 7.58756,
9.83852,
8.17682,
2.50268,
5.23172,
- 7.92024,
5.01744,
5.09797,
7.46597,
9.34297,
- 0.468231,
- 3.80247,
0.211944,
- 9.41715,
- 5.22827,
- 3.02846,
- 0.304284,
1.65036
	};

	

	//for (int i = 0; i < 101; i++) {
	//	x0(0) = (O(0, 1) - O(0, 0)) * R() / R.max() + O(0, 0);
	//	x0(1) = (O(1, 1) - O(1, 0)) * R() / R.max() + O(1, 0);
	//	x0_array[i] = x0(0);
	//	x1_array[i] = x0(1);
	//}
	
	ofstream S("a1.csv");
	for (int i = 0; i < 101; i++) {

		x0(0) = x0_array[i];
		x0(1) = x1_array[i];

		O(0, 2) = w;
		solution opt = Powell(x0, epsilon, Nmax, O);
		S << x0(0) << ";" << x0(1) << ";";
		S << opt.x(0) << ";" << opt.x(1) << ";" << opt.y(0) << ";" << opt.y(1) << ";" << opt.f_calls << ";" << endl;
		solution::clear_calls();

		w += 0.01;
	}
	S.close();


	//cout << "x0 = " << x0 << endl;

	//O(0, 2) = w;
	//solution opt = Powell(x0, epsilon, Nmax, O);
	//cout <<"opt = "<< opt << endl << endl;
	//solution::clear_calls();

#elif LAB_PART == 2
	matrix x0(2, 1);
	matrix O(2, 3);
	int Nmax = 5000;
	double epsilon = 1e-3;
	random_device R;
	double w = 0;
	//w = 0,5;
	//ograniczenia
	O(0, 0) = 0.2;
	O(0, 1) = 1;
	O(1, 0) = 0.01;
	O(1, 1) = 0.05;
	ofstream S("a_real.csv");
	
	S << "x0(0);" << "x0(1);" << "x(0);"<<"x(1);"<<"y(0);"<<"y(1);"<< endl;

	for (int i = 0; i < 101; i++) {
		x0(0) = (O(0, 1) - O(0, 0)) * R() / R.max() + O(0, 0);
		x0(1) = (O(1, 1) - O(1, 0)) * R() / R.max() + O(1, 0);

		O(0, 2) = w;

		solution opt = Powell(x0, epsilon, Nmax, O);

		cout << "opt = " << opt << endl << endl;

		S << x0(0) << ";" << x0(1) << ";" << opt.x(0) << ";" << opt.x(1) << ";" << opt.y(0) << ";" << opt.y(1) <<";"<< opt.f_calls << ";" << endl;

		solution::clear_calls();
		w += 0.01;
	}

	S.close();


#endif

#elif LAB_NO==7
	matrix O(1,1);
	matrix Ox(2,2);
	int Nmax = 5000;
	double epsilon = 1e-3;
	random_device R;
	matrix x0(2, 1);

	//ograniczenia
	Ox(0, 0) = Ox(1, 0) = -5;
	Ox(0, 1) = Ox(1, 1) = 5;

	O(0, 0) = 100; // sigma
	ofstream S("ea_100.csv");
	S << "x*(0);x*(1);y*;f_calls;" << endl;

	for (int i = 0; i < 100; i++) {
			x0(0) = (Ox(0, 1) - Ox(0, 0)) * R() / R.max() + Ox(0, 0);
			x0(1) = (Ox(1, 1) - Ox(1, 0)) * R() / R.max() + Ox(1, 0);


		solution opt = EA(2, Ox, epsilon, Nmax, O);
		S << opt.x(0) << ";" << opt.x(1) << ";" << opt.y(0) << ";" << solution::f_calls << ";" << endl;
	
		solution::clear_calls();
	}
	




#endif
	}
	catch (char * EX_INFO)
	{
		cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}

//lab3
//w przypadku funkcji testowej, 3 dlugosci kroki, po 100 razy, wyniki w tabeli 1
//w tabeli 2 tylko te ktore byly sukcesem, minimum globalne, jest wartosc 0, 

//liczyc ile razy skonczylo sie minimum globalne
// drugie minimum jest w odlgeosci przynajmniej 0.5
//dla jednego wybranego przypadku, narysowac wykres, naniesc rozwiazanie w kazdej iteracji, w HJ punkty bazowe, W rosen modyfikacja kroku, kierunku
//przed zmiana kierunku wypisywac x

//ramie robota
//tylko raz
//narysowac wykres dla k1 i k2, (P[0],P 1 )
//jak zmienia sie alfa i omega



//lab2
//olac NAN, liczyc srednie z pozostalych
//jak sa ujemne to tez uwzglednic. ( objetosc A), 10e-6


//zakladamy ze moment sily jest staly
//zakladamy ze przez 1s stan wody jest staly


// ok 24 cm2 zbiornik


