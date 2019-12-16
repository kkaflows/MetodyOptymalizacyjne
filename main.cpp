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

#elif LAB_NO==7

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


