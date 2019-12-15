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
		matrix Y0 = matrix(new double[3]{ 5,1,10 }, 3); //5 m3 w A, 1 m3 w B , 10 stopni celsjusza w B
		matrix *Y = solve_ode(0, 1, 1000, Y0); // t0, krok_czasu, czas_symulacji
		ofstream S("sim_t.csv");
		S << Y[0];
		S.close();
		S.open("sim_y.csv");
		S << Y[1];
		S.close();

#elif LAB_NO==2
#if LAB_PART == 1
		double x0, d = 1, alfa = 1.410, epsilon = 1e-5, gamma = 1e-200;
		int Nmax = 1000;
		random_device R;
		/*ofstream file("outputFile.txt");
		for (int i = 0; i < 100; i++) {
			x0 = 200.0*R() / R.max() - 100;	
			file << x0 << ";";
			double *ab = expansion(x0, d, alfa, Nmax);
			file << ab[0] << ';' << ab[1] << ";" << ab[2]<< ";";
			solution::clear_calls();

			solution opt_F = fib(ab[0], ab[1], epsilon);
			file << opt_F.x << opt_F.y << solution::f_calls << ";" << ";";
			solution::clear_calls();

			solution opt_L = lag(ab[0], ab[1], epsilon, gamma, Nmax);
			file << opt_L.x << opt_L.y << solution::f_calls << ";" << ";"<<endl;
			solution::clear_calls;
		}
		file.close();*/
		/*x0 = 200.0*R() / R.max() - 100;
		cout << x0 << endl << endl;
		double *ab = expansion(x0, d, alfa, Nmax);
		cout << ab[0] << '\t' << ab[1] << endl << endl;
		solution::clear_calls();*/

		double ab[2]{ -100,100 };
		solution opt_F = fib(ab[0], ab[1], epsilon);
		cout << opt_F << endl << endl;
		solution::clear_calls();

		ab[0] = -100;
		ab[1] = 100;
		solution opt_L = lag(ab[0], ab[1], epsilon, gamma, Nmax);
		cout << opt_L << endl << endl;
		solution::clear_calls;
#elif LAB_PART == 2
		double x0, d = 1e-4, alfa = 1.5, epsilon = 1e-10, gamma = 1e-200;
		int Nmax = 1000;
		random_device R;
		x0 = (1e-2 - 1e-4)*R() / R.max() + 1e-4;
		//x0 = 0.0016;
		cout << x0 << endl << endl;
		double *ab = expansion(x0, d, alfa, Nmax);
		cout << ab[0] << '\t' << ab[1] << endl << endl;
		solution::clear_calls();

		solution opt_F = fib(ab[0], ab[1], epsilon);
		cout << opt_F << endl << endl;
		solution::clear_calls();

		solution opt_L = lag(ab[0], ab[1], epsilon, gamma, Nmax);
		cout << opt_L << endl << endl;
		solution::clear_calls;
#endif

#elif LAB_NO==3
#if LAB_PART == 1

		ofstream File("opt34_test.csv");
		locale loc("");
		File.imbue(loc);
		File << "x0;x1;MHJ;;;;;MR;;;;\n";
		File << "x0;x1;x1*;x2*;y*;lwywolan;min;x1*;x2*;y*;lwywolan;min\n";

		matrix x0(2,1);
		int Nmax = 1000;
		double alfa, beta, epsilon = 1e-3, s = 0.136558;
		random_device R;

		int fcalls1, fcalls2;

		for (int i = 0; i < 1; i++)
		{
			x0(0) = 2.0*R() / R.max() - 1;
			x0(1) = 2.0*R() / R.max() - 1;

			x0(0) = -0.413967;
			x0(1) = -0.994215;



			cout << x0 << endl << endl;
			alfa = 0.5;
			solution opt_HJ = HJ(x0, s, alfa, epsilon, Nmax);
			fcalls1 = solution::f_calls;
			File << x0(0) << ";" << x0(1) << ";" << opt_HJ.x(0) << ";" << opt_HJ.x(1) << ";" << opt_HJ.y << fcalls1 << ";;";

			cout << opt_HJ << endl << endl;
			solution::clear_calls();
			cout << endl;

			alfa = 2;
			beta = 0.5;
			matrix s0(2, 1);
			s0(0) = s;
			s0(1) = s;
			solution opt_R = Rosen(x0, s0, alfa, beta, epsilon, Nmax);
			fcalls2 = solution::f_calls;
			File << opt_R.x(0) << ";" << opt_R.x(1) << ";" << opt_R.y << fcalls2 << ";\n";

			cout << opt_R << endl << endl;
			solution::clear_calls();
			cout << endl;

			cout << endl << endl;
			/*cout << x0 << opt_HJ;
			cout << fit1 << "  ";
			cout << opt_R;
			cout << fit2 << endl;*/
		}

#elif LAB_PART == 2
		matrix x0(2, 1);
		int Nmax = 1000;
		double alfa, beta, epsilon = 1e-3, s = 0.32658;
		random_device R;
		x0(0) = 10.0*R() / R.max();
		x0(1) = 10.0*R() / R.max();
		cout << x0 << endl << endl;
		alfa = 0.5;
		solution opt_HJ = HJ(x0, s, alfa, epsilon, Nmax);
		cout << opt_HJ << endl << endl;
		
		solution::clear_calls();

		alfa = 2;
		beta = 0.5;
		matrix s0(2, 1);
		s0(0) = s;
		s0(1) = s;
		solution opt_R = Rosen(x0, s0, alfa, beta, epsilon, Nmax);
		cout << opt_R << endl << endl;
		solution::clear_calls();

#elif LAB_PART == 3
		matrix Y0 = matrix(new double[2]{ 0,0 }, 2); //polozenie 0, predkoscramienia 0
		matrix P = matrix(new double[2]{ 2.7688, 3.18129 }, 2); //k1, k2
		matrix* Y = solve_ode(0, 0.1, 100, Y0, P); // t0, krok_czasu, czas_symulacji
		ofstream S("sim_y.csv");
		S << Y[1];
		S.close();

		int Nmax = 1000;
		double alfa, beta, epsilon = 1e-3, s = 0.1;
		random_device R;

		S.open("Hpart3.csv");
		alfa = 0.5;
		solution opt_HJ = HJ(Y0, s, alfa, epsilon, Nmax);
		cout << opt_HJ << endl << endl;
		S << opt_HJ.x(0) << ";" << opt_HJ.x(1) << ";" << opt_HJ.y << ";;";

		S.close();


#endif

#elif LAB_NO==4

	#if LAB_PART == 1

	//kara zewnetrzna
	matrix x0(2);

	double a, c = 1, dc = 2, epsilon = 1e-4;
	double a_wew, c_wew = 10, dc_wew = 0.5, epsilon_wew = 1e-4;
	int Nmax = 5000;
	random_device R;

	//a przyjmuje 3 wartosci
	//wyznaczamy punkt poczatkowy

	//w przypadku kary zewnetrznej, caly simpleks moze lezec poza obszarem
	//w karze wewnetrznej nie moze byc poza, (przynajmniej jeden musi byc w srodku?)

	//wartosc funkcji zmodyfikowanej wynosi nieskonczonosc, dlatego nie moze (? xD)
	ofstream S("kara_zewnetrzna_a_4.csv");
	S << "x0(0);x0(1);x_opt(0);x_opt(1);r;y;f_calls;" << endl;

	ofstream S2("kara_wewnetrzna_a_4.csv");
	S2 << "x0(0);x0(1);x_opt(0);x_opt(1);r;y;f_calls;" << endl;

	for (int i = 0; i <= 100; i++) {

		a = 4;
		do {
			x0(0) = 4.0 * R() / R.max() + 1;
			x0(1) = 4.0 * R() / R.max() + 1;
		} while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) > a);

		//kara zewnetrzna
		solution opt = pen(x0, c, dc, epsilon, Nmax, a);

		double r_opt = sqrt(pow(opt.x(0), 2) + pow(opt.x(1), 2));

		//wpisywanie do pliku
		S << x0(0) << ";" << x0(1) << ";" << opt.x(0) << ";" << opt.x(1) << ";" << r_opt << ";" << opt.y <<";"<<solution::f_calls<<";"<< endl;

		solution::clear_calls();

		//kara wewnetrzna
		solution opt_wew = pen(x0, c_wew, dc_wew, epsilon_wew, Nmax, a);
		double r_opt_wew = sqrt(pow(opt_wew.x(0), 2) + pow(opt_wew.x(1), 2));
		S2 << x0(0) << ";" << x0(1) << ";" << opt_wew.x(0) << ";" << opt_wew.x(1) << ";" << r_opt_wew << ";" << opt_wew.y << ";" << solution::f_calls << ";" << endl;


		solution::clear_calls();

	}




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


	
	
#endif

#elif LAB_NO==5
		
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