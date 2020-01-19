//Do not edit the code below (unless you know what you are doing)

#include"solution.h"

int solution::f_calls = 0;
int solution::g_calls = 0;
int solution::H_calls = 0;

solution::solution(double L)
{
	x = matrix(L);
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(const matrix &A)
{
	x = A;
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(double *A, int n)
{
	x = matrix(A, n);
	g = NAN;
	H = NAN;
	y = NAN;
}

void solution::clear_calls()
{
	f_calls = 0;
	g_calls = 0;
	H_calls = 0;
}

ostream &operator<<(ostream &S, const solution &A)
{
	S << "x = " << A.x << endl;
	S << "y = " << A.y << endl;
	S << "f_calls = " << solution::f_calls << endl;
	S << "g_calls = " << solution::g_calls << endl;
	S << "H_calls = " << solution::H_calls << endl;
	return S;
}

//You can edit the following code

void solution::fit_fun(matrix O)
{
#if LAB_NO == 2
	#if LAB_PART == 1
		y = -cos(0.1*x(0))*exp(-pow(0.1*x(0) - 2 * 3.14, 2)) + 0.002*pow(0.1*x(0), 2);


	#elif LAB_PART == 2
		matrix Y0 = matrix(new double[3]{ 5,1,10 }, 3);
		matrix *Y = solve_ode(0, 1, 1000, Y0, x);
		int *w = get_size(Y[1]);
		double max = Y[1](0, 2);
		for (int i = 1; i < w[0]; ++i) {
			if (max < Y[1](i, 2))
				max = Y[1](i, 2);
		}
		y = abs(50 - max)
	#endif

#elif LAB_NO==3
	#if LAB_PART==1
		y = x(0) * x(0) + x(1) * x(1) - cos(2.5 *3.14 * x(0)) - cos(2.5 *3.14 * x(1)) + 2; // testowa funkcja celu
	

		//rozwiazanie: w pierwszej kolumnie polozenie (alfa), w drugiej predkosc(omega)



	#elif LAB_PART == 2
	double a_ref = 3.14, o_ref = 0;
	matrix Y0(2); //alfa i omega sa 0
	matrix *Y = solve_ode(0, 0.1, 100, Y0, x);

	//obliczanie calki
	//Y(0) zawiera czas
	//Y(1) dwie kolumny: alfa i omega
	int *n = get_size(Y[1]);
	y(0) = 0;
	for (int i = 0; i < n[0]; ++i) {
		y = y + 10 * pow(a_ref - Y[1](i, 1), 2) + pow(o_ref - Y[1](i, 1), 2) + pow(x(0)*(a_ref - Y[1](i, 0)) + x(1)*(o_ref - Y[1](i, 1)), 2);

	}
	y = y * 0.1; // pomnozenie przez dt (krok czasowy)



	#endif

#elif LAB_NO == 4
	#if LAB_PART == 1
	//metoda kary zewnetrznej
	//bierzemy kwadrat funkcji jako kara
	//pierwsza funkcja celu z konspektu

	//dodajemy kare jak jest na zewnatrz
	double arg = 3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2));
	y = sin(arg) / arg;


	//ogeraniczenia z konspektu
	//jak przekracza funkcje kary to wtedy dopiero dodajemy kare
	if (-x(0) + 1 > 0)
		y = y + O(0)*pow(-x(0) + 1, 2); // y + c*funkcja kary

	if (-x(1) + 1 > 0)
		y = y + O(0)*pow(-x(1) + 1, 2);

	if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - O(1) > 0)
		y = y + O(0)*pow(sqrt(pow(x(0), 2) + pow(x(1), 2)) - O(1), 2);

	#elif LAB_PART == 2
	//dodajemy kare jak jest wewnatrz
	
	//metoda kary wewnetrznej
	// jako kare bierzemy odwrotnosc funkcji kary
		double arg = 3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2));
		y = sin(arg) / arg;

		if (-x(0) + 1 > 0) {
			y = 1e10;
			return;
		}
		else {
			y = y - O(0) / (-x(0) + 1);
		}
		if (-x(1) + 1 > 0) {
			y = 1e10;
			return;
		}
		else {
			y = y - O(0) / (-x(1) + 1);
		}


		if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - O(1) > 0) {
			y = 1e10;
			return;
		}
		else {
			y = y - O(0) / (sqrt(pow(x(0), 2) + pow(x(1), 2)) - O(1));
		}

	#elif LAB_PART == 3
	//problem rzeczywisty

		//y = 100, bo bedzie spadac,
		//Y0(4) bedzie na minusie bo spada
		
		//predkosc pilki do koszykowki jak spada okolo 20

	

	matrix Y0(new double[4]{ 0,x(0), 100,0 }, 4); // poczatkowe polozenie x,predkosc x ,poczatkowe polozenie y, predkosc y
	matrix *Y = solve_ode(0, 0.01, 7, Y0, x(1));

	int *n = get_size(Y[1]);

	//polozenie pilki to x przy zerowym y
	//sprawdzic czy przelecial przez otwor
	//x50 musi nalezec do rpzedzialu 4,6, jesli nie to kara

	int i50 = 0;
	int i0 = 0;

	for (int i = 0; i < n[0]; i++) {
		if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50)) {
			i50 = i;
		}
		if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2))) {
			i0 = i;
		}
	}

	y = -1*Y[1](i0, 0); // szukamy maximum dlatego przemnazamy przez -1, jesli bym nie przemnozyl to znalazlby minimum

	//ograniczenia, na x, na omege, przejscie
	//jakie jest x w (50,0)

	//index i50 pokazuje w jakim miejscu jest w x = 50, sprawdzic czy nalezy do [4 ,6]

	//ogranizenia





#endif

#elif LAB_NO == 5
	int *n = get_size(O);
	if (n[1] == 1)
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2); // tu jest oryginalna funkcja celu
	else
	{
		solution temp;
		temp.x = O(0) + x*O[1];
		temp.fit_fun();
		y = temp.y;
		--f_calls; // zmniejszamy, bo wywolujemy jeszcze raz; nalezy je zwiekszac tylko obliczajac oryginalna funkcje celu

	}
#elif LAB_NO == 6
#if LAB_PART == 1
int *n = get_size(O);
if (n[1] == 1) {
	int  a = 1;
		y = matrix(2, 1);
		y(0) = a*pow(x(0) - 5, 2) + pow(x(1) - 5, 2);
		y(1) = 1.0 / a*pow(x(0) + 5, 2) + pow(x(1) + 5, 2);
}
else {
	solution temp;
	temp.x = O[0] + x*O[1];
	temp.fit_fun();
	y = O(0, 2) * temp.y(0) + (1 - O(0, 2))*temp.y(1); // srednia wazona
	--f_calls;
}

#elif LAB_PART == 2
int *n = get_size(O);
if (n[1] == 1) {
	double ro = 7800;
	double P = 1e3; //Newton
	double E = 207e9; //Pascal
	y = matrix(3, 1);

	y(0) = ro*x(0)*3.14*pow(x(1),2)/4;
	y(1) = 64 * P*pow(x(0), 3) / (3 * E*3.14*pow(x(1), 4));
	y(2) = 32 * P*x(0) / (3.14*pow(x(1), 3));

}
else {
	solution temp;
	temp.x = O[0] + x*O[1];
	temp.fit_fun();
	y = O(0, 2)*temp.y(0) + (1 - O(0, 2))*temp.y(1);
	if (temp.y(1) > 0.005)
		y = y + 1e6*pow(temp.y(1) - 0.005, 2);

	if (temp.y(2) > 300e6)
		y = y + 1e6*pow(temp.y(2) - 300e6, 2);

	--f_calls;
}


#endif
#elif LAB_NO == 7
	y = matrix(1, 1);
	y(0) = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;

#endif
	++f_calls;
}

void solution::grad(matrix O)
{
	g = matrix(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 8 * x(0) + 10 * x(1) - 38;
	++g_calls;
}

void solution::hess(matrix O)
{
	H = matrix(2, 2);
	H(0, 0) = H(1, 1) = 10;
	H(0, 1) = H(1, 0) = 8;
	++H_calls;
}
