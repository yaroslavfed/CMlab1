#include<iostream>
#include<vector>
#include<fstream>
#include<iomanip>
#include<cmath>


using namespace std;

typedef float type;
typedef double sum;
const type admissible_precision{ (type)1e-15 };

#pragma region FunctionInitialization

void input(string filename, vector<vector<type>>& a);

template <typename T>
void input(string filename, vector<T>& v);

inline void input(string filename, size_t& N);

type operator*(vector<type>& a, vector<type>& b);

void LUS(size_t& N, vector<type>& di, vector<size_t>& ai, vector<type>& au, vector<type>& al);

void forwardStep(size_t& N, vector<size_t>& ai, vector<type>& al, vector<type>& y, vector<type>& b);

void backStep(size_t& N, vector<size_t>& ai, vector<type>& au, vector<type>& di, vector<type>& x, vector<type>& y);

void Gauss(size_t& N, vector<vector<type>>& a, vector<type>& b, vector<type>& x);

void Gilbert(size_t& k, vector<size_t>& ai, vector<type>& al, vector<type>& au, vector<type>& di, vector<type>& b);

template <typename T>
void write(string filename, const vector<T>& v);
#pragma endregion


int main()
{
	int c;
	cout << "1 - LU, 2 - Gilbert, 3 - Gauss" << "\n";
	cin >> c;

	switch (c)
	{
	case(1):	//LU*
		try
		{
			size_t N{ 0 };
			input("size.txt", N);
			vector<type> di(N), al, au, b(N), x(N), y(N);
			vector<size_t> ai(N + 1);
			input("di.txt", di);
			input("ai.txt", ai);
			al.resize(ai[N]);
			au.resize(ai[N]);
			input("au.txt", au);
			input("al.txt", al);
			input("b.txt", b);
			type k{ 1e-15 };
			di[0] += k;
			b[0] += k;
			LUS(N, di, ai, au, al);
			forwardStep(N, ai, al, y, b);
			backStep(N, ai, au, di, x, y);
			write("result.txt", x);
		}
		catch (string error) { cout << error << endl; }
		break;


	case(2):	//Gilbert
		try
		{
			size_t N{ 15 };
			vector<type> di(N), al, au, b(N), x(N), y(N);
			vector<size_t> ai(N + 1);
			Gilbert(N, ai, al, au, di, b);
			LUS(N, di, ai, au, al);
			forwardStep(N, ai, al, y, b);
			backStep(N, ai, au, di, x, y);
			write("result.txt", x);
		}
		catch (string error) { cout << error << endl; }
		break;


	case(3):	//Gauss
		try
		{
			size_t N{ 0 };
			input("size.txt", N);
			vector<type> b(N), x(N);
			input("b.txt", b);
			vector<vector<type>> a(N);
			input("a.txt", a);
			type k{ 1e-15 };
			Gauss(N, a, b, x);
			write("result.txt", x);
		}
		catch (string error) { cout << error << endl; }
		break;


	default:
		cout << "Unknown command";
		break;
	}
}




void input(string filename, vector<vector<type>>& a)
{
	ifstream in(filename);
	if (in.is_open())
	{
		size_t N = a.size();
		for (size_t i{ 0 }; i < N; i++)
			for (size_t j{ 0 }; j < N; j++)
			{
				a[i].resize(N);
				in >> a[i][j];
			}
		in.close();
		return;
	}
	throw string("InputError");
}


template <typename T>
void input(string filename, vector<T>& v)
{
	ifstream in(filename);

	if (in.is_open())
	{
		size_t N = v.size();
		for (size_t i{ 0 }; i < N; i++)
			in >> v[i];
		in.close();
		return;
	}

	throw string("InputError");
}

inline void input(string filename, size_t& N)
{
	ifstream in(filename);

	if (in.is_open())
	{
		in >> N;
		in.close();
		return;
	}

	throw string("InputError");
}

type operator*(vector<type>& a, vector<type>& b)
{
	if (a.size() != b.size()) throw string("size a != size b");

	size_t size = a.size();
	type res{ 0 };

	for (size_t i{ 0 }; i < size; i++)
		res += a[i] * b[i];

	return res;
}


void LUS(size_t& N, vector<type>& di, vector<size_t>& ai, vector<type>& au, vector<type>& al)
{
	if (!N) throw string("Empty!");

	for (size_t i{ 0 }; i < N; i++)
	{
		if (abs(di[i]) < admissible_precision) throw string("DataError");

		sum suml{ 0 }, sumu{ 0 }, sumd{ 0 };										 // Переменные суммирования
		size_t j0{ i - (ai[i + 1] - ai[i]) };										 // Начало i профиля нижнего треугольника

		for (size_t j{ j0 }; j < i; j++)
		{
			suml = 0.0, sumu = 0.0;
			size_t i0{ j - (ai[j + 1] - ai[j]) };									 // Начало j профиля верхнего треугольника
			for (size_t shift{ i0 > j0 ? i0 : j0 }; shift < j; shift++)     // Выбираем максимально отдаленный профиль от 0 индекса
			{
				suml += au[ai[j] + shift - i0] * al[ai[i] + shift - j0];    // Ходим по столбцам в матрице
				sumu += au[ai[i] + shift - j0] * al[ai[j] + shift - i0];	 // Ходим транспонированно
			}
			al[ai[i] + j - j0] -= suml;                                     // (j - j0) - номер элемента в профиле
			al[ai[i] + j - j0] /= di[j];
			au[ai[i] + j - j0] -= sumu;
		}

		for (size_t p{ j0 }; p < i; p++)
			sumd += au[ai[i] + p - j0] * al[ai[i] + p - j0];                // Для диагонали считаем отдельно
		di[i] -= sumd;
	}
}


void forwardStep(size_t& N, vector<size_t>& ai, vector<type>& al, vector<type>& y, vector<type>& b)
{
	if (!N) throw string("Empty!");

	for (size_t i{ 0 }; i < N; i++)                                         // Ходим по строкам
	{
		sum suml{ 0 };
		size_t j0{ i - (ai[i + 1] - ai[i]) };

		for (size_t k{ j0 }; k < i; k++)
			suml += y[k] * al[ai[i] + k - j0]; 		// Берем элементы i профиля

		y[i] = b[i] - suml;
	}
}


void backStep(size_t& N, vector<size_t>& ai, vector<type>& au, vector<type>& di, vector<type>& x, vector<type>& y)
{
	if (!N) throw string("Empty!");

	for (size_t i{ 0 }; i < N; i++)
	{
		if (abs(di[i]) < admissible_precision) throw string("DataError");
		x[i] = y[i];
	}

	for (int j{ (int)N - 1 }; j >= 0; j--)                                  // Номер столбца
	{
		x[j] /= di[j];
		size_t i0{ j - (ai[j + 1] - ai[j]) };                               // Индекс 1го элемента в столбце
		for (int i{ j - 1 }; i >= (int)i0; i--)
			x[i] -= x[j] * au[ai[j + 1] + i - j];
	}
}


void Gauss(size_t& N, vector<vector<type>>& a, vector<type>& b, vector<type>& x)
{
	type max{ 0 };
	size_t index{ 0 };

	//Прямой ход
	for (size_t j{ 0 }; j < N - 1; j++)
	{
		max = a[j][j];
		index = j;

		for (size_t i{ j }; i < N; i++)
			if (abs(max) < abs(a[i][j]))
			{
				index = i;
				max = a[i][j];
			}

		if (abs(max) < admissible_precision) throw string("Data Error");

		if (index > j)
		{
			type tmp{ 0 };
			for (size_t k{ 0 }; k < N; k++)
			{
				tmp = a[j][k];
				a[j][k] = a[index][k];
				a[index][k] = tmp;
			}
			tmp = b[j];
			b[j] = b[index];
			b[index] = tmp;
		}

		for (size_t k{ j }; k < N; k++)
			a[j][k] /= max;

		b[j] /= max;

		for (size_t i{ j + 1 }; i < N; i++)
		{
			type dev{ a[i][j] };
			if (dev)
			{
				for (size_t k{ 0 }; k < N; k++)
					a[i][k] -= dev * a[j][k];
				b[i] -= dev * b[j];
			}
		}
	}

	//Обратный ход
	for (size_t i{ 0 }; i < N; i++)
		x[i] = 0;

	for (int i{ (int)N - 1 }; i >= 0; i--)
	{
		if (abs(a[i][i]) < admissible_precision) throw string("Data Error");

		for (size_t j{ 0 }; j < N; j++)
		{
			x[i] -= a[i][j] * x[j];
		}
		x[i] += b[i];
		x[i] /= a[i][i];
	}
}


void Gilbert(size_t& k, vector<size_t>& ai, vector<type>& al, vector<type>& au, vector<type>& di, vector<type>& b)
{
	if (!k) throw string("Empty!");
	vector<size_t> x_real(k);
	for (size_t i{ 0 }; i < k; i++)
		x_real[i] = i + 1;
	type tmp{ 0 };
	ai[0] = 0;
	for (size_t i{ 0 }; i < k; i++)
	{
		ai[i + 1] = i > 0 ? i + ai[i] : 0;
		for (size_t j{ 0 }; j < k; j++)
		{
			size_t t{ i + j + 1 };
			tmp = (type)1 / t;
			if (i == j) di[i] = tmp;
			if (j < i) al.push_back(tmp);
		}
	}

	size_t size{ al.size() };
	au.resize(size);

	for (size_t i{ 0 }; i < size; i++)
		au[i] = al[i];

	for (size_t i{ 0 }; i < k; i++)
		b[i] = 0;

	for (size_t i{ 0 }; i < k; i++)
	{
		for (size_t j{ 0 }; j < i; j++)
			b[i] += al[ai[i] + j] * x_real[j];
		b[i] += di[i] * x_real[i];
		for (size_t j{ i + 1 }; j < k; j++)
			b[i] += au[ai[j] + i] * x_real[j];
	}
}


template <typename T>
void write(string filename, const vector<T>& v)
{
	size_t N = v.size();
	ofstream out(filename);

	if (out.is_open())
		for (size_t i{ 0 }; i < N; i++)
			out << scientific << setprecision(15) << v[i] << endl;

	out.close();
}