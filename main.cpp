#include <iostream>
#include "ArbInt.h"

using namespace std;

inline int eq(const char*s1, const char*s2)
{
	return (strcmp(s1, s2) == 0);
}

void testadd(const arbint& A, const arbint& B)
{
	arbint sum = A + B;
	cout << A << " + " << B << " = " << sum << "\n";
}

void testsub(const arbint& A, const arbint& B)
{
	arbint diff = A - B;
	cout << A << " - " << B << " = " << diff << "\n";
}

void testneg(const arbint& A)
{
	arbint neg = -A;
	cout << " 0, 0 - " << A << " = " << neg << "\n";
}

void testmul(const arbint& A, const arbint& B)
{
	arbint prod = A * B;
	cout << A << " * " << B << " = " << prod << "\n";
}

void testdiv(const arbint& A, const arbint& B)
{
	arbint quot = A / B;
	cout << A << " / " << B << " = " << quot << "\n";
}

void testmod(const arbint& A, const arbint& B)
{
	arbint rem = A % B;
	cout << A << " % " << B << " = " << rem << "\n";
}

void testdivmod(const arbint& A, const arbint& B)
{
	arbint quot, rem;
	dodivmod(A, B, quot, rem);
	cout << A << " / " << B << " = " << quot << "\n";
	cout << A << " % " << B << " = " << rem << "\n";
}


int main()
{
	cout << "Testing..." << endl;
	/*
	cout << "sizeof(unsigned short) = " << sizeof(unsigned short) << endl;
	cout << "sizeof(unsigned int) = " << sizeof(unsigned int) << endl;
	cout << "sizeof(unsigned long) = " << sizeof(unsigned long) << endl;
	*/

	arbint x, y;
	char op[100];

	//cin >> x;
	//cout << x << endl;

	while (cin >> x >> op >> y)
	{
		if (eq(op, "+"))
			testadd(x, y);
		else if (eq(op, "-"))
			testsub(x, y);
		else if (eq(op, "--"))
			testneg(x);
		else if (eq(op, "*"))
			testmul(x, y);
		else if (eq(op, "/")) 
			testdiv(x, y);
		else if (eq(op, "%")) 
			testmod(x, y);
		else if (eq(op, "/%")) 
			testdivmod(x, y);
	}
	

	return 0;
}
