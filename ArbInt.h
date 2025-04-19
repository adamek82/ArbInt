/*
	arbitrary precision arithmetic

	All numbers are stored in an array of unsigned shorts in two's complement format. Its length is all maintained separately.
*/
#ifndef ARBINT_H
#define ARBINT_H

#include <iostream>
#include <minmax.h>
#include <string.>
//#include <error.h>

using namespace std;

typedef unsigned short ARB_type;
typedef unsigned long ARB_Ltype;
const ARB_Ltype ARB_base = 0x10000;

class arbint
{
	struct arep
	{
		ARB_type *value;	// ptr to value
		int length;			// length of data
		int refcnt;			// reference count
	};

	arep *p;				// ptr to data

	arbint(ARB_type*, int);
	friend void dodivmod(const arbint&, const arbint&, arbint&, arbint&);
	int isneg() const
	{
		return (p->value[0] & (ARB_base >> 1)) != 0;
	}
	friend int arb_cmp(const ARB_type*, const ARB_type*, int);

public:
	arbint();			     // arbint x; or
							 // new arbint;
	~arbint();			     // delete arbint;

	arbint(const arbint&);	 // arbint x = y;
	arbint(long);		     // arbint x = 35L;
	arbint(unsigned long);	 // arbint x = 35LU;
	arbint(double);		     // arbint x = 35.0;

	arbint& operator=(const arbint&);	// x = y;
	arbint& operator=(long);			// x = 35L;
	arbint& operator=(unsigned long);	// x = 35LU;
	arbint& operator=(double);			// x = 35.0;

	static void dassign(arep *p, double n);

	friend arbint operator+(const arbint&, const arbint&);
	friend arbint operator-(const arbint&, const arbint&);
	friend arbint operator*(const arbint&, const arbint&);
	friend arbint operator/(const arbint&, const arbint&);
	friend arbint operator%(const arbint&, const arbint&);
	friend arbint operator+(const arbint&);
	friend arbint operator-(const arbint&);

	friend int operator==(const arbint&, const arbint&);
	friend int operator!=(const arbint&, const arbint&);
	friend int operator<(const arbint&, const arbint&);
	friend int operator>=(const arbint&, const arbint&);
	friend int operator>(const arbint&, const arbint&);
	friend int operator<=(const arbint&, const arbint&);

	friend ostream& operator<<(ostream&, const arbint&);
	friend istream& operator>>(istream&, arbint&);
};

#endif
