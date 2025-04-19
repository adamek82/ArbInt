#include <iomanip>
#include "ArbInt.h"
#include "error.h"

arbint::arbint()	    // arbint x; or new arbint;
{
	p = new arep;
	p->length = 1;
	p->refcnt = 1;
	p->value = new ARB_type[1];
	p->value[0] = 0;
}

arbint::arbint(const arbint &x)   // arbint x = y;
{
	x.p->refcnt++;
	p = x.p;
}

arbint::~arbint()	    // destructor or delete x;
{
	if (--p->refcnt == 0)
	{
		delete p->value;
		delete p;
	}
}

const int arbperlong = sizeof(long) / sizeof(ARB_type);

arbint::arbint(long n)		// arbint x = 36L;
{
	// allocate data area
	p = new arep;
	p->refcnt = 1;
	p->length = arbperlong;
	p->value = new ARB_type[arbperlong];

	// assign the new value
	for (int i = arbperlong - 1; i >= 0; i--)
	{
		p->value[i] = n % ARB_base;
		n /= ARB_base;
	}
}

arbint::arbint(unsigned long n)	    // arbint x = 3UL;
{
	int nlen = arbperlong;
	if ((long)n < 0)
		nlen++;

	// allocate data area
	p = new arep;
	p->refcnt = 1;
	p->length = nlen;
	p->value = new ARB_type[nlen];

	// assign the new value
	for (int i = nlen - 1; i >= 0; i--)
	{
		p->value[i] = n % ARB_base;
		n /= ARB_base;
	}
}

arbint& arbint::operator=(long n)	// x = 36L;
{
	// get rid of old data
	if (p->refcnt > 1)
	{
		p->refcnt--;
		p = new arep;
		p->refcnt = 1;
		p->length = arbperlong;
		p->value = new ARB_type[arbperlong];
	}

	// delete old value
	else if (p->length != arbperlong)
	{
		delete p->value;
		p->length = arbperlong;
		p->value = new ARB_type[arbperlong];
	}

	// assign the new value
	for (int i = arbperlong - 1; i >= 0; i--)
	{
		p->value[i] = n % ARB_base;
		n /= ARB_base;
	}
	return *this;
}

arbint& arbint::operator=(unsigned long n)  // x = 3UL;
{
	int nlen = arbperlong;
	if (((long)n) < 0)
		nlen++;

	// get rid of old data
	if (p->refcnt > 1)
	{
		p->refcnt--;
		p = new arep;
		p->refcnt = 1;
		p->length = arbperlong;
		p->value = new ARB_type[arbperlong];
	}

	// reallocate old value
	else if (p->length != nlen)
	{
		delete p->value;
		p->length = nlen;
		p->value = new ARB_type[nlen];
	}

	// assign the new value
	for (int i = nlen - 1; i >= 0; i--)
	{
		p->value[i] = n % ARB_base;
		n /= ARB_base;
	}
	return *this;
}

/*
	Assignment operator for type arbint.
*/
arbint& arbint::operator=(const arbint &x)	// x = y;
{
	x.p->refcnt++;
	if (--p->refcnt == 0)
	{
		delete p->value;
		delete p;
	}

	p = x.p;
	return *this;
}

void arbint::dassign(arep * p, double n)
{
	p->refcnt = 1;

	/* declare working storage for dtoarb */
#   define BITS(type)	(CHAR_BIT * (int)sizeof(type))
	const int dlen = DBL_MAX_EXP / BITS(ARB_type);
	ARB_type s[dlen + 1];

	/* convert the double to positive */
	int neg = 0;
	if (n < 0)
	{
		n = -n;
		neg = 1;
	}

	// convert the double
	double i = n;
	ARB_type *s1 = (s + neg);	// run forward on digits

	// this leaves the digits in reverse order
	do {
		*s1++ = (ARB_type)fmod(i, ARB_base);
		i /= ARB_base;
	} while (i >= 1.0);

	// run forward and backward on the digits reversing them
	p->length = (int)(s1-- - (s + neg));
	ARB_type *s2 = s;
	while (s2 < s1)
		swap(*s2++, *s1--);

	// if necessary, convert to negative
	if (neg)
	{
		ARB_Ltype k = 1;
		s[0] = 0;
		p->length++;
		for (int j = p->length - 1; ; )
		{
			ARB_type l1 = ~s[j];
			ARB_Ltype l = l1 + k;
			s[j] = (ARB_type)l;
			if (--j <= 0)
				break;
			k = (l / ARB_base) ? 1 : 0;
		}
	}

	// copy to allocated space
	p->value = new ARB_type[p->length];
	memcpy((char *)p->value, (char *)s, p->length * sizeof(ARB_type));
}

arbint::arbint(ARB_type *w, int wlen)
{
	p = new arep;
	const ARB_type hibit = (ARB_base >> 1);

	// Leave i pointing to first normal digit.
	int i = 0, j = 1;
	if (w[0] & hibit)	// negative
	{
		// set up a sentinel of non-0, non-~0, but negative
		ARB_type save_w_n = w[wlen - 1];
		w[wlen - 1] = hibit;

		// Increment while digit is all sign bit and the next digit has a sign bit...
		for (; (w[i] == ~0) && (w[j] & hibit); i++, j++)
			;

		// restore sentinel digit
		w[wlen - 1] = save_w_n;
	}

	else		// positive
	{
		// set up a sentinel of non-0, non-~0, but positive
		ARB_type save_w_n = w[wlen - 1];
		if (save_w_n == 0)
			w[wlen - 1] = 1;

		// Increment while digit is all 0 and the next digit is not
		for (; (w[i] == 0) && !(w[j] & hibit); i++, j++)
			;

		// restore sentinel digit
		w[wlen - 1] = save_w_n;
	}

	p->length = wlen - i;
	p->refcnt = 1;

	// copy the data
	if (i == 0)
		p->value = w;

	else
	{
		p->value = new ARB_type[p->length];
		memcpy((char*)p->value, (char*)&w[i], p->length * sizeof(ARB_type));
		delete w;
	}
}

/*
	Add u[1..n] + v[1..n] to form w[0..n]

	Based on:

	The Art of Computer Programming, volume 2
	D. Knuth, Section 4.3.1, Algorithm A
*/
arbint operator+(const arbint& u, const arbint& v)
{
	int ulen = u.p->length;
	int vlen = v.p->length;
	int wlen = max(ulen, vlen) + 1;
	ARB_type *wv = new ARB_type[wlen];
	ARB_type *uv = u.p->value;
	ARB_type *vv = v.p->value;

	/*
		A1 [Initialize]
			set j <- n
			k <- 0
			{ modification: uj for u, vj for v }

		A3(a) [Loop on j]
			decrease j by one
			{ modification: decrease uj and vj }
	*/
	ARB_Ltype k = 0;
	int uj, vj, wj;
	for (uj = ulen - 1, vj = vlen - 1, wj = wlen - 1;
		uj >= 0 && vj >= 0;
		uj--, vj--, wj--)
	{
		/*
			A2(a) [Add digits]
				set w[j] <- (u[j] + v[j] + k) mod b
				{ modification: w[wj] <- (u[uj] + v[vj] + k) mod b }
		*/
		ARB_Ltype t = uv[uj] + vv[vj] + k;
		wv[wj] = (ARB_type)t;	// % ARB_base;

		/*
			A2(b)
				k <- (u[j] + v[j] + k) / b
		*/
		k = (t / ARB_base) ? 1 : 0;
	}

	const ARB_type hibit = (ARB_base >> 1);
	const ARB_type u_sign = uv[0] & hibit;
	const ARB_type v_sign = vv[0] & hibit;

	/*
		NOTE: either (uj >= 0) or (vj >= 0), but not both. As long as the carry is
		non-zero, it must be added to the next digit.
	*/

	if (uj >= 0)
	{
		/*
			Take care of u[1..uj].
		*/
		const ARB_type v_rest = v_sign ? ~0 : 0;
		for (; uj >= 0; uj--, wj--)
		{
			ARB_Ltype t = uv[uj] + v_rest + k;
			wv[wj] = (ARB_type)t;		// % ARB_base;
			k = (t / ARB_base) ? 1 : 0;
		}

		/*
			If there are any digits left to be filled in (there should be at most 1)
			they must be filled with either 0 or ~0 depending on the sign of wv[wj+1]
		*/
		if (wj >= 0)
		{
			ARB_type rest = 0;
			// overflow never occurs when adding two numbers with different signs
			if (u_sign != v_sign)
				rest = (wv[wj + 1] & hibit) ? ~0 : 0;
			// if both both numbers have the same sign, then overflow occurs if and only if the result has the opposite sign
			else
				rest = k ? ~0 : 0;
			do {
				wv[wj--] = rest;
			} while (wj >= 0);
		}
	}

	else if (vj >= 0)
	{
		/*
			Take care of v[1..vj]
		*/
		const ARB_type u_rest = u_sign ? ~0 : 0;
		for (; vj >= 0; vj--, wj--)
		{
			ARB_Ltype t = vv[vj] + u_rest + k;
			wv[wj] = (ARB_type)t;		// % ARB_base;
			k = (t / ARB_base) ? 1 : 0;
		}

		/*
			If there are any digits left (there should be at most 1)
			they must be filled with either 0 or ~0 depending on the sign of wv[wj+1]
		*/
		if (wj >= 0)
		{
			ARB_type rest = 0;
			// overflow never occurs when adding two numbers with different signs
			if (u_sign != v_sign)
				rest = (wv[wj + 1] & hibit) ? ~0 : 0;
			// if both both numbers have the same sign, then overflow occurs if and only if the result has the opposite sign
			else
				rest = k ? ~0 : 0;
			do {
				wv[wj--] = rest;
			} while (wj >= 0);
		}
	}

	else
	{
		/*
			If there is a digit left (there should be at most 1), it should be set to the sign of w[wj]
			when adding two numbers with different signs.
			If both both numbers have the same sign, a carry here means that the number must be
			sign-extended;

			A3(b)
				Set w[0] <- k
				{ modification, w[0] <- k ? sign(w[wj]) : 0	if sign(u) <> sign(v) }
				{               w[0] <- sign (w[wj]) otherwise					  }
		*/
		ARB_type rest = 0;
		// overflow never occurs when adding two numbers with different signs
		if (u_sign != v_sign)
			rest = (wv[wj + 1] & hibit) ? ~0 : 0;
		// if both both numbers have the same sign, then overflow occurs if and only if the result has the opposite sign
		else
			rest = k ? ~0 : 0;
		while (wj >= 0)
			wv[wj--] = rest;
	}

	/* Normalize and return */
	arbint w(wv, wlen);
	return w;
}

arbint operator+(const arbint& j)
{
	int nlen = j.p->length;
	ARB_type *nv = new ARB_type[nlen];
	memcpy((char*)nv, (char*)j.p->value, nlen * sizeof(ARB_type));
	arbint k(nv, nlen);
	return k;
}

arbint operator-(const arbint& u)
{
	int ulen = u.p->length;
	int wlen = ulen + u.isneg();
	ARB_type *wv = new ARB_type[wlen];
	ARB_type *uv = u.p->value;

	/*
		A1 [Initialize]
			set j <- n
			k <- 0	becomes	k <- 1

		A3(a) [Loop on j]
			decrease j by one
	*/
	ARB_Ltype k = 1;
	int uj, wj;
	for (uj = ulen - 1, wj = wlen - 1; uj >= 0; uj--, wj--)
	{
		/*
			A2(a) [Add digits]
				set w[j] <- (u[j] + v[j] + k) mod b
				becomes
				set w[j] <- (~u[j] + k) mod b
		*/
		ARB_type l1 = ~uv[uj];
		ARB_Ltype l = l1 + k;
		wv[wj] = (ARB_type)l;	// % ARB_base

		/*
			A2(b)
				k <- (u[j] + v[j] + k) / b
		*/
		k = (l / ARB_base) ? 1 : 0;
	}

	/*
		A3(b)
			Set w[0] <- 0
	*/
	if (wj == 0)
		wv[0] = 0;

	/* Normalize and return */
	arbint w(wv, wlen);
	return w;
}

arbint operator-(const arbint& u, const arbint& v)
{
	int ulen = u.p->length;
	int vlen = v.p->length;
	int wlen = max(ulen, vlen) + 1;
	ARB_type *wv = new ARB_type[wlen];
	ARB_type *uv = u.p->value;
	ARB_type *vv = v.p->value;

	/*
		A1 [Initialize]
			set j <- n
			k <- 0
			{ modification: uj for u, vj for v
			  k <- 1 }

		A3(a) [Loop on j]
			decrease j by one
			{ modification: decrease uj and vj }
	*/
	ARB_Ltype k = 1;
	int uj, vj, wj;
	for (uj = ulen - 1, vj = vlen - 1,
		wj = wlen - 1;
		uj >= 0 && vj >= 0;
		uj--, vj--, wj--)
	{
		/*
			A2(a) [Add digits]
				set w[j] <- (u[j] + v[j] + k) mod b
				{ modification: w[wj] <- (u[uj] + ~v[vj] + k) mod b }
		*/
		ARB_Ltype t = uv[uj];
		ARB_type t1 = ~vv[vj];
		t += t1;
		t += k;
		wv[wj] = (ARB_type)t;	// % ARB_base

		/*
			A2(b)
				k <- (u[j] + v[j] + k) / b
		*/
		k = (t / ARB_base) ? 1 : 0;
	}

	const ARB_type hibit = (ARB_base >> 1);
	const ARB_type u_sign = uv[0] & hibit;
	const ARB_type v_sign = vv[0] & hibit;

	/*
		NOTE: either (uj >= 0) or (vj >= 0), but not both. As long as the carry is
		non-zero, it must be added to the next digit.
	*/

	if (uj >= 0)
	{
		const ARB_type v_rest = v_sign ? 0 : ~0;
		/*
			Take care of u[1..uj].
		*/
		for (; uj >= 0; uj--, wj--)
		{
			ARB_Ltype t = uv[uj] + v_rest + k;
			wv[wj] = (ARB_type)t;		// % ARB_mask;
			k = (t / ARB_base) ? 1 : 0;
		}

		/*
			If there are any digits left (there should be at most 1)
			they must be filled with either 0 or ~0 depending on the sign of wv[wj+1]
		*/
		if (wj >= 0)
		{
			ARB_type rest = 0;
			if (u_sign == v_sign)
				rest = (wv[wj + 1] & hibit) ? ~0 : 0;
			else
				rest = k ? ~0 : 0;
			do {
				wv[wj--] = rest;
			} while (wj >= 0);
		}
	}

	else if (vj >= 0)
	{
		const ARB_type u_rest = u_sign ? ~0 : 0;
		/*
			Take care of v[1..vj]
		*/
		for (; vj >= 0; vj--, wj--)
		{
			ARB_Ltype t = u_rest + (ARB_type)~vv[vj] + k;
			wv[wj] = (ARB_type)t;		// % ARB_base
			k = (t / ARB_base) ? 1 : 0;
		}

		/*
			If there are any digits left (there should be at most 1)
			they must be filled with either 0 or ~0 depending on the sign of wv[wj+1]
		*/
		if (wj >= 0)
		{
			ARB_type rest = 0;
			if (u_sign == v_sign)
				rest = (wv[wj + 1] & hibit) ? ~0 : 0;
			else
				rest = k ? ~0 : 0;
			do {
				wv[wj--] = rest;
			} while (wj >= 0);
		}
	}

	else
	{
		/*
			If there is a digit left (there should be at most 1), it should be set to k.
			Because we are dealing with 2's complement, a carry here means that the number must be
			sign-extended; otherwise the digit is set to 0.

			A3(b)
				Set w[0] <- k
				{ modification, w[0] <- k ? sign(w[wj]) : 0	if sign(u) == sign(v) }
				{               w[0] <- sign (w[wj]) otherwise					  }
		*/
		ARB_type rest = 0;
		if (u_sign == v_sign)
			rest = (wv[wj + 1] & hibit) ? ~0 : 0;
		else
			rest = k ? ~0 : 0;
		while (wj >= 0)
			wv[wj--] = rest;
	}

	/* Normalize and return */
	arbint w(wv, wlen);
	return w;
}

arbint operator*(const arbint& tmp_u, const arbint& tmp_v)
{
	// Make u (the multiplicand) and v (the multiplier) positive.
	int negu = tmp_u.isneg();
	int negv = tmp_v.isneg();
	arbint u(+tmp_u);
	if (negu) u = -tmp_u;
	arbint v(+tmp_v);
	if (negv) v = -tmp_v;

	// The product will be negative if the signs of the multiplier 
	// and multiplicand do not match.
	int negans = (negu != negv);

	int n = u.p->length;
	int m = v.p->length;
	int wlen = n + m;
	ARB_type *wv = new ARB_type[wlen];
	ARB_type *uv = u.p->value;
	ARB_type *vv = v.p->value;

	/*
		M1(a) [Initialize]
			Set w[m+1..m+n] <- 0
	*/
	memset((char*)&wv[m], 0, n * sizeof(ARB_type));

	/*
		M1(b) [Initialize]
			Set j <- m
		M6 [Loop on j]
			decrease j by one
			if j > 0, goto M2
	*/
	for (int j = m - 1; j >= 0; j--)
	{
		/*
			M2 [zero multiplier?]
				if v[j] == 0
					w[j] <- 0
					goto M6
		*/
		if (vv[j] == 0)
		{
			wv[j] = 0;
			continue;
		}

		/*
			M3 [Initialize i]
				set i <- n,
					k <- 0
			M5(a) [Loop on i]
				decrease i by one
				if i > 0, goto M4
		*/
		ARB_Ltype v_j = vv[j];
		ARB_Ltype k = 0;
		for (int i = n - 1, iplusj = i + j + 1;
			i >= 0;
			i--, iplusj--)

		{
			/*
				M4 [multiply and add]
					set t <- u[i] * v[j] + w[i+j] + k
					w[i+j] <- t % b
					k <- t / b
			*/
			ARB_Ltype t = uv[i] * v_j +	wv[iplusj] + k;
			wv[iplusj] = (ARB_type)t;		// % ARB_base;
			k = t / ARB_base;
		}

		/*
			M5(b) [Loop on i]
				if i <= 0,
					set w[j] <- k
		*/
		wv[j] = (ARB_type)k;
	}

	/* Normalize */
	arbint w(wv, wlen);

	/* Restore sign and return */
	if (negans)
		return -w;
	else
		return w;
}

int arb_cmp(const ARB_type *l1, const ARB_type *l2, int n)
{
	if (l1 != l2)
		while (--n >= 0)
			if (*l1++ != *l2++)
				return (l1[-1] > l2[-1] ? 1 : -1);
	return (0);
}

void dodivmod(const arbint &tmp_u, const arbint &tmp_v,	arbint &quotient, arbint &remainder)
{
	//	Make u (the dividend) and v (the divisor) positive.
	int negu = tmp_u.isneg();
	int negv = tmp_v.isneg();
	arbint u(+tmp_u);
	if (negu) u = -tmp_u;
	arbint v(+tmp_v);
	if (negv) v = -tmp_v;

	/*
		The sign of the remainder will match the sign of the dividend.
		The sign of the quotient will be negative if the sign of the divisor
		and dividend do not match, else	positive.
	*/
	int negremainder = negu;
	int negquotient = (negu != negv);

	// Set local variables.
	ARB_type *uv = u.p->value;
	ARB_type *vv = v.p->value;
	int m_n = u.p->length;	// m + n
	int n = v.p->length;
	int m = m_n - n;

	/*
		For n == 1, use simpler algorithm
	*/
	if (n == 1)
	{
		if (vv[0] == 0)
		{
			error("division by zero!");
			quotient = remainder = +u;
		}

		else
		{
			ARB_type *qv = new ARB_type[m_n];
			ARB_Ltype prevu = 0;
			ARB_type v1 = vv[0];

			for (int r = 0; r < m_n; r++)
			{
				ARB_Ltype t = uv[r] + prevu * ARB_base;
				ARB_type tmpq = ARB_type(t / v1);
				qv[r] = tmpq;	// % ARB_base
				prevu = t - v1 * tmpq;
			}

			arbint q(qv, m_n);
			if (negquotient)
				quotient = -q;
			else
				quotient = q;

			if (negremainder)
				remainder = -(arbint)prevu;
			else
				remainder = prevu;
		}
		return;
	}

	/*
		Degenerate case of length(u) < length(v) i.e., m < 0, implying that u < v
	*/
	else if (m_n < n)
	{
		quotient = 0L;
		if (negremainder)
			remainder = -u;
		else
			remainder = +u;
		return;
	}

	/*
		Degenerate case of length(u) == length(v) i.e., m == 0, possibly implying that
		u < v or u == v
	*/
	else if (m_n == n)
	{
		int cmp = arb_cmp(uv, vv, m_n);
		if (cmp < 0)	// u < v
		{
			quotient = 0L;
			if (negremainder)
				remainder = -u;
			else
				remainder = +u;
			return;
		}

		else if (cmp == 0)	// u == v
		{
			if (negquotient)
				quotient = -1L;
			else
				quotient = 1L;
			remainder = 0L;
			return;
		}
	}

	/*
		Now call out all of the guns from Knuth
	*/

	// In rare circumstances, the first digit of u or v may be 0. This digit must now be discarded.
	while ((*uv == 0) && (m_n > 1))
	{
		uv++; m_n--;
	}
	while ((*vv == 0) && (n > 1))
	{
		vv++; n--;
	}
	m = m_n - n;

	/*
		D1(a) [Normalize.]
			Set d <- b/(v1+1).
	*/
	ARB_type d = (ARB_type)(ARB_base / (vv[0] + 1));

	/*
		D1(b) [Normalize.]
			Set (u[0]u[1]...u[m+n]) to (u[1]u[2]...u[m+n]) * d
	*/
	ARB_type *Ouv = uv;
	ARB_type k;
	uv = new ARB_type[m_n + 1];
	if (d == 1)
	{
		// copy old value
		uv[0] = 0;
		memcpy((char*)&uv[1], (char*)&Ouv[0], m_n * sizeof(ARB_type));
	}

	else
	{
		// multiply u by d
		k = 0;
		int Oi, i;
		for (Oi = m_n - 1, i = m_n;
			i > 0;
			Oi--, i--)
		{
			ARB_Ltype t = Ouv[Oi] * d + k;
			uv[i] = (ARB_type)t;		// % ARB_base;
			k = (ARB_type)(t / ARB_base);
		}
		uv[i] = k;
	}

	/*
		D1(c) [Normalize.]
			Set (v[1]v[2]...v[n]) to (v[1]v[2]...v[n]) * d
	*/
	ARB_type *Ovv = vv;
	vv = new ARB_type[n];
	if (d == 1)
	{
		// copy old value
		memcpy((char*)&vv[0], (char*)&Ovv[0], n * sizeof(ARB_type));
	}

	else
	{
		// multiply v by d
		k = 0;
		int Oi, i;
		for (Oi = n - 1, i = n - 1;
			i >= 0;
			Oi--, i--)
		{
			ARB_Ltype t = Ovv[Oi] * d + k;
			vv[i] = (ARB_type)t;		// % ARB_base;
			k = (ARB_type)(t / ARB_base);
		}
	}

	/*
		D2 [Initialize j]
			Set j <- 0

		D7 [Loop on j]
			Set j <- j+1
			Loop if j <= m
	*/
	ARB_type v1 = vv[0];
	ARB_type v2 = vv[1];
	ARB_type *qv = new ARB_type[m + 1];
	ARB_type *nv = new ARB_type[n + 1];

	for (int j = 0; j <= m; j++)
	{
		/*
			D3 [Calculate q^]
				If u[j] == v[1]
					Set q^ <- base - 1
				Else
					Set q^ <- (u[j] * base + u[j+1]) / v[1]
					While v[2] * q^ > (u[j] * base + u[j+1] - q^ * v[1]) * base + u[j+2]
						Set q^ <- q^ - 1
		*/

		ARB_Ltype q_hat;

		if (uv[j] == v1)
			q_hat = ARB_base - 1;
		else
			q_hat = (uv[j] * ARB_base + uv[j + 1]) / v1;

		ARB_Ltype u_j = uv[j];
		ARB_Ltype u_j1 = uv[j + 1];
		ARB_Ltype u_j2 = uv[j + 2];

		for (; ; q_hat--)
		{
			/*
				if ((v2 * q_hat) <=	((u_j * ARB_base + u_j1 - v1 * q_hat) * ARB_base + u_j2))
					q^--;
			*/
			ARB_Ltype u_j_q_hat = u_j * ARB_base + u_j1;
			u_j_q_hat -= v1 * q_hat;
			if (u_j_q_hat / ARB_base != 0)
				break;

			u_j_q_hat *= ARB_base;
			u_j_q_hat += u_j2;
			if ((v2 * q_hat) <= u_j_q_hat)
				break;
		}

		/*
			D4 [Multiply and subtract.]
				Replace u[j..j+n] by u[j..j+n] - q^ * v[1..n]
		*/
		// set nv <- q^ * (v[1..n])
		k = 0;
		int dl = n, vl = n - 1;
		for (; vl >= 0; vl--, dl--)
		{
			ARB_Ltype t = vv[vl] * q_hat + k;
			nv[dl] = (ARB_type)t;		// % ARB_base;
			k = (ARB_type)(t / ARB_base);
		}
		nv[0] = k;

		// subtract nv[0..n] from u[j..j+n]
		int borrow = 0;
		int ul = j + n;
		for (dl = n; dl >= 0; dl--, ul--)
		{
			ARB_Ltype t = uv[ul] - nv[dl] - borrow;
			uv[ul] = (ARB_type)t;		// % ARB_base
			borrow = (t / ARB_base) ? 1 : 0;
		}

		qv[j] = (ARB_type)q_hat;
		if (borrow != 0)
		{
			for (k = 0, ul = j + n, vl = n - 1;
				vl >= 0;
				vl--, ul--)
			{
				ARB_Ltype t = uv[ul] + vv[vl] + k;
				uv[ul] = (ARB_type)t;	// % ARB_base
				k = (ARB_type)(t / ARB_base);
			}

			uv[j] += k;
			qv[j]--;
		}
	}

	/*
		D8 [Unnormalize]
			q[0..m] is quotient
			u[m+1..m+n] / d is remainder
	*/
	arbint q(qv, m + 1);
	if (negquotient)
		quotient = -q;
	else
		quotient = q;

	// divide u[m+1..m+n] by d
	ARB_type *rem = new ARB_type[n];

	if (d == 1)		// nothing special to do
		memcpy((char*)rem, (char*)&uv[m + 1], n * sizeof(ARB_type));
	else
	{
		ARB_Ltype prevu = 0;

		// do division by single digit
		for (int rl = 0, ul = m + 1;
			ul <= m_n;
			ul++, rl++)
		{
			ARB_Ltype t = uv[ul] + prevu * ARB_base;
			ARB_type tmpq = (ARB_type)(t / d);
			rem[rl] = tmpq;	// % ARB_base
			prevu = t - d * tmpq;
		}
	}

	arbint r(rem, n);
	if (negremainder)
		remainder = -r;
	else
		remainder = r;

	delete uv;
	delete vv;
	delete nv;
}

arbint operator%(const arbint &u, const arbint &v)
{
	arbint quot, rem;
	dodivmod(u, v, quot, rem);
	return rem;
}

arbint operator/(const arbint &u, const arbint &v)
{
	arbint quot, rem;
	dodivmod(u, v, quot, rem);
	return quot;
}


ostream& operator<<(ostream& out, const arbint& j)
{
	arbint val = +j;

	// reverse the sign of negative numbers
	if (j.isneg())
	{
		val = -j;
		out << "-";
	}

	// output can be stored in 6*val.length hyper-decimal digits, each digit holding a value 0..9999.
	const ARB_type dec_base = 10000;
	ARB_type *digs = new ARB_type[val.p->length * ARB_base / dec_base];
	ARB_type *svdigs = digs;
	ARB_type *val_val = val.p->value;
	int val_len = val.p->length;

	// store "digits" in reverse order
	for (;;)
	{
		// *digs++ = val % 10000
		// val /= 10000;
		ARB_Ltype prevu = 0;

		for (int r = 0; r < val_len; r++)
		{
			ARB_Ltype tmp = val_val[r] + prevu * ARB_base;
			ARB_type tmpq = (ARB_type)(tmp / dec_base);
			val_val[r] = tmpq;		// % ARB_base
			prevu = tmp - dec_base * tmpq;
		}
		*digs++ = (ARB_type)prevu;
		
		// if (val == 0)
		//    break;
		int allzero = 1;
		for (int r = 0; r < val_len; r++)
			if (val_val[r] != 0)
			{
				allzero = 0;
				break;
			}

		if (allzero)
			break;
	}

	// Output digits in forward order. All but first digit must be 
	// 4 decimal digits in lengths, including leading zeros.
	out << *--digs;
	while (digs-- > svdigs)
		out << setw(4) << setprecision(4) << setfill('0') << (unsigned short)*digs;

	delete svdigs;
	return out;
}

const int chunksize = 5;

// convert hexadecimal 0-9, a-f, A-F to its numeric value
inline int hexvalue(char c)
{
	if (isdigit(c))
		return c - '0';
	else if (c >= 'a' && c <= 'f')
		return c - 'a' + 10;
	else
		return c - 'A' + 10;
}

// return true if c is an octal digit
inline int isodigit(int c)
{
	return isdigit(c) && c < '8';
}

static ARB_type *i_s;
static int i_len;

/*
	Add to the length of i_s by allocated new space, copying the old data to the new space, and deleting the
	old space. The new space is initialized to all zeros, except for the last digit whose value is passed in.
*/
static void extend_s(int k)
{
	int olen = i_len;
	i_len += chunksize;

	ARB_type *new_i_s = new ARB_type[i_len];
	memcpy((char*)&new_i_s[chunksize], (char*)i_s, olen * sizeof(ARB_type));

	memset((char*)new_i_s, 0, 4 * sizeof(ARB_type));
	new_i_s[4] = (ARB_type)k;

	delete i_s;
	i_s = new_i_s;
}

// perform: i = i * base - newdigit
static void addonedigit(int multiplier, int addend)
{
	// i *= base
	// almost identical to that used for operator*
	ARB_Ltype k = 0;
	for (int l = i_len - 1; l >= 0; l--)
	{
		ARB_Ltype t = i_s[l] * multiplier + k;
		i_s[l] = (ARB_type)t;	// % ARB_base
		k = t / ARB_base;
	}

	// add digits, if necessary
	if ((k != 0) || (i_s[0] & 0x8000))
		extend_s(k);

	// i += addend
	// almost identical to that used for operator+
	k = addend;
	for (int l = i_len - 1; l >= 0; l--)
	{
		ARB_Ltype t = i_s[l] + k;
		i_s[l] = (ARB_type)t;
		k = (t / ARB_base) ? 1 : 0;
	}

	// add digits, if necessary
	if ((k != 0) || (i_s[0] & 0x8000))
		extend_s(k);
}

istream& operator>> (istream& in, arbint& i)
{
	char c, sign = 0;
	in >> ws;		// skip white space
	if (!in.get(c))
	{
		i = 0L;
		return in;
	}

	switch (c)		// remember the sign
	{
	case '+': case '-':
		sign = c;
		if (!in.get(c))
		{
			i = 0L;
			return in;
		}
		break;
	}

	// allocate the initial data space
	i_s = new ARB_type[chunksize];
	i_len = chunksize;
	memset((char*)i_s, 0, chunksize * sizeof(ARB_type));

	switch (c)
	{
	case '0':		// hex or octal
		if (!in.get(c))
			break;

		switch (c)
		{
		case 'x':		// hex number
		case 'X':
			for (in.get(c);
				in && isxdigit(c);
				in.get(c))
				addonedigit(16, hexvalue(c));
			break;

		default:		// octal number
			for (;
				in && isodigit(c);
				in.get(c))
				addonedigit(8, (c - '0'));
			break;
		}
		break;

	default:		// decimal number
		for (; in && isdigit(c); in.get(c))
			addonedigit(10, (c - '0'));
		break;
	}

	if (in)
		in.putback(c);

	// change sign, if necessary, and save data
	arbint ret(i_s, i_len);
	if (sign == '-')
		i = -ret;
	else
		i = ret;

	return in;
}
