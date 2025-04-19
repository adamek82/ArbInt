#ifndef SWAP_H
#define SWAP_H

#define defswap(typex)						\
    inline void swap(typex &a, typex &b)	\
    {										\
		typex tmp = a;						\
		a = b;								\
		b = tmp;							\
    }

// all basic types
defswap(char)
defswap(short)
defswap(int)
defswap(long)

#ifndef NO_UNSIGNED_OVERLOADING
defswap(unsigned char)
defswap(unsigned short)
defswap(unsigned int)
defswap(unsigned long)
#endif /* NO_UNSIGNED_OVERLOADING */

defswap(float)
defswap(double)

// pointers to each basic type
defswap(char *)
defswap(short *)
defswap(int *)
defswap(long *)

defswap(unsigned char *)
defswap(unsigned short *)
defswap(unsigned int *)
defswap(unsigned long *)

defswap(float *)
defswap(double *)

#endif /* SWAP_H */