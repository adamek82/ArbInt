# Arbitrary Precision Integer Arithmetic

A C++ library and CLI tool providing arbitrary-precision signed integer arithmetic using two’s-complement arrays of 16-bit “digits.” Supports addition, subtraction, multiplication, division, modulo, I/O, and conversion from built-in types.

---

## Features

- **`arbint` class**
  - Reference-counted storage of big integers in two’s-compliment  
  - Constructors from `long`, `unsigned long`, `double`  
  - Copy-on-write assignment and destructor manage memory automatically  
  - Overloaded operators: `+`, `-`, `*`, `/`, `%`, unary `+`/`-`, comparisons  
  - Stream I/O (`<<`, `>>`) for decimal text  

- **Knuth-style algorithms**
  - Multi-digit addition/subtraction (TAOCP Vol. 2, Alg. A)  
  - Multiplication (Alg. M)  
  - Long division (Alg. D) with `dodivmod` helper  

- **Utility headers**
  - `swap.h`: type-safe inline `swap()` macro  
  - `error.h` / `error.cpp`: printf-style `error()` and `exit(1)`  

- **Test CLI**
  - `main.cpp`: reads expressions like `123 + 456` or `123 /% 10` and prints results  

---

## Project Structure

```text
/ArbInt
│
├── .gitattributes        # enforce LF eol, UTF-8 encoding
├── .gitignore            # ignore build artifacts, binaries
├── LICENSE.txt           # GNU GPL v3
│
├── ArbInt.sln            # Visual Studio solution
├── ArbInt.vcxproj        # Project file & filters
├── ArbInt.vcxproj.user   # User settings (ignored)
│
├── ArbInt.h              # arbint interface & type definitions
├── ArbInt.cpp            # arbint implementation (constructors, operators)
│
├── error.h               # error(...) prototype
├── error.cpp             # error(...) implementation
│
├── swap.h                # inline defswap(type) macro
│
└── main.cpp              # CLI test harness (parses, computes, prints)
```

---

## Building

1. **Windows / Visual Studio**  
   - Open `ArbInt.sln` in Visual Studio (2015+)  
   - Build in Debug or Release mode  

2. **Cross-platform**  
   ```bash
   g++ -std=c++11 ArbInt.cpp error.cpp main.cpp -o arbint
   ```

---

## Usage

Run the test CLI, then enter operations:

```console
$ arbint
Testing...
12345678901234567890 + 98765432109876543210 = 111111111011111111100
1000000000000000000000 - 1 = 999999999999999999999
12345 * 6789 = 83810205
100 / 3 = 33
100 % 3 = 1
100 /% 3
100 / 3 = 33
100 % 3 = 1
```

Supported operators:

- `+`   addition  
- `-`   subtraction  
- `--`  unary negate (e.g. `-- 123 = -123`)  
- `*`   multiplication  
- `/`   division  
- `%`   modulo  
- `/%`  quotient & remainder  

---

## License

This project is released under the **GNU General Public License v3.0**. See `LICENSE.txt` for full terms.
