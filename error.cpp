#include <iostream>
#include <cstdarg>
#include <cstdlib>

using namespace std;

void error(const char *fmt ...)
{
	va_list ap;
	va_start(ap, fmt);
	char ch;

	// loop across format string
	while (ch = *fmt++)
	{
		// output normal chars
		if (ch != '%')
			cerr.put(ch);
		else
		{
			// found % sequence
			switch (ch = *fmt++)
			{
			// %% becomes a single %
			case '%':
				cerr.put('%');
				break;

			// output string
			case 's':
			{
				char *s = va_arg(ap, char*);
				cerr << s;
			}
			break;

			// output decimal integer
			case 'd':
			{
				int s = va_arg(ap, int);
				cerr << s;
			}
			break;

			// output character
			case 'c':
			{
				int s = va_arg(ap, int);
				cerr.put((char)s);
			}
			break;

			default:
				cerr << "\nunknown % sequence: %" << char(ch) << "\n";
				break;
			}
		}
	}
	// all done
	va_end(ap);
	exit(1);
}
