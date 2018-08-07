// bigInt.h
#ifndef BIGINT_H_
#define BIGINT_H_

#include <iostream>

namespace BIGINT
{
	class BigInt
	{
		static const long long MAX_DIGITS = 1000;
		enum SIGN { POSITIVE, NEGATIVE };

	private:
		char digits[MAX_DIGITS];
		SIGN sign = POSITIVE;
	public:

		// constructors
		BigInt(long long someInt = 0);
		BigInt(const BigInt & bigInt);

		// operators
		BigInt & operator=(const BigInt & a);
		friend BigInt power(const BigInt & a, const BigInt & b);
		friend BigInt factorial(const BigInt & a);
		friend BigInt operator+(const BigInt & a, const BigInt & b);
		friend BigInt operator+(const BigInt & a, long long b) { return a + BigInt(b); }
		friend BigInt operator+(long long a, const BigInt & b) { return BigInt(a) + b; }
		friend BigInt operator-(const BigInt & a, const BigInt & b);
		friend BigInt operator*(const BigInt & a, const BigInt & b);
		friend BigInt operator/(const BigInt & a, const BigInt & b);
		friend BigInt operator/(const BigInt & a, long long b) { return a / BigInt(b); }
		friend BigInt operator/(long long a, const BigInt & b) { return b / a; }
		friend BigInt operator%(const BigInt & a, const BigInt & b) { return a - (a / b) * b; }
		friend std::ostream & operator<<(std::ostream & os, const BigInt & a);

		// comparators
		friend bool operator>=(const BigInt & a, const BigInt & b);
		friend bool operator<=(const BigInt & a, const BigInt & b) { return b >= a; }
		friend bool operator==(const BigInt & a, const BigInt & b);
		friend bool operator>(const BigInt & a, const BigInt & b) { return (a >= b && !(a == b)); }
		friend bool operator<(const BigInt & a, const BigInt & b) { return (b > a); }

	};

	BigInt::BigInt(long long someInt)
	{
		long long temp = someInt;

		if (temp < 0)
		{
			sign = NEGATIVE;
			temp *= -1;
		}
		else
			sign = POSITIVE;

		for ( long long i = 0; i < MAX_DIGITS; ++i)
			digits[i] = '0';

		 long long j = 0;
		while (temp > 0 && j < MAX_DIGITS)
		{
			digits[j] = temp % 10 + '0';
			temp /= 10;
			++j;
		}
	}

	BigInt::BigInt(const BigInt & bigInt)
	{
		for ( long long i = 0; i < MAX_DIGITS; ++i)
			digits[i] = bigInt.digits[i];
		this->sign = bigInt.sign;
	}

	BigInt & BigInt::operator=(const BigInt & a)
	{
		for ( long long i = 0; i < MAX_DIGITS; ++i)
			digits[i] = a.digits[i];
		this->sign = a.sign;
		return *this;
	}

	BigInt operator+(const BigInt & a, const BigInt & b)
	{
		if (a == 0)
			return BigInt(b);
		if (b == 0)
			return BigInt(a);

		if (a.sign == BigInt::POSITIVE && b.sign == BigInt::POSITIVE)
		{
			BigInt c;
			bool carry = false;
			for ( long long i = 0; i < BigInt::MAX_DIGITS; ++i)
			{
				int x = (a.digits[i] - '0') + (b.digits[i] - '0');
				if (carry)
				{
					++x;
					carry = false;
				}
				if (x >= 10)
				{
					carry = true;
					x -= 10;
				}
				c.digits[i] = x + '0';
			}
			return c;
		}
		else if (a.sign == BigInt::POSITIVE && b.sign == BigInt::NEGATIVE)
		{
			BigInt check = b;
			check.sign = BigInt::POSITIVE;
			if (check > a)
			{
				// compute b-a, make it negative, and return it
				BigInt c;
				bool carry = false;
				for ( long long i = 0; i < BigInt::MAX_DIGITS; ++i)
				{
					int x = b.digits[i] - a.digits[i];
					if (carry)
					{
						--x;
						carry = false;
					}
					if (x < 0)
					{
						carry = true;
						x += 10;
					}
					c.digits[i] = x + '0';
				}
				c.sign = BigInt::NEGATIVE;
				return c;
			}
			else if (check == a)
				return BigInt(0);
			else
			{
				// compute a-b, make it positive, and return it
				BigInt c;
				bool carry = false;
				for ( long long i = 0; i < BigInt::MAX_DIGITS; ++i)
				{
					int x = a.digits[i] - b.digits[i];
					if (carry)
					{
						--x;
						carry = false;
					}
					if (x < 0)
					{
						carry = true;
						x += 10;
					}
					c.digits[i] = x + '0';
				}
				c.sign = BigInt::POSITIVE;
				return c;
			}
		}
		else if (a.sign == BigInt::NEGATIVE && b.sign == BigInt::POSITIVE)
		{
			return b + a;
		}
		else if (a.sign == BigInt::NEGATIVE && b.sign == BigInt::NEGATIVE)
		{
			BigInt c = a, d = b;
			c.sign = BigInt::POSITIVE;
			d.sign = BigInt::POSITIVE;
			BigInt e = c + d;
			e.sign = BigInt::NEGATIVE;
			return e;
		}
	}

	BigInt operator-(const BigInt & a, const BigInt & b)
	{
		BigInt c = b;
		if (b.sign == BigInt::POSITIVE)
		{
			c.sign = BigInt::NEGATIVE;
			return a + c;
		}
		else
		{
			c.sign = BigInt::POSITIVE;
			return a + c;
		}
	}

	BigInt operator*(const BigInt & a, const BigInt & b)
	{

		BigInt ret = 0;

		if ((a.sign == BigInt::POSITIVE) && (b.sign == BigInt::POSITIVE))
			ret.sign = BigInt::POSITIVE;
		else if ((a.sign == BigInt::NEGATIVE) && (b.sign == BigInt::POSITIVE))
			ret.sign = BigInt::NEGATIVE;
		else if ((a.sign == BigInt::POSITIVE) && (b.sign == BigInt::NEGATIVE))
			ret.sign = BigInt::NEGATIVE;
		else if ((a.sign == BigInt::NEGATIVE) && (b.sign == BigInt::NEGATIVE))
			ret.sign = BigInt::POSITIVE;

		// if either a or b are zero, then we return 0
		if (a == 0 || b == 0)
			return BigInt(0);

		// if either a or b are one, then we return the other
		if (a == 1)
			return (ret = b);
		if (b == 1)
			return (ret = a);

		// integer multiplication

		// take a digit x from B
		// multiply each digit of A by x
		// and store the result plus any carries
		// into a temporary variable C
		// Note that you don't want to lose the carry over
		// when the loop terminates!
		// Add C to a running sum
		// take another digit, but this time
		// store the digits in C shifted one
		// digit
		// Repeat this procedure until all the
		// digits of B have been used

		// get the start of a
		long long astart = BigInt::MAX_DIGITS - 1;
		while (a.digits[astart] == '0') --astart;

		// get the start of b
		long long bstart = BigInt::MAX_DIGITS - 1;
		while (b.digits[bstart] == '0') --bstart;

		BigInt C = 0;
		for (long long i = 0; i <= bstart; ++i)
		{
			int x = b.digits[i] - '0', carry = 0;
			long long j = 0;
			for (j = 0; j <= astart; ++j)
			{
				int t = (a.digits[j] - '0') * x;
				if (carry)
				{
					t += carry;
					carry = 0;
				}

				if (t >= 10)
				{
					carry = (t / 10);
				}
				C.digits[i + j] = t % 10 + '0';
			}

			// don't lose the last carry factor
			// j is incremented after the loop ends
			if (carry && i + j < BigInt::MAX_DIGITS)
				C.digits[i + j] = carry + '0';

			ret = (ret + C);
			C = 0;
		}

		return ret;
	}

	/*
	BigInt operator*(const BigInt & a, const BigInt & b)
	{
		// slow integer multiplication
		BigInt temp = a, ret = 0;

		if ((a.sign == BigInt::POSITIVE) && (b.sign == BigInt::POSITIVE))
			ret.sign = BigInt::POSITIVE;
		else if ((a.sign == BigInt::NEGATIVE) && (b.sign == BigInt::POSITIVE))
			ret.sign = BigInt::NEGATIVE;
		else if ((a.sign == BigInt::POSITIVE) && (b.sign == BigInt::NEGATIVE))
			ret.sign = BigInt::NEGATIVE;
		else if ((a.sign == BigInt::NEGATIVE) && (b.sign == BigInt::NEGATIVE))
			ret.sign = BigInt::POSITIVE;

		while (temp > 0)
		{
			ret = ret + b;
			temp = temp - 1;
		}

		return ret;
	}*/

	BigInt operator/(const BigInt & a, const BigInt & b)
	{

		BigInt::SIGN pos;

		// set return sign value
		if ((a.sign == BigInt::POSITIVE) && (b.sign == BigInt::POSITIVE))
			pos = BigInt::POSITIVE;
		else if ((a.sign == BigInt::NEGATIVE) && (b.sign == BigInt::POSITIVE))
			pos = BigInt::NEGATIVE;
		else if ((a.sign == BigInt::POSITIVE) && (b.sign == BigInt::NEGATIVE))
			pos = BigInt::NEGATIVE;
		else if ((a.sign == BigInt::NEGATIVE) && (b.sign == BigInt::NEGATIVE))
			pos = BigInt::POSITIVE;

		if ((b > a) || (b == 0))
			return 0;

		BigInt ret;
		if (b == 1)
		{
			ret = a;
			ret.sign = pos;
			return ret;
		}

		// long division

		// get start of a
		long long i = BigInt::MAX_DIGITS - 1, astart = 0;
		while (a.digits[i] == '0')
			--i;
		astart = i;

		// get start of b
		long long j = BigInt::MAX_DIGITS - 1, bstart = 0;
		while (b.digits[j] == '0')
			--j;
		bstart = j;


		// get a sizable chunk of A starting from i
		BigInt x;
		x.digits[0] = a.digits[i];

		while (i >= 0)
		{
			while (b > x)
			{
				if (i == 0)
				{
					// we're finished with the division process here
					--i;
					break;
				}

				// get a larger chunk until b <= x
				x = x * 10;
				x.digits[0] = a.digits[--i];
			}

			// check to see if we're finished
			if (i < 0)
				break;

			// once we have obtained a large enough chunk,
			// we can divide b into x and record the result
			int count = 0;
			while (x >= b)
			{
				x = (x - b);
				++count;
			}
			ret.digits[i] = count + '0';

			// we can continue to get another chunk and
			// repeat the process
		}

		ret.sign = pos;
		return ret;
	}

	/*BigInt operator/(const BigInt & a, const BigInt & b)
	{
		if ((b > a) || (b == 0))
			return BigInt(0);

		// slow integer division
		BigInt c = a, ret = 0;

		if ((a.sign == BigInt::POSITIVE) && (b.sign == BigInt::POSITIVE))
			ret.sign = BigInt::POSITIVE;
		else if ((a.sign == BigInt::NEGATIVE) && (b.sign == BigInt::POSITIVE))
			ret.sign = BigInt::NEGATIVE;
		else if ((a.sign == BigInt::POSITIVE) && (b.sign == BigInt::NEGATIVE))
			ret.sign = BigInt::NEGATIVE;
		else if ((a.sign == BigInt::NEGATIVE) && (b.sign == BigInt::NEGATIVE))
			ret.sign = BigInt::POSITIVE;

		while (c >= b)
		{
			c = c - b;
			ret = ret + 1;
		}
			
		return ret;
	}*/

	BigInt power(const BigInt & a, const BigInt & b)
	{

		// slow power function

		if (b.sign == BigInt::NEGATIVE)
			return 0;
		
		BigInt c = 1;
		for (BigInt x = b; x > 0; x = x - 1)
			c = c * a;

		return c;
	}


	BigInt factorial(const BigInt & a)
	{
		if (a < 1)
			return BigInt(-1);

		// slow factorial function
		BigInt ret = 1, temp = a;
		while (temp > 1)
		{
			ret = ret * temp;
			temp = temp - 1;
		}
		return ret;
	}

	std::ostream & operator<<(std::ostream & os, const BigInt & a)
	{

		if (a == 0)
		{
			os << '0';
			return os;
		}

		if (a.sign == BigInt::NEGATIVE)
			os << "-";

		// get the beginning of the integer
		long long k = BigInt::MAX_DIGITS - 1;
		while (a.digits[k] == '0' && k >= 0)
			--k;

		for (long long i = k; i >= 0; --i)
			os << a.digits[i];

		return os;
	}

	bool operator>=(const BigInt & a, const BigInt & b)
	{
		// get largest digit of a and b
		long long i = BigInt::MAX_DIGITS - 1, j = BigInt::MAX_DIGITS - 1;

		while (a.digits[i] == '0' && i >= 0)
			--i;

		while (b.digits[j] == '0' && j >= 0)
			--j;

		if (i > j)
			return true;
		else if (i == j)
		{
			while (i >= 0)
			{
				if (a.digits[i] > b.digits[i])
					return true;
				else if (a.digits[i] < b.digits[i])
					return false;
				--i;
			}
			return true;
		}

		return false;
	}

	bool operator==(const BigInt & a, const BigInt & b)
	{
		long long i = 0;
		while (i < BigInt::MAX_DIGITS)
		{
			if (a.digits[i] != b.digits[i])
				return false;
			++i;
		}
		return true;
	}

}


#endif
