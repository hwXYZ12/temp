// Triple.h - euler488
#ifndef TRIPLE_H_
#define TRIPLE_H_

namespace TRIPLE
{
	template <class T1, class T2, class T3>
	class Triple
	{
	private:
		T1 a;
		T2 b;
		T3 c;
	public:
		T1 & first();
		T2 & second();
		T3 & third();
		T1 first() const		{ return a; }
		T2 second() const		{ return b; }
		T3 third() const		{ return c; }
		Triple(const T1 & aval, const T2 & bval, const T3 & cval) : a(aval), b(bval), c(cval)	{}
		Triple()	{}

		//template <class T1, class T2>
		//friend bool operator==(const Pair &, const Pair &);
	};

	template <class T1, class T2, class T3>
	T1 & Triple<T1, T2, T3>::first()
	{
		return a;
	}

	template <class T1, class T2, class T3>
	T2 & Triple<T1, T2, T3>::second()
	{
		return b;
	}

	template <class T1, class T2, class T3>
	T3 & Triple<T1, T2, T3>::third()
	{
		return c;
	}

}


#endif
