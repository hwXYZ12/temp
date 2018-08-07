// project euler problem 17 - euler17.cpp
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <tuple>
#include "Triple.h"
#include "RBTree.h"
#include <string>
#include "BigInt.h"

typedef long long INT;
typedef TRIPLE::Triple<INT, INT, INT> INTTriple;
typedef std::vector<INTTriple> VectorOfTriples;
typedef RBTREE::RBTree<INTTriple> TreeOfINTTriples;
typedef TRIPLE::Triple<INT, TreeOfINTTriples *, TreeOfINTTriples *> INTandTwoTrees;
typedef RBTREE::RBTree<INTandTwoTrees> TreeOfINTandTwoTrees;

static TreeOfINTandTwoTrees losingPositions[3] = {
		TreeOfINTandTwoTrees(INTandTwoTrees(0,
											new TreeOfINTTriples(INTTriple(0, 1, 2), [](const INTTriple & x, const INTTriple & y){return x.second() < y.second(); },
																					[](const INTTriple & x, const INTTriple & y){return x.second() == y.second(); }),
											new TreeOfINTTriples(INTTriple(0, 1, 2), [](const INTTriple & x, const INTTriple & y){return x.third() < y.third(); },
																					[](const INTTriple & x, const INTTriple & y){return x.third() == y.third(); })),
							[](const INTandTwoTrees & x, const INTandTwoTrees & y) { return x.first() < y.first(); },
							[](const INTandTwoTrees & x, const INTandTwoTrees & y) { return x.first() == y.first(); } ),
		TreeOfINTandTwoTrees(INTandTwoTrees(1,
											new TreeOfINTTriples(INTTriple(0, 1, 2), [](const INTTriple & x, const INTTriple & y){return x.first() < y.first(); },
																					[](const INTTriple & x, const INTTriple & y){return x.first() == y.first(); }),
											new TreeOfINTTriples(INTTriple(0, 1, 2), [](const INTTriple & x, const INTTriple & y){return x.third() < y.third(); },
																					[](const INTTriple & x, const INTTriple & y){return x.third() == y.third(); })),
							[](const INTandTwoTrees & x, const INTandTwoTrees & y) { return x.first() < y.first(); },
							[](const INTandTwoTrees & x, const INTandTwoTrees & y) { return x.first() == y.first(); } ),
		TreeOfINTandTwoTrees(INTandTwoTrees(2,
											new TreeOfINTTriples(INTTriple(0, 1, 2), [](const INTTriple & x, const INTTriple & y){return x.first() < y.first(); },
																					[](const INTTriple & x, const INTTriple & y){return x.first() == y.first(); }),
											new TreeOfINTTriples(INTTriple(0, 1, 2), [](const INTTriple & x, const INTTriple & y){return x.second() < y.second(); },
																					[](const INTTriple & x, const INTTriple & y){return x.second() == y.second(); })),
							[](const INTandTwoTrees & x, const INTandTwoTrees & y) { return x.first() < y.first(); },
							[](const INTandTwoTrees & x, const INTandTwoTrees & y) { return x.first() == y.first(); } ),
												};

static INTandTwoTrees dummy(0, new TreeOfINTTriples(INTTriple(0, 0, 0), 
												[](const INTTriple & x, const INTTriple & y){return x.second() < y.second(); }, // don't even use these, just for the compiler
												[](const INTTriple & x, const INTTriple & y){return x.second() == y.second(); }),
							   new TreeOfINTTriples(INTTriple(0, 0, 0), 
												[](const INTTriple & x, const INTTriple & y){return x.third() < y.third(); },
												[](const INTTriple & x, const INTTriple & y){return x.third() == y.third(); }));
static INTTriple dummyINTTriple(0, 0, 0);

// helper method to efficiently search losing positions
// that have already been determined
INTTriple * searchLosingPositions(INT key1, INT key2);

// helper to determine losing positions
void generateLosingPositions(INT n);

// helper predicate to check a triple against a triple
bool sameTriple(INTTriple & t1, INTTriple & t2);

// helper predicate to check a triple against a triple in terms of relative size
inline bool lessThanTriple(INTTriple & t1, INTTriple & t2);

/*	Checks whether a move to a losing position
	is possible, if it is then this is NOT a
	losing position but otherwise it IS a
	losing position.
*/
bool isLosingPosition(INTTriple &);

// helper that reorders an INT triple based on value
void reorderTriple(INTTriple & t);

// subroutine to compute the desired F(n) without
// computing the losing positions of the nim game
INT computeFN1(INT n);

// computes part of the final solution
INT computeFN2();

// get first 9 digits
INT get9Digits(INT x);

// computes F(n) using a big integer class
BIGINT::BigInt computeFN3(long long n);

// prints some of the symbols used to compute the final
// solution
void printSymbols();

// myriad functions used to compute F(n)
BIGINT::BigInt beta(INT r);
BIGINT::BigInt theta(INT r);
BIGINT::BigInt alpha(INT n, INT r);
BIGINT::BigInt m(INT n, INT r);
BIGINT::BigInt size(INT r);
BIGINT::BigInt useful(INT n, INT r);
BIGINT::BigInt sumXr(INT r);
BIGINT::BigInt s1(INT n, INT r);
BIGINT::BigInt s2(INT n, INT r);
BIGINT::BigInt s3(INT n, INT r);
BIGINT::BigInt sum(INT n, INT r);
BIGINT::BigInt total(INT n, INT r);

// computes f(n)
BIGINT::BigInt computeFN4(long long n);

static std::ofstream out;

int main()
{

	using std::cout;
	using std::endl;
	using RBTREE::RBTree;
	//static const INT n = 100;

	/*	call helper method to generate
		all losing positions for which every value
		in the position is strictly less than n.
		Note that this method is more or less the heart
		of this program and the main algorithm that we
		are INTerested in here.
	*/
	//generateLosingPositions(n);

	out.open("thePattern.txt");

	
	// print losing positions & determine the sum
	/*
	static INT sum = 0;
	static INT numberOfEach[n] = { 0 };
	cout << "Losing Positions: \n";
	losingPositions[0].inOrderTraversal([](INTandTwoTrees & data)
	{
		data.second()->inOrderTraversal([](INTTriple & dat)
										{
											std::cout << "("
													<< dat.first()
													<< ", " << dat.second() 
													<< ", " << dat.third()
													<< ")" << std::endl;
											out << "("
												<< dat.first()
												<< ", " << dat.second()
												<< ", " << dat.third()
												<< ")" << std::endl;
											if (dat.first() != 0)
											{
												sum += dat.first();
												sum += dat.second();
												sum += dat.third();
											}
											//++numberOfEach[dat.first()];
										});
	});
	*/
	/*
	cout << "LPs that start with..." << endl;
	out << "LPs that start with..." << endl;
	for (int i = 0; i < n; ++i)
	{
		//cout << i << " : " << numberOfEach[i] << endl;
		//out << i << " : " << numberOfEach[i] << endl;
	}
	*/

	/*
	int q = 0;
	for (int i = 0; i < n; ++i)
	{
		q = std::max<int>(floor(n - pow(2, floor(log(i + 1) / log(2)) + 1) + 1) / 2, 0);
		cout << i << " : " << q << endl;
		out << i << " : " << q << endl;
	}
	*/

	/*
	cout << "The sum F(" << n << ") as requested is: " << sum << endl;
	out << "The sum F(" << n << ") as requested is: " << sum << endl;
	*/

	/*
	INT fn = computeFN1(n);
	cout << "The guess we have for F(" << n << ") is: " << fn << endl;
	out << "The guess we have for F(" << n << ") is: " << fn << endl;
	
	INT part = computeFN2();
	cout << "Part of the solution: " << part << endl;
	out << "Part of the solution: " << part << endl;

	printSymbols();
	*/

	/*
	INT test = computeFN1(n);
	cout << "Solution: " << test << endl;
	out << "Solution: " << test << endl;*/
	
	
	using BIGINT::BigInt;
	BigInt otherSum = computeFN4((long long)pow(10,18));
	cout << "Solution: " << otherSum << endl;
	out << "Solution: " << otherSum << endl;
	

	std::cout << "Press ENTER to exit.";
	std::cin.get();
	return 0;
}

// helper that handles placing a losing position into our
// custom data structure
void insertLosingPosition(const INTTriple & x)
{
	dummy.first() = x.first();
	INTandTwoTrees * a = losingPositions[0].search(dummy);
	if (a)
	{
		// insert into 2nd tree as well as 3rd tree
		a->second()->insert(x);
		a->third()->insert(x);
	}
	else
	{
		/*	create a node in the losing positions tree
			which also creates the trees that are attached
			to that node.
			Note that each of the trees have their own rules for organizing the INTTriple nodes!
		*/
		losingPositions[0].insert(INTandTwoTrees(x.first(),
												new TreeOfINTTriples(x,
																[](const INTTriple & x, const INTTriple & y){return x.second() < y.second(); },
																[](const INTTriple & x, const INTTriple & y){return x.second() == y.second(); }),
												new TreeOfINTTriples(x,
																[](const INTTriple & x, const INTTriple & y){return x.third() < y.third(); },
																[](const INTTriple & x, const INTTriple & y){return x.third() == y.third(); }))
								 );
	}

	dummy.first() = x.second();
	a = losingPositions[1].search(dummy);
	if (a)
	{
		// insert into 2nd tree as well as 3rd tree
		a->second()->insert(x);
		a->third()->insert(x);
	}
	else
	{
		/*	create a node in the losing positions tree
		which also creates the trees that are attached
		to that node.
		Note that each of the trees have their own rules for organizing the INTTriple nodes!
		*/
		losingPositions[1].insert(INTandTwoTrees(x.second(),
												new TreeOfINTTriples(x,
																	[](const INTTriple & x, const INTTriple & y){return x.first() < y.first(); },
																	[](const INTTriple & x, const INTTriple & y){return x.first() == y.first(); }),
												new TreeOfINTTriples(x,
																	[](const INTTriple & x, const INTTriple & y){return x.third() < y.third(); },
																	[](const INTTriple & x, const INTTriple & y){return x.third() == y.third(); }))
								);
	}

}

// helper to determine losing positions
void generateLosingPositions(INT n)
{
	using TRIPLE::Triple;
	/*	The data structure we are going to use is a construct
		built out of RBTrees. The idea is that we use 3 different
		trees - the first is sorted by the 1st element of the
		losing positions, the 2nd is sorted by the 2nd element of
		the losing positions, and the 3rd is sorted by the 3rd
		element of the losing positions.
		Since losing positions can share 1 element (although,
		by nature, cannot share 2 elements) we require a method
		to quickly search through these elements. Thus additional
		RBTrees are used for this purpose.
	*/
	for (INT k = 0; 2 * k + 2 < n; ++k)
		insertLosingPosition(INTTriple(0, 2 * k + 1, 2 * k + 2));


	/*	we test all losing positions starting
		with those that begin with 1.
	*/
	for (INT i = 1; i < n - 2; ++i) // O(n)
	{
		INT j = i + 1;
		INT k = 0;
		while (j < n - 1) // O(n)
		{
			if (!searchLosingPositions(i, j)) // O(logn)
			{
				k = j + 1;
				INTTriple t(i, j, k);

				while (k < n) // O(n) as well???
				{

					if (isLosingPosition(t)) // O(logn)
					{
						insertLosingPosition(t); // O(logn)
						break;
					}
					else
					{
						++k;
						t.third() = k;
					}
				}
			}
			++j;
		}
	} // I think this is something like O(n^3lognlogn)

}

// helper that reorders an INT triple based on value
void reorderTriple(INTTriple & t)
{
	// put smallest first
	INT min = t.first();
	if (t.second() < min)
	{
		min = t.second();
		t.second() = t.first();
		t.first() = min;
	}
	if (t.third() < min)
	{
		min = t.third();
		t.third() = t.first();
		t.first() = min;
	}
	if (t.second() > t.third())
	{
		min = t.third();
		t.third() = t.second();
		t.second() = min;
	}
}

/*	Checks whether a move to a losing position
	is possible, if it is then this is NOT a
	losing position but otherwise it IS a
	losing position.
*/
bool isLosingPosition(INTTriple & t)
{
	INTTriple temp(t);
	for (int perm = 0; perm < 3; ++perm)
	{
		// cycle (123), 3 times
		INT var = temp.first();
		temp.first() = temp.second();
		temp.second() = temp.third();
		temp.third() = var;

		for (int i = 0; i < 3; ++i)
		{
			dummy.first() = temp.first();
			INTandTwoTrees * x = losingPositions[i].search(dummy);
			if (x)
			{
				INTTriple * y = nullptr;
				y = x->second()->search(temp);
				if (y)
				{
					// in the case that we searched for a losing position
					// that has already been discovered, we return true
					if (sameTriple(*y, t))
						return true;

					// otherwise we are one move from a losing position
					// and return false
					return false;
				}
				else
				{
					// now search the next tree of triples
					y = x->third()->search(temp);
					if (y)
					{
						if (sameTriple(*y, t))
							return true;
						return false;
					}

				}
			}
		}
	}

	temp.first() = t.second();
	temp.second() = t.first();
	temp.third() = t.third();
	for (int perm = 0; perm < 3; ++perm)
	{
		// cycle (123), 3 times from a disjoint starting point
		INT var = temp.first();
		temp.first() = temp.second();
		temp.second() = temp.third();
		temp.third() = var;

		for (int i = 0; i < 3; ++i)
		{
			dummy.first() = temp.first();
			INTandTwoTrees * x = losingPositions[i].search(dummy);
			if (x)
			{
				INTTriple * y = nullptr;
				y = x->second()->search(temp);
				if (y)
				{
					// in the case that we searched for a losing position
					// that has already been discovered, we return true
					if (sameTriple(*y, t))
						return true;

					// otherwise we are one move from a losing position
					// and return false
					return false;
				}
				else
				{
					// now search the next tree of triples
					y = x->third()->search(temp);
					if (y)
					{
						if (sameTriple(*y, t))
							return true;
						return false;
					}

				}
			}
		}
	}

	// we haven't found a corresponding losing position
	// that we could make one move to
	// which implies that this IS a losing position
	return true;
}

	// loop through all possible moves on the 
	// first pile
	/*INTTriple theMove;
	for (INT i = t.first() - 1; i >= 0; --i)
	{

		theMove = INTTriple(i, t.second(), t.third());

		// in the case that an invalid move is made, skip ahead
		if (i == t.second() || i == t.third())
			continue;

		// order the move in terms of size
		reorderTriple(theMove);

		// check that the move isn't a winning one
		dummy.first() = theMove.first();
		INTandTwoTrees * x = losingPositions[0].search(dummy);
		if (x)
		{
			INTTriple * y = x->second()->search(theMove);
			if (y && (sameTriple(*y, theMove)))
					return false;
		}
		
		// if the move isn't a winning one then we
		// assume it's a losing one and continue
	}

	// loop through all possible moves on the 
	// second pile
	for (INT i = t.second() - 1; i >= 0; --i)
	{
		theMove = INTTriple(t.first(), i, t.third());

		// in the case that an invalid move is made, skip ahead
		if (i == t.first() || i == t.third())
			continue;

		// order the move in terms of size
		reorderTriple(theMove);

		// check that the move isn't a winning one
		dummy.first() = theMove.first();
		INTandTwoTrees * x = losingPositions[0].search(dummy);
		if (x)
		{
			INTTriple * y = x->second()->search(theMove);
			if (y && (sameTriple(*y, theMove)))
				return false;
		}

		// if the move isn't a winning one then we
		// assume it's a losing one and continue
	}


	// loop through all possible moves on the 
	// third pile
	for (INT i = t.third(); i >= 0; --i)
	{

		theMove = INTTriple(t.first(), t.second(), i);

		// in the case that an invalid move is made, skip ahead
		if (i == t.first() || i == t.second())
			continue;

		// order the move in terms of size
		reorderTriple(theMove);

		// check that the move isn't a winning one
		dummy.first() = theMove.first();
		INTandTwoTrees * x = losingPositions[0].search(dummy);
		if (x)
		{
			INTTriple * y = x->second()->search(theMove);
			if (y && (sameTriple(*y, theMove)))
				return false;
		}

		// if the move isn't a winning one then we
		// assume it's a losing one and continue
	}*/

	
	/*INTTriple * pos = nullptr;
	INTTriple temp;
	for (INT i = t.first() - 1; i >= 0; --i)
	{
		pos = searchLosingPositions(i, t.second());
		temp = INTTriple(i, t.second(), t.third());
		if (pos && lessThanTriple(*pos, temp))
		{
			// found a losing position that shares 2 values
			if (sameTriple(*pos, temp))
				// NOT a losing position, since it is actually a winning position!
				return false;
			// since we found a losing position, we can assume that
			// this move can be moved back to a losing position
			// if this holds for all moves, then we return true
		}
		// What happens if we haven't found a losing position???
		// We're going to assume that we'll find one when we search the
		// first and third values of the move but ignore the case in which
		// no losing position is found for either key searches.
		pos = searchLosingPositions(i, t.third());
		if (pos && lessThanTriple(*pos, temp))
		{
			if (sameTriple(*pos, INTTriple(i, t.second(), t.third())))
				return false;
		}
	}

	// check all moves made on the 2nd pile
	pos = searchLosingPositions(t.first(), t.third());
	if (pos)
	{
		if (sameTriple(*pos, t))
			return false;
	}
	for (INT i = t.second() - 1; i >= 0; --i)
	{
		pos = searchLosingPositions(t.first(), i);
		if (pos)
		{
			if (sameTriple(*pos, t))
				return false;
		}
	}


	pos = searchLosingPositions(t.first(), t.second());
	if (pos)
	{
		if (sameTriple(*pos, t))
			return false;
	}
	pos = nullptr;
	for (INT i = t.third() - 1; i >= 0; --i)
	{
		pos = searchLosingPositions(t.first(), i);
		if (pos)
		{
			if (sameTriple(*pos, t))
				return false;
		}
	}

	return true;
	*/

// helper method to efficiently search losing positions
// that have already been determined
INTTriple * searchLosingPositions(INT key1, INT key2)
{
	/*	We've already got all of the losing positions
		stored in the first RBTree. We need only search
		the first two keys as these keys represent the first
		and the second values of the position.
		*/
	dummy.first() = key1;
	INTandTwoTrees * x = losingPositions[0].search(dummy);
	if (x)
	{
		dummyINTTriple.first() = key1;
		dummyINTTriple.second() = key2;
		INTTriple * y = x->second()->search(dummyINTTriple);
		if (y)
			return y;
		else
			return x->third()->search(dummyINTTriple);
	}

	return nullptr;
}

// helper predicate to check a triple against a triple in terms of equality
inline bool sameTriple(INTTriple & t1, INTTriple & t2)
{
	return (t1.first() == t2.first()
		&& t1.second() == t2.second()
		&& t1.third() == t2.third());
}

// helper predicate to check a triple against a triple in terms of relative size
inline bool lessThanTriple(INTTriple & t1, INTTriple & t2)
{
	return (t1.first() <= t2.first()
		&& t1.second() <= t2.second()
		&& t1.third() <= t2.third());
}


// subroutine to compute the desired F(n) without
// computing the losing positions of the nim game
INT computeFN1(INT n)
{

	using std::cout;
	using std::endl;

	auto beta = [](INT r){return (INT) pow(2, r+1) - 1; };
	auto theta = [](INT r){return (INT) pow(2, r + 1); };
	auto alpha = [&beta, &theta, &n](INT r){return (INT)ceil((double)(n - beta(r)) / theta(r)); };
	auto m = [&alpha, &beta, &theta, &n](INT r)
								{ return (n >= beta(r) + (alpha(r) - 1)*theta(r) + theta(r) / 2) ?
										 n
										 :
										 beta(r) + (alpha(r) - 1)*theta(r) + theta(r) / 2; };
	auto size = [](INT r){return (INT) pow(2, r); };
	auto useful = [&alpha, &beta, &theta, &m](INT r) {return (INT)(beta(r) + alpha(r)*theta(r) - m(r)); };
	auto sumXr = [](INT r){ return (INT) (3 * pow(2, 2 * r - 1) - 3 * pow(2, r - 1)); };
	auto s1 = [&useful, &size, &m](INT r){ return (INT) size(r) * ((useful(r) - 1)*useful(r) / 2 + m(r)*useful(r)); };
	auto s2 = [&useful, &alpha, &beta, &theta](INT r){ return (INT) useful(r)*((beta(r)+(alpha(r) - 1)*theta(r))*theta(r)/2 + ((theta(r) / 2 - 1)*theta(r) / 4)); };
	auto s3 = [&useful, &sumXr](INT r){ return (INT) useful(r)*sumXr(r); };
	auto Sr = [&s1, &s2, &s3](INT r){ return (INT) s1(r) + s2(r) + s3(r); };
	auto Totalr = [&alpha, &beta, &theta, &size, &sumXr](INT r)
				  {return (INT)((beta(r)*alpha(r)*theta(r) + (alpha(r)*theta(r) - 1)*alpha(r)*theta(r) / 2)*size(r) + sumXr(r)*alpha(r)*theta(r) / 2); };
	INT Sn = 0, Tn = 0;
	INT R = (INT) floor(log(((INT)n) / 3) / log(2));

	// we aren't including the triples that start with zero!
	for (INT r = 1; r <= R; ++r)
	{
		Sn += Sr(r);
		Tn += Totalr(r);
		cout << "alpha(" << r << ") : " << alpha(r) << endl
			<< "beta(" << r << ") : " << beta(r) << endl
			<< "theta(" << r << ") : " << theta(r) << endl
			<< "m(" << r << ") : " << m(r) << endl
			<< "size(" << r << ") : " << size(r) << endl
			<< "sumXr(" << r << ") : " << sumXr(r) << endl
			<< "useful(" << r << ") : " << useful(r) << endl
			<< "s1(" << r << ") : " << s1(r) << endl
			<< "s2(" << r << ") : " << s2(r) << endl
			<< "s3(" << r << ") : " << s3(r) << endl
			<< "Sr(" << r << ") : " << Sr(r) << endl
			<< "Totalr(" << r << ") : " << Totalr(r) << endl
			<< "Sn from r = 1 to r = " << r << " : " << Sn << endl
			<< "Totaln from r = 1 to r = " << r << " : " << Tn << endl;
		out << "alpha(" << r << ") : " << alpha(r) << endl
			<< "beta(" << r << ") : " << beta(r) << endl
			<< "theta(" << r << ") : " << theta(r) << endl
			<< "m(" << r << ") : " << m(r) << endl
			<< "size(" << r << ") : " << size(r) << endl
			<< "sumXr(" << r << ") : " << sumXr(r) << endl
			<< "useful(" << r << ") : " << useful(r) << endl
			<< "s1(" << r << ") : " << s1(r) << endl
			<< "s2(" << r << ") : " << s2(r) << endl
			<< "s3(" << r << ") : " << s3(r) << endl
			<< "Sr(" << r << ") : " << Sr(r) << endl
			<< "Totalr(" << r << ") : " << Totalr(r) << endl
			<< "Sn from r = 1 to r = " << r << " : " << Sn << endl
			<< "Totaln from r = 1 to r = " << r << " : " << Tn << endl;
	}

	return (Tn - Sn);
}

// computes part of the final solution
INT computeFN2()
{

	using std::cout;
	using std::endl;

	long long sum = 0, temp = 0;
	cout << "First 17 sums of F(10^18): " << endl;
	out << "First 17 sums of F(10^18): " << endl;
	for (int r = 1; r < 18; ++r)
	{
		temp = - ( 7 * pow(2, 3 * r - 1) ) + ( 9  * pow(2, 2 * r - 1));
		sum += temp;
		out << r << ": " << temp << endl;
		out << "tally: " << sum << endl;
		cout << r << ": " << temp << endl;
		cout << "tally: " << sum << endl;
	}
	sum = get9Digits(pow(10, 18) + sum);

	cout << "Total: " << sum << endl;
	out << "Total: " << sum << endl;

	return sum;
}

BIGINT::BigInt beta(INT r)
{	
	return BIGINT::power(2, r + 1) - 1; 
}

BIGINT::BigInt theta(INT r)
{
	return BIGINT::power(2, r + 1);
}

BIGINT::BigInt alpha(INT n, INT r)
{
	using BIGINT::BigInt;

	BigInt x = n - BIGINT::power(2, r + 1) + 1, y = BIGINT::power(2, r + 1);

	if ((x % y) == 0)
		return (x / y);
	return (x / y) + 1;
}

BIGINT::BigInt m(INT n, INT r)
{
	return (n >= beta(r) + (alpha(n, r) - 1)*theta(r) + theta(r) / 2) ?
			n
			:
			beta(r) + (alpha(n, r) - 1)*theta(r) + theta(r) / 2;
}

BIGINT::BigInt size(INT r)
{
	return BIGINT::power(2, r);
}

BIGINT::BigInt useful(INT n, INT r)
{
	return beta(r) + alpha(n, r)*theta(r) - m(n, r);
}

BIGINT::BigInt sumXr(INT r)
{
	return 3 * BIGINT::power(2, 2 * r - 1) - 3 * BIGINT::power(2, r - 1);
}

BIGINT::BigInt s1(INT n, INT r)
{
	return size(r) * ((useful(n, r) - 1)*useful(n, r) / 2 + m(n, r)*useful(n, r));
}

BIGINT::BigInt s2(INT n, INT r)
{
	return useful(n, r)*((beta(r) + (alpha(n, r) - 1)*theta(r))*theta(r) / 2 + ((theta(r) / 2 - 1)*theta(r) / 4));
}

BIGINT::BigInt s3(INT n, INT r)
{
	return useful(n, r)*sumXr(r);
}

BIGINT::BigInt sum(INT n, INT r)
{
	return  s1(n, r) + s2(n, r) + s3(n, r);
}

BIGINT::BigInt total(INT n, INT r)
{
	return ((beta(r)*alpha(n, r)*theta(r) + (alpha(n, r)*theta(r) - 1)*alpha(n, r)*theta(r) / 2)*size(r) + sumXr(r)*alpha(n, r)*theta(r) / 2);
}

// computes F(n) using a big integer class
BIGINT::BigInt computeFN4(long long n)
{

	using std::cout;
	using std::endl;
	using BIGINT::BigInt;

	BigInt Sn = 0, Tn = 0;
	INT R = (INT)floor(log(((INT)n) / 3) / log(2));

	// we aren't including the triples that start with zero!
	for (INT r = 1; r <= R; ++r)
	{
		Sn = Sn + sum(n,r);
		Tn = Tn + total(n,r);
		cout << "alpha(" << r << ") : " << alpha(n, r) << endl
			<< "beta(" << r << ") : " << beta(r) << endl
			<< "theta(" << r << ") : " << theta(r) << endl
			<< "m(" << r << ") : " << m(n, r) << endl
			<< "size(" << r << ") : " << size(r) << endl
			<< "sumXr(" << r << ") : " << sumXr(r) << endl
			<< "useful(" << r << ") : " << useful(n, r) << endl
			<< "s1(" << r << ") : " << s1(n, r) << endl
			<< "s2(" << r << ") : " << s2(n, r) << endl
			<< "s3(" << r << ") : " << s3(n, r) << endl
			<< "Sr(" << r << ") : " << sum(n, r) << endl
			<< "Totalr(" << r << ") : " << total(n, r) << endl
			<< "Sn from r = 1 to r = " << r << " : " << Sn << endl
			<< "Totaln from r = 1 to r = " << r << " : " << Tn << endl;
		out << "alpha(" << r << ") : " << alpha(n, r) << endl
			<< "beta(" << r << ") : " << beta(r) << endl
			<< "theta(" << r << ") : " << theta(r) << endl
			<< "m(" << r << ") : " << m(n, r) << endl
			<< "size(" << r << ") : " << size(r) << endl
			<< "sumXr(" << r << ") : " << sumXr(r) << endl
			<< "useful(" << r << ") : " << useful(n, r) << endl
			<< "s1(" << r << ") : " << s1(n, r) << endl
			<< "s2(" << r << ") : " << s2(n, r) << endl
			<< "s3(" << r << ") : " << s3(n, r) << endl
			<< "Sr(" << r << ") : " << sum(n, r) << endl
			<< "Totalr(" << r << ") : " << total(n, r) << endl
			<< "Sn from r = 1 to r = " << r << " : " << Sn << endl
			<< "Totaln from r = 1 to r = " << r << " : " << Tn << endl;
	}

	return (Tn - Sn);
}

// computes F(n) using a big integer class
BIGINT::BigInt computeFN3(long long n)
{

	using std::cout;
	using std::endl;
	using BIGINT::BigInt;

	auto beta = [](INT r){return BigInt((INT)pow(2, r + 1) - 1); };
	auto theta = [](INT r){return BigInt((INT)pow(2, r + 1)); };
	auto alpha = [&beta, &theta, &n](INT r){return BigInt((INT)ceil((double)(n - pow(2, r+1) + 1) / pow(2, r+1))); };
	auto m = [&alpha, &beta, &theta, &n](INT r)
	{ return (n >= beta(r) + (alpha(r) - 1)*theta(r) + theta(r) / 2) ?
	n
	:
	beta(r) + (alpha(r) - 1)*theta(r) + theta(r) / 2; };
	auto size = [](INT r){return BigInt((INT)pow(2, r)); };
	auto useful = [&alpha, &beta, &theta, &m](INT r) {return (beta(r) + alpha(r)*theta(r) - m(r)); };
	auto sumXr = [](INT r){ return BigInt((INT)(3 * pow(2, 2 * r - 1) - 3 * pow(2, r - 1))); };
	auto s1 = [&useful, &size, &m](INT r){ return size(r) * ((useful(r) - 1)*useful(r) / 2 + m(r)*useful(r)); };
	auto s2 = [&useful, &alpha, &beta, &theta](INT r){ return useful(r)*((beta(r) + (alpha(r) - 1)*theta(r))*theta(r) / 2 + ((theta(r) / 2 - 1)*theta(r) / 4)); };
	auto s3 = [&useful, &sumXr](INT r){ return useful(r)*sumXr(r); };
	auto Sr = [&s1, &s2, &s3](INT r){ return s1(r) + s2(r) + s3(r); };
	auto Totalr = [&alpha, &beta, &theta, &size, &sumXr](INT r)
	{return ((beta(r)*alpha(r)*theta(r) + (alpha(r)*theta(r) - 1)*alpha(r)*theta(r) / 2)*size(r) + sumXr(r)*alpha(r)*theta(r) / 2); };
	BigInt Sn = 0, Tn = 0;
	INT R = (INT)floor(log(((INT)n) / 3) / log(2));

	// we aren't including the triples that start with zero!
	for (INT r = 1; r <= R; ++r)
	{
		Sn = Sn + Sr(r);
		Tn = Tn + Totalr(r);
		cout << "alpha(" << r << ") : " << alpha(r) << endl
			<< "beta(" << r << ") : " << beta(r) << endl
			<< "theta(" << r << ") : " << theta(r) << endl
			<< "m(" << r << ") : " << m(r) << endl
			<< "size(" << r << ") : " << size(r) << endl
			<< "sumXr(" << r << ") : " << sumXr(r) << endl
			<< "useful(" << r << ") : " << useful(r) << endl
			<< "s1(" << r << ") : " << s1(r) << endl
			<< "s2(" << r << ") : " << s2(r) << endl
			<< "s3(" << r << ") : " << s3(r) << endl
			<< "Sr(" << r << ") : " << Sr(r) << endl
			<< "Totalr(" << r << ") : " << Totalr(r) << endl
			<< "Sn from r = 1 to r = " << r << " : " << Sn << endl
			<< "Totaln from r = 1 to r = " << r << " : " << Tn << endl;
		out << "alpha(" << r << ") : " << alpha(r) << endl
			<< "beta(" << r << ") : " << beta(r) << endl
			<< "theta(" << r << ") : " << theta(r) << endl
			<< "m(" << r << ") : " << m(r) << endl
			<< "size(" << r << ") : " << size(r) << endl
			<< "sumXr(" << r << ") : " << sumXr(r) << endl
			<< "useful(" << r << ") : " << useful(r) << endl
			<< "s1(" << r << ") : " << s1(r) << endl
			<< "s2(" << r << ") : " << s2(r) << endl
			<< "s3(" << r << ") : " << s3(r) << endl
			<< "Sr(" << r << ") : " << Sr(r) << endl
			<< "Totalr(" << r << ") : " << Totalr(r) << endl
			<< "Sn from r = 1 to r = " << r << " : " << Sn << endl
			<< "Totaln from r = 1 to r = " << r << " : " << Tn << endl;
	}

	return (Tn - Sn);
}

// subroutine to print some of the symbols required to compute
// the final solution
void printSymbols()
{

	using std::cout;
	using std::endl;
	using std::string;
	using std::to_string;
	using std::stoll;

	long long alphas[7] = 
	{	((pow(5,18)-1)/2),
	((pow(5, 18) - 1) / 4),
	((pow(5, 18) - 1) / 8),
	((pow(5, 18) - 1) / 16 - 1 / 2),
	((pow(5, 18) - 1) / 32 - 3 / 4),
	((pow(5, 18) - 1) / 64 - 7 / 8),
	((pow(5, 18) - 1) / 128 - 15 / 16)
	};
	long long ms[7] = 
	{
		(pow(10, 18)),
		(pow(10,18)+pow(2,18)-1),
		(pow(10,18)+3*pow(2,18)-1),
		(pow(10,18)),
		(pow(10,18)),
		(pow(10,18)),
		(pow(10,18)),
	};

	auto beta = [](int r){return (string) to_string((long long)pow(2,r+1) - 1); };
	auto theta = [](int r){return (string) to_string((long long)pow(2, r + 1)); };
	auto alpha = [&alphas](int r){return (string) to_string(alphas[r - 18]); };
	auto m = [&ms](int r){ return (string)to_string(ms[r - 18]); };
	auto size = [](int r){return (string) to_string((long long)pow(2, r));  };
	auto useful = [&alpha, &beta, &theta, &m](int r) {return (string) to_string(stoll(beta(r))+stoll(alpha(r))*stoll(theta(r))-stoll(m(r))); };
	auto sumXr = [](int r){ return (string) to_string((long long)(3*(pow(2, 2 * r - 1) - pow(2, r - 1)))); };
	auto s1 = [&useful, &size, &m](int r){ return (string) "(" + size(r) + " * " + "((" + useful(r) + " - 1)*" + useful(r) + " / 2 + " + m(r) + "*" + useful(r) + "))"; };
	auto s2 = [&useful, &alpha, &beta, &theta](int r)
	{ return (string) "(" + useful(r) + "*((" + beta(r) + "+" + "(" + alpha(r) + " - 1)*" + theta(r) + ")*" + theta(r) + " / 2 + (( " + theta(r) + " / 2 - 1)*" + theta(r) + " / 4)) )"; };
	auto s3 = [&useful, &sumXr](int r){ return (string) "("+useful(r) + "*" + sumXr(r)+")"; };
	auto Sr = [&s1, &s2, &s3](INT r){ return (string) s1(r) + "+" + s2(r) + "+" + s3(r); };
	auto Totalr = [&alpha, &beta, &theta, &size, &sumXr](int r)
	{return (string)"((" + beta(r) + "*" + alpha(r) + "*" + theta(r) + "+ (" + alpha(r) + "*" + theta(r) + "- 1)*" + alpha(r) + "*" + theta(r) + " / 2)*" + size(r) + "+" + sumXr(r)
	+ "*" + alpha(r) + "*" + theta(r) + " / 2)"; };

	string Sn = "";
	string Tn = "";
	for (int r = 18; r <= 24; ++r)
	{
		Sn = Sn + Sr(r);
		Tn = Tn + Totalr(r);
		cout << "alpha(" << r << ") : " << alpha(r) << endl
			<< "beta(" << r << ") : " << beta(r) << endl
			<< "theta(" << r << ") : " << theta(r) << endl
			<< "m(" << r << ") : " << m(r) << endl
			<< "size(" << r << ") : " << size(r) << endl
			<< "sumXr(" << r << ") : " << sumXr(r) << endl
			<< "useful(" << r << ") : " << useful(r) << endl
			<< "s1(" << r << ") : " << s1(r) << endl
			<< "s2(" << r << ") : " << s2(r) << endl
			<< "s3(" << r << ") : " << s3(r) << endl
			<< "Sr(" << r << ") : " << Sr(r) << endl
			<< "Totalr(" << r << ") : " << Totalr(r) << endl
			<< "Sn from r = 18 to r = " << r << " : " << Sn << endl
			<< "Totaln from r = 18 to r = " << r << " : " << Tn << endl;
		out << "alpha(" << r << ") : " << alpha(r) << endl
			<< "beta(" << r << ") : " << beta(r) << endl
			<< "theta(" << r << ") : " << theta(r) << endl
			<< "m(" << r << ") : " << m(r) << endl
			<< "size(" << r << ") : " << size(r) << endl
			<< "sumXr(" << r << ") : " << sumXr(r) << endl
			<< "useful(" << r << ") : " << useful(r) << endl
			<< "s1(" << r << ") : " << s1(r) << endl
			<< "s2(" << r << ") : " << s2(r) << endl
			<< "s3(" << r << ") : " << s3(r) << endl
			<< "Sr(" << r << ") : " << Sr(r) << endl
			<< "Totalr(" << r << ") : " << Totalr(r) << endl
			<< "Sn from r = 18 to r = " << r << " : " << Sn << endl
			<< "Totaln from r = 18 to r = " << r << " : " << Tn << endl;
	}

}

// get first 9 digits
INT get9Digits(INT x)
{
	int digits[9] = { 0 };
	INT temp = x;
	for (int i = 0; i < 9; ++i)
	{
		digits[i] = temp % 10;
		temp /= 10;
	}
	
	INT ret = 0;
	for (int i = 0; i < 9; ++i)
		ret += (digits[i] * (int) pow(10, i));
	
	return ret;
}