
#include "common.h"

class TermPredicate;
class RealPolynomialTermPredicate;
class RealPolynomial;
typedef std::pair<int, Int2IntPtr> RealPolynomialTerm;
typedef std::map<RealPolynomialTerm, Number, RealPolynomialTermPredicate> RealPolynomialTerms;
typedef std::shared_ptr<RealPolynomialTerms> RealPolynomialTermsPtr;

std::ostream & operator<<(std::ostream & stream, Int2Int const & rhs);
std::ostream & operator<<(std::ostream & stream, RealPolynomialTerms const & rhs);
std::ostream & operator<<(std::ostream & stream, RealPolynomial const & rhs);

RealPolynomial operator+(RealPolynomial const &);
RealPolynomial operator-(RealPolynomial const &);

RealPolynomial operator+(RealPolynomial const &, RealPolynomial const &);
RealPolynomial operator-(RealPolynomial const &, RealPolynomial const &);
RealPolynomial operator*(RealPolynomial const &, RealPolynomial const &);
RealPolynomial operator/(RealPolynomial const &, RealPolynomial const &);


class RealPolynomialTermPredicate {
public:
	bool operator()(RealPolynomialTerm const & p1, RealPolynomialTerm const & p2) const{
		if (p1.first == p2.first)
			return *p1.second < *p2.second;
		else
			return p1.first < p2.first;
	}
};

class RealPolynomial{
public:
	enum Operation {
		SUM,DIFF,PROD
	};
public:
	friend class Problem;
	friend class FunctionComplex;
	friend RealPolynomial operator+(RealPolynomial const &);
	friend RealPolynomial operator-(RealPolynomial const &);

	friend RealPolynomial operator+(RealPolynomial const &, RealPolynomial const &);
	friend RealPolynomial operator-(RealPolynomial const &, RealPolynomial const &);
	friend RealPolynomial operator*(RealPolynomial const &, RealPolynomial const &);
	friend RealPolynomial operator/(RealPolynomial const &, RealPolynomial const &);
public:
	RealPolynomial();
	RealPolynomial(Number const &);
	RealPolynomial(NumberPtr const &);
	~RealPolynomial();
public:
	RealPolynomial clone()const;
public:
	void clear();

	std::ostream & print(std::ostream & stream, Problem const &)const;
	std::ostream & print(std::ostream & stream)const;
public:
	RealPolynomialTerms & terms();
	RealPolynomialTerms const & terms()const;

	int & degree();
	int const & degree()const;

	void set(int, Number);
	void set(int);
	void set(Number);

	Number constant()const;
private:
	void add(RealPolynomial const & key, Number value, Number factor);
	void add(RealPolynomial const & key, Number factor);
	void add(RealPolynomialTerm const & key, Number value);
	void add(typename RealPolynomialTerms::value_type const & kvp, Number factor);
private:
	RealPolynomialTermsPtr _terms;
	IntPtr _degree;
};