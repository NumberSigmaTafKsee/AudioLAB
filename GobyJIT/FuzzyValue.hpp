#pragma once

#include <vector>
#include <string>

class MembershipFunction
{
 public:
  virtual long double getMembership(long double) = 0;
  long double getTypicalValue();
  
 protected:
  long double typical_value;
};

class StandardMBF_Z : public MembershipFunction
{
 public:
  StandardMBF_Z(long double, long double);
  long double getMembership(long double);

 protected:
  long double min;
  long double max;
  long double slope;
  long double y_intercept;
};

class StandardMBF_S : public MembershipFunction
{
 public:
  StandardMBF_S(long double, long double);
  long double getMembership(long double);

 protected:
  long double min;
  long double max;
  long double slope;
  long double y_intercept;
};

class StandardMBF_Lambda : public MembershipFunction
{
 public:
  StandardMBF_Lambda(long double, long double, long double);
  long double getMembership(long double);

 protected:
  long double min;
  long double mid;
  long double max;
  long double slope_up;
  long double y_intercept_up;
  long double slope_down;
  long double y_intercept_down;
};

class StandardMBF_Pi : public MembershipFunction
{
 public:
  StandardMBF_Pi(long double, long double, long double, long double);
  long double getMembership(long double);

 protected:
  long double min;
  long double lower_mid;
  long double higher_mid;
  long double max;
  long double slope_up;
  long double y_intercept_up;
  long double slope_down;
  long double y_intercept_down;
};

class LinguisticSet
{
public:
  LinguisticSet(MembershipFunction*, std::string);
  MembershipFunction* getMembershipFunction();
  std::string getName();

protected:
  std::string name;
  MembershipFunction* membership_function;
};

class LinguisticDomain
{
 public:
  LinguisticDomain(std::string);
  void addLinguisticSet(LinguisticSet*);
  std::vector<long double> getMembership(long double);
  long double getMembership(long double, std::string);
  int getNumberOfLinguisticSets();
  std::vector<LinguisticSet*> getLinguisticSetList();

 protected:
  std::vector<LinguisticSet*> linguistic_set_list;
  std::string name;
};

class FuzzyValue
{
 public:
  FuzzyValue(LinguisticDomain*);
  void setCrispValue(long double);
  long double getCrispValue();
  void setSetMembership(long double, std::string);
  long double getSetMembership(std::string);
  void clear();
  long double AND(std::string, FuzzyValue*, std::string);
  void OR_setSetMembership(long double, std::string);

 protected:
  LinguisticDomain* linguistic_domain;
  std::vector<long double> values;
};

long double MembershipFunction::getTypicalValue()
{
  return typical_value;
}

StandardMBF_Z::StandardMBF_Z(long double mn, long double mx)
{
  min = mn;
  max = mx;
  typical_value = min;
  slope = -1 * ( (1-0) / (max - min) );
  y_intercept = -1 * slope * min + 1;
}

long double StandardMBF_Z::getMembership(long double x)
{
  if (x <= min)
    {
      return 1;
    }
  if ( (x >= min) && (x <= max) )
    {
      return slope * x + y_intercept;
    }
  if (x >= max)
    {
      return 0;
    }
    return -1;
}
 
StandardMBF_S::StandardMBF_S(long double mn, long double mx)
{
  min = mn;
  max = mx;
  typical_value = max;
  slope = (1-0) / (max - min);
  y_intercept = slope * min * -1;
}

long double StandardMBF_S::getMembership(long double x)
{
  if (x <= min)
    {
      return 0;
    }
  if ( (x >= min) && (x <= max) )
    {
      return slope * x + y_intercept;
    }
  if (x >= max)
    {
      return 1;
    }
  return -1;
}

StandardMBF_Lambda::StandardMBF_Lambda(long double mn, long double md, long double mx)
{
  min = mn;
  mid = md;
  max = mx;
  typical_value = mid;
  slope_up = (1-0) / (mid - min);
  y_intercept_up = slope_up * min * -1;
  slope_down = -1 * ( (1-0) / (max - mid) );
  y_intercept_down = -1 * slope_down * mid + 1;
}

long double StandardMBF_Lambda::getMembership(long double x)
{
  if (x <= min)
    {
      return 0;
    }
  if ( (x >= min) && (x <= mid) )
    {
      return slope_up * x + y_intercept_up;
    }
  if ( (x >= mid) && (x <= max) )
    {
      return slope_down * x + y_intercept_down;
    }
  if (x >= max)
    {
      return 0;
    }
    return -1;
}

StandardMBF_Pi::StandardMBF_Pi(long double mn, long double lmd, long double hmd, long double mx)
{
  min = mn;
  lower_mid = lmd;
  higher_mid = hmd;
  max = mx;
  typical_value = (lower_mid + higher_mid) / 2;
  slope_up = (1-0) / (lower_mid - min);
  y_intercept_up = slope_up * min * -1;
  slope_down = -1 * ( (1-0) / (max - higher_mid) );
  y_intercept_down = -1 * slope_down * higher_mid + 1;
}

long double StandardMBF_Pi::getMembership(long double x)
{
  if (x <= min)
    {
      return 0;
    }
  if ( (x >= min) && (x <= lower_mid) )
    {
      return slope_up * x + y_intercept_up;
    }
  if ( (x >= lower_mid) && (x <= higher_mid) )
    {
      return 1;
    }

  if ( (x >= higher_mid) && (x <= max) )
    {
      return slope_down * x + y_intercept_down;
    }
  if (x >= max)
    {
      return 0;
    }
    return -1;
}

LinguisticSet::LinguisticSet(MembershipFunction* mf, std::string nm)
{
  membership_function = mf;
  name = nm;
}

MembershipFunction* LinguisticSet::getMembershipFunction()
{
  return membership_function;
}

std::string LinguisticSet::getName()
{
  return name;
}

LinguisticDomain::LinguisticDomain(std::string nm)
{
  name = nm;
}

void LinguisticDomain::addLinguisticSet(LinguisticSet* l)
{
  linguistic_set_list.push_back(l);
}

std::vector<long double> LinguisticDomain::getMembership(long double x)
{
  std::vector<long double> membership_vector;
  for (int i = 0; i < linguistic_set_list.size(); i++)
    {
     membership_vector.push_back(linguistic_set_list[i]->getMembershipFunction()->getMembership(x));
    }
  return membership_vector;
}

long double LinguisticDomain::getMembership(long double x, std::string s)
{
  for (int i = 0; i < linguistic_set_list.size(); i++)
    {
      if (linguistic_set_list[i]->getName() == s)
	{
	  return linguistic_set_list[i]->getMembershipFunction()->getMembership(x);
	}
    }
  return 0.;
}

int LinguisticDomain::getNumberOfLinguisticSets()
{
  return linguistic_set_list.size();
}

std::vector<LinguisticSet*> LinguisticDomain::getLinguisticSetList()
{
  return linguistic_set_list;
}

FuzzyValue::FuzzyValue(LinguisticDomain* ld)
{
  linguistic_domain = ld;
  long double dummy = 0;
  int i;
  for (i=0; i<ld->getNumberOfLinguisticSets(); i++)
    {
      values.push_back(dummy);
    }
}

void FuzzyValue::setCrispValue(long double x)
{
  int i;
  std::vector<long double> membership_vector = linguistic_domain->getMembership(x);
  for (i=0; i<linguistic_domain->getNumberOfLinguisticSets(); i++)
    {
      values[i] = membership_vector[i];
    }
}

long double FuzzyValue::getCrispValue()
{

  // center-of-maximum method

  int i;
  std::vector<long double> typical_value_vector;
  for (i=0; i<linguistic_domain->getNumberOfLinguisticSets(); i++)
    {
      typical_value_vector.push_back( linguistic_domain->getLinguisticSetList()[i]->getMembershipFunction()->getTypicalValue());
    }

  long double answer = 0;
  for (i=0; i<linguistic_domain->getNumberOfLinguisticSets(); i++)
    {
      answer = answer + typical_value_vector[i] * values[i];
    }
  return answer;
}

void FuzzyValue::clear()
{
  int i;
  for (i=0; i<linguistic_domain->getNumberOfLinguisticSets(); i++)
    {
      values[i] = 0;
    }
}

void FuzzyValue::setSetMembership(long double x, std::string nm)
{
  int i;
  for (i=0; i<linguistic_domain->getNumberOfLinguisticSets(); i++)
    {
      if (linguistic_domain->getLinguisticSetList()[i]->getName() == nm)
	{
	  values[i] = x;
	}
    }
}

long double FuzzyValue::getSetMembership(std::string nm)
{
  int i;
  for (i=0; i<linguistic_domain->getNumberOfLinguisticSets(); i++)
    {
      if (linguistic_domain->getLinguisticSetList()[i]->getName() == nm)
	{
	  return values[i];
	}
    }
    return -1;
}

long double FuzzyValue::AND(std::string sa, FuzzyValue* fv, std::string sb)
{
  long double a = getSetMembership(sa);
  long double b = fv->getSetMembership(sb);
  if (a <= b)
    {
      return a;
    }
else
  {
    return b;
  }
}

void FuzzyValue::OR_setSetMembership(long double x, std::string nm)
{
  int i;
  for (i=0; i<linguistic_domain->getNumberOfLinguisticSets(); i++)
    {
      if (linguistic_domain->getLinguisticSetList()[i]->getName() == nm)
	{
	  if (x > values[i])
	    {
	      values[i] = x;
	    }
	}
    }
}