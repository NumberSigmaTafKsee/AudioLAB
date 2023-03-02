#include "FuzzyValue.hpp"
#include <iostream>

int main()
{
  std::cout << std::endl;
  std::cout << "This program demonstrates Dan Williams' yet untitled" << std::endl;
  std::cout << "fuzzy logic toolkit." << std::endl;
  std::cout << std::endl;
  std::cout << "It implements the crane example in Constantin von Altrock's book " << std::endl;
  std::cout << "'Fuzzy Logic & Neurofuzzy Applications Explained'" << std::endl;

  std::cout << std::endl;
  std::cout << "Distance sensor reads 12 yards." << std::endl;
  std::cout << "Angle sensor reads 4 degrees." << std::endl;

  /*
  create and name a linguistic domain for distance
  */
  LinguisticDomain* distance_domain = new LinguisticDomain("distance_domain");

  /*
  define some linguistic values and membership functions for the distance domain
  */

  // too_far
  StandardMBF_Z* distance_domain_too_far_membership = new StandardMBF_Z(-5, 0);
  LinguisticSet* distance_domain_too_far = new LinguisticSet(distance_domain_too_far_membership, "too_far");

  // zero
  StandardMBF_Lambda* distance_domain_zero_membership = new StandardMBF_Lambda(-5, 0, 5);
  LinguisticSet* distance_domain_zero = new LinguisticSet(distance_domain_zero_membership, "zero");

  // close
  StandardMBF_Lambda* distance_domain_close_membership = new StandardMBF_Lambda(0, 5, 10);
  LinguisticSet* distance_domain_close = new LinguisticSet(distance_domain_close_membership, "close");

  // medium
  StandardMBF_Lambda* distance_domain_medium_membership = new StandardMBF_Lambda(5, 10, 30);
  LinguisticSet* distance_domain_medium = new LinguisticSet(distance_domain_medium_membership, "medium");

  // far
  StandardMBF_S* distance_domain_far_membership = new StandardMBF_S(10, 30);
  LinguisticSet* distance_domain_far = new LinguisticSet(distance_domain_far_membership, "far");

  /*
  Add the linguistic values to the distance domain
  */
  distance_domain->addLinguisticSet(distance_domain_too_far);
  distance_domain->addLinguisticSet(distance_domain_zero);
  distance_domain->addLinguisticSet(distance_domain_close);
  distance_domain->addLinguisticSet(distance_domain_medium);
  distance_domain->addLinguisticSet(distance_domain_far);

  /*
  create and name a linguistic domain for angle
  */
  LinguisticDomain* angle_domain = new LinguisticDomain("angle_domain");

  /*
  define some linguistic values and membership functions for the angle domain
  */

  // neg_big
  StandardMBF_Z* angle_domain_neg_big_membership = new StandardMBF_Z(-45, -5);
  LinguisticSet* angle_domain_neg_big = new LinguisticSet(angle_domain_neg_big_membership, "neg_big");

  // neg_small
  StandardMBF_Lambda* angle_domain_neg_small_membership = new StandardMBF_Lambda(-45, -5, 0);
  LinguisticSet* angle_domain_neg_small = new LinguisticSet(angle_domain_neg_small_membership, "neg_small");

  // zero
  StandardMBF_Lambda* angle_domain_zero_membership = new StandardMBF_Lambda(-5, 0, 5);
  LinguisticSet* angle_domain_zero = new LinguisticSet(angle_domain_zero_membership, "zero");

  // pos_small
  StandardMBF_Lambda* angle_domain_pos_small_membership = new StandardMBF_Lambda(0, 5, 45);
  LinguisticSet* angle_domain_pos_small = new LinguisticSet(angle_domain_pos_small_membership, "pos_small");

  // pos_big
  StandardMBF_S* angle_domain_pos_big_membership = new StandardMBF_S(5, 45);
  LinguisticSet* angle_domain_pos_big = new LinguisticSet(angle_domain_pos_big_membership, "pos_big");

  /*
    Add the linguistic values to the angle domain
  */
  angle_domain->addLinguisticSet(angle_domain_neg_big);
  angle_domain->addLinguisticSet(angle_domain_neg_small);
  angle_domain->addLinguisticSet(angle_domain_zero);
  angle_domain->addLinguisticSet(angle_domain_pos_small);
  angle_domain->addLinguisticSet(angle_domain_pos_big);

  /*
  create and name a linguistic domain for power
  */
  LinguisticDomain* power_domain = new LinguisticDomain("power_domain");

  /*
  define some linguistic values and membership functions for the power domain
  */

  // neg_high
  StandardMBF_Lambda* power_domain_neg_high_membership = new StandardMBF_Lambda(-30, -25, -8);
  LinguisticSet* power_domain_neg_high = new LinguisticSet(power_domain_neg_high_membership, "neg_high");

  // neg_medium
  StandardMBF_Lambda* power_domain_neg_medium_membership = new StandardMBF_Lambda(-25, -8, 0);
  LinguisticSet* power_domain_neg_medium = new LinguisticSet(power_domain_neg_medium_membership, "neg_medium");

  // zero
  StandardMBF_Lambda* power_domain_zero_membership = new StandardMBF_Lambda(-8, 0, 8);
  LinguisticSet* power_domain_zero = new LinguisticSet(power_domain_zero_membership, "zero");

  // pos_medium
  StandardMBF_Lambda* power_domain_pos_medium_membership = new StandardMBF_Lambda(0, 8, 25);
  LinguisticSet* power_domain_pos_medium = new LinguisticSet(power_domain_pos_medium_membership, "pos_medium");

  // pos_high
  StandardMBF_Lambda* power_domain_pos_high_membership = new StandardMBF_Lambda(8, 25, 20);
  LinguisticSet* power_domain_pos_high = new LinguisticSet(power_domain_pos_high_membership, "pos");

  /*
  add the linguistic values to the power domain
  */
  power_domain->addLinguisticSet(power_domain_neg_high);
  power_domain->addLinguisticSet(power_domain_neg_medium);
  power_domain->addLinguisticSet(power_domain_zero);
  power_domain->addLinguisticSet(power_domain_pos_medium);
  power_domain->addLinguisticSet(power_domain_pos_high);

  /* 
  "Fuzzify" sensor readings
  */
  FuzzyValue* distance = new FuzzyValue(distance_domain);
  distance->setCrispValue(12);
  FuzzyValue* angle = new FuzzyValue(angle_domain);
  angle->setCrispValue(4);

  /*
  Create a fuzzy variable to store power inference calculations
  */
  FuzzyValue* power = new FuzzyValue(power_domain);

  /*
  Fuzzy inference of power value
  */
  power->OR_setSetMembership( distance->AND("medium", angle, "pos_small"), "pos_medium" );
  power->OR_setSetMembership( distance->AND("medium", angle, "zero"), "zero" );
  power->OR_setSetMembership( distance->AND("far", angle, "zero"), "pos_medium" );

  /* 
  "Defuzzify" infered power value
  */
  long double power_setting;
  power_setting = power->getCrispValue();
  std::cout << "Set power to " << power_setting << " kW." << std::endl;

  std::cout << std::endl;
  return 0;
}
