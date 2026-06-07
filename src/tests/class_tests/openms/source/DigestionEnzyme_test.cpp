// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Alen Šarić   $
// $Authors: Alen Šarić$
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

//////////////////////////////////////////

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h> //Needed for Empty initialization
#include <sstream>

using namespace OpenMS;
using namespace std;

//////////////////////////////////////////

START_TEST(DigestionEnzyme,"$ID")

//////////////////////////////////////////

START_SECTION(bool setValueFromFile(const String& key, const String&  value))
  DigestionEnzymeProtein enzyme;

  // Test the Name Setting.
  TEST_TRUE(enzyme.setValueFromFile("test:Name","Trypsin"))
  TEST_EQUAL(enzyme.getName(),"Trypsin")

  // Test the RegEx Setting.
  TEST_TRUE(enzyme.setValueFromFile("test:RegEx","Reg"))
  TEST_EQUAL(enzyme.getRegEx(), "Reg")

  // Test the RegExDescription Setting.
  TEST_TRUE(enzyme.setValueFromFile("test:RegExDescription","Desc"))
  TEST_EQUAL(enzyme.getRegExDescription(),"Desc")

  // Test Synonym Setting
  TEST_TRUE(enzyme.setValueFromFile("syn:Synonyms:","Trypsin"))
  TEST_TRUE(enzyme.setValueFromFile("test:Synonyms:","TrypsinI"))

  // Since Synonyms are a set, test using set functions.
  TEST_EQUAL(enzyme.getSynonyms().count("Trypsin"),1)
  TEST_EQUAL(enzyme.getSynonyms().size(),2)

  // Test incorrect keys.
  TEST_FALSE(enzyme.setValueFromFile("test","Tryp-Like"))
END_SECTION

START_SECTION(bool operator==(const String& cleavage_regex) const)
  DigestionEnzymeProtein enzyme;
  enzyme.setRegEx("Verify");

  TEST_TRUE(enzyme == "Verify")
  TEST_FALSE(enzyme == "Accept")
END_SECTION

START_SECTION(bool operator!=(const String& cleavage_regex) const)
  DigestionEnzymeProtein enzyme;
  enzyme.setRegEx("Verify");

  TEST_TRUE(enzyme != "Accept")
  TEST_FALSE(enzyme != "Verify")
END_SECTION

// < compares the names of the enzymes.
START_SECTION(bool operator<(const DigestionEnzyme& enzyme) const)
  DigestionEnzymeProtein e1,e2;

  e1.setName("A_Enzyme");
  e2.setName("B_ENZYME");

  TEST_TRUE(e1 < e2)
  TEST_FALSE(e2 < e1)

  // Safety test: when names are same, neither greater nor smaller, whatever the regex.
  DigestionEnzymeProtein e3;
  e3.setName("A_Enzyme");
  e3.setRegEx("Greater");

  TEST_FALSE(e1 < e3)
  TEST_FALSE(e3 < e1)
END_SECTION

START_SECTION(std::ostream& operator<<(std::ostream& os, const DigestionEnzyme& enzyme))
  DigestionEnzymeProtein enzyme;
  enzyme.setName("TestEnzyme");
  enzyme.setRegEx("[K]");
  enzyme.setRegExDescription("cuts at K");

  stringstream ss;
  ss << enzyme;
  String output = ss.str();

  TEST_TRUE(output.hasSubstring("digestion enzyme:TestEnzyme"))
  TEST_TRUE(output.hasSubstring("(cleavage: [K] - cuts at K)"))
END_SECTION

END_TEST
