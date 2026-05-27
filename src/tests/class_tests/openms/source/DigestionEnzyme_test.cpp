// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Alen Saric   $
// $Authors: Alen Saric$
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

//////////////////////////////////////////

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <sstream>

using namespace OpenMS;
using namespace std;

//////////////////////////////////////////

START_TEST(DigestionEnzyme,"$ID")

//////////////////////////////////////////

DigestionEnzyme* e_ptr = nullptr;
DigestionEnzyme* e_null = nullptr;


START_SECTION((DigestionEnzyme()))
  e_ptr = new DigestionEnzyme();
  TEST_NOT_EQUAL(e_ptr,e_null)
END_SECTION

START_SECTION((virtual ~DigestionEnzyme()))
  delete e_ptr;
END_SECTION

START_SECTION(bool setValueFromFile(const String& key, const String&  value))
  DigestionEnzyme enzyme;

  // Test the Name Setting.
  TEST_EQUAL(enzyme.setValueFromFile("test:Name","Trypsin"),true)
  TEST_EQUAL(enzyme.getName(),"Trypsin")

  // Test the RegEx Setting.
  TEST_EQUAL(enzyme.setValueFromFile("test:RegEx","Reg"),true)
  TEST_EQUAL(enzyme.getRegEx(), "Reg")

  // Test the RegExDescription Setting.
  TEST_EQUAL(enzyme.setValueFromFile("test:RegExDescription","Desc"),true)
  TEST_EQUAL(enzyme.getRegExDescription(),"Desc")

  // Test Synonym Setting
  TEST_EQUAL(enzyme.setValueFromFile("syn:Synonyms:","Trypsin"),true)
  TEST_EQUAL(enzyme.setValueFromFile("test:Synonyms:","TrypsinI"),true)

  // Since Synonyms are a set, test using set functions.
  TEST_EQUAL(enzyme.getSynonyms().count("Trypsin"),1)
  TEST_EQUAL(enzyme.getSynonyms().size(),2)

  // Test incorrect keys.
  TEST_NOT_EQUAL(enzyme.setValueFromFile("test","Tryp-Like"),true)
END_SECTION

START_SECTION(bool operator==(const String& cleavage_regex) const)
  DigestionEnzyme enzyme;
  enzyme.setRegEx("Verify");

  TEST_EQUAL(enzyme == "Verify",true)
  TEST_NOT_EQUAL(enzyme == "Accept",true)
END_SECTION

START_SECTION(bool operator!=(const String& cleavage_regex) const)
  DigestionEnzyme enzyme;
  enzyme.setRegEx("Verify");

  TEST_EQUAL(enzyme != "Accept",true)
  TEST_NOT_EQUAL(enzyme != "Verify",true)
END_SECTION

// < compares the names of the enzymes.
START_SECTION(bool operator<(const DigestionEnzyme& enzyme) const)
  DigestionEnzyme e1,e2;

  e1.setName("A_Enzyme");
  e2.setName("B_ENZYME");

  TEST_EQUAL(e1 < e2, true)
  TEST_EQUAL(e2 < e1, false)

  // Safety test: when names are same, neither greater nor smaller, whatever the regex.
  DigestionEnzyme e3;
  e3.setName("A_Enzyme");
  e3.setRegEx("Greater");

  TEST_NOT_EQUAL(e1 < e3, true)
  TEST_NOT_EQUAL(e3 < e1, true)
END_SECTION

START_SECTION(std::ostream& operator<<(std::ostream& os, const DigestionEnzyme& enzyme))
  DigestionEnzyme enzyme;
  enzyme.setName("TestEnzyme");
  enzyme.setRegEx("[K]");
  enzyme.setRegExDescription("cuts at K");

  stringstream ss;
  ss << enzyme;
  String output = ss.str();

  TEST_EQUAL(output.hasSubstring("digestion enzyme:TestEnzyme"),true)
  TEST_EQUAL(output.hasSubstring("(cleavage: [K] - cuts at K)"),true)
END_SECTION

END_TEST
