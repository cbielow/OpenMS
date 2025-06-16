// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FEATUREFINDER/IsotopeModel.h>
#include <sstream>


///////////////////////////

START_TEST(IsotopeModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

// default ctor
IsotopeModel* ptr = nullptr;
IsotopeModel* nullPointer = nullptr;
START_SECTION((IsotopeModel()))
ptr = new IsotopeModel();
TEST_EQUAL(ptr->getName(), "IsotopeModel")
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~IsotopeModel()))
delete ptr;
END_SECTION

// assignment operator
START_SECTION((virtual IsotopeModel & operator=(const IsotopeModel& source)))
IsotopeModel im1;

Param tmp;
tmp.setValue("charge", 3);
tmp.setValue("isotope:mode:GaussianSD", 0.8);
tmp.setValue("statistics:mean", 670.5);
im1.setParameters(tmp);

IsotopeModel im2;
im2 = im1;

IsotopeModel im3;
im3.setParameters(tmp);

im1 = IsotopeModel();
TEST_EQUAL(im3.getParameters(), im2.getParameters())
END_SECTION

// copy ctor
START_SECTION((IsotopeModel(const IsotopeModel& source)))
IsotopeModel im1;

Param tmp;
tmp.setValue("charge", 3);
tmp.setValue("isotope:mode:GaussianSD", 0.8);
tmp.setValue("statistics:mean", 670.5);
im1.setParameters(tmp);

IsotopeModel im2(im1);
IsotopeModel im3;
im3.setParameters(tmp);

im1 = IsotopeModel();
TEST_EQUAL(im3.getParameters(), im2.getParameters())
END_SECTION

START_SECTION([EXTRA] DefaultParamHandler::setParameters(...))
TOLERANCE_ABSOLUTE(0.001)
IsotopeModel im1;
Param tmp;
tmp.setValue("charge", 3);
tmp.setValue("isotope:mode:GaussianSD", 0.8);
tmp.setValue("statistics:mean", 670.5);
im1.setParameters(tmp);

IsotopeModel im2;
im2.setParameters(im1.getParameters());

std::vector<Peak1D> dpa1;
std::vector<Peak1D> dpa2;
im1.getSamples(dpa1);
im2.getSamples(dpa2);

TOLERANCE_ABSOLUTE(0.00001)
TEST_EQUAL(dpa1.size(), dpa2.size())
ABORT_IF(dpa1.size() != dpa2.size());
for (Size i = 0; i < dpa1.size(); ++i)
{
  TEST_REAL_SIMILAR(dpa1[i].getPosition()[0], dpa2[i].getPosition()[0])
  TEST_REAL_SIMILAR(dpa1[i].getIntensity(), dpa2[i].getIntensity())
}
END_SECTION

START_SECTION(UInt getCharge())
// can only reliably be tested after fitting, only sanity check here
IsotopeModel im1;
TEST_EQUAL(im1.getCharge() == 1, true) // default charge is 1
END_SECTION

START_SECTION(CoordinateType getCenter() const)
// can only reliably be tested after fitting, only sanity check here
IsotopeModel im1;
TEST_EQUAL(im1.getCenter() == 0, true)
END_SECTION

START_SECTION(void setSamples(const EmpiricalFormula& formula))
IsotopeModel im1;
Param tmp;
EmpiricalFormula ef("C66H129O3");
tmp.setValue("statistics:mean", ef.getAverageWeight() / 1);
tmp.setValue("interpolation_step", 0.01); // .003
tmp.setValue("charge", 1);
tmp.setValue("isotope:maximum", 100);
tmp.setValue("isotope:mode:mode", "Gaussian");
tmp.setValue("isotope:mode:GaussianSD", 0.15);
tmp.setValue("isotope:pattern_mode", "coarse");

im1.setParameters(tmp);
im1.setSamples(ef);
/*
{
  double data[] = {
    0.000429512, 0.00093697,  0.00196383,  0.00395466,  0.00765145,  0.0142235,    0.0254037,    0.043593,     0.0718726,   0.113852,    0.173278,
    0.253381,    0.355987,    0.480533,    0.623218,    0.776577,    0.929731,     1.06945,      1.18192,      1.25501,     1.28036,     1.25501,
    1.18192,     1.06945,     0.929731,    0.776577,    0.623218,    0.480533,     0.355987,     0.253381,     0.173278,    0.113852,    0.0718726,
    0.0439064,   0.0260875,   0.0156567,   0.0105376,   0.00953883,  0.0123444,    0.019477,     0.0318149,    0.0524539,   0.0830909,   0.126461,
    0.184922,    0.259806,    0.350701,    0.454835,    0.566759,    0.678534,     0.7805,       0.862586,     0.915925,    0.934428,    0.915925,
    0.862586,    0.7805,      0.678534,    0.566759,    0.454835,    0.350701,     0.259806,     0.184922,     0.126461,    0.0830909,   0.0524539,
    0.0318149,   0.0186555,   0.0106322,   0.00611169,  0.00394849,  0.00348857,   0.00450454,   0.00682396,   0.01171,     0.0193065,   0.0305829,
    0.0465459,   0.0680634,   0.0956255,   0.129081,    0.167409,    0.208604,     0.249745,     0.287275,     0.317488,    0.33712,     0.34393,
    0.33712,     0.317488,    0.287275,    0.249745,    0.208604,    0.167409,     0.129081,     0.0956255,    0.0680634,   0.0465459,   0.0305829,
    0.0193065,   0.0117385,   0.00688626,  0.00395131,  0.0023183,   0.00157108,   0.00147331,   0.00194089,   0.00289868,  0.00477912,  0.00757047,
    0.011522,    0.0168484,   0.0236711,   0.0319527,   0.0414404,   0.0516379,    0.0618218,    0.071112,     0.0785909,   0.0834507,   0.0851365,
    0.0834507,   0.0785909,   0.071112,    0.0618218,   0.0516379,   0.0414404,    0.0319527,    0.0236711,    0.0168484,   0.011522,    0.00757047,
    0.00477912,  0.00289868,  0.00169455,  0.000957445, 0.000533225, 0.000312194,  0.000225836,  0.000239371,  0.00031625,  0.000542686, 0.000894739,
    0.00141733,  0.00215713,  0.00315433,  0.00443167,  0.00598213,  0.0077584,    0.00966757,   0.0115742,    0.0133135,   0.0147137,   0.0156235,
    0.0159391,   0.0156235,   0.0147137,   0.0133135,   0.0115742,   0.00966757,   0.0077584,    0.00598213,   0.00443167,  0.00315433,  0.00215713,
    0.00141733,  0.000894739, 0.000542686, 0.00031625,  0.000177068, 9.52526e-005, 4.92314e-005, 2.44476e-005, 1.16643e-005};
  int size = sizeof(data) / sizeof(data[0]);
  std::vector<double> dpa2(data, &data[size]);

  std::vector<Peak1D> dpa1;
  im1.getSamples(dpa1);

  TEST_EQUAL(dpa1.size(), dpa2.size())
  ABORT_IF(dpa1.size() != dpa2.size());
  for (Size i = 0; i < dpa1.size(); ++i)
  {
    TEST_REAL_SIMILAR(dpa1[i].getIntensity(), dpa2[i])
  }
}

{
  // lorentzian
  tmp.setValue("isotope:mode:mode", "Lorentzian");
  tmp.setValue("isotope:mode:LorentzFWHM", 0.05);

  im1.setParameters(tmp);
  im1.setSamples(ef);

  double data[]
    = {0.0249619,  0.0291547,   0.0344977,   0.0414526,   0.0507371,   0.0635168,   0.0817848,   0.109176,    0.152888,   0.228787,   0.377365,
       0.725701,   1.80202,     5.53034,     3.91171,     1.28304,     0.570747,    0.315088,    0.198245,    0.1358,     0.0986954,  0.0749089,
       0.0587688,  0.0473237,   0.0389177,   0.0325645,   0.027647,    0,           0,           0,           0,          0,          0,
       0.0182176,  0.0212776,   0.0251771,   0.0302528,   0.0370288,   0.0463557,   0.059688,    0.0796788,   0.111581,   0.166973,   0.275407,
       0.52963,    1.31515,     4.03614,     2.85483,     0.936385,    0.416541,    0.229957,    0.144682,    0.0991093,  0.0720296,  0.0546698,
       0.0428905,  0.0345377,   0.0284028,   0.0237661,   0.0201772,   0,           0,           0,           0,          0,          0,
       0,          0.00670527,  0.00783155,  0.0092668,   0.011135,    0.013629,    0.0170619,   0.0219691,   0.029327,   0.041069,   0.061457,
       0.101368,   0.194938,    0.48406,     1.48556,     1.05076,     0.344651,    0.153314,    0.0846392,   0.0532526,  0.0364787,  0.0265116,
       0.0201221,  0.0157865,   0.0127121,   0.0104541,   0.00874748,  0.00742654,  0,           0,           0,          0,          0,
       0,          0.00165982,  0.00193862,  0.0022939,   0.00275636,  0.00337373,  0.0042235,   0.00543822,  0.0072596,  0.0101662,  0.015213,
       0.0250926,  0.048255,    0.119824,    0.367736,    0.260106,    0.0853148,   0.0379514,   0.0209516,   0.0131821,  0.00902993, 0.00656267,
       0.00498101, 0.00390778,  0.00314675,  0.00258781,  0.00216535,  0.00183836,  0,           0,           0,          0,          0,
       0,          0,           0.000310749, 0.000362945, 0.000429461, 0.000516041, 0.000631624, 0.000790718, 0.00101813, 0.00135913, 0.0019033,
       0.00284816, 0.00469779,  0.00903422,  0.0224333,   0.068847,    0.0486966,   0.0159725,   0.0071052,   0.00392252, 0.00246794, 0.00169057,
       0.00122865, 0.000932537, 0.00073161,  0.00058913,  0.000484485, 0.000405393, 0.000344176};
  int size = sizeof(data) / sizeof(data[0]);
  std::vector<double> dpa2(data, &data[size]);

  std::vector<Peak1D> dpa1;
  im1.getSamples(dpa1);

  TEST_EQUAL(dpa1.size(), dpa2.size())
  ABORT_IF(dpa1.size() != dpa2.size());
  for (Size i = 0; i < dpa1.size(); ++i)
  {
    TEST_REAL_SIMILAR(dpa1[i].getIntensity(), dpa2[i])
  }
}
*/
IsotopeModel im2;
Param tmp2;
tmp2.setValue("statistics:mean", ef.getAverageWeight() / 1);
tmp2.setValue("interpolation_step", 0.01);
tmp2.setValue("charge", 1);
tmp2.setValue("isotope:maximum", 100);
tmp2.setValue("isotope:mode:mode", "Gaussian");
tmp2.setValue("isotope:mode:GaussianSD", 0.15);
tmp2.setValue("isotope:pattern_mode", "fine");
// ef = EmpiricalFormula("C100H159N31O31S2");
//AASequence aa = AASequence::fromString("ACDEFGHIKLMNPQRSTVWY");
//ef = aa.getFormula(Residue::Full);
// ef = EmpiricalFormula("C100H202");

im2.setParameters(tmp2);
im1.setSamples(ef);
im2.setSamples(ef);

{
  std::vector<Peak1D> coarse_peaks;
  std::vector<Peak1D> fine_peaks;

  im1.getSamples(coarse_peaks);
  im2.getSamples(fine_peaks);

  int coarse_mono_idx = 0;
  int fine_mono_idx = 0;

  for (Size i = 0; coarse_peaks[i].getIntensity() < coarse_peaks[i+1].getIntensity(); ++i)
  {
    coarse_mono_idx = i;
  }

   for (Size i = 0; fine_peaks[i].getIntensity() < fine_peaks[i+1].getIntensity(); ++i)
  {
    fine_mono_idx = i;
  }

  OPENMS_LOG_INFO << "Monoisotopic coarse: POS: " << coarse_peaks[coarse_mono_idx].getMZ() << " INT: " << coarse_peaks[coarse_mono_idx].getIntensity() << std::endl;
  OPENMS_LOG_INFO << "Monoisotopic fine:   POS: " << fine_peaks[fine_mono_idx].getMZ() << " INT: " << fine_peaks[fine_mono_idx].getIntensity() << std::endl;

  // Log zur Kontrolle
  OPENMS_LOG_INFO << "First coarse: POS: " << coarse_peaks[0].getMZ() << " INT: " << coarse_peaks[0].getIntensity() << std::endl;
  OPENMS_LOG_INFO << "First fine:   POS: " << fine_peaks[0].getMZ() << " INT: " << fine_peaks[0].getIntensity() << std::endl;

  // Vergleich über Toleranz – Position und Intensität
  TEST_REAL_SIMILAR(coarse_peaks[0].getMZ(), fine_peaks[0].getMZ());

  OPENMS_LOG_INFO << "---- COARSE Peaks ----" << std::endl;
  for (Size i = 0; i < 20; ++i)
  {
    OPENMS_LOG_INFO << "COARSE " << i << ": MZ = " << coarse_peaks[i].getMZ() << "  INT = " << coarse_peaks[i].getIntensity() << std::endl;
  }

  OPENMS_LOG_INFO << "---- FINE Peaks ----" << std::endl;
  for (Size i = 0; i < 20; ++i)
  {
    OPENMS_LOG_INFO << "FINE " << i << ": MZ = " << fine_peaks[i].getMZ() << "  INT = " << fine_peaks[i].getIntensity() << std::endl;
  }

  double sum_coarse = 0.0;
  for (Size i = 0; i < coarse_peaks.size(); ++i)
  {
    sum_coarse += coarse_peaks[i].getIntensity();
  }

  double sum_fine = 0.0;
  for (Size i = 0; i < fine_peaks.size(); ++i)
  {
    sum_fine += fine_peaks[i].getIntensity();
  }

  OPENMS_LOG_INFO << "Sum coarse: " << sum_coarse << std::endl;
  OPENMS_LOG_INFO << "Sum fine:   " << sum_fine << std::endl;
}

// getcenter vergleichen
OPENMS_LOG_INFO << "Coarse center: " << im1.getCenter() << std::endl;
OPENMS_LOG_INFO << "Fine center:   " << im2.getCenter() << std::endl;


END_SECTION

START_SECTION(void setOffset(CoordinateType offset))
TOLERANCE_ABSOLUTE(0.1)
IsotopeModel im1;
Param tmp;
tmp.setValue("charge", 3);
tmp.setValue("isotope:mode:GaussianSD", 0.8);
tmp.setValue("statistics:mean", 670.5);
im1.setParameters(tmp);
im1.setOffset(673.5);

IsotopeModel im2;
im2.setParameters(im1.getParameters());
im2.setOffset(673.5);

std::vector<Peak1D> dpa1;
std::vector<Peak1D> dpa2;
im1.getSamples(dpa1);
im2.getSamples(dpa2);

TEST_EQUAL(dpa1.size(), dpa2.size())
ABORT_IF(dpa1.size() != dpa2.size());
for (Size i = 0; i < dpa1.size(); ++i)
{
  TEST_REAL_SIMILAR(dpa1[i].getPosition()[0], dpa2[i].getPosition()[0])
  TEST_REAL_SIMILAR(dpa1[i].getIntensity(), dpa2[i].getIntensity())
}
END_SECTION

START_SECTION(CoordinateType getOffset())
TOLERANCE_ABSOLUTE(0.1)
IsotopeModel im1;
Param tmp;
tmp.setValue("charge", 3);
tmp.setValue("isotope:mode:GaussianSD", 0.8);
tmp.setValue("statistics:mean", 670.5);
im1.setParameters(tmp);
im1.setOffset(673.5);

IsotopeModel im2;
im2.setParameters(im1.getParameters());
im2.setOffset(im1.getOffset());

std::vector<Peak1D> dpa1;
std::vector<Peak1D> dpa2;
im1.getSamples(dpa1);
im2.getSamples(dpa2);

TEST_EQUAL(dpa1.size(), dpa2.size())
ABORT_IF(dpa1.size() != dpa2.size());
for (Size i = 0; i < dpa1.size(); ++i)
{
  TEST_REAL_SIMILAR(dpa1[i].getPosition()[0], dpa2[i].getPosition()[0])
  TEST_REAL_SIMILAR(dpa1[i].getIntensity(), dpa2[i].getIntensity())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
