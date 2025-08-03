// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/RawMSSignalSimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(RawMSSignalSimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RawMSSignalSimulation* ptr = nullptr;
RawMSSignalSimulation* nullPointer = nullptr;
SimTypes::MutableSimRandomNumberGeneratorPtr empty_rnd_gen(new SimTypes::SimRandomNumberGenerator);
// const unsigned long rnd_gen_seed = 1;

START_SECTION((RawMSSignalSimulation(SimRandomNumberGeneratorPtr rng)))
{
  ptr = new RawMSSignalSimulation(empty_rnd_gen);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~RawMSSignalSimulation())
{
  delete ptr;
}
END_SECTION

START_SECTION((RawMSSignalSimulation(const RawMSSignalSimulation& source)))
{
  RawMSSignalSimulation source(empty_rnd_gen);
  Param p = source.getParameters();
  p.setValue("peak_fwhm", 0.3);
  source.setParameters(p);

  RawMSSignalSimulation target(source);
  TEST_EQUAL(source.getParameters(), target.getParameters())
}
END_SECTION

START_SECTION((RawMSSignalSimulation & operator=(const RawMSSignalSimulation& source)))
{
  RawMSSignalSimulation source(empty_rnd_gen);
  RawMSSignalSimulation target(source);

  Param p = source.getParameters();
  p.setValue("peak_fwhm", 0.3);
  source.setParameters(p);
  TEST_NOT_EQUAL(source.getParameters(), target.getParameters())

  target = source;

  TEST_EQUAL(source.getParameters(), target.getParameters())
}
END_SECTION

START_SECTION((void generateRawSignals(SimTypes::FeatureMapSim& features,
                                       SimTypes::MSSimExperiment& experiment,
                                       SimTypes::MSSimExperiment& experiment_ct,
                                       SimTypes::FeatureMapSim& contaminants)))
{
  // TODO
}
END_SECTION


START_SECTION((void loadContaminants()))
{
  // TODO
}
END_SECTION

START_SECTION((void compressSignalsIonMobility_(SimTypes::MSSimExperiment& experiment)))
{
  using namespace OpenMS;
  using SimTypes::MSSimExperiment;

  MSSimExperiment exp;
  MSSpectrum spectrum;

  std::vector<double> ims = {1.06, 1.061, 1.062, 1.12, 1.12, 1.086, 1.063, 1.09, 1.04, 1.02, 1.06, 1.06};
  std::vector<double> mzs = {400, 400, 402, 123, 500, 495, 401, 12, 403, 404, 400, 400};
  std::vector<double> ints = {100, 120, 110, 80, 90, 70, 130, 60, 85, 95, 110, 69};

  // add peaks
  for (Size i = 0; i < ims.size(); ++i)
  {
    Peak1D p;
    p.setMZ(mzs[i]);
    p.setIntensity(ints[i]);
    spectrum.push_back(p);
  }

  // add ionmobility data
  MSSpectrum::FloatDataArrays& fda = spectrum.getFloatDataArrays();
  fda.resize(1);
  IMDataConverter::setIMUnit(fda[0], DriftTimeUnit::VSSC);
  fda[0].setMetaValue("cv accession", "MS:1003008");
  fda[0].setMetaValue("unit_accession", "MS:1002814");
  fda[0].setMetaValue("unit_name", "volt-second per square centimeter");
  fda[0].setMetaValue("unit_cv_ref", "MS");
  fda[0].insert(fda[0].begin(), ims.begin(), ims.end());


  InstrumentSettings is;
  is.getScanWindows().resize(1);
  is.getScanWindows()[0].begin = 10.0;
  is.getScanWindows()[0].end = 510.0;
  spectrum.setInstrumentSettings(is);

  exp.addSpectrum(spectrum);


  std::cout << "-------- BEFORE COMPRESSION --------" << std::endl;
  const auto& im_data_before = exp[0].getFloatDataArrays()[0];
  for (Size i = 0; i < exp[0].size(); ++i)
  {
    std::cout << "m/z: " << exp[0][i].getMZ() << " | Intensity: " << exp[0][i].getIntensity() << " | IM: " << im_data_before[i] << std::endl;
  }


  exp[0].MSSpectrum::initializeIMFloatDataArray("vssc");
  // Compression aufrufen
  RawMSSignalSimulation sim;
  sim.setIMGridWidth_(0.01);
  sim.compressSignalsIonMobility_(exp);


  std::cout << "-------- AFTER COMPRESSION --------" << std::endl;
  const auto& im_data_after = exp[0].getFloatDataArrays()[0];
  for (Size i = 0; i < exp[0].size(); ++i)
  {
    std::cout << "m/z: " << exp[0][i].getMZ() << " | Intensity: " << exp[0][i].getIntensity() << " | IM: " << im_data_after[i] << std::endl;
  }


  TEST_EQUAL(exp[0][2].getIntensity(), 399)

  const auto& im_array = exp[0].getFloatDataArrays()[0];

  for (Size i = 1; i < im_array.size(); ++i)
  {
    TEST_TRUE(im_array[i - 1] <= im_array[i])
  }

  TEST_EQUAL(im_array.size(), 9)
}
END_SECTION

START_SECTION((Test Coarse and FineIsotopePatternGeneration))
{
  EmpiricalFormula formula = EmpiricalFormula("C200H299N51O62S5");
  IsotopeDistribution fine_dist = formula.getIsotopeDistribution(FineIsotopePatternGenerator(0.00005,ProbabilityMode::relative));
  IsotopeDistribution coarse_dist = formula.getIsotopeDistribution(CoarseIsotopePatternGenerator(100, false));

  // Log zur Kontrolle
  for (Size i = 0; i < coarse_dist.size(); ++i)
  {
    OPENMS_LOG_INFO << "COARSE " << i << ": MZ = " << coarse_dist[i].getMZ() << "  INT = " << coarse_dist[i].getIntensity() << std::endl;
  }
  for (Size i = 0; i < fine_dist.size(); ++i)
  {
    OPENMS_LOG_INFO << "FINE " << i << ": MZ = " << fine_dist[i].getMZ() << "  INT = " << fine_dist[i].getIntensity() << std::endl;
  }
  // Vergleich über Toleranz – Position und Intensität
  TEST_REAL_SIMILAR(coarse_dist[0].getMZ(), fine_dist[0].getMZ());
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
