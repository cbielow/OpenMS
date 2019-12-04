#include <OpenMS/config.h>

#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

#include <omp.h>
#include <numeric>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  std::cout << "usage: <mzMLfile> <#threads>\n";  
  omp_set_num_threads(atoi(argv[2]));

  //input file names
  String in = argv[1];
  std::cout << "Read method: indexed (parallel)" << std::endl;

  IndexedMzMLFileLoader imzml;

  // load data from an indexed MzML file
  OnDiscPeakMap map2;
  map2.openFile(in, true);
  map2.setSkipXMLChecks(true);

  double TIC = 0.0;
  long int nr_peaks = 0;
  auto nr_spec = (SignedSize)map2.getNrSpectra();

  // firstprivate means that each thread has its own instance of the
  // variable, each copy initialized with the initial value 
#pragma omp parallel for firstprivate(map2) 
  for (SignedSize i = 0; i < nr_spec; i++)
  {
    long nr_peaks_l;
    double TIC_l;
    if (1) // activate this line to cause a bug
    {
      auto dd = map2.getSpectrumById(i)->getIntensityArray()->data;
      nr_peaks_l = dd.size();
      TIC_l = std::accumulate(dd.begin(), dd.end(), 0.0);
    }

    //#pragma omp flush(TIC, nr_peaks)
#pragma omp critical
    {
      //std::cout << i << " " << nr_peaks_l << "\n";
      TIC = TIC + TIC_l;
      nr_peaks = nr_peaks + nr_peaks_l;
    }
  }

  std::cout << "There are " << map2.getNrSpectra() << " spectra and " << nr_peaks << " peaks in the input file." << std::endl;
  std::cout << "The total ion current is " << TIC << std::endl;


  return 0;
}

/// @endcond
