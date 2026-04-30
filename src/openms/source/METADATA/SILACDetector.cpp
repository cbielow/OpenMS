// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SILACDetector.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/KERNEL/MSExperiment.h>
//#include <map>
//#include <cmath>

namespace OpenMS
{

  // def process_spectrum(self, spectrum, binned_spectrum){...}
  /*
  if spectrum.charge < 1:        // ACHTUNG!!! in MSExperiment -> set to 1 if charge is 0, otherwise div/0 below => müssten schauen, ob charge = 1 ist statt < 1!
            # can't use spectra of unknown charge
            return
        self.scan_numbers.append(spectrum.scan_number)
        precursor_mass = calc_mplush_from_mz_charge(spectrum.precursor_mz, spectrum.charge)
        binidx = calc_binidx_for_mass_precursor(precursor_mass)
        self.precursor_mass_bins.append(binidx)
  */

  /*
  void process_spectrum(MSExperiment experiment, MSSpectrum spectrum, MSSpectrum binned_spectrum) //??????????????
  {
    //experiment.getCharge() // gibt es anscheinend nicht in MSExperiment, aber wird dort verwendet??? Wo bekommen wir das her?
    //this ??????

    double charge = 0.5; //(wäre hier testweise charge = 2), von wo auch immer wir das her bekommen :/
    if (charge == 1) // ACHTUNG!!! in MSExperiment -> set to 1 if charge is 0, otherwise div/0 below => müssten schauen, ob charge = 1 ist statt < 1!
    {
      return;
    }
    //this.scan_numbers.append(spectrum.scan_number) ?????????????

  }
  
  bool detectSILAC(MSExperiment experiment)
  {
    
    /*
    CONTROL_BIN_DISTANCES = [11, 14, 15, 21, 23, 27]
    # maximum separation in scans over which to count a pair as present
    MAX_SCAN_SEPARATION = 50

    SILAC_ZSCORE_CUTOFF = 4.0

     maps giving exact mass difference for nominal mass differences for K and R labels.
    # One suggestion is made at random from the options in Unimod

    // Kommt alles in den Konstruktor rein?

    std::map<int, double> silac_mod_K_exactmass_map =
    {
      {4, 4.025107},  //# http://www.unimod.org/modifications_view.php?editid1=481
      {6, 6.020129},  //# http://www.unimod.org/modifications_view.php?editid1=188
      {8, 8.014199},  //# http://www.unimod.org/modifications_view.php?editid1=259
      {10, 10.008269}  //# http://www.unimod.org/modifications_view.php?editid1=267
    };

    std::map<int, double> silac_mod_R_exactmass_map =
    {
      {6, 6.020129},  //# http://www.unimod.org/modifications_view.php?editid1=188
      {10, 10.008269}  //# http://www.unimod.org/modifications_view.php?editid1=267
    };

    // In Crux wurde das etwas anders implementiert -> mit pushBack
    
    // Sets sind bereits in C++ vorhanden -> brauchen kein Vector!!!
    std::set<int> control_dist_set = {11, 14, 15, 21, 23, 27};;
    std::set<int> silac_dist_set = {4, 6, 8, 10};

    std::set<int> separations_evaluate_set = silac_dist_set;
    set_union(separations_evaluate_set.begin(), separations_evaluate_set.end(), control_dist_set.begin(), control_dist_set.end(), separations_evaluate_set);
    set_union(set1.begin(), set1.end(), set2.begin(),
              set2.end(), inserter(result, result.begin()));

    
    std::vector<int> control_distances = {11, 14, 15, 21, 23, 27};
    std::vector<int> silac_distances = {4, 6, 8, 10};
    //separations_to_evaluate = set(SILAC_MOD_BIN_DISTANCES + SILACDetector.CONTROL_BIN_DISTANCES) -> Set von allen Massenabständen (4, 6, 8, 10, Kontrollwerte)

    // Wollen wir dem Nutzer erlauben, eigene Werte einzufügen? Vermutlich nicht.
    std::vector<int> separations_to_evaluate = {};

    // concatenate both distance_vectors in the separations_to_evaluate
    separations_to_evaluate.reserve( silac_distances.size() + control_distances.size() );
    separations_to_evaluate.insert( separations_to_evaluate.end(), silac_distances.begin(), silac_distances.end() );
    separations_to_evaluate.insert( separations_to_evaluate.end(), control_distances.begin(), control_distances.end() );

    // remove all duplicated by ordering the vector and putting these at the end
    sort(separations_to_evaluate.begin(), separations_to_evaluate.end());
    auto it = unique(separations_to_evaluate.begin(), separations_to_evaluate.end());
    separations_to_evaluate.erase(it, separations_to_evaluate.end());

    //# paranoia
    //    if len(separations_to_evaluate) < len(SILAC_MOD_BIN_DISTANCES) + len(SILACDetector.CONTROL_BIN_DISTANCES):
    //        logger.warn("A specified separation is also a control separation! Specified: %s" % str(SILAC_MOD_BIN_DISTANCES))

    //# initialize a map from separation distances to counts of pairs with that separation
    //    counts_with_separations = {}
    //    for separation in separations_to_evaluate:
    //        counts_with_separations[separation] = 0

    std::map<int, int> counts_with_separations = {};
    for (int i = 0; i < separations_to_evaluate.size(); i++)
    {
      counts_with_separations.insert({separations_to_evaluate[i], 0});
    }

    //# keep track of the scan window defined by the minimum and maximum scan index to consider
    //    minidx = 0
    //    maxidx = 0

    int min_idx = 0;
    int max_idx = 0;

    /*for i in xrange(0, len(self.scan_numbers)):
            # determine the minimum and maximum scan number currently in range
            scan_number = self.scan_numbers[i]
            min_scan_number = scan_number - SILACDetector.MAX_SCAN_SEPARATION
            max_scan_number = scan_number + SILACDetector.MAX_SCAN_SEPARATION
            while self.scan_numbers[minidx] < min_scan_number:
                minidx += 1
            while self.scan_numbers[maxidx] < max_scan_number and maxidx < len(self.scan_numbers) - 1:
                maxidx += 1

            # within the scan window, increment the separations that we care about with any that involve this scan
            for j in xrange(minidx, maxidx):
                separation = abs(self.precursor_mass_bins[i] - self.precursor_mass_bins[j])
                if separation in separations_to_evaluate:
                    counts_with_separations[separation] += 1

    // ich mache hier einfach mal irgendwas rein als Test, wird dann richtig später eingesetzt wenn ich weiß wie...
    // Sollte irgendwie intern berechnet und als Variable angelegt werden, aber WTF?
    std::vector<int> scan_numbers = {2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 16, 20, 23, 27, 29, 35, 39, 44, 46, 48, 50}; // ????????????
    int scan_number = 0;
    double min_scan_number = 0;
    double max_scan_number = 0;
    for (int i = 0; i < scan_numbers.size(); i++)
    {
      scan_number = scan_numbers[i];
      min_scan_number = scan_number - 0.5; //SILACDetector.MAX_SCAN_SEPARATION ??
      max_scan_number = scan_number + 0.5; //SILACDetector.MAX_SCAN_SEPARATION ??

      while (scan_numbers[min_idx] < min_scan_number){
        min_idx++;
      }
      while (scan_numbers[max_idx] < max_scan_number){
        max_idx++;
      }
      /* # within the scan window, increment the separations that we care about with any that involve this scan
      for j in xrange(minidx, maxidx):
                 separation = abs(self.precursor_mass_bins[i] - self.precursor_mass_bins[j]) -> for separation in separations_to_evaluate
                 if separation in separations_to_evaluate:
                     counts_with_separations[separation] += 1
      int separation = 0;
      for (int j = min_idx; j < max_idx; j++)
      {
        // experiment.getPrecursorSpectrum() gibt es
        //separation = abs(self.precursor_mass_bins[i] - self.precursor_mass_bins[j])
        separation = abs(experiment.getPrecursorSpectrum(i) - experiment.getPrecursorSpectrum(j)); // ??? MSExperiment.h Zeile 1173 könnte helfen?
        // if separation in separations_to_evaluate:
        auto iterator = find(separations_to_evaluate.begin(), separations_to_evaluate.end(), separation);
        if (iterator != separations_to_evaluate.end())
        {
          counts_with_separations[separation] += 1;
        }
      } // for j
    } // for i

    /*# summarize the control separations -> Brauchte ChatGPT um die Syntax zu verstehen wegen dem "x for y in z" !??!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!
        mean_control_count = (float(sum([counts_with_separations[separation] for separation in
                                         SILACDetector.CONTROL_BIN_DISTANCES])) /
                              len(SILACDetector.CONTROL_BIN_DISTANCES))
    // könnten auch double nehmen?
    double mean_control_count = 0.0;
    int sum = 0;
    for(int i = 0; i < control_distances.size(); i++)
    {
      sum = sum + counts_with_separations[i];
      mean_control_count = static_cast<double>(sum) / static_cast<double>(control_distances.size());
      //mean_control_group = mean_control_group + static_cast<double>(counts_with_separations[i]);
    }
    /* Ist nur Log, können wir erst einmal weglassen
    for separation in SILACDetector.CONTROL_BIN_DISTANCES:
            logger.debug("  %d: %d" % (separation, counts_with_separations[separation]))
        logger.debug("SILAC: Mean control separation count: %.05f" % mean_control_count)*/

    /*if mean_control_count > 0:
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug("SILAC: Counts for each separation:")
                for separation in SILAC_MOD_BIN_DISTANCES:
                    proportion_to_control = float(counts_with_separations[separation]) / mean_control_count
                    logger.debug("  %d: %d (proportion=%.05f)" % (separation, counts_with_separations[separation], proportion_to_control))
    
    if (mean_control_count > 0)
    {
      //if logger.isEnabledFor(logging.DEBUG):
      //          logger.debug("SILAC: Counts for each separation:") ???

      //          for separation in SILAC_MOD_BIN_DISTANCES:
      //              proportion_to_control = float(counts_with_separations[separation]) / mean_control_count
      //              logger.debug("  %d: %d (proportion=%.05f)" % (separation, counts_with_separations[separation], proportion_to_control))
      double proportion_to_control = 0.0;
      for (int i = 0; i < silac_distances.size(); i++)
      {
        proportion_to_control = static_cast<double>(counts_with_separations[i]) / mean_control_count;
      }
    }
    else // TODO!!!!!!!!!!
    {
      //logger.warn("SILAC: No counts for any control separation pairs! Cannot estimate prevalence of SILAC separations.")
      //      # make a dummy result with no significant inferences

      // TODO!!!!!!!!!!
      // RunAttributeResult() -> "Holds a result of a RunAttributeDetector analysis" ???
      /*for separation in SILAC_MOD_BIN_DISTANCES:
                result.name_value_pairs['SILAC_%dDa_present' % separation] = 'ERROR'
                result.name_value_pairs['SILAC_%dDa_statistic' % separation] = 'ERROR'
      
      //          return result
      

    }
    // control_sd = np.std([counts_with_separations[separation] for separation in SILACDetector.CONTROL_BIN_DISTANCES])

    // CONGRATULATIONS!!! I COULD NOT FIND A FUNCTION WITH STANDARD DEVIATION!!! I found something in the Documentation about stddev(), but where is it?

    // TODO

    std::vector<int> sum_vec = {};
    for(int i = 0; i < control_distances.size(); i++)
    {
      sum_vec.push_back(counts_with_separations[i]);
    }

    //double control_sd = stddev(sum_vec)
    double control_sd = 3.1; // nur damit es "läuft"...?
    
    //result = util.RunAttributeResult() // ???
    int result = 0;
    std::vector<int>significant_separations = {};

    /*for separation in SILAC_MOD_BIN_DISTANCES:
            # z-score is checked against a cutoff to determine significance
            zscore_to_control = float(counts_with_separations[separation] - mean_control_count) / control_sd
    
    double zscore_to_control = 0.0;
    for (int i = 0; i < control_distances.size(); i++)
    {
      zscore_to_control = static_cast<double>(counts_with_separations[i] - mean_control_count) / control_sd;

      /*if zscore_to_control > SILAC_ZSCORE_CUTOFF:
                significant_separations.append(separation)
                logger.info("SILAC: %dDa separation detected." % separation)
      double zscore_cutoff = 4.0; // müsste dann noch weg
      if (zscore_to_control > zscore_cutoff)
      {
        significant_separations.push_back(i);
        //# paranoia
        //        if separation not in SILAC_MOD_K_EXACTMASS_MAP and separation not in SILAC_MOD_R_EXACTMASS_MAP:
        //            raise ValueError('Unknown SILAC separation %d' % separation)
        
      /*# figure out the exact appropriate mass for search.
                if separation in SILAC_MOD_K_EXACTMASS_MAP:
                    result.search_modifications.append(util.Modification("K", SILAC_MOD_K_EXACTMASS_MAP[separation], True))
                if separation in SILAC_MOD_R_EXACTMASS_MAP:
                    result.search_modifications.append(util.Modification("R", SILAC_MOD_R_EXACTMASS_MAP[separation], True))
                result.name_value_pairs['SILAC_%dDa_present' % separation] = 'T'

        // Diese sind jetzt fast ganz oben in der Funktion gelandet -> TODO: nochmal besseren Ort finden!

        auto iterator = find(silac_mod_K_exactmass_map.begin(), silac_mod_K_exactmass_map.end(), i);
        if (iterator != silac_mod_K_exactmass_map.end())
        {
          //result.search_modifications.append(util.Modification("K", SILAC_MOD_K_EXACTMASS_MAP[separation], True))
        }
        iterator = find(silac_mod_R_exactmass_map.begin(), silac_mod_R_exactmass_map.end(), i);
        if (iterator != silac_mod_R_exactmass_map.end())
        {
          //result.search_modifications.append(util.Modification("R", SILAC_MOD_R_EXACTMASS_MAP[separation], True))
        }
        //result.name_value_pairs['SILAC_%dDa_present' % separation] = 'T'

      }
    }

    return false;
  }*/

  bool SILACDetector::detectSILAC(MSExperiment experiment)
  {
    int spectrum_number = 0;
    std::vector<int> control_distances = {11, 14, 15, 21, 23, 27};
    std::vector<int> silac_distances = {4, 6, 8, 10};

    std::map<int,int> distance_count = {{4,0},{6,0},{8,0},{10,0},{11,0},{14,0},{15,0},{21,0},{23,0},{27,0}};
    double RT_window = 5;
    int min_index = 0;
    int max_index = 0;
    
    
    // vielleicht von anfang an nur MS2 spektren nehmen

     for (const auto& spectrum : experiment)
    {
       
      std::cout << "Spectrumnumber: " << spectrum_number << std::endl;
      spectrum_number++;
    
      if (2 == spectrum.getMSLevel())
      {

       
        double spectrum_rt = spectrum.getRT();
        while (experiment[min_index].getRT() < spectrum.getRT() - RT_window)
        {
          min_index++;
        } 
        while ((experiment[max_index].getRT() < spectrum.getRT() + RT_window) && max_index < experiment.size())
        {
          max_index++;
        } 
        for (int i = min_index; i < max_index; i++)
        {
          if (2 == experiment[i].getMSLevel())
          {
            if (spectrum.getPrecursors()[0].getCharge() == experiment[i].getPrecursors()[0].getCharge())
            {
              double distance = spectrum.getPrecursors()[0].getMZ() - experiment[i].getPrecursors()[0].getMZ();
              distance *= spectrum.getPrecursors()[0].getCharge();
              distance = std::abs(distance);
              distance += 0.5;
              int a = distance;
              //std::cout << distance << std::endl;
              if (distance_count.find(a) != distance_count.end())
              {
                distance_count[a]++;
              }
            }
            
          }
          
        }
      }

      
      
      
    }
    double n = 6;

    double control_mean = (distance_count[11] + distance_count[14] + distance_count[15] + distance_count[21] + distance_count[23] + distance_count[27]) / 6;
    double control_sd = 0;
    double s;
    double sd_sum=0;
    for (int i = 0; i < 6; i++)
    {
      s = distance_count[control_distances[i]]-control_mean;
      s *= s;
      sd_sum += s;

    }
    control_sd = std::sqrt(sd_sum/(n-1));
 
    
    


   /*  for (long unsigned int i = 0; i < experiment.size();i++)
  {
    
    if(experiment[i].getMSLevel()==2)
    {
      experiment[i].getPrecursors()[0].getCharge();
      experiment[i].getPrecursors()[0].getMZ();
      experiment[i].getPrecursors()[0].getUnchargedMass();
    }

  } */
    std::cout << "Distanz 4: " << distance_count[4] << " Zscore: " << (distance_count[4]-control_mean)/control_sd<< std::endl;
    std::cout << "Distanz 6: " << distance_count[6] << " Zscore: " << (distance_count[6]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 8: " << distance_count[8] << " Zscore: " << (distance_count[8]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 10: " << distance_count[10] << " Zscore: " << (distance_count[10]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 11: " << distance_count[11] << " Zscore: " << (distance_count[11]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 14: " << distance_count[14] << " Zscore: " << (distance_count[14]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 15: " << distance_count[15] << " Zscore: " << (distance_count[15]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 21: " << distance_count[21] << " Zscore: " << (distance_count[21]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 23: " << distance_count[23] << " Zscore: " << (distance_count[23]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 27: " << distance_count[27] << " Zscore: " << (distance_count[27]-control_mean)/control_sd<<std::endl;

    
    return true;

  }

  int SILACDetector::testfunktion()
  {
    return 1;
  }
} // namespace OpenMS