#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/SpanningGraph.h>
#include <boost/regex.hpp>
#include <algorithm>
#include <numeric>
#include <math.h>  

using namespace OpenMS;
using namespace std;
