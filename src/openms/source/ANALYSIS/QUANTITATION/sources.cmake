### the directory name
set(directory source/ANALYSIS/QUANTITATION)

### list all filenames of the directory here
set(sources_list
AbsoluteQuantitation.cpp
AbsoluteQuantitationMethod.cpp
DDAWorkflowCommons.cpp
IsobaricChannelExtractor.cpp
IsobaricIsotopeCorrector.cpp
IsobaricKitDetection.cpp
IsobaricNormalizer.cpp
IsobaricQuantifier.cpp
IsobaricQuantifierStatistics.cpp
IsobaricQuantitationMethod.cpp
IsotopeLabelingMDVs.cpp
ItraqConstants.cpp
ItraqEightPlexQuantitationMethod.cpp
ItraqFourPlexQuantitationMethod.cpp
KDTreeFeatureMaps.cpp
KDTreeFeatureNode.cpp
LabellingDetector.cpp
PeptideAndProteinQuant.cpp
TMTEighteenPlexQuantitationMethod.cpp
TMTElevenPlexQuantitationMethod.cpp
TMTSixPlexQuantitationMethod.cpp
TMTSixteenPlexQuantitationMethod.cpp
TMTTenPlexQuantitationMethod.cpp
TMTThirtyFivePlexQuantitationMethod.cpp
TMTThirtyTwoPlexQuantitationMethod.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\QUANTITATION" FILES ${sources})
