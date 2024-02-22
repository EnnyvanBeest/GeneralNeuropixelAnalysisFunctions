# GeneralNeuropixelAnalysisFunctions
Bunch of scripts that might be useful for other people doing Neuropixel analysis in Matlab,
who prefer to integrate everything in one Matlab pipeline rather than converting data back and forth

This has been developed using PyKS2 (Mouseland/Kush) and phy2 spikesorting

Functions in Histology might be useful for general public
Functions in SyncImecAndTimeline are most useful when you record timestamps on a separate system and want to match those to 
Imec time, either directly or via NIDQ.

ExampleMPEPAnalysis.m included as a helping hand to put these functions in a logical position within your analysis. 
This script assumes data is at least processed by Kilosort, and preferably also spikesorted with Phy.
It also requires some form of histology output (i.e. coordinates with area names, see the the histology scripts for details)

Toolboxes used:
https://github.com/cortex-lab/spikes
https://github.com/petersaj/AP_histology
https://github.com/kwikteam/npy-matlab

Pre-processing toolboxes:
https://github.com/cortex-lab/KiloSort
https://github.com/cortex-lab/phy
https://github.com/cortex-lab/allenCCF
https://github.com/brainglobe

If there are issues I'll do my best to help, provided that these are easy fixes. However this toolbox is not officially maintained to be working for everyone's individual needs and my priority is doing science. On that note, if others have made improvements, please feel free to submit pull requests to make them more useful for everyone.

