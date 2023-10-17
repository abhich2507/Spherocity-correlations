# Spherocity-correlations
Forward-backward Spherocity correlations in pp  collisions at sqrt(s)=13 TeV
# Description
1.)tut4b.cc is for generating events(10K).
2.) SpherocityAnalysis1x.C files contain the code for analysis.
tut4b.cc will provide output tut4.root that contains a tree having multiple branches. Using this output, SpherocityAnalysis1x.C gives the Spherocity for each event, and It made two trees, one for Spherocity in forward pseudorapidity and another for backward. Then I add both trees using the addingtree1() function in macro and perform calculations for correlation coefficients. In the end, I used it multiple times for different eta widths.
