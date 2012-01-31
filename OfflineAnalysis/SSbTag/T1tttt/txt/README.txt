Format of text file:

At the top, can have comment lines preceeded by #


Column1: gluino mass

Column2: LSP mass

Column3: efficiency = numb events accepted/numb events generated

Column4: 1 + error due to JES and btag (the error from lepton eff is
         not included since it is constant)

Column5: Observed upper limit --- if it is not calculated, put 0.0
         (it can be calculated on the fly)

Column6: Expected upper limit --- ditto: 0.0 if not calculated



Column3 should eventually include the correction for btag scale 
factors, trigger efficiencies, and be vtx reweighted.
