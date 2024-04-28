# project tasks
  
  - to do (optimization)
    - change Newton_Raphson to Fisher_Scoring (name only)
    - make wrapper less awkward (somewhat done)
    - (maybe) `c++` implementation of the optimization (doesn't really seem worth it,
      but go follow your dreams Manu ;) ) 
    
  - to do (plots)
  
    - (maybe) grid scatter plot
    - (maybe) fully integrate shiny plotting
    
  - to do (simulation study) ###expand this point###
      - include scenarios with violations of the error assumptions
      - compare different optimizers wrt their bias
      - make some nice plots for everything
  
  - to do (general):
    - (maybe) split up R scripts (class/optimization/plotting)
    - documentation
    - expand README 
      - finish writing up derivation
      - add some plots
    - add `verbose` messages
    - include more user feedback (mostly meaningful error messages)
    
  - to do (testing):
    - figure out why expect_equivalent throws an error even though
      derivation is only e-7 
    - test additions to the class (residuals, other active fields)
    - some more scenarios for the optimizers
    
  - done:
    - read into theory
    - establish git framework, basic description etc.
    - basic `plot` method functionality
    - WLS algorithm
    - gradient descent
    - newton raphson
    - midterm report
    - discuss which plots to offer
    - more test routines
    - second midterm report
    - implement first set of diagnosis plots
    - cook's distance and influence plot
    - some more plots regarding the scale/std (look at gamlss package)
    - make plots prettier
    - spread level plot
    - allow for output of single plots 
    - split up plot function into influence_plot, etc.
    - include option to label interesting observations 
    - improve wrapper / user interaction with plotting feature
    
    