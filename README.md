# R implementation of JICO

All relevant helper functions have been organized into [continumm/function/JICO.R](https://github.com/peiyaow/continuum/blob/master/function/JICO.R).
- The main function that conducts the multigroup JICO algorithm is [continuum.multigroup.iter](https://github.com/peiyaow/continuum/blob/master/function/JICO.R#L160); we also implement [cv.continuum.iter](https://github.com/peiyaow/continuum/blob/master/function/JICO.R#L447) to tune gamma and ranks using cross-validation.
- The function that conducts general continuum algorithm is [continuum](https://github.com/peiyaow/continuum/blob/master/function/JICO.R#L41).

We provide a [toy example script](https://github.com/peiyaow/continuum/blob/master/sim/toy.R) to get you familarize with the code.

We also include all the scripts we used for simulations in [continumm/sim](https://github.com/peiyaow/continuum/tree/master/sim). 

For further inquiry please reach out to the repository owner: peiyaow76@gmail.com. We will address your questions at our earliest convenience. 
