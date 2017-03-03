Wednesday 1st March 2017

<< Changes from "matchsticks.sh" (version 1) to "makeMatchsticks_v2.sh" (version 2) >>

The number of matchsticks being made has been increased to better ascertain any non-linearities in the ADCs we used.

Additionally, the matchsticks corresponding to 1V pulser voltage has double the statistics than the other matchsticks. This is to ensure we can easily pin/identify it in the matchsticks spectra.

The previous matchsticks.sh script took pulser voltages of:
0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0

In the new makeMatchsticks_v2.sh script, this has been expanded to:
0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0

As such, this new script will take a lot longer to run, but should provide much better matchsticks data.

Thursday 2nd March 2017

<< Changes from "makeMatchsticks_v2.sh" (version 2) to "makeMatchsticks_v3.sh" (version 3) >>

Removed the matchsticks at 0.01V and 0.05V and all those at 0.X5V (X=1-9), to avoid too many matchsticks cluttered up in a relatively small region of channels. Added matchsticks at 0.03V and 0.07V instead.
