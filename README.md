DRIVER.m - run this file to find fixed points of the tick-host system and determine local stability.

tickMap.m - this is called by the "mySystem.m" and "mySystem2.m" files during the system solving process; describes the tick-host model.

mySystem.m - this is called by the "DRIVER.m" and "colorgraphs.m" files during the system solving process; it iterates the tickMap function once to find fixed points.

mySystem2.m - this is called by the "DRIVER.m" and "colorgraphs.m" files during the system solving process; it iterates the tickMap function twice to find 2-cycle points.

makeJacobian.m - this is called by the "DRIVER.m" and "colorgraphs.m" files during to assess local stability; it numerically finds the Jacobian of the tick-host system.

plotter.m - this script produces Figures 2, 3, and 4 from our manuscript, depicting tick populations with respect to the demographic reproductive number and the questing behavior parameter c. This script should be run after running DRIVER.m

colorgraphs.m - this script produces Figures 5, 6, and 7 from our manuscript, which depicts the change in stability when key demographic parameters are varied.
