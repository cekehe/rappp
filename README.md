# rappp

This is the read me for the package rappp

ALWAYS UPDATE DEV NUMBER +1 BEFORE DOING A PULL REQUEST!

See further information about version numbering below. 

Ideas for the future:
---------------------

- colorful correlation values to pairs2 (Cecilia)


Some important commands for the terminal:
-----------------------------------------

git pull upstream master --ff-only

git pull upstream master

Some useful roxygen2 commands:
-----------------------------------------
\\cr - line break without empty line in between in neither code nor output.

This is how you should name versions: 
-------------------------------------

Version numbering of packages. 

0.0.0.9000

major.minor.patch.dev

Major: Large changes, not always backwards compatible. Usually 1 upon first release out of dev

Minor: Bug fixes & new features. Most common

Patch: Small bugfixes, no new features.

Dev: Only used while under development. Always starts at 9000

Install current master-version from GitHub: 
-------------------------------------------
install.packages("devtools")

library(devtools)

install_github("cekehe/rappp")

library(rappp)
