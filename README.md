# rappp

This is the read me for the package rappp

FOR NEW RELEASES, IE PULL REQUESTS:
-----------------------------------
- UPDATE DEV NUMBER +1 IN DESCRIPTION
- UPDATE DATE IN DESCRIPTION
- DO PULL REQUEST
- CREATE A .tar.gz-file AND PUT ON BOX (PAPP>Data>R codes>R_packages>rappp (under development))

See further information about version numbering and creating tar.gz-files below. 

Ideas for the future:
---------------------




Some important commands for the terminal:
-----------------------------------------

git pull upstream master --ff-only # Fast forward merge, try this first

git pull upstream master # If Fast forward is not possible and you have to solve a merge conflict.

git reset HEAD~1 # Undo the last commit (for Windows), but keep all changes (the changes in the last commit will be unstaged).
Some good info about the reset command at https://stackoverflow.com/questions/2530060/in-plain-english-what-does-git-reset-do

Some useful roxygen2 commands:
------------------------------
\\cr - line break without empty line in between in neither code nor output.

Some useful general package tips:
---------------------------------
https://kbroman.org/pkg_primer/pages/depends.html - Good info on dependencies.

https://www.w3schools.com/cssref/default.asp - Good reference list for css-syntax (eg useful in Rmarkdown).

Create .tar.gz-file - In RStudio terminal (go one folder up and create the file there):<br/>
```
cd..
R CMD build rappp
```

Load RData-file and save loaded variable as .rda instead:<br/>
```
load("data/Scrambled_data_for_rappp.RData")
usethis::use_data(MockSBA, overwrite=T)
```

Console commands useful when updating package locally:<br/>
```
devtools::build_vignettes() # Create vignette so it will be included in package.
devtools::document() # Update any documentation, like help or vignette.
# Cmd/Ctrl + Shift + B  - bulid, ie. install current package content.
browseVignettes("packagename") # Read a package vignette.
```


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
```
install.packages("devtools")<br/>
library(devtools)<br/>
install_github("cekehe/rappp")<br/>
library(rappp)
```
