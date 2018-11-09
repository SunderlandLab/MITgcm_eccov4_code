README -- describes contents of MITgcm_code
14 Jan 2016
GEOS-Chem Support Team
geos-chem.as.harvard.edu


Overview:
==========

This directory is a local Git repository of the MITgcm code for use by
members of the Harvard Biogeochemistry  of Global Pollutants Group,
(aka Harvard BGC group), which is led by Prof. Elsie Sunderland and
Dr. Katherine von Stackelberg.

The remote repository is stored on Bitbucket.org at:

  https://bitbucket.org/harvard-bgc/mitgcm_code


To check out a fresh version of this repository, use this command:

  git clone https://USERNAME@bitbucket.org/harvard-bgc/mitgcm_code MITgcm_code

  where USERNAME is your Bitbucket user name. 
  You will then be prompted to supply your Bitbucket password as well.


To pull further updates from this repository, use this command:

 
  git pull origin BRANCH-NAME

  where BRANCH-NAME is the name of the branch that you want to get.
  Most of the time you will be pulling from the master branch, i.e.

  git pull origin master

 
NOTE: The initial import from the MITgcm CVS repository was done on 
12 Jan 2016.  We were not able to preserve the MITgcm version history.



Documentation:
==============

To view the list of changes (aka "commits") that were made, you can use
the gitk browser.

  cd MITgcm_code
  gitk --all &

For more information about the MITgcm, please view the following resources:

  Harvard BGC Group home:
  http://bgc.seas.harvard.edu
   

  Harvard BGC Group wiki:
  https://wiki.harvard.edu/confluence/display/BGC/Home

  MITgcm User Manual
  http://mitgcm.org/public/r2_manual/latest/online_documents/manual.html