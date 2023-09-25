This build system is based on GNU auto tools: https://www.gnu.org/software/automake/manual/html_node/GNU-Build-System.html
mkDepends and mkSrcfiles are from NCAR models. 

If one wants to customize some setup or debug something, only files that are needed to modify are: 

-build.sh
-configure.ac
-Makefile.am
-path/Filepath.XXX
