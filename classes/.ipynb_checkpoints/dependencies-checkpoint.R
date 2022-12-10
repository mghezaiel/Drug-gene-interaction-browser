#brew install libgit2 (to install devtools)
packages = c("BiocManager","devtools","biomaRt","optparse","UniProt.ws")
for (package in packages){
      if (package=="UniProt.ws")
          {
          BiocManager::install("UniProt.ws")
          }
      else
          {
          install.packages(package)
      }
   
    }