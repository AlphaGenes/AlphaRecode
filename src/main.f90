program AlphaRecode
use pedigreeModule
use AlphaRecodeSpecFileModule
implicit none
type(alphaRecodeSpecFileHolder) :: spec
type(PedigreeHolder) :: pedigree
spec = alphaRecodeSpecFileHolder()  
    
pedigree = PedigreeHolder(spec%pedigreeFile)
call  pedigree%outputSortedPedigreeInAlphaImputeFormat(spec%outputFile)
call pedigree%destroyPedigree

end program AlphaRecode