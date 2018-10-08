program AlphaRecode
use pedigreeModule
use AlphaRecodeSpecFileModule
use PedigreeExtractionModule, only : outputSortedPedigreeRecoded
implicit none
type(alphaRecodeSpecFileHolder) :: spec
type(PedigreeHolder) :: pedigree
spec = alphaRecodeSpecFileHolder()  
    
pedigree = PedigreeHolder(spec%pedigreeFile)
! call  pedigree%outputSortedPedigreeInAlphaImputeFormat(spec%outputFile)
call outputSortedPedigreeRecoded(pedigree, "out.txt", .false., spec%startGeneration, spec%endGeneration)
call pedigree%destroyPedigree

end program AlphaRecode