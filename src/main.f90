program AlphaRecode
use pedigreeModule
use AlphaRecodeSpecFileModule
use PedigreeExtractionModule, only : outputSortedPedigreeRecoded
implicit none
type(alphaRecodeSpecFileHolder) :: spec
type(PedigreeHolder) :: pedigree
spec = alphaRecodeSpecFileHolder()  
    
call initPedigree(pedigree,spec%pedigreeFile)
! call  pedigree%outputSortedPedigreeInAlphaImputeFormat(spec%outputFile)
call outputSortedPedigreeRecoded(pedigree, "out.txt", .false., spec%startGeneration, spec%endGeneration)

end program AlphaRecode