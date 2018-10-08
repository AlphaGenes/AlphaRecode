module PedigreeExtractionModule


contains

subroutine outputSortedPedigreeRecoded(pedigree, file, oldID, startGenerationIn, endGenerationIn)
    use PedigreeModule
    use iso_fortran_env, only : output_unit
    type(PedigreeHolder) :: pedigree
    character(len=*), intent(in), optional :: file
    logical, intent(in), optional :: oldID
    logical :: useOldId
    integer, intent(in), optional :: startGenerationIn, endGenerationIn
    character(len=IDLENGTH) :: tmpsireId,tmpDamId
    integer :: unit, i,h,startGeneration,endGeneration 
    type(IndividualLinkedListNode), pointer :: tmpIndNode
        
    if (.not. allocated(pedigree%generations)) then
        call pedigree%setPedigreeGenerationsAndBuildArrays
    endif
    print *,"finished building generation arrays"

    if (present(file)) then
        open(newUnit=unit, file=file, status="unknown")
    else
        unit = output_unit
    endif
    if (present(oldID)) then
        useOldId = oldID
    else 
        useOldId = .false.
    endif

    if (present(startGenerationIn)) then
        startGeneration = startGenerationIn
    else
        startGeneration = 0
    endif
    if (present(endGenerationIn)) then
        endGeneration = endGenerationIn
    else
        endGeneration = pedigree%maxGeneration
    endif  
    block 
        integer :: sireId, damId
        if (useOldId) then
            write (unit,'(6a20)') "originalID","Sire ID","DamID","generation","NumberOffsprings", "gender"
        else
            write (unit,'(7a20)') "New ID","Sire ID","DamID", "originalID","generation","NumberOffsprings", "gender"
        endif
        do i=startGeneration, endGeneration
            tmpIndNode => pedigree%generations(i)%first
            do h=1, pedigree%generations(i)%length
                    if (associated(tmpIndNode%item%damPointer)) then
                        damId = tmpIndNode%item%damPointer%id
                        tmpdamId = tmpIndNode%item%damPointer%originalID
                    else
                        damId = 0
                        tmpdamId = "0"
                    endif
                    if (associated(tmpIndNode%item%sirePointer)) then
                        sireId = tmpIndNode%item%sirePointer%id
                        tmpsireId = tmpIndNode%item%sirePointer%originalID
                    else
                        sireId = 0
                        tmpdamId=""
                    endif
                    if (useOldId) then
                        write (unit,'(3a20,3i20)') tmpIndNode%item%originalID,sireId,damId, tmpIndNode%item%originalID,tmpIndNode%item%generation,tmpIndNode%item%nOffs, tmpIndNode%item%gender
                    else
                        write (unit,'(3i20,a20,3i20)') tmpIndNode%item%id,sireId,damId, tmpIndNode%item%originalID,tmpIndNode%item%generation,tmpIndNode%item%nOffs, tmpIndNode%item%gender
                    endif
                    ! write(*,'(a,",",a,",",a,",",i8)') tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId,tmpIndNode%item%generation
                    tmpIndNode => tmpIndNode%next
            end do
        enddo
    endblock
    if (present(file)) then !avoids closing stdout
        close(unit)
    endif


end subroutine outputSortedPedigreeRecoded


end module PedigreeExtractionModule