!#############################################################################################################################################################################################################################

module PedigreeTable
    implicit none
    integer,parameter :: lengan=20
    public :: pedigreeLine
    
    type pedigreeLine
        character*(lengan) :: ID
        character*(lengan) :: Sire
        character*(lengan) :: Dam
        integer :: generation

    end type pedigreeLine
    type(pedigreeLine), parameter :: DICT_NULL = pedigreeLine('','','', 0)

    interface operator ( == )
      module procedure comparepedigreeline
   end interface operator ( == )
    contains

        logical function comparePedigreeLine(l1,l2)
            class(pedigreeLine), intent(in) :: l1,l2

            if (l1%ID == l2%ID .and. l1%Sire == l2%Sire .and. l1%Dam == l2%Dam) then
                comparePedigreeLine=.true.
            else 
                comparePedigreeLine=.false.
            endif


            return
        end function comparePedigreeLine




end module PedigreeTable



module dictModule
use PedigreeTable, DICT_DATA => pedigreeLine

implicit none 
integer, parameter :: DICT_KEY_LENGTH = lengan


include "dictionary.f90"


end module dictModule

!#############################################################################################################################################################################################################################

module GlobalPedigree
use PedigreeTable

implicit none


real(kind=4),allocatable :: xnumrelmatHold(:)
integer :: NRMmem, shell, shellmax, shellWarning,nAnisP,nAnisRawPedigree
integer,allocatable:: seqid(:),seqsire(:),seqdam(:),RecodeGenotypeId(:),passedorder(:),RecPed(:,:)
character*(lengan),allocatable :: Id(:),sire(:),dam(:)
type(pedigreeLine), allocatable :: ped(:)
integer :: GlobalExtraAnimals   !Change John Hickey


end module GlobalPedigree



!#############################################################################################################################################################################################################################

program AlphaRecode
use GlobalPedigree
implicit none

character (len=1000) :: dumC,FilIn
integer :: i

open (unit=1,file="AlphaRecodeSpec.txt",status="unknown")
read (1,*) dumC,FilIn

open (unit=2,file=trim(FilIn),status="unknown")
open (unit=3,file="AlphaRecodedPedigree.txt",status="unknown")

call CountInData
call ReadInData
call sortPedigree
! call PVseq(nAnisRawPedigree,nAnisP)

! allocate(RecPed(0:nAnisP,3))
  
! RecPed(0,:)=0
! do i=1,nAnisP
!   RecPed(i,1)=i
! enddo
! RecPed(1:nAnisP,2)=seqsire(1:nAnisP)
! RecPed(1:nAnisP,3)=seqdam(1:nAnisP)

! do i=1,nAnisP
!   write (3,'(3i20,a2,a20)') RecPed(i,:),'  ',Id(i)
! enddo


! deallocate(seqid)
! deallocate(seqsire)
! deallocate(seqdam)

end program AlphaRecode

!#############################################################################################################################################################################################################################

subroutine CountInData
use GlobalPedigree
implicit none

integer :: k
character (len=300) :: dumC

nAnisRawPedigree=0
do
  read (2,*,iostat=k) dumC
  nAnisRawPedigree=nAnisRawPedigree+1
    if (k/=0) then
        nAnisRawPedigree=nAnisRawPedigree-1
      exit
      endif
enddo
rewind(2)

print*, " ",nAnisRawPedigree," individuals in the pedigree file"

end subroutine CountInData

!#############################################################################################################################################################################################################################

subroutine ReadInData
use GlobalPedigree
implicit none

integer :: i,j,k
character(len=300) :: dumC

allocate(Ped(nAnisRawPedigree))

do i=1,nAnisRawPedigree
    read(2,*) ped(i)%id, ped(i)%sire, ped(i)%dam
    ! ped(i)%generation = 0
enddo

end subroutine ReadInData




!#############################################################################################################################################################################################################################
subroutine sortPedigree
use GlobalPedigree
use dictModule

implicit none

type(DICT_STRUCT), pointer     :: dict



type(pedigreeLine), allocatable :: SortedPed(:)
type(pedigreeLine) :: temp
integer :: i = 1, maxNumberToLookAt
integer(kind=1) :: switch = 0



maxNumberToLookAt = size(ped)

do while (maxNumberToLookAt /= 0)
  

  ! If founder
  if(ped(i)%sire == '0' .and. ped(i)%dam == '0') then
    if (maxNumberToLookAt== size(ped)) then
      call dict_create( dict, ped(i)%id, ped(i) )
    else 
      call dict_add_key( dict, ped(i)%id, ped(i) )
    endif
    write(3,*) ped(i)%ID,ped(i)%sire,ped(i)%dam
    temp = ped(maxNumberToLookAt)
    ped(maxNumberToLookAt) = ped(i)
    ped(i) = temp
    maxNumberToLookAt = maxNumberToLookAt - 1 

  endif 

  ! If both parents have been defined
  if (dict_has_key(dict,ped(i)%sire) .and. dict_has_key(dict,ped(i)%dam)) then
    call dict_add_key( dict, ped(i)%id, ped(i))
    write(3, *) ped(i)%ID,ped(i)%sire,ped(i)%dam
    temp = ped(maxNumberToLookAt)
    ped(maxNumberToLookAt) = ped(i)
    ped(i) = temp
    maxNumberToLookAt = maxNumberToLookAt - 1 
  endif 

  ! if one pedigree is unknown and loop has already executed once
  if (switch == 2 .and. (ped(i)%sire == '0' .or. ped(i)%dam == '0')) then
    call dict_add_key( dict, ped(i)%id, ped(i) )
    write(3,*) ped(i)%ID,ped(i)%sire,ped(i)%dam
    temp = ped(maxNumberToLookAt)
    ped(maxNumberToLookAt) = ped(i)
    ped(i) = temp
    maxNumberToLookAt = maxNumberToLookAt - 1 
    switch = 0
  endif

  i = i + 1 
  if (i > maxNumberToLookAt) then
    switch = switch +1
    i = 1 ! reset counter to avoid decuring
  endif
enddo

call dict_destroy( dict )
close(3)
end subroutine sortPedigree


!#############################################################################################################################################################################################################################

subroutine PVseq(nObs,nAnisPedigree)
USE GlobalPedigree
implicit none

character (LEN=lengan), ALLOCATABLE :: holdsireid(:), holddamid(:)
character (LEN=lengan), ALLOCATABLE :: holdid(:), SortedId(:), SortedSire(:), SortedDam(:)
character (LEN=lengan)              :: IDhold
integer, ALLOCATABLE                :: SortedIdIndex(:), SortedSireIndex(:), SortedDamIndex(:)
integer, ALLOCATABLE                :: OldN(:), NewN(:), holdsire(:), holddam(:)
INTEGER :: mode    ! mode=1 to generate dummy ids where one parent known.  Geneprob->1  Matesel->0
INTEGER :: i, j, k, kk, newid, itth, itho, ihun, iten, iunit
integer :: nsires, ndams, newsires, newdams, nbisexuals, flag
INTEGER :: ns, nd, iextra, oldnobs, kn, kb, oldkn, ks, kd
INTEGER :: Noffset, Limit, Switch, ihold, ipoint
integer :: nObs,nAnisPedigree,verbose
character (LEN=lengan) :: path

mode=1 

allocate(id(0:nobs),sire(nobs),dam(nobs),seqid(nobs),seqsire(nobs),seqdam(nobs))

do i=1,nobs

        id(i)=                ped(i)%id

        sire(i)=ped(i)%sire

        dam(i)=ped(i)%dam

end do



nAnisPedigree=nObs
path=".\"

Verbose=1


do j = 1, nobs

  If (dam(j) == ''.or. dam(j) == '0'.or. dam(j) == '#'.or. dam(j) == '*' .or. dam(j) == '.') Then

    dam(j) = '0'

    seqdam(j)=0

  endif

  If (sire(j) == ''.or.sire(j) == '0'.or.sire(j) == '#'.or.sire(j) == '*'.or.sire(j) == '.') Then

    sire(j) = '0'

    seqsire(j)=0

  endif

enddo !j



if(mode.eq.1) then

 !PRINT*,  ' Inserting dummy IDs ... '

 newid=0

 do j = 1, nobs

   if(((sire(j) == '0').and.(dam(j).ne.'0'))  .or. ((sire(j).ne.'0').and.(dam(j) == '0'))) then

         newid=newid+1

         if(newid.gt.99999) then

  !         PRINT*, newid, ' ...'

           stop 'too many dummy single parent IDs'

         endif

         itth=int(newid/10000)

         itho=int(newid/1000)-10*itth

         ihun=int(newid/100)-10*itho-100*itth

         iten=int(newid/10)-10*ihun-100*itho-1000*itth

         iunit=newid-10*iten-100*ihun-1000*itho-10000*itth

         if(sire(j) == '0') sire(j)='dum'//achar(48+itth)//achar(48+itho)//achar(48+ihun)//achar(48+iten)//achar(48+iunit)

         if( dam(j) == '0')  dam(j)='dum'//achar(48+itth)//achar(48+itho)//achar(48+ihun)//achar(48+iten)//achar(48+iunit)

   endif

 enddo

endif



!PRINT*,  ' Sorting Sires ... '



ALLOCATE  (SortedId(nobs), SortedIdIndex(nobs))



SortedId(1:nobs) = Sire(1:nobs)



  Noffset = INT(nobs/2)

  DO WHILE (Noffset>0)

      Limit = nobs - Noffset

      switch=1

    DO WHILE (Switch.ne.0)

       Switch = 0

       do i = 1, Limit

          IF (SortedId(i).gt.SortedId(i + Noffset)) THEN

               IDhold=SortedId(i)

               SortedId(i)=SortedId(i + Noffset)

               SortedId(i + Noffset)=IDhold



               Switch = i

          endif

       enddo

       Limit = Switch - Noffset

    enddo

    Noffset = INT(Noffset/2)

  enddo



nsires=0

IF(SortedId(1) /= '0') nsires=1

do i=2,nobs

  IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nsires=nsires+1

end do



ALLOCATE  (SortedSire(0:nsires), SortedSireIndex(nsires))

SortedSire(0) = '0'



nsires=0

IF(SortedId(1) /= '0') THEN

 nsires=1

 SortedSire(1) = SortedId(1)

ENDIF

do i=2,nobs

  IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') then

   nsires=nsires+1

   SortedSire(nsires) = SortedId(i)

  ENDIF

end do



!PRINT*,  ' Sorting Dams ... '



SortedId(1:nobs) = Dam(1:nobs)



  Noffset = INT(nobs/2)

  DO WHILE (Noffset>0)

      Limit = nobs - Noffset

      switch=1

    DO WHILE (Switch.ne.0)

       Switch = 0

       do i = 1, Limit

          IF (SortedId(i).gt.SortedId(i + Noffset)) THEN

               IDhold=SortedId(i)

               SortedId(i)=SortedId(i + Noffset)

               SortedId(i + Noffset)=IDhold



               Switch = i

          endif

       enddo

       Limit = Switch - Noffset

    enddo

    Noffset = INT(Noffset/2)

  enddo



nDams=0

IF(SortedId(1) /= '0') nDams=1

do i=2,nobs

  IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nDams=nDams+1

end do



ALLOCATE  (SortedDam(0:nDams), SortedDamIndex(ndams))

SortedDam(0)='0'



nDams=0

IF(SortedId(1) /= '0') THEN

 nDams=1

 SortedDam(1) = SortedId(1)

ENDIF

do i=2,nobs

  IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') then

   nDams=nDams+1

   SortedDam(nDams) = SortedId(i)

  ENDIF

end do





!PRINT*,  ' Sorting IDs ... '



SortedId(1:nobs) = ID(1:nobs)

do i=1,nobs

 SortedIdIndex(i) = i

end do



  Noffset = INT(nobs/2)

  DO WHILE (Noffset>0)

      Limit = nobs - Noffset

      switch=1

    DO WHILE (Switch.ne.0)

       Switch = 0

       do i = 1, Limit

          IF (SortedId(i).gt.SortedId(i + Noffset)) THEN

               IDhold=SortedId(i)

               SortedId(i)=SortedId(i + Noffset)

               SortedId(i + Noffset)=IDhold



               ihold=SortedIdIndex(i)

               SortedIdIndex(i)=SortedIdIndex(i + Noffset)

               SortedIdIndex(i + Noffset)=ihold



               Switch = i

          endif

       enddo

       Limit = Switch - Noffset

    enddo

    Noffset = INT(Noffset/2)

  enddo



!PRINT*,  ' Check for duplicate IDs ... '

flag = -1

Do i = 2, nobs

  If (SortedID(i) == SortedID(i - 1)) Then

   If (flag == -1) Then

     open (1,FILE='ID_err.txt',STATUS = 'unknown')

     WRITE(1,*) 'Duplicated IDs ...'

     flag = 0

   End If

   WRITE(1,*) SortedID(i)

   flag = flag + 1

  End If

enddo



 If (flag > -1) Then

  Close (1)

 ! PRINT*, flag,' case(s) of duplicate ID. See ID_ERR.TXT        <------------ WARNING !!!'

 End If



!PRINT*,  ' Males ... '

!PRINT*,  '  Find or set sire indices ... '



newsires = 0



do j=1,nsires



! check if already listed as an individual

   ipoint=INT(nobs/2)

   Noffset = INT(ipoint/2)

   do while (Noffset>1)

    IF (SortedSire(j).lt.SortedId(ipoint)) THEN

     ipoint = ipoint - Noffset

     Noffset = INT(Noffset/2)

    else

     ipoint = ipoint + Noffset

     Noffset = INT(Noffset/2)

    endif

   enddo



    kn=0

    if (SortedSire(j)==SortedId(ipoint)) kn=1

    do while (ipoint<nobs .and. kn==0 .and. SortedSire(j) > SortedId(ipoint))

     ipoint=ipoint+1

    enddo

    if (SortedSire(j)==SortedId(ipoint)) kn=1

    do while (ipoint>1 .and. kn==0 .and. SortedSire(j) < SortedId(ipoint))

     ipoint=ipoint-1

    enddo

    if (SortedSire(j)==SortedId(ipoint)) kn=1



    IF(kn==1) then

     SortedSireIndex(j) = SortedIdIndex(ipoint)

    else    ! sire is unlisted base sire

     newsires = newsires + 1

     SortedSireIndex(j) = nobs + newsires ! for now

    endif

end do !j



 ALLOCATE  (holdsireid(newsires))

 kn=0

 do j=1,nsires

  if (SortedSireIndex(j) > nobs) then

   kn=kn+1

   holdsireid(SortedSireIndex(j)-nobs) = SortedSire(j)

  end if

 enddo

 IF(kn /= newsires) stop'newsires error'





!PRINT*,  '  Find seqsire ... '



do j = 1, nobs



  If (sire(j) == '0') Then

    seqsire(j)=0

  else



    ipoint=INT(nsires/2)

    Noffset = INT(ipoint/2)

   do while (Noffset>1)

    IF (Sire(j).lt.SortedSire(ipoint)) THEN

     ipoint = ipoint - Noffset

     Noffset = INT(Noffset/2)

    else

     ipoint = ipoint + Noffset

     Noffset = INT(Noffset/2)

    endif

   enddo



    kn=0

    if (Sire(j)==SortedSire(ipoint)) kn=1

    do while (ipoint<nsires .and. kn==0 .and. Sire(j) > SortedSire(ipoint))

     ipoint=ipoint+1

    enddo

    if (Sire(j)==SortedSire(ipoint)) kn=1

    do while (ipoint>1 .and. kn==0 .and. Sire(j) < SortedSire(ipoint))

     ipoint=ipoint-1

    enddo

    if (Sire(j)==SortedSire(ipoint)) kn=1

    IF(kn==1) then

     seqsire(j) = SortedSireIndex(ipoint)

    else

     !PRINT*, ' Error: Sire missing: ', Sire(j)

     stop

    endif



  endif



ENDDO !j



!PRINT*,  '  Sires: ',newsires,' unlisted, ',nsires,' in total'



!PRINT*,  ' Females ... '

!PRINT*,  '  Find or set dam indices ... '



newdams = 0

nbisexuals = 0



do j=1,ndams



! check if already listed as an individual

   ipoint=INT(nobs/2)

   Noffset = INT(ipoint/2)

   do while (Noffset>1)

    IF (Sorteddam(j).lt.SortedId(ipoint)) THEN

     ipoint = ipoint - Noffset

     Noffset = INT(Noffset/2)

    else

     ipoint = ipoint + Noffset

     Noffset = INT(Noffset/2)

    endif

   enddo



    kn=0

    if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint  ! store ipoint here as ipoint can change with bisexuals

    do while (ipoint<nobs .and. kn==0 .and. Sorteddam(j) > SortedId(ipoint))

     ipoint=ipoint+1

    enddo

    if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint

    do while (ipoint>1 .and. kn==0 .and. Sorteddam(j) < SortedId(ipoint))

     ipoint=ipoint-1

    enddo

    if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint

! check if already listed as a sire (and therefore bisexual)

   ipoint=INT(nsires/2)

   Noffset = INT(ipoint/2)

   do while (Noffset>1)

    IF (SortedDam(j).lt.SortedSire(ipoint)) THEN

     ipoint = ipoint - Noffset

     Noffset = INT(Noffset/2)

    else

     ipoint = ipoint + Noffset

     Noffset = INT(Noffset/2)

    endif

   enddo



    kb=0

    if (SortedDam(j)==SortedSire(ipoint)) kb=1

    do while (ipoint<nsires .and. kb==0 .and. SortedDam(j) > SortedSire(ipoint))

     ipoint=ipoint+1

    enddo

    if (SortedDam(j)==SortedSire(ipoint)) kb=1

    do while (ipoint>1 .and. kb==0 .and. SortedDam(j) < SortedSire(ipoint))

     ipoint=ipoint-1

    enddo

    if (SortedDam(j)==SortedSire(ipoint)) kb=1



    IF(kb==1) then

      nbisexuals = nbisexuals + 1

      open (1,FILE='bisex.txt',position = 'append')

       WRITE(1,*) SortedDam(j)

      close(1)

    endif



    if (kb==1) then

     SorteddamIndex(j) = SortedSireIndex(ipoint)

    elseif (kn>=1) then

     SorteddamIndex(j) = SortedIdIndex(kn)

    else    ! dam is unlisted base dam

     newdams = newdams + 1

     SorteddamIndex(j) = nobs + newsires + newdams ! for now

    endif



end do !j



If (nbisexuals > 0)  PRINT*, nbisexuals,' bisexual parent(s) found. See file bisex.txt.  <------------ WARNING !!!'



 ALLOCATE  (holddamid(newdams))

 kn=0

 do j=1,ndams

  if (SortedDamIndex(j) > nobs+newsires) then

   kn=kn+1

   holddamid(SortedDamIndex(j)-nobs-newsires) = SortedDam(j)

  end if

 enddo

 IF(kn /= newdams) stop'newdams error'







!PRINT*,  '  Find seqdam ... '



do j = 1, nobs



  If (dam(j) == '0') Then

    seqdam(j)=0

  else



    ipoint=INT(ndams/2)

    Noffset = INT(ipoint/2)

   do while (Noffset>1)

    IF (dam(j).lt.Sorteddam(ipoint)) THEN

     ipoint = ipoint - Noffset

     Noffset = INT(Noffset/2)

    else

     ipoint = ipoint + Noffset

     Noffset = INT(Noffset/2)

    endif

   enddo



    kn=0

    if (dam(j)==Sorteddam(ipoint)) kn=1

    do while (ipoint<ndams .and. kn==0 .and. dam(j) > Sorteddam(ipoint))

     ipoint=ipoint+1

    enddo

    if (dam(j)==Sorteddam(ipoint)) kn=1

    do while (ipoint>1 .and. kn==0 .and. dam(j) < Sorteddam(ipoint))

     ipoint=ipoint-1

    enddo

    if (dam(j)==Sorteddam(ipoint)) kn=1

    IF(kn==1) then

     seqdam(j) = SorteddamIndex(ipoint)

    else

    ! PRINT*, ' Error: dam missing: ', dam(j)

     stop

    endif



  endif



ENDDO !j



!PRINT*,  '  Dams: ',newdams,' unlisted, ',ndams,' in total'



!PRINT*,  ' Arranging unlisted base parents ... '





iextra = newsires + newdams



If (iextra > 0) then

     ! PRINT*, ' ', iextra, ' unlisted base parents found.'



 ! SortedId and SortedIdIndex just used as a holder while redimensioning



 SortedId(1:nobs)=id(1:nobs)

 deallocate (id)

 ALLOCATE(id(nobs+iextra))

 id(1+iextra:nobs+iextra)=SortedId(1:nobs)



 SortedId(1:nobs)=sire(1:nobs)

 deallocate (sire)

 ALLOCATE(sire(nobs+iextra))

 sire(1+iextra:nobs+iextra)=SortedId(1:nobs)



 SortedId(1:nobs)=dam(1:nobs)

 deallocate (dam)

 ALLOCATE(dam(nobs+iextra))

 dam(1+iextra:nobs+iextra)=SortedId(1:nobs)



 SortedIdIndex(1:nobs)=seqsire(1:nobs)

 deallocate (seqsire)

 ALLOCATE(seqsire(nobs+iextra))

 seqsire(1+iextra:nobs+iextra)=SortedIdIndex(1:nobs)



 SortedIdIndex(1:nobs)=seqdam(1:nobs)

 deallocate (seqdam)

 ALLOCATE(seqdam(nobs+iextra))

 seqdam(1+iextra:nobs+iextra)=SortedIdIndex(1:nobs)



endif



!PRINT*, ' Inserting unlisted base parents ...'



oldnobs = nobs

nobs = nobs + iextra

!PRINT*, ' Total number of animals = ',nobs



ALLOCATE (passedorder(nobs))

passedorder=0



do i = 1+iextra, nobs

 passedorder(i)= i-iextra



 If (sire(i) == '0')then

   seqsire(i) = 0

 Else

   seqsire(i) = iextra + seqsire(i)

   If (seqsire(i) > nobs)  seqsire(i) = seqsire(i) - nobs  ! for unlisted sires

 End If



 If (dam(i) == '0') Then

   seqdam(i) = 0

  Else

   seqdam(i) = iextra + seqdam(i)

   If (seqdam(i) > nobs)  seqdam(i) = seqdam(i) - nobs

  End If

ENDDO !i



do i = 1, newsires

 ID(i) = holdsireid(i)

 passedorder(i)=0

 seqsire(i) = 0

 seqdam(i) = 0

ENDDO !i

do i = newsires + 1, newsires + newdams

 ID(i) = holddamid(i - newsires)

 passedorder(i)=0

 seqsire(i) = 0

 seqdam(i) = 0

ENDDO !i



DEALLOCATE(holdsireid, holddamid, SortedIdIndex, SortedId)



flag = 0

Do i = 1, nobs

If (i <= seqsire(i) .Or. i <= seqdam(i) ) flag = 1

enddo !i

!If (flag == 0) !PRINT*, 'not needed'!return



!PRINT*, ' Re-Ordering pedigree ...'





Allocate ( OldN(0:nobs), NewN(0:nobs) )

ALLOCATE ( holdid(0:nobs), holdsire(nobs), holddam(nobs) )



OldN(0) = 0

NewN=0

!seqsire(0) = 0 !not needed !

!seqdam(0) = 0



holdid(1:nobs) = ID(1:nobs)

holdsire = seqsire

holddam = seqdam



!Find base ancestors ...

kn = 0

do i = 1, nobs

 If (seqsire(i) == 0 .And. seqdam(i) == 0) Then

      kn = kn + 1

      NewN(i) = kn

      OldN(kn) = i

 End If

ENDDO !i



!Re-order pedigree ...

NewN(0) = nobs + 1

flag = 0

Do While (kn < nobs)

 oldkn = kn

 do i = 1, nobs

  If (NewN(i) == 0) Then !And ID(i) <> 'UniqueNULL' Then

    Ks = seqsire(i)

    Kd = seqdam(i)

    If (NewN(Ks) > 0 .And. NewN(Kd) > 0) Then

      kn = kn + 1

      NewN(i) = kn

      OldN(kn) = i

    End If

  End If

 enddo !i



 ! to avoid hang on unexpected problem ...

 If (kn == oldkn) Then

  flag = flag + 1

 Else

  flag = 0

 endif



 If (flag > 10) Then

   open(1,file='ped_err.txt',status='unknown')

   write(1,*) 'Pedigree errors found involving two or more of the following relationships ...'

   write(1,*)

   write(1,*) '       Index numbers are followed by names.'

   write(1,*) '       Index number 0 means unknown, whence name is blank.'

   write(1,*)

   do i = 1, nobs

    If (NewN(i) == 0) Then

     write(1,*) 'Individual:',          i, ':  ', ID(i)

     write(1,*) '    Father:', seqsire(i), ':  ', ID(seqsire(i))

     write(1,*) '    Mother:',  seqdam(i), ':  ', ID(seqdam(i))

     write(1,*)

    End If

   ENDDO !i

   Close (1)

   PRINT*,  'Logical error when re-ordering pedigree - see details in file PED_ERR.TXT'

   stop

 End If

ENDDO !while



NewN(0) = 0



do i = 1, nobs

 ID(i) = holdid(OldN(i))

enddo



do i = 1, nobs

seqsire(i) = NewN(holdsire(OldN(i)))

seqdam(i) = NewN(holddam(OldN(i)))

If (i <= NewN(holdsire(OldN(i))) .Or. i <= NewN(holddam(OldN(i)))) then

   !PRINT*,  'out of order'

   stop

endif

ENDDO !i



DO i = 1, nobs

  holdsire(i) = passedorder(i)  ! holdsire just because it is free

enddo

DO i = 1, nobs

  passedorder(i) = holdsire(OldN(i))

enddo



deallocate ( OldN, NewN, holdid, holdsire, holddam) ! holdrec)



!do i = 1, nobs

! PRINT'(3i5,2x,3a4,i5)', i, seqsire(i), seqdam(i), id(i), sire(i), dam(i), passedorder(i)

!enddo

nAnisPedigree=nObs

GlobalExtraAnimals=iextra   !Change John Hickey


end subroutine PVseq

!#############################################################################################################################################################################################################################
