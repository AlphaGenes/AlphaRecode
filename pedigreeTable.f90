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
    type(pedigreeLine), parameter :: DICT_NULL = pedigreeLine('','','', huge(lengan))

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
