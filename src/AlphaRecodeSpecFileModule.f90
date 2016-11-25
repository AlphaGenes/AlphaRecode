module AlphaRecodeSpecFileModule



    type alphaRecodeSpecFileHolder
        character(len=256) :: pedigreeFile = "Pedigree.txt" 
        character(len=256) :: genderfile = "" 
        character(len=256) :: outputFile = "recoded.txt"


        !for pedigree extraction
        integer :: startGeneration
        integer :: endGeneration

    end type alphaRecodeSpecFileHolder


    interface alphaRecodeSpecFileHolder
        module procedure readInSpecFile
    end interface alphaRecodeSpecFileHolder

    contains

    function readInSpecFile(fileIn) result(input)
        use AlphaHouseMod, only: parseToFirstWhitespace,splitLineIntoTwoParts,toLower

        character(len=*), optional :: fileIn
        integer :: unit

        character(len=300) :: first, line
        character(len=:), allocatable::dumStr
        character(len=300),dimension(:),allocatable :: second
        character(len=256):: specFile
        integer :: stat

        type(alphaRecodeSpecFileHolder) :: input

        if (.not. present(fileIn)) then
            write(specFile, "(A)") "AlphaRecodeSpec.txt"
        else
            write(specFile, "(A)") fileIn
        end if

        open(newUnit=unit, file=specFile, action="read", status="old")

        READFILE: do while (stat==0)
            read(unit,"(A)", IOStat=stat)  line
            if (len_trim(line)==0) then
                CYCLE
            end if

            call splitLineIntoTwoParts(trim(line), first, second)

            dumStr = parseToFirstWhitespace(first)
            if (first(1:1)=="=" .or. len(trim(line))==0) then
              cycle
            else

                select case(trim(dumStr))
                    case("pedigreefilepath")
                        if (.not. allocated(second)) then
                          write(*, "(A,A)") "No pedigreeFilePath  found. Using default filename: ", input%pedigreeFile
                        else
                            print *, second(1)
                          write(input%pedigreeFile, *) trim(second(1))
                        end if

                    case("genderfilepath")
                        if (.not. allocated(second)) then
                          write(*, "(A,A)") "No genderfile  found. Using default filename: ", input%genderfile
                        else
                          write(input%genderfile, *) trim(second(1))
                        end if

                    case("outputfilepath")
                         if (.not. allocated(second)) then
                          write(*, "(A,A)") "No genderfile  found. Using default filename: ", input%outputFile
                        else
                          write(input%outputFile, *) trim(second(1))
                        end if

                    case("startgeneration")
                        read(second(1),*) input%startGeneration


                    case("endgeneration")
                        read(second(1),*) input%endGeneration

                    case default
                        write(*,"(A,A)") trim(dumStr), " is not a valid input for macs"
                    cycle
                end select
            end if
            end do READFILE



end function readInSpecFile




end module AlphaRecodeSpecFileModule