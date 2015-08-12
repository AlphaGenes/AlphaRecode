module priority_queue_mod
use PedigreeTable
implicit none
 

 
type queue
  type(pedigreeLine), allocatable :: buf(:)
  integer :: n = 0 !Length of queue
contains
  procedure :: top
  procedure :: enqueue
  procedure :: shift
end type
 
contains
 
subroutine shift(this, a)
  class (queue)           :: this
  integer                 :: a, parent, child
  associate (x => this%buf)
  parent = a
  do while(parent*2 <= this%n)
    child = parent*2
    if (child + 1 <= this%n) then 
      if (x(child+1)%generation < x(child)%generation ) then
        child = child +1 
      end if
    end if
    if (x(parent)%generation > x(child)%generation) then
      x([child, parent]) = x([parent, child])
      parent = child
    else
      exit
    end if  
  end do      
  end associate
end subroutine
 
function top(this) result (res)
  class(queue) :: this
  type(pedigreeLine)   :: res
  res = this%buf(1)
  this%buf(1) = this%buf(this%n)
  this%n = this%n - 1
  call this%shift(1)
end function
 
subroutine enqueue(this, anim)
  class(queue), intent(inout) :: this
  type(pedigreeLine)          :: anim
  type(pedigreeLine), allocatable     :: tmp(:)
  integer                     :: i
  this%n = this%n +1  
  if (.not.allocated(this%buf)) allocate(this%buf(1))
  if (size(this%buf)<this%n) then
    allocate(tmp(2*size(this%buf)))
    tmp(1:this%n-1) = this%buf
    call move_alloc(tmp, this%buf)
  end if
  this%buf(this%n) = anim
  i = this%n
  do 
    i = i / 2
    if (i==0) exit
    call this%shift(i)
  end do
end subroutine
end module 