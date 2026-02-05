module fftpack_kinds
   implicit none(type, external)
   private
   public :: sp, dp

   !> Single precision real numbers.
   integer, parameter :: sp = selected_real_kind(6)

   !> Double precision real numbers.
   integer, parameter :: dp = selected_real_kind(8)
end module fftpack_kinds
