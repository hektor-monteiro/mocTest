code = """
do i = 1, grid%nx
   do j = 1, yTop
      do k = 1, grid%nz
         read(77, *) grid%xAxis(i), grid%yAxis(j), grid%zAxis(k), HdenTemp(i,j,k)
      end do
   end do
end do
"""
print("Analyzed.")
