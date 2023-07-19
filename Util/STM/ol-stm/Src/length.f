      real*8 function length(a)
*
*     This function calculates the length of a vector defined in
*    cartesian coordinates.
*
      real*8 a(3)
*
*     a(3) : cartesian coordinates of the vector.  
*
      length = sqrt(a(1)**2 + a(2)**2 + a(3)**2)

      return
      end

