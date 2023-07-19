      double precision function force_underflow(x)
      double precision x
C
C     Avoid very small numbers that might need a three-character
C     exponent field in formatted output
C      
      if (abs(x) .lt. 1.0d-99) then
         force_underflow = 0.0d0
      else
         force_underflow = x
      endif

      end
