      SUBROUTINE amux_msr(n, x, y, va, da, ja, ia)
c     AMUX_MSR - matrix-vector product
c        computes y <- A * x, where a is a matrix stored in MSR format
      IMPLICIT NONE

      INTEGER n
      REAL*8 x(*), y(*), va(*), da(*)
      INTEGER ja(*), ia(*)

      REAL*8 s, v, xi
      INTEGER i, j, k

      DO k = 1, n
         y(k) = da(k) * x(k)
      END DO

      DO i = 1, n
         xi = x(i)
         s = 0d0
         DO k = ia(i)+1, ia(i+1)
            j = ja(k) + 1
            v = va(k)
            s = s + v * x(j)
            y(j) = y(j) + v  * xi
         ENDDO
         y(i) = y(i) + s
      ENDDO

      END

      SUBROUTINE amux_csr(n, x, y, va, ja, ia)
c     AMUX_CSR - matrix-vector product
c        computes y <- A * x, where a is a matrix stored in CSR format
      IMPLICIT NONE

      INTEGER n
      REAL*8 x(*), y(*), va(*)
      INTEGER ja(*), ia(*)

      REAL*8 s
      INTEGER i, k

      DO  i = 1,n
         s = 0d0
         DO k = ia(i), ia(i+1)-1
            s = s + va(k)*x(ja(k))
         ENDDO
         y(i) = s
      ENDDO
      END

      SUBROUTINE atmux_csr(m, n, x, y, va, ja, ia)
c     ATMUX_CSR - matrix-vector product with transposed matrix
c        computes y <- A^T * x, where a is a m-by-n matrix stored in CSR format
      IMPLICIT NONE

      INTEGER m, n
      REAL*8 x(*), y(*), va(*)
      INTEGER ia(*), ja(*)

      INTEGER i, k

      DO i = 1, n
         y(i) = 0d0
      ENDDO

      DO i = 1, m
         DO k = ia(i), ia(i+1)-1
            y(ja(k)) = y(ja(k)) + x(i)*va(k)
         ENDDO
      ENDDO

      RETURN
      END
