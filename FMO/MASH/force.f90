module force
    use params
    use functions
    implicit none

    contains

    subroutine pot(q, V)
        !! Diabatic potential matrix
        real(kind=wp1), intent(in) :: q(:)
        real(kind=wp1), intent(inout) :: V(nsys, nsys)
        V = Hamil_e

        call pot_add_bath(q, V)
        
    end subroutine pot

    subroutine pot_add_bath(q, V)
        !! reshapes q and avoids complications
        real(kind=wp1), intent(in) :: q(nbath, nsys)
        real(kind=wp1), intent(inout) :: V(nsys, nsys)
        integer :: ig, gerr

        real(kind=wp1) :: V0
        real(kind=wp1), allocatable :: V1(:)
        allocate(V1(nsys), stat=gerr)


        V0 = sum(matmul(mw2, q**2))*half
        V1 = matmul(ci, q)
        do ig = 1,nsys
            V(ig, ig) = V(ig, ig) + V0 + V1(ig)
        end do

        deallocate(V1, stat=gerr)

    end subroutine pot_add_bath

    !! get the gradients
    subroutine grad(q, G)
        real(kind=wp1), intent(in) :: q(:)
        real(kind=wp1), intent(out) :: G(nbath_total, nsys, nsys)

        call get_grad(q, G)
    end subroutine grad

    subroutine get_grad(q, G)
        !! reshape and compute gradient
        real(kind=wp1), intent(in) :: q(nbath, nsys)
        real(kind=wp1), intent(out) :: G(nbath, nsys, nsys, nsys)
        integer :: i,j

        G = rzero
        do i = 1,nsys
            do j = 1,nsys
                G(:,i,j,j) = mw2*q(:,i)
            end do
            G(:,i,i,i) = G(:,i,i,i) + ci
        end do
    end subroutine get_grad

    subroutine grad_adiabatic(q, U, Gad)
        !! use gradient , only : grad, get_grad
        real(kind=wp1), intent(in) :: q(:), U(:,:)
        real(kind=wp1), intent(out) :: Gad(nbath_total, nsys, nsys)
        real(kind=wp1), allocatable :: Gdia(:,:,:)
        integer :: i, ierror

        allocate(Gdia(nbath_total, nsys, nsys), stat=ierror)

        call grad(q, Gdia)
        do i=1,nbath_total
            Gad(i,:,:) = matmul(transpose(U), matmul(Gdia(i,:,:), U))
        end do
        deallocate(Gdia, stat=ierror)

    end subroutine grad_adiabatic

    subroutine grad_adiabatic_diagonal(q, U, a, dvdq)
        real(kind=wp1), intent(in) :: q(:), U(:,:)
        integer, intent(in) :: a
        real(kind=wp1), intent(out) :: dvdq(nbath_total)
        real(kind=wp1), allocatable :: Gdia(:,:,:)
        integer :: ii, ierror

        allocate(Gdia(nbath_total, nsys, nsys), stat=ierror)

        call grad(q, Gdia)

        !! mydot can be replaced with dot_product if blas is slow
        do ii=1,nbath_total
            dvdq(ii) = mydot(U(:,a), matmul(Gdia(ii,:,:), U(:,a)))
        end do
        deallocate(Gdia, stat=ierror)
    end subroutine grad_adiabatic_diagonal

    !! non adiabatic coupling
    subroutine nonAdCoup(q,Vad,U,d)
        !! with Hellman-Feynman theorem
        real(kind=wp1), intent(in) :: q(:), Vad(:), U(:,:)
        real(kind=wp1), intent(out) :: d(nbath_total, nsys, nsys)
        integer :: gi, i, j

        real(kind=wp1), allocatable :: Gad(:,:,:)
        allocate(Gad(nbath_total, nsys, nsys), stat=gi)

        call grad_adiabatic(q,U,Gad)
        
        do i=1,nsys
            d(:,i,i) = rzero
            do j=i+1,nsys
                d(:,i,j) = Gad(:,i,j)/(Vad(j)-Vad(i))
                d(:,j,i) = -d(:,i,j)
            end do
        end do

        deallocate(Gad, stat=gi)

    end subroutine nonAdCoup

    subroutine nonAdCoupDirection(q, cad, Vad, U, n, m, dj)
        !! rescale momentum or change directtion
        real(kind=wp1), intent(in) :: q(:), Vad(:), U(:,:)
        complex(kind=wp1), intent(in) :: cad(:)
        real(kind=wp1), intent(out) :: dj(nbath_total)
        integer, intent(in) :: n, m
        integer :: i, gi

        real(kind=wp1), allocatable :: d(:,:,:)
        allocate(d(nbath_total, nsys, nsys), stat=gi)

        call nonAdCoup(q, Vad, U, d)

        dj = rzero
        ! print*,n,m
        ! Eq E3 in ref 1
        do i=1,nsys
            dj = dj + d(:,i,n)*real(conjg(cad(i))*cad(n))
            dj = dj - d(:,i,m)*real(conjg(cad(i))*cad(m))
        end do

        deallocate(d, stat=gi)

    end subroutine nonAdCoupDirection
end module force