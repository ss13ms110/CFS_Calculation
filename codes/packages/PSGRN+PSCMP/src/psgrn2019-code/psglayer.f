      subroutine psglayer(ierr)
      use psgalloc
      implicit none
      integer*4 ierr
c
      integer*4 l,n,li,lp0
      real*8 x1
c
      real*8, allocatable:: x(:),x0(:)
c
      lp0=n0+2
c
      allocate(x(lp0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psglayer: x not allocated!'
      allocate(x0(lp0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psglayer: x0 not allocated!'
c
      lp0=1
      x0(lp0)=0.d0
      do n=1,n0-1
        lp0=lp0+1
        x0(lp0)=x0(lp0-1)+h(n)
      enddo
      lp0=lp0+1
      x0(lp0)=zrec
      lp0=lp0+1
      x0(lp0)=zs
c
c     sort the z0-profile
c
      do l=1,lp0-1
        do li=l+1,lp0
          if(x0(li).lt.x0(l))then
            x1=x0(l)
            x0(l)=x0(li)
            x0(li)=x1
          endif
        enddo
      enddo
c
c     delete duplicates
c
      lp=1
      x(lp)=0.d0
      do l=2,lp0
        if(x0(l).gt.x(lp))then
          hp(lp)=x0(l)-x(lp)
          lp=lp+1
          x(lp)=x0(l)
        endif
      enddo
      hp(lp)=0.d0
c
c     determine ls,lzrec
c
      do l=1,lp
        if(x(l).eq.zs)ls=l
        if(x(l).eq.zrec)lzrec=l
      enddo
c
c     determine layer no of each depth
c
      li=1
      x1=h(1)
      nno(1)=1
      do l=2,lp
        if(x(l).ge.x1.and.li.lt.n0)then
          li=li+1
          x1=x1+h(li)
        endif
        nno(l)=li
      enddo
c
      ierr=0
c
      deallocate(x,x0)
c
      return
      end
